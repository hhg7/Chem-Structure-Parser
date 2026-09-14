#ifndef _GNU_SOURCE
#define _GNU_SOURCE // glibc / Linux
#endif
#ifndef __EXTENSIONS__
#define __EXTENSIONS__ 1 // Solaris/illumos
#endif
#define PERL_NO_GET_CONTEXT
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "ppport.h"
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
/*Chem::Structure::Parser -- the parts of reading a PDB file worth doing in C.

A PDB file is one record per line with every field at a fixed column range
(the format is punched-card old).  Two things about that make XS pay for
itself here, and only two: the coordinate section is enormous -- the largest
file in PDBbind v2020 is 33 MB, ~400,000 ATOM/HETATM lines -- and every one
of those lines wants a dozen substr()s and three numeric conversions, which
in Perl is a dozen SVs per line born and buried; and residue names are looked
up constantly (three letters -> one letter, "amino acid, nucleotide or
water?"), a switch on three bytes in C against a hash lookup plus a sub call
in Perl.

So C does exactly that: one pass over the bytes, ATOM/HETATM split into
parallel column arrays (one array per field, not one hash per atom, which
keeps peak memory down by roughly an order of magnitude), residue boundaries
marked as it goes, every other record handed back as raw lines grouped by
record name.  The header records (HEADER, COMPND, REMARK, SEQRES, ...) are
parsed in Perl: a few dozen of them, irregular, and the part most likely to
need fixing later.*/

/*CSP_RESTRICT -- "no other pointer here reaches this object".  Worth spelling
out because perl's ccflags carry -fno-strict-aliasing, which switches off the
type-based aliasing analysis a C compiler would otherwise use; restrict is one
of the few aliasing facts still available to it, and the parse below is full of
the pattern it wants -- a const char * into the file buffer read alongside a
char * being written, both the same type, neither able to touch the other.

Not spelled `restrict' directly: that is C99, and the perls this builds against
are configured with whatever -std their Configure chose (5.10 on an old gcc gets
gnu89).  GCC and clang take __restrict__ in any mode, MSVC takes __restrict,
and a compiler with neither loses the hint and nothing else.  Predefine it empty
-- OPTIMIZE='-O2 -DCSP_RESTRICT=' -- to build without the hint, which is how the
before/after timings were taken and the escape hatch if a compiler mishandles
it.*/
#ifndef CSP_RESTRICT
#if defined(__cplusplus)
#  if defined(__GNUC__) || defined(_MSC_VER)
#    define CSP_RESTRICT __restrict
#  else
#    define CSP_RESTRICT
#  endif
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
#  define CSP_RESTRICT restrict
#elif defined(__GNUC__) && __GNUC__ >= 3
#  define CSP_RESTRICT __restrict__
#elif defined(_MSC_VER) && _MSC_VER >= 1400
#  define CSP_RESTRICT __restrict
#else
#  define CSP_RESTRICT
#endif
#endif

/*Width-correct libm for NV.

C has no type-generic <math.h>: sqrt(), cos() and acos() are declared to take a
double, so calling one with an NV on a long-double or __float128 build converts
the argument down, computes at 53 bits of mantissa and converts back.  Nothing
warns and nothing fails to compile; the answer is simply less accurate than the
perl running it.  So every NV-valued libm call in this file goes through the
nv_* macros below, which paste on the suffix for the width NV actually is: none
for double, `l' for long double, `q' for __float128 (perl.h has already included
<quadmath.h> by this point, and the quadmath perl's $Config{perllibs} already
carries -lquadmath, so the q functions need nothing at link time).

The long-double row is conditional because the `l' variants are C99 but not
universally present -- some BSD libms are thin on them.  Makefile.PL link-tests
the set and defines CSP_HAVE_LONG_DOUBLE_MATH only when every one resolves;
without it this falls back to the double functions, which costs accuracy on a
long-double perl and nothing else.

This layer is the same one Stats::LikeR carries, cut down to the six functions
the feature calculations below need.  Adding a sixth means adding it to the
Makefile.PL probe in the same edit, or the probe passes on a libm that does not
have it and the build breaks on a machine nobody here owns.*/
#if defined(USE_QUADMATH)
#  define CSP_NVFN(base) base ## q
#elif defined(USE_LONG_DOUBLE) && defined(CSP_HAVE_LONG_DOUBLE_MATH)
#  define CSP_NVFN(base) base ## l
#else
#  define CSP_NVFN(base) base
#endif
#define nv_fabs(x) CSP_NVFN(fabs)(x)
#define nv_sqrt(x) CSP_NVFN(sqrt)(x)
#define nv_acos(x) CSP_NVFN(acos)(x)
#define nv_cos(x)  CSP_NVFN(cos)(x)
#define nv_sin(x)  CSP_NVFN(sin)(x)
#define nv_atan2(y,x) CSP_NVFN(atan2)((y),(x))

/*pi at the full width of an NV.  A literal cannot supply this: an unsuffixed C
floating constant is a double, so spelling the digits out would round pi to 53
bits on a long-double or __float128 build, and the `L' / `Q' suffixes that would
avoid that are per-build (and `Q' is a GNU extension).  acos(-1) is pi correctly
rounded to the working width on every libm, and nv_acos() picks the right
width's acos.  M_PI is not used: it is not in C89 or C99, MSVC hides it behind
_USE_MATH_DEFINES, and where it exists it is a double.*/
#define CSP_PI nv_acos((NV)-1.0)

/*residue table.  One table serves aa3to1(), res1() and res_type(), so the
three can never disagree about what a residue is.  Keyed on the three-byte
name packed into a U32, left-padded with blanks, so that "ALA", " DA" and
"  A" are all one switch on an integer.*/
#define RT_OTHER 0 //a ligand, an ion, a sugar: not in the table
#define RT_AA    1 //amino acid, standard or modified
#define RT_NUC   2 //nucleotide, DNA or RNA
#define RT_WATER 3 //HOH, DOD, WAT, and the rest of the spellings

/*U32 and not unsigned: the key is three bytes wide and the top one is shifted
16 places, which plain unsigned is only guaranteed 16 bits for.*/
#define K3(a, b, c) (((U32)(a) << 16) | ((U32)(b) << 8) | (U32)(c))

typedef struct {
	char one;           //single-letter code; '\0' when there isn't one
	unsigned char type; //RT_*
} res_info;

//pack a (possibly blank-padded, possibly short) residue name into a key
static U32 res_key(const char *CSP_RESTRICT s, STRLEN len)
{
	char b[3];
	STRLEN i;
	while (len && (*s == ' ' || *s == '\t')) { s++; len--; }
	while (len && (s[len - 1] == ' ' || s[len - 1] == '\t')) len--;
	if (len == 0 || len > 3) return 0;
	b[0] = b[1] = b[2] = ' ';
	for (i = 0; i < len; i++) b[3 - len + i] = (char)toupper((unsigned char)s[i]);
	return K3((unsigned char)b[0], (unsigned char)b[1], (unsigned char)b[2]);
}

/*res_lookup() -- classify a residue name.  True when the name is known, false
when it is not (a ligand, an ion, a sugar: anything the caller has to work out
for itself from the atoms).*/
static bool res_lookup(const char *CSP_RESTRICT s, STRLEN len, res_info *CSP_RESTRICT out)
{
	char one;
	unsigned char type = RT_AA;
	switch (res_key(s, len)) {
		//the twenty, plus the ambiguity codes the format allows
		case K3('A','L','A'): one = 'A'; break;
		case K3('A','R','G'): one = 'R'; break;
		case K3('A','S','N'): one = 'N'; break;
		case K3('A','S','P'): one = 'D'; break;
		case K3('C','Y','S'): one = 'C'; break;
		case K3('G','L','N'): one = 'Q'; break;
		case K3('G','L','U'): one = 'E'; break;
		case K3('G','L','Y'): one = 'G'; break;
		case K3('H','I','S'): one = 'H'; break;
		case K3('I','L','E'): one = 'I'; break;
		case K3('L','E','U'): one = 'L'; break;
		case K3('L','Y','S'): one = 'K'; break;
		case K3('M','E','T'): one = 'M'; break;
		case K3('P','H','E'): one = 'F'; break;
		case K3('P','R','O'): one = 'P'; break;
		case K3('S','E','R'): one = 'S'; break;
		case K3('T','H','R'): one = 'T'; break;
		case K3('T','R','P'): one = 'W'; break;
		case K3('T','Y','R'): one = 'Y'; break;
		case K3('V','A','L'): one = 'V'; break;
		case K3('A','S','X'): one = 'B'; break; //ASX ASP or ASN
		case K3('G','L','X'): one = 'Z'; break; //GLX GLU or GLN
		case K3('X','L','E'): one = 'J'; break; //XLE LEU or ILE
		case K3('S','E','C'): one = 'U'; break; //SEC selenocysteine
		case K3('P','Y','L'): one = 'O'; break; //PYL pyrrolysine
		case K3('U','N','K'): case K3('X','A','A'): one = 'X'; break;

		/*modified residues, mapped to the parent they were made from.  A
		structure that soaked in selenomethionine is still the same sequence,
		and a sequence with an 'X' every seventh position is no use to anyone.*/
		case K3('M','S','E'): case K3('M','H','O'): case K3('F','M','E'): case K3('C','X','M'):
		case K3('S','M','E'): case K3('M','E','D'): one = 'M'; break; //MSE selenomethionine
		case K3('C','S','O'): case K3('C','S','D'): case K3('C','S','S'): case K3('C','S','X'):
		case K3('C','S','W'): case K3('C','M','E'): case K3('C','M','T'): case K3('C','Y','X'):
		case K3('C','A','S'): case K3('C','A','F'): case K3('O','C','S'): case K3('S','M','C'):
		case K3('S','N','C'): case K3('Y','C','M'): one = 'C'; break;
		case K3('S','E','P'): case K3('S','A','C'): one = 'S'; break; //SEP phosphoserine
		case K3('T','P','O'): one = 'T'; break; //TPO phosphothreonine
		case K3('P','T','R'): case K3('T','Y','S'): case K3('T','Y','I'): case K3('T','Y','Q'):
		case K3('T','P','Q'): case K3('P','A','Q'): case K3('S','T','Y'): case K3('I','Y','R'):
			one = 'Y'; break; //PTR phosphotyrosine
		case K3('K','C','X'): case K3('L','L','P'): case K3('M','L','Y'): case K3('M','L','Z'):
		case K3('M','3','L'): case K3('A','L','Y'): case K3('L','Y','Z'): one = 'K'; break;
		case K3('H','Y','P'): case K3('H','Y','3'): one = 'P'; break; //HYP hydroxyproline
		//PCA pyroglutamate
		case K3('P','C','A'): case K3('C','G','U'): case K3('G','M','A'): one = 'E'; break;
		case K3('H','I','C'): case K3('H','I','D'): case K3('H','I','E'): case K3('H','I','P'):
		case K3('H','S','D'): case K3('H','S','E'): case K3('H','S','P'): case K3('M','H','S'):
		case K3('N','E','P'): one = 'H'; break;
		case K3('A','I','B'): case K3('A','B','A'): case K3('A','L','M'): case K3('A','Y','A'):
		case K3('B','A','L'): one = 'A'; break;
		case K3('S','A','R'): one = 'G'; break; //SAR sarcosine
		case K3('M','L','E'): case K3('N','L','E'): one = 'L'; break;
		case K3('M','V','A'): one = 'V'; break;
		case K3('C','I','R'): one = 'R'; break; //CIR citrulline
		case K3('O','R','N'): case K3('D','A','B'): one = 'X'; break; //ORN ornithine
		case K3('T','R','O'): one = 'W'; break;
		case K3('P','H','I'): case K3('P','H','L'): case K3('M','E','A'): one = 'F'; break;

		/*D-amino acids -- common in the peptide ligands this was written for,
		and the same letter as their L partner*/
		case K3('D','A','L'): one = 'A'; break;
		case K3('D','A','R'): one = 'R'; break;
		case K3('D','S','G'): one = 'N'; break;
		case K3('D','A','S'): one = 'D'; break;
		case K3('D','C','Y'): one = 'C'; break;
		case K3('D','G','N'): one = 'Q'; break;
		case K3('D','G','L'): one = 'E'; break;
		case K3('D','H','I'): one = 'H'; break;
		case K3('D','I','L'): one = 'I'; break;
		case K3('D','L','E'): one = 'L'; break;
		case K3('D','L','Y'): one = 'K'; break;
		case K3('D','P','N'): one = 'F'; break;
		case K3('D','P','R'): one = 'P'; break;
		case K3('D','S','N'): one = 'S'; break;
		case K3('D','T','H'): one = 'T'; break;
		case K3('D','T','R'): one = 'W'; break;
		case K3('D','T','Y'): one = 'Y'; break;
		case K3('D','V','A'): case K3('D','I','V'): one = 'V'; break;

		//nucleotides: DNA (" DA"), RNA ("  A"), and the pre-v3 spellings
		case K3(' ','D','A'): one = 'A'; type = RT_NUC; break;
		case K3(' ','D','C'): one = 'C'; type = RT_NUC; break;
		case K3(' ','D','G'): one = 'G'; type = RT_NUC; break;
		case K3(' ','D','T'): one = 'T'; type = RT_NUC; break;
		case K3(' ','D','U'): one = 'U'; type = RT_NUC; break;
		case K3(' ','D','I'): one = 'I'; type = RT_NUC; break;
		case K3(' ',' ','A'): one = 'A'; type = RT_NUC; break;
		case K3(' ',' ','C'): one = 'C'; type = RT_NUC; break;
		case K3(' ',' ','G'): one = 'G'; type = RT_NUC; break;
		case K3(' ',' ','T'): one = 'T'; type = RT_NUC; break;
		case K3(' ',' ','U'): one = 'U'; type = RT_NUC; break;
		case K3(' ',' ','I'): one = 'I'; type = RT_NUC; break;
		case K3(' ',' ','N'): one = 'N'; type = RT_NUC; break;
		case K3('A','D','E'): one = 'A'; type = RT_NUC; break;
		case K3('C','Y','T'): one = 'C'; type = RT_NUC; break;
		case K3('G','U','A'): one = 'G'; type = RT_NUC; break;
		case K3('T','H','Y'): one = 'T'; type = RT_NUC; break;
		case K3('U','R','I'): case K3('P','S','U'): case K3('H','2','U'): case K3('4','S','U'):
			one = 'U'; type = RT_NUC; break; //PSU pseudouridine
		case K3('5','M','U'): one = 'T'; type = RT_NUC; break;
		case K3('5','M','C'): case K3('O','M','C'): one = 'C'; type = RT_NUC; break;
		case K3('1','M','A'): one = 'A'; type = RT_NUC; break;
		case K3('2','M','G'): case K3('7','M','G'): case K3('1','M','G'): case K3('M','2','G'):
		case K3('O','M','G'): one = 'G'; type = RT_NUC; break;

		//water, under every name it gets written with
		case K3('H','O','H'): case K3('W','A','T'): case K3('D','O','D'): case K3('H','2','O'):
		case K3('S','O','L'): case K3('T','I','P'): one = '\0'; type = RT_WATER; break;

		default: return FALSE;
	}
	out->one  = one;
	out->type = type;
	return TRUE;
}

/*the reverse direction.  res_lookup() is many-to-one: sixty-odd names come
back as 'C'.  Going the other way there is one answer worth giving, the parent
residue the letter stands for, so this is a plain table indexed by the letter
-- 104 bytes of read-only memory, one bounds check and one index per call, and
no hash to build at load time.

Every letter A-Z is spoken for, because the ambiguity codes have names of
their own.  There is nothing here for the nucleotides -- 'A' is ALA on the way
back, and a caller who wants " DA" knows already that the chain is DNA and
does not need this to guess.

t/residue_names.t rounds every letter back through aa3to1() so that this table
and the switch above cannot drift apart.*/
static const char aa1to3_name[26][4] = {
	"ALA", "ASX", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", //A-I
	"XLE", "LYS", "LEU", "MET", "ASN", "PYL", "PRO", "GLN", "ARG", //J-R
	"SER", "THR", "SEC", "VAL", "TRP", "UNK", "TYR", "GLX"         //S-Z
};//B ASP or ASN, J LEU or ILE, O pyrrolysine, U selenocysteine, Z GLU or GLN

/*aa1to3_lookup() -- the name for a single-letter code, or NULL.  Blanks and
case are forgiven, as they are in res_key(), because the letter usually
arrives out of a sequence string that has been through something else.*/
static const char *aa1to3_lookup(const char *CSP_RESTRICT s, STRLEN len)
{
	int c;
	while (len && (*s == ' ' || *s == '\t')) { s++; len--; }
	while (len && (s[len - 1] == ' ' || s[len - 1] == '\t')) len--;
	if (len != 1) return NULL;
	c = toupper((unsigned char)*s);
	if (c < 'A' || c > 'Z') return NULL;
	return aa1to3_name[c - 'A'];
}

/*column fields.  Every accessor below clips to the line's real length first:
PDB files in the wild are right-trimmed, so a line can stop before the element
column (or before the B-factor) and that is not an error.*/

static STRLEN trim_ptr(const char **CSP_RESTRICT sp, STRLEN len)
{
	const char *s = *sp;
	while (len && (*s == ' ' || *s == '\t')) { s++; len--; }
	while (len && (s[len - 1] == ' ' || s[len - 1] == '\t')) len--;
	*sp = s;
	return len;
}

//fld() -- the trimmed contents of columns [from, to], zero-based inclusive
static STRLEN fld(const char *CSP_RESTRICT line, STRLEN llen, STRLEN from, STRLEN to,
                  const char **CSP_RESTRICT out)
{
	const char *s;
	STRLEN n;
	*out = line;
	if (from >= llen) return 0;
	if (to >= llen) to = llen - 1;
	s = line + from;
	n = to - from + 1;
	n = trim_ptr(&s, n);
	*out = s;
	return n;
}

/*str2iv() -- integers, exactly, without going through a double.  '*****'
(what a serial number becomes once it overflows five columns) and anything
else non-numeric report failure rather than a wrong number.*/
static bool str2iv(const char *CSP_RESTRICT s, STRLEN n, IV *CSP_RESTRICT out)
{
	bool neg = FALSE, seen = FALSE;
	IV v = 0;
	STRLEN i = 0;
	if (i < n && (s[i] == '+' || s[i] == '-')) { neg = (s[i] == '-'); i++; }
	for (; i < n; i++) {
		if (!isdigit((unsigned char)s[i])) return FALSE;
		v = v * 10 + (s[i] - '0');
		seen = TRUE;
	}
	if (!seen) return FALSE;
	*out = neg ? -v : v;
	return TRUE;
}

/*STR2NV_STRTOD() -- the strtod() that reads a whole NV.  A perl built with
-Duselongdouble or -Dusequadmath has an NV wider than a double, and parsing
into a double first loses the bits that the wider type would have kept: 0.60
becomes 0.599999999999999978 once perl prints it back at full NV precision.
Perl's own Strtod() picks the right one, but it only exists from 5.22 on, so
older perls get the choice made here.*/
#ifdef Strtod
#  define STR2NV_STRTOD(s, e) Strtod((s), (e))
#elif defined(USE_QUADMATH)
#  define STR2NV_STRTOD(s, e) strtoflt128((s), (e))
#elif defined(USE_LONG_DOUBLE) && defined(HAS_STRTOLD)
#  define STR2NV_STRTOD(s, e) strtold((s), (e))
#else
#  define STR2NV_STRTOD(s, e) strtod((s), (e))
#endif

/*str2nv_slow() -- whatever the fixed-point reader below declines: an exponent,
a hex float, inf/nan, a field with rubbish after the number.  strtod() on a
stack copy, because the fields are adjacent in the record (x ends where y
begins, with no separator when a coordinate is wide) and so strtod() cannot be
pointed at the record buffer itself.*/
static bool str2nv_slow(const char *CSP_RESTRICT s, STRLEN n, NV *CSP_RESTRICT out)
{
	/*64 is far above the eight columns a coordinate field has, and a longer
	string is declined below rather than truncated into the buffer.*/
	char buf[64], *end;
	NV v;
	if (n == 0 || n >= sizeof(buf)) return FALSE;
	memcpy(buf, s, n);
	buf[n] = '\0';
	errno = 0;
	v = STR2NV_STRTOD(buf, &end);
	if (end == buf) return FALSE;
	*out = v;
	return TRUE;
}

/*str2nv_fixed() -- [+-]?digits[.digits], which is every number a coordinate
section actually contains.  Worth a reader of its own because glibc's strtod()
was 22% of the instructions in a whole-file parse (____strtod_l_internal,
str_to_mpn, round_and_return): it is a fully general conversion -- locale, hex
floats, exponents, an arbitrary-precision slow path -- and none of that can
appear in a field eight columns wide.

It returns 1 only on consuming the whole field, so anything it does not
understand falls through to str2nv_slow() and keeps its old meaning.

What it produces is not an approximation of strtod(), it is the same NV bit for
bit, and the reasoning is worth keeping because a change here could silently
cost a digit.  Both operands of the division are exact: the mantissa is
accumulated as an integer in an NV, and every integer up to 2**53 is exact even
on the narrowest NV perl configures (a plain double), which the fifteen-digit
cap stays well inside; 10**frac is exact for every power of ten to 10**22.  An
IEEE 754 division of two exact operands is correctly rounded, and correctly
rounding the true decimal value is exactly what strtod() is required to return
-- so they agree, and agree the same way on a long double or __float128 perl.
Dividing by 10**frac rather than multiplying by 10**-frac is the part that has
to be this way round: 1e-3 is not representable, so multiplying would round
once into the reciprocal and again into the product.

Leading zeros count against the cap, which costs a little generality that no
eight-column field could use, and buys the invariant that keeps the table index
in range: frac never exceeds the digit count, so it never exceeds the cap.*/
#define STR2NV_MAX_DIGITS 15
static const NV str2nv_pow10[STR2NV_MAX_DIGITS + 1] = {
	1.0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7,
	1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15
};

static bool str2nv_fixed(const char *CSP_RESTRICT s, STRLEN n, NV *CSP_RESTRICT out)
{
	NV m = 0;
	STRLEN i = 0, nd = 0, frac = 0;
	bool neg = FALSE, dot = FALSE;

	if (n == 0) return FALSE;
	if (s[0] == '+' || s[0] == '-') { neg = (s[0] == '-'); i = 1; }
	for (; i < n; i++) {
		const char c = s[i];
		if (c >= '0' && c <= '9') {
			if (++nd > STR2NV_MAX_DIGITS) return FALSE;
			m = m * 10 + (NV)(c - '0');
			if (dot) frac++;
		} else if (c == '.' && !dot) {
			dot = TRUE;
		} else {
			return FALSE; //an exponent, a stray letter, the '*****' of an overflowed field
		}
	}
	if (nd == 0) return FALSE; //"", "-", ".", "+."
	if (frac) m /= str2nv_pow10[frac];
	*out = neg ? -m : m;
	return TRUE;
}

//str2nv() -- coordinates and the like
static bool str2nv(const char *CSP_RESTRICT s, STRLEN n, NV *CSP_RESTRICT out)
{
	return str2nv_fixed(s, n, out) || str2nv_slow(s, n, out);
}

/*The two above are reachable individually from Perl as _str2nv_paths(), which
exists for t/numbers.t and asserts that they return the same NV bit for bit.
There is no other way to compare them: only one of the two is ever reached for
any given field, and Perl's own string-to-NV conversion cannot referee the
question because on an older -Duselongdouble perl Perl_my_atof is a hand-rolled
decimal accumulator rather than strtod().  The XSUB is at the foot of the file
with the rest of them.*/

static bool fld_iv(const char *CSP_RESTRICT line, STRLEN llen, STRLEN from, STRLEN to,
                  IV *CSP_RESTRICT out)
{
	const char *s;
	STRLEN n = fld(line, llen, from, to, &s);
	return n ? str2iv(s, n, out) : 0;
}

static bool fld_nv(const char *CSP_RESTRICT line, STRLEN llen, STRLEN from, STRLEN to,
                  NV *CSP_RESTRICT out)
{
	const char *s;
	STRLEN n = fld(line, llen, from, to, &s);
	return n ? str2nv(s, n, out) : 0;
}

/*charge_ok() -- columns 79-80 hold a digit and a sign, '1+' as the format
writes it and '+1' as some programs do.  A file old enough to keep the entry id
in columns 73-80 has the tail of it here instead, and 'DR' is not a charge, nor
is '09'.  A field that is not one reads as the empty string a blank field would
have given, which is the answer for "the file does not say".*/
static bool charge_ok(const char *CSP_RESTRICT s, STRLEN n)
{
	STRLEN i, digits = 0, signs = 0;
	if (n == 0 || n > 2) return FALSE;
	for (i = 0; i < n; i++) {
		if (s[i] >= '0' && s[i] <= '9') digits++;
		else if (s[i] == '+' || s[i] == '-') signs++;
		else return FALSE;
	}
	return digits <= 1 && signs <= 1;
}

/*guess_element() -- only when columns 77-78 are absent, which happens in
files written before the element column existed and in files written by
programs that should know better.

The atom name is right-justified from column 14 for one-letter elements and
from column 13 for two-letter ones, so " CA " is a carbon alpha and "CA  " is
a calcium.  That rule is worth following exactly, because the alternative --
taking the first two letters -- turns every "HG11" hydrogen into mercury.*/
static STRLEN guess_element(const char *CSP_RESTRICT raw, STRLEN rawlen,
                            char *CSP_RESTRICT buf)
{
	if (rawlen == 0) return 0;
	/*A hydrogen with a long name -- HG11, HD22, HE21 -- fills all four columns
	and starts in column 13, exactly where a two-letter element starts.  Column
	position alone cannot separate the two, so take the full field as the
	tie-breaker: no two-letter element has a four character atom name, while
	nearly every hydrogen past the first does.  Without this, every HG11 in the
	file becomes mercury.*/
	if (rawlen >= 4 && raw[0] != ' ' && raw[1] != ' ' && raw[2] != ' ' && raw[3] != ' '
	    && (raw[0] == 'H' || raw[0] == 'D' || raw[0] == 'h' || raw[0] == 'd')) {
		buf[0] = (char)toupper((unsigned char)raw[0]);
		return 1;
	}
	/*a name that starts in column 13 with two letters is a two-letter element
	(FE, ZN, CL, MG); a digit there is a hydrogen count (1HB)*/
	if (raw[0] != ' ' && isalpha((unsigned char)raw[0])
	    && rawlen >= 2 && isalpha((unsigned char)raw[1])) {
		buf[0] = (char)toupper((unsigned char)raw[0]);
		buf[1] = (char)toupper((unsigned char)raw[1]);
		return 2;
	}
	{	//otherwise the first letter in the field is the element
		STRLEN i;
		for (i = 0; i < rawlen; i++) {
			if (isalpha((unsigned char)raw[i])) {
				buf[0] = (char)toupper((unsigned char)raw[i]);
				return 1;
			}
		}
	}
	return 0;
}

/*elem_case() -- the official spelling of an element symbol.

Columns 77-78 are uppercase in the PDB format, mmCIF's type_symbol is uppercased
on the way in to match, and guess_element() uppercases what it takes from the
atom name, so a symbol arrives here as 'ZN' whichever of the three read it.
IUPAC writes it Zn, and a caller tallying elements wants the symbol rather than
the file's shouted version of it.  Only the second letter is ever in question:
every symbol begins with a capital, and every one-letter symbol is that capital
by itself.

A field that spells no element is left exactly as it came, so the 'XX' of a
broken file stays 'XX' instead of turning into a plausible-looking 'Xx'.  D for
deuterium is left alone for the same reason -- the PDB writes it, IUPAC has no
symbol to correct it to, and the readers below test for it by name.

The list is the 118 named elements of IUPAC's 2021 periodic table (the red book
table, https://iupac.org/what-we-do/periodic-table-of-elements/), less the
fourteen whose symbol is a single letter -- H B C N O F P S K V Y I W U -- which
need nothing done to them.*/
#define K2(a, b) (((U32)(a) << 8) | (U32)(b))

static void elem_case(char *CSP_RESTRICT buf, STRLEN n)
{
	if (n != 2) return;
	switch (K2((unsigned char)buf[0], (unsigned char)buf[1])) {
		case K2('H','E'): case K2('L','I'): case K2('B','E'): case K2('N','E'):
		case K2('N','A'): case K2('M','G'): case K2('A','L'): case K2('S','I'):
		case K2('C','L'): case K2('A','R'): case K2('C','A'): case K2('S','C'):
		case K2('T','I'): case K2('C','R'): case K2('M','N'): case K2('F','E'):
		case K2('C','O'): case K2('N','I'): case K2('C','U'): case K2('Z','N'):
		case K2('G','A'): case K2('G','E'): case K2('A','S'): case K2('S','E'):
		case K2('B','R'): case K2('K','R'): case K2('R','B'): case K2('S','R'):
		case K2('Z','R'): case K2('N','B'): case K2('M','O'): case K2('T','C'):
		case K2('R','U'): case K2('R','H'): case K2('P','D'): case K2('A','G'):
		case K2('C','D'): case K2('I','N'): case K2('S','N'): case K2('S','B'):
		case K2('T','E'): case K2('X','E'): case K2('C','S'): case K2('B','A'):
		case K2('L','A'): case K2('C','E'): case K2('P','R'): case K2('N','D'):
		case K2('P','M'): case K2('S','M'): case K2('E','U'): case K2('G','D'):
		case K2('T','B'): case K2('D','Y'): case K2('H','O'): case K2('E','R'):
		case K2('T','M'): case K2('Y','B'): case K2('L','U'): case K2('H','F'):
		case K2('T','A'): case K2('R','E'): case K2('O','S'): case K2('I','R'):
		case K2('P','T'): case K2('A','U'): case K2('H','G'): case K2('T','L'):
		case K2('P','B'): case K2('B','I'): case K2('P','O'): case K2('A','T'):
		case K2('R','N'): case K2('F','R'): case K2('R','A'): case K2('A','C'):
		case K2('T','H'): case K2('P','A'): case K2('N','P'): case K2('P','U'):
		case K2('A','M'): case K2('C','M'): case K2('B','K'): case K2('C','F'):
		case K2('E','S'): case K2('F','M'): case K2('M','D'): case K2('N','O'):
		case K2('L','R'): case K2('R','F'): case K2('D','B'): case K2('S','G'):
		case K2('B','H'): case K2('H','S'): case K2('M','T'): case K2('D','S'):
		case K2('R','G'): case K2('C','N'): case K2('N','H'): case K2('F','L'):
		case K2('M','C'): case K2('L','V'): case K2('T','S'): case K2('O','G'):
			buf[1] = (char)tolower((unsigned char)buf[1]);
			break;
		default: break; //not a symbol: leave it as the file wrote it
	}
}

/*elem_bump() -- add one to a tally kept in a hash.

Written out rather than left to sv_inc() so that the count is a UV: these are
counted up from nothing and never down, an IV would promise a sign the tally
cannot have, and a Perl caller printing one gets the same digits either way.*/
static void elem_bump(pTHX_ HV *CSP_RESTRICT tally, const char *CSP_RESTRICT sym, STRLEN len)
{
	SV **cnt = hv_fetch(tally, sym, (I32)len, 1);
	if (cnt && *cnt) sv_setuv(*cnt, SvOK(*cnt) ? SvUV(*cnt) + 1 : 1);
}

/*sub_hv() -- the hash under a key, made on first use.  Two levels of these are
what the per-chain tallies hang from; the same shape as the meta table below,
which grows an array under a record name the same way.*/
static HV *sub_hv(pTHX_ HV *CSP_RESTRICT parent, const char *CSP_RESTRICT key, STRLEN klen)
{
	SV **slot = hv_fetch(parent, key, (I32)klen, 0);
	HV *h;
	if (slot && *slot && SvROK(*slot) && SvTYPE(SvRV(*slot)) == SVt_PVHV)
		return (HV *)SvRV(*slot);
	h = newHV();
	(void)hv_store(parent, key, (I32)klen, newRV_noinc((SV *)h), 0);
	return h;
}

/*chain_tally() -- the element tally of one chain of one model.

Kept per model as well as per chain because with model => 'all' the same chain
id is a different chain in each model and Perl builds them separately; with the
default model => 1 there is only ever the one key.  The model number is spelled
here the way Perl spells it when it looks the chain up, which is what a decimal
IV stringifies to on both sides.

Reached once per residue rather than once per atom -- chain and model are both
part of a residue's identity, so neither can change without the readers below
noticing a new residue -- which leaves the per-atom cost at the one hash
increment the whole-structure tally already costs.*/
static HV *chain_tally(pTHX_ HV *CSP_RESTRICT root, IV model,
                       const char *CSP_RESTRICT chain, STRLEN chain_len)
{
	char mbuf[32]; //an IV prints in at most 21 characters, sign included
	STRLEN mlen = (STRLEN)my_snprintf(mbuf, sizeof(mbuf), "%" IVdf, model);
	return sub_hv(aTHX_ sub_hv(aTHX_ root, mbuf, mlen), chain, chain_len);
}

//options
static SV *opt_get(pTHX_ HV *CSP_RESTRICT o, const char *CSP_RESTRICT k)
{
	SV **p;
	if (!o) return NULL;
	p = hv_fetch(o, k, (I32)strlen(k), 0);
	return (p && *p && SvOK(*p)) ? *p : NULL;
}

static IV opt_iv(pTHX_ HV *CSP_RESTRICT o, const char *CSP_RESTRICT k, IV dflt)
{
	SV *v = opt_get(aTHX_ o, k);
	return v ? SvIV(v) : dflt;
}

static bool opt_bool(pTHX_ HV *CSP_RESTRICT o, const char *CSP_RESTRICT k, bool dflt)
{
	SV *v = opt_get(aTHX_ o, k);
	return v ? (SvTRUE(v) ? TRUE : FALSE) : dflt;
}

//the parse
#define NCOL 17 //parallel per-atom columns
static const char *const col_name[NCOL] = {
	"serial", "name", "altloc", "resname", "chain", "resseq", "icode",
	"x", "y", "z", "occupancy", "bfactor", "element", "charge",
	"het", "model", "line"
};
enum {
	C_SERIAL, C_NAME, C_ALTLOC, C_RESNAME, C_CHAIN, C_RESSEQ, C_ICODE,
	C_X, C_Y, C_Z, C_OCC, C_B, C_ELEMENT, C_CHARGE, C_HET, C_MODEL, C_LINENO
};

/*The per-residue sums, kept as sums rather than means so that Perl can add
two of them together: a residue whose atoms are written in two runs -- an
alternate conformer group split apart, a ligand interleaved with the polymer
-- arrives here as two entries and is one residue by the time Perl is done
with it.  Means computed here could not be recombined without the counts, and
with the counts they may as well be sums.*/
#define NRSUM 6
static const char *const res_sum_name[NRSUM] = { "sx", "sy", "sz", "n_xyz", "sb", "n_b" };
enum { R_SX, R_SY, R_SZ, R_NXYZ, R_SB, R_NB };

static void flush_residue(pTHX_ AV **CSP_RESTRICT rs, NV sx, NV sy, NV sz, UV nc, NV sb, UV nb)
{
	av_push(rs[R_SX],   nc ? newSVnv(sx) : newSVsv(&PL_sv_undef));
	av_push(rs[R_SY],   nc ? newSVnv(sy) : newSVsv(&PL_sv_undef));
	av_push(rs[R_SZ],   nc ? newSVnv(sz) : newSVsv(&PL_sv_undef));
	av_push(rs[R_NXYZ], newSVuv(nc));
	av_push(rs[R_SB],   nb ? newSVnv(sb) : newSVsv(&PL_sv_undef));
	av_push(rs[R_NB],   newSVuv(nb));
}

/*buf is the whole file, and restrict here is the one that earns the most: every
pointer the loop reads a field through -- line, s, nm_raw -- is derived from it,
so one qualifier on the parameter covers all of them, and the compiler no longer
has to assume that av_push() and the newSV*() calls between two reads of a
record might have rewritten the record.  It holds because the buffer belongs to
the caller for the duration: slurp()'s Newx() block in _parse_file, and an SV we
never touch again in _parse_string.  Nothing in here writes through it.*/
static HV *parse_buf(pTHX_ const char *CSP_RESTRICT buf, STRLEN len, HV *CSP_RESTRICT opts)
{
	HV *out  = newHV();
	HV *meta = newHV();
	HV *elements = newHV();
	HV *chain_elements = newHV();
	/*the tally of the chain the current residue belongs to, refreshed at each
	residue boundary; see chain_tally()*/
	HV *cur_chain = NULL;
	AV *atom_hv = newAV();
	AV *col[NCOL];
	AV *res_sum[NRSUM];
	AV *res_first = newAV(), *res_last = newAV();
	AV *ter = newAV(), *model_nums = newAV();
	HV *want_chain = NULL;
	//whole-structure statistics over the atoms that were kept
	UV n_hydrogen = 0, n_water_atom = 0, bn = 0;
	NV bmin = 0, bmax = 0, bsum = 0;
	NV xmin = 0, ymin = 0, zmin = 0, xmax = 0, ymax = 0, zmax = 0;
	bool have_bbox = FALSE, have_res = FALSE;
	//the residue being accumulated
	NV rsx = 0, rsy = 0, rsz = 0, rsb = 0;
	UV rnc = 0, rnb = 0;
	/*counts of what was read, all of them only ever counted up.  The model and
	residue numbers next to them stay signed: those are read out of the file's
	own columns, where a negative number is legal (residue numbering routinely
	starts before 1) or at worst is damage, and either way is better carried as
	the negative it was written as than wrapped to a huge unsigned one.*/
	UV n_models = 0, n_anisou = 0, n_skipped = 0, n_atom = 0;
	UV n_atom_rec = 0, n_het_rec = 0, lineno = 0;
	IV want_model, cur_model = 1;
	bool keep_h, keep_water, keep_het, keep_meta, keep_anisou, keep_lineno, build_atoms;
	//previous kept atom's residue identity, for boundary detection
	char p_chain[8], p_icode[4], p_resname[8];
	IV p_resseq = 0, p_model = 0;
	bool p_het = FALSE;
	STRLEN pos = 0;

	/*a negative model number means every model.  Zero cannot be the sentinel:
	an ensemble whose models are numbered from 0 is unusual but legal, and
	"model 0" would then quietly mean "all of them".*/
	want_model  = opt_iv(aTHX_ opts, "model", 1);
	keep_h      = opt_bool(aTHX_ opts, "hydrogens", TRUE);
	keep_water  = opt_bool(aTHX_ opts, "waters", TRUE);
	keep_het    = opt_bool(aTHX_ opts, "hetatm", TRUE);
	keep_meta   = opt_bool(aTHX_ opts, "meta", TRUE);
	keep_anisou = opt_bool(aTHX_ opts, "anisou", FALSE);
	keep_lineno = opt_bool(aTHX_ opts, "lineno", FALSE);
	build_atoms = opt_bool(aTHX_ opts, "atom_hashes", FALSE);
	{
		SV *c = opt_get(aTHX_ opts, "chains");
		if (c && SvROK(c) && SvTYPE(SvRV(c)) == SVt_PVHV) want_chain = (HV *)SvRV(c);
	}
	p_chain[0] = p_icode[0] = p_resname[0] = '\0';

	for (unsigned short int i = 0; i < NCOL; i++)  col[i] = newAV();
	for (unsigned short int i = 0; i < NRSUM; i++) res_sum[i] = newAV();

	while (pos < len) {
		const char *CSP_RESTRICT line = buf + pos;
		const char *nl = (const char *)memchr(line, '\n', len - pos);
		STRLEN llen = nl ? (STRLEN)(nl - line) : (len - pos);
		pos += llen + (nl ? 1 : 0);
		lineno++;
		while (llen && (line[llen - 1] == '\r' || line[llen - 1] == ' ' || line[llen - 1] == '\t'))
			llen--; //DOS line ends, and the trailing blanks of a padded record
		if (llen == 0) continue;

		if (llen >= 6 && (memcmp(line, "ATOM  ", 6) == 0 || memcmp(line, "HETATM", 6) == 0)) {
			const char *s;
			const char *nm_raw;
			STRLEN nm_rawlen, n;
			char elbuf[4];
			STRLEN ellen;
			bool het = (line[0] == 'H');
			res_info ri;
			bool known;
			IV iv;
			char chain[4], icode[2], resname[8];
			STRLEN resname_len;

			if (het) n_het_rec++; else n_atom_rec++;
			if (want_model >= 0 && cur_model != want_model) { n_skipped++; continue; }
			if (het && !keep_het) { n_skipped++; continue; }

			//residue name, columns 18-20
			n = fld(line, llen, 17, 19, &s);
			resname_len = n > sizeof(resname) - 1 ? sizeof(resname) - 1 : n;
			memcpy(resname, s, resname_len);
			resname[resname_len] = '\0';
			known = res_lookup(resname, resname_len, &ri);
			if (!keep_water && known && ri.type == RT_WATER) { n_skipped++; continue; }

			/*chain, column 22 -- with column 21 as a fallback, because a
			two-character chain id (some large assemblies) spills left*/
			chain[0] = '\0';
			if (llen > 21 && line[21] != ' ') { chain[0] = line[21]; chain[1] = '\0'; }
			else if (llen > 20 && line[20] != ' ') { chain[0] = line[20]; chain[1] = '\0'; }
			if (want_chain) {
				if (!hv_exists(want_chain, chain, (I32)strlen(chain))) { n_skipped++; continue; }
			}

			//atom name, columns 13-16, kept untrimmed for the element rule
			nm_rawlen = 0;
			nm_raw = line;
			if (llen > 12) {
				STRLEN to = 15 >= llen ? llen - 1 : 15;
				nm_raw = line + 12;
				nm_rawlen = to - 12 + 1;
			}

			//element, columns 77-78, guessed from the name when absent
			ellen = fld(line, llen, 76, 77, &s);
			if (ellen) {
				STRLEN k;
				if (ellen > 2) ellen = 2;
				for (k = 0; k < ellen; k++) elbuf[k] = (char)toupper((unsigned char)s[k]);
				/*and not believed unless it spells an element.  A file written
				before the element column existed keeps the entry id in columns
				73-80 instead, so what sits where the element goes is '1G' and
				every atom in pdb1gdr reads as element '1' -- which also stops
				hydrogens => 0 from finding the hydrogens.  A field that is not
				letters is not an element, and the atom name knows better.*/
				for (k = 0; k < ellen; k++) {
					if (!isALPHA(elbuf[k])) { ellen = 0; break; }
				}
			}
			if (!ellen) ellen = guess_element(nm_raw, nm_rawlen, elbuf);
			/*settled: the one spelling everything downstream sees, so that the
			atom's own element, the whole-structure tally and the per-chain
			tally cannot disagree about what to call a zinc*/
			elem_case(elbuf, ellen);
			if (!keep_h && ellen == 1 && (elbuf[0] == 'H' || elbuf[0] == 'D')) { n_skipped++; continue; }

			/*kept.  The residue boundary is settled first, before this atom has
			been added to anything, so that closing the previous residue is a
			matter of flushing what is already there.  Everything the test needs
			is read from the record's own columns.*/
			{
				IV rs = 0;
				bool changed;
				fld_iv(line, llen, 22, 25, &rs);
				icode[0] = (llen > 26 && line[26] != ' ') ? line[26] : '\0';
				icode[1] = '\0';
				changed = !have_res
				        || rs != p_resseq
				        || p_model != cur_model
				        || p_het != het
				        || strcmp(p_chain, chain) != 0
				        || strcmp(p_icode, icode) != 0
				        || strcmp(p_resname, resname) != 0;
				if (changed) {
					if (have_res) {
						av_push(res_last, newSVuv(n_atom - 1));
						flush_residue(aTHX_ res_sum, rsx, rsy, rsz, rnc, rsb, rnb);
					}
					rsx = rsy = rsz = rsb = 0; rnc = rnb = 0;
					av_push(res_first, newSVuv(n_atom));
					have_res = TRUE;
					p_resseq = rs;
					p_model  = cur_model;
					p_het    = het;
					my_strlcpy(p_chain, chain, sizeof(p_chain));
					my_strlcpy(p_icode, icode, sizeof(p_icode));
					my_strlcpy(p_resname, resname, sizeof(p_resname));
					cur_chain = chain_tally(aTHX_ chain_elements, cur_model,
					                        chain, strlen(chain));
				}
			}

			//the atom's own fields
			{
				const char *nm_s, *alt_s, *chg_s;
				STRLEN nm_n, alt_n, chg_n;
				IV serial = 0;
				bool have_serial, have_occ;
				NV xv, yv, zv, bv, ov;
				bool have_xyz = fld_nv(line, llen, 30, 37, &xv)
				             & fld_nv(line, llen, 38, 45, &yv)
				             & fld_nv(line, llen, 46, 53, &zv);
				bool have_b = fld_nv(line, llen, 60, 65, &bv);

				have_occ    = fld_nv(line, llen, 54, 59, &ov);
				have_serial = fld_iv(line, llen, 6, 10, &serial);
				nm_n  = fld(line, llen, 12, 15, &nm_s);
				alt_n = fld(line, llen, 16, 16, &alt_s);
				chg_n = fld(line, llen, 78, 79, &chg_s);
				if (!charge_ok(chg_s, chg_n)) chg_n = 0;

				/*The sums and extremes are gathered here rather than in Perl
				because they have to touch every atom whether or not the caller
				wanted per-atom hashes.  Doing them in the loop that is already
				reading the numbers costs nothing; doing them again in Perl
				costs more than the whole parse.*/
				if (have_xyz) {
					rsx += xv; rsy += yv; rsz += zv; rnc++;
					if (!have_bbox) {
						xmin = xmax = xv; ymin = ymax = yv; zmin = zmax = zv;
						have_bbox = TRUE;
					} else {
						if (xv < xmin) xmin = xv; else if (xv > xmax) xmax = xv;
						if (yv < ymin) ymin = yv; else if (yv > ymax) ymax = yv;
						if (zv < zmin) zmin = zv; else if (zv > zmax) zmax = zv;
					}
				}
				if (have_b) {
					rsb += bv; rnb++;
					if (!bn++) { bmin = bmax = bv; }
					else if (bv < bmin) bmin = bv;
					else if (bv > bmax) bmax = bv;
					bsum += bv;
				}
				if (ellen) {
					elem_bump(aTHX_ elements, elbuf, ellen);
					if (cur_chain) elem_bump(aTHX_ cur_chain, elbuf, ellen);
					if (ellen == 1 && (elbuf[0] == 'H' || elbuf[0] == 'D')) n_hydrogen++;
				}
				if (known && ri.type == RT_WATER) n_water_atom++;

/*Two shapes to hand back, and never both, because building every
atom twice -- once as columns here, once as a hash in Perl
afterwards -- was costing more than everything else in the read
put together.  The hash is what Chem::Structure::Parser wants, built here
where the fields already are; the columns are what the low-level
parse hands to anyone calling it directly.*/
				if (build_atoms) {
					HV *a = newHV();
					(void)hv_stores(a, "name",   newSVpvn(nm_s, nm_n));
					(void)hv_stores(a, "serial", have_serial ? newSViv(serial) : newSVsv(&PL_sv_undef));
					(void)hv_stores(a, "altloc", newSVpvn(alt_s, alt_n));
					(void)hv_stores(a, "x", have_xyz ? newSVnv(xv) : newSVsv(&PL_sv_undef));
					(void)hv_stores(a, "y", have_xyz ? newSVnv(yv) : newSVsv(&PL_sv_undef));
					(void)hv_stores(a, "z", have_xyz ? newSVnv(zv) : newSVsv(&PL_sv_undef));
					(void)hv_stores(a, "occupancy", have_occ ? newSVnv(ov) : newSVsv(&PL_sv_undef));
					(void)hv_stores(a, "bfactor",   have_b   ? newSVnv(bv) : newSVsv(&PL_sv_undef));
					(void)hv_stores(a, "element", newSVpvn(elbuf, ellen));
					(void)hv_stores(a, "charge",  newSVpvn(chg_s, chg_n));
					(void)hv_stores(a, "hetero",  newSViv(het));
					av_push(atom_hv, newRV_noinc((SV *)a));
				} else {
					av_push(col[C_SERIAL], have_serial ? newSViv(serial) : newSVsv(&PL_sv_undef));
					av_push(col[C_NAME],   newSVpvn(nm_s, nm_n));
					av_push(col[C_ALTLOC], newSVpvn(alt_s, alt_n));
					av_push(col[C_X], have_xyz ? newSVnv(xv) : newSVsv(&PL_sv_undef));
					av_push(col[C_Y], have_xyz ? newSVnv(yv) : newSVsv(&PL_sv_undef));
					av_push(col[C_Z], have_xyz ? newSVnv(zv) : newSVsv(&PL_sv_undef));
					av_push(col[C_OCC], have_occ ? newSVnv(ov) : newSVsv(&PL_sv_undef));
					av_push(col[C_B],   have_b   ? newSVnv(bv) : newSVsv(&PL_sv_undef));
					av_push(col[C_ELEMENT], newSVpvn(elbuf, ellen));
					av_push(col[C_CHARGE],  newSVpvn(chg_s, chg_n));
				}
			}

/*The residue's identity.  Perl reads these at the index a residue
starts on, once per residue rather than once per atom, but they are
emitted per atom because that is where they are read from.*/
			av_push(col[C_RESNAME], newSVpvn(resname, resname_len));
			av_push(col[C_CHAIN], newSVpvn(chain, strlen(chain)));
			if (fld_iv(line, llen, 22, 25, &iv)) av_push(col[C_RESSEQ], newSViv(iv));
			else                                 av_push(col[C_RESSEQ], newSVsv(&PL_sv_undef));
			av_push(col[C_ICODE], newSVpvn(icode, icode[0] ? 1 : 0));
			av_push(col[C_HET], newSViv(het));
			av_push(col[C_MODEL], newSViv(cur_model));
			if (keep_lineno) av_push(col[C_LINENO], newSVuv(lineno));

			n_atom++;
			continue;
		}

		if (llen >= 6 && memcmp(line, "ANISOU", 6) == 0) {
			n_anisou++;
			if (!keep_anisou) continue;
			//fall through to meta
		}

		if (llen >= 5 && memcmp(line, "MODEL", 5) == 0) {
			IV m;
			n_models++;
			cur_model = fld_iv(line, llen, 10, 13, &m) ? m : (IV)n_models;
			av_push(model_nums, newSViv(cur_model));
			continue;
		}
		if (llen >= 6 && memcmp(line, "ENDMDL", 6) == 0) continue;

		if (llen >= 3 && memcmp(line, "TER", 3) == 0) {
			HV *t = newHV();
			const char *s;
			STRLEN n;
			IV iv;
			n = fld(line, llen, 17, 19, &s); (void)hv_stores(t, "resname", newSVpvn(s, n));
			n = fld(line, llen, 21, 21, &s); (void)hv_stores(t, "chain",   newSVpvn(s, n));
			n = fld(line, llen, 26, 26, &s); (void)hv_stores(t, "icode",   newSVpvn(s, n));
			(void)hv_stores(t, "resseq", fld_iv(line, llen, 22, 25, &iv) ? newSViv(iv) : newSVsv(&PL_sv_undef));
			(void)hv_stores(t, "model", newSViv(cur_model));
			av_push(ter, newRV_noinc((SV *)t));
			continue;
		}

		/*everything else is a header/annotation record: keep the line whole and
		let Perl take it apart by column*/
		if (keep_meta) {
			const char *s;
			STRLEN n = fld(line, llen, 0, 5, &s);
			if (n) {
				SV **slot = hv_fetch(meta, s, (I32)n, 0);
				AV *av;
				if (slot && *slot && SvROK(*slot)) {
					av = (AV *)SvRV(*slot);
				} else {
					av = newAV();
					(void)hv_store(meta, s, (I32)n, newRV_noinc((SV *)av), 0);
				}
				av_push(av, newSVpvn(line, llen));
			}
		}
	}
	if (have_res) {
		av_push(res_last, newSVuv(n_atom - 1));
		flush_residue(aTHX_ res_sum, rsx, rsy, rsz, rnc, rsb, rnb);
	}

	(void)hv_stores(out, "atoms", newRV_noinc((SV *)atom_hv));
	for (unsigned short int i = 0; i < NRSUM; i++)
		(void)hv_store(out, res_sum_name[i], (I32)strlen(res_sum_name[i]),
		               newRV_noinc((SV *)res_sum[i]), 0);
	(void)hv_stores(out, "elements",     newRV_noinc((SV *)elements));
	(void)hv_stores(out, "chain_elements", newRV_noinc((SV *)chain_elements));
	(void)hv_stores(out, "n_hydrogens",  newSVuv(n_hydrogen));
	(void)hv_stores(out, "n_water_atoms", newSVuv(n_water_atom));
	if (bn) {
		HV *b = newHV();
		(void)hv_stores(b, "min",  newSVnv(bmin));
		(void)hv_stores(b, "max",  newSVnv(bmax));
		(void)hv_stores(b, "mean", newSVnv(bsum / (NV)bn));
		(void)hv_stores(b, "n",    newSVuv(bn));
		(void)hv_stores(out, "bfactor_stats", newRV_noinc((SV *)b));
	} else {
		(void)hv_stores(out, "bfactor_stats", newSVsv(&PL_sv_undef));
	}
	if (have_bbox) {
		HV *bb = newHV();
		AV *ctr = newAV();
		(void)hv_stores(bb, "xmin", newSVnv(xmin));
		(void)hv_stores(bb, "ymin", newSVnv(ymin));
		(void)hv_stores(bb, "zmin", newSVnv(zmin));
		(void)hv_stores(bb, "xmax", newSVnv(xmax));
		(void)hv_stores(bb, "ymax", newSVnv(ymax));
		(void)hv_stores(bb, "zmax", newSVnv(zmax));
		av_push(ctr, newSVnv((xmin + xmax) / 2));
		av_push(ctr, newSVnv((ymin + ymax) / 2));
		av_push(ctr, newSVnv((zmin + zmax) / 2));
		(void)hv_stores(out, "bbox",   newRV_noinc((SV *)bb));
		(void)hv_stores(out, "center", newRV_noinc((SV *)ctr));
	} else {
		(void)hv_stores(out, "bbox",   newSVsv(&PL_sv_undef));
		(void)hv_stores(out, "center", newSVsv(&PL_sv_undef));
	}

	for (unsigned short int i = 0; i < NCOL; i++) {
		if (i == C_LINENO && !keep_lineno) { SvREFCNT_dec((SV *)col[i]); continue; }
		(void)hv_store(out, col_name[i], (I32)strlen(col_name[i]), newRV_noinc((SV *)col[i]), 0);
	}
	(void)hv_stores(out, "res_first",    newRV_noinc((SV *)res_first));
	(void)hv_stores(out, "res_last",     newRV_noinc((SV *)res_last));
	(void)hv_stores(out, "ter",          newRV_noinc((SV *)ter));
	(void)hv_stores(out, "meta",         newRV_noinc((SV *)meta));
	(void)hv_stores(out, "model_numbers", newRV_noinc((SV *)model_nums));
	(void)hv_stores(out, "n_atoms",      newSVuv(n_atom));
	(void)hv_stores(out, "n_residues",   newSVuv((UV)(av_len(res_first) + 1)));
	(void)hv_stores(out, "n_models",     newSVuv(n_models ? n_models : 1));
	(void)hv_stores(out, "n_anisou",     newSVuv(n_anisou));
	(void)hv_stores(out, "n_skipped",    newSVuv(n_skipped));
	(void)hv_stores(out, "n_atom_records",   newSVuv(n_atom_rec));
	(void)hv_stores(out, "n_hetatm_records", newSVuv(n_het_rec));
	(void)hv_stores(out, "n_lines",      newSVuv(lineno));
	(void)hv_stores(out, "format",       newSVpvs("pdb"));
	return out;
}

/*mmCIF -- the same structures, written down differently.

The same structures, written down differently.  A PDB file is fixed columns; an
mmCIF file is tag/value pairs and `loop_' tables, free-form, with quoting.  So
none of the column arithmetic above applies and the whole reader is separate --
but it fills in exactly the same output hash, atom for atom and field for field,
because everything above this point in the module and everything in Perl below
it is written against that hash and not against either file format.  A caller
who reads 1cka.pdb and 1cka.cif gets two structures that compare equal.

Where the two formats disagree about how to write the same fact, this reader
converts to what the PDB reader would have produced: a formal charge of -1 comes
back as "1-", auth_* identifiers are preferred over label_* ones because auth_*
is what the PDB record carried, and a `.' or `?' -- mmCIF for "not applicable"
and "unknown" -- comes back as the empty field the PDB record would have had.*/

//the token kinds the lexer below returns
#define CT_EOF   0 //end of the buffer
#define CT_VALUE 1 //a value, quoted or bare, including '.' and '?'
#define CT_TAG   2 //_category.item
#define CT_LOOP  3 //loop_
#define CT_DATA  4 //data_<block>, with t->s the block name
#define CT_SAVE  5 //save_<frame>, which this reader steps over
#define CT_STOP  6 //stop_ or global_, likewise

typedef struct {
	const char *buf;
	const char *end;
	const char *p;
} cif_lex;

typedef struct {
	const char *s;
	STRLEN n;
	unsigned short int kind; //CT_*: seven token kinds and no sentinel
	/*a quoted value is text and nothing else: '.' in quotes is a full stop, and
	only a bare '.' is the null the format means*/
	bool quoted;
} cif_tok;

static bool cif_iskw(const char *CSP_RESTRICT s, STRLEN n, const char *CSP_RESTRICT kw,
                    STRLEN kwn)
{
	STRLEN i;
	if (n < kwn) return FALSE;
	for (i = 0; i < kwn; i++)
		if (tolower((unsigned char)s[i]) != kw[i]) return FALSE;
	return TRUE;
}

static bool cif_space(char c)
{
	return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}

//a bare '.' (not applicable) or '?' (unknown): the field was not filled in
static bool cif_null(const cif_tok *CSP_RESTRICT t)
{
	return !t->quoted && t->n == 1 && (t->s[0] == '.' || t->s[0] == '?');
}

/*cif_next() -- one token.  The four shapes a value can arrive in are a bare
word, a quoted string, and a semicolon-delimited text field, which is the only
one that can span lines and the only one whose delimiter is position-sensitive:
a ';' is a delimiter only as the first character of a line, and anywhere else it
is an ordinary character in a value.*/
static void cif_next(cif_lex *CSP_RESTRICT lx, cif_tok *CSP_RESTRICT t)
{
	const char *p = lx->p, *end = lx->end;
	t->quoted = FALSE;
	for (;;) {
		while (p < end && cif_space(*p)) p++;
		if (p >= end) { lx->p = end; t->kind = CT_EOF; t->s = end; t->n = 0; return; }
		if (*p == '#') {                       //a comment runs to the end of its line
			while (p < end && *p != '\n') p++;
			continue;
		}
		break;
	}

	if (*p == ';' && (p == lx->buf || p[-1] == '\n')) {
		const char *s = p + 1, *q = s;
		for (;;) {
			const char *CSP_RESTRICT nl = (const char *)memchr(q, '\n', (STRLEN)(end - q));
			if (!nl) { q = end; lx->p = end; break; }
			if (nl + 1 < end && nl[1] == ';') { q = nl; lx->p = nl + 2; break; }
			if (nl + 1 >= end) { q = nl; lx->p = end; break; }
			q = nl + 1;
		}
		t->s = s;
		t->n = (STRLEN)(q - s);
		t->kind = CT_VALUE;
		t->quoted = TRUE;
		return;
	}

	if (*p == '\'' || *p == '"') {
		/*the closing quote is one followed by whitespace or the end of the file.
		Anything else is a quote inside the value, which is how a CIF writes
		O5' without escaping it.*/
		char qc = *p;
		const char *CSP_RESTRICT s = p + 1, *q = s;
		for (;;) {
			while (q < end && *q != qc) q++;
			if (q >= end) break;
			if (q + 1 >= end || cif_space(q[1])) break;
			q++;
		}
		t->s = s;
		t->n = (STRLEN)(q - s);
		t->kind = CT_VALUE;
		t->quoted = TRUE;
		lx->p = q < end ? q + 1 : end;
		return;
	}

	{
		const char *s = p;
		STRLEN n;
		while (p < end && !cif_space(*p)) p++;
		n = (STRLEN)(p - s);
		lx->p = p;
		t->s = s;
		t->n = n;
		if (s[0] == '_')                             t->kind = CT_TAG;
		else if (n == 5 && cif_iskw(s, n, "loop_", 5))  t->kind = CT_LOOP;
		else if (n == 5 && cif_iskw(s, n, "stop_", 5))  t->kind = CT_STOP;
		else if (n == 7 && cif_iskw(s, n, "global_", 7)) t->kind = CT_STOP;
		else if (n >= 5 && cif_iskw(s, n, "data_", 5)) { t->kind = CT_DATA; t->s = s + 5; t->n = n - 5; }
		else if (n >= 5 && cif_iskw(s, n, "save_", 5)) { t->kind = CT_SAVE; t->s = s + 5; t->n = n - 5; }
		else t->kind = CT_VALUE;
	}
}

/*a tag is _category.item.  Both halves are wanted -- the category says which
table a loop is, the item says which column -- and both are matched without
regard to case, which the format allows and real files use.*/
static STRLEN cif_split_tag(const char *CSP_RESTRICT s, STRLEN n,
                            const char **CSP_RESTRICT item, STRLEN *CSP_RESTRICT itemn)
{
	const char *dot = (const char *)memchr(s, '.', n);
	if (!dot) { *item = s; *itemn = n; return n; }   //core CIF: no category half
	*item  = dot + 1;
	*itemn = n - (STRLEN)(dot - s) - 1;
	return (STRLEN)(dot - s);
}

//the _atom_site columns this reader knows what to do with
static const char *const cif_atom_item[] = {
	"group_pdb", "id", "type_symbol", "label_atom_id", "label_alt_id",
	"label_comp_id", "label_asym_id", "label_seq_id", "pdbx_pdb_ins_code",
	"cartn_x", "cartn_y", "cartn_z", "occupancy", "b_iso_or_equiv",
	"pdbx_formal_charge", "auth_seq_id", "auth_comp_id", "auth_asym_id",
	"auth_atom_id", "pdbx_pdb_model_num"
};
enum {
	A_GROUP, A_ID, A_SYMBOL, A_LATOM, A_ALT, A_LCOMP, A_LASYM, A_LSEQ, A_ICODE,
	A_X, A_Y, A_Z, A_OCC, A_B, A_CHARGE, A_ASEQ, A_ACOMP, A_AASYM, A_AATOM,
	A_MODEL, A_NFIELD
};

static short int cif_atom_field(const char *CSP_RESTRICT item, STRLEN n)
{
	short int i;
	for (i = 0; i < A_NFIELD; i++) {
		STRLEN kn = strlen(cif_atom_item[i]);
		if (kn == n && cif_iskw(item, n, cif_atom_item[i], kn)) return i;
	}
	return -1;
}

/*mmCIF writes a formal charge as a signed integer and PDB as a magnitude
followed by its sign, so -1 there is "1-" here.  A value that is not an integer
is passed through: some writers already put the PDB spelling in this field.

A charge of zero comes back as "0" rather than as nothing, because the two are
different answers and the PDB reader keeps them apart: "the field was blank"
against "the field said zero".  mmCIF spells the first '?', which arrives here
as a field that was never set, so nothing is lost by reading the second
literally -- and 6cc9, which writes a 0 in columns 79-80, reads the same from
either format because of it.*/
static STRLEN cif_charge(const char *CSP_RESTRICT s, STRLEN n, char *CSP_RESTRICT buf)
{
	IV c;
	if (n == 0) return 0;
	if (!str2iv(s, n, &c)) {
		if (n > 2) n = 2;
		memcpy(buf, s, n);
		return n;
	}
	if (c == 0) { buf[0] = '0'; return 1; }
	if (c < 0) { c = -c; buf[1] = '-'; } else buf[1] = '+';
	if (c > 9) return 0;
	buf[0] = (char)('0' + c);
	return 2;
}

/*guess_element() reads an atom name out of the four columns a PDB record gives
it, and its rule is where in those four the name starts.  mmCIF has no columns,
so the name goes back where the PDB convention would have put it before the
question is asked.  Only reached when _atom_site.type_symbol is missing, which
is rare -- but a file without it should still not read CA as calcium.*/
static STRLEN cif_guess_element(const char *CSP_RESTRICT s, STRLEN n, char *CSP_RESTRICT buf)
{
	char pad[4];
	STRLEN i;
	if (n >= 4) return guess_element(s, 4, buf);
	pad[0] = ' ';
	for (i = 0; i < n; i++)     pad[1 + i] = s[i];
	for (i = n + 1; i < 4; i++) pad[i] = ' ';
	return guess_element(pad, 4, buf);
}

/*Everything one _atom_site row needs to add itself to, gathered up so that a
row can be handed to one function.  A loop_ is the ordinary way to write the
table and a run of plain tags is the legal way to write a table with one row in
it; both end up in cif_atom_row() rather than in two copies of it.*/
typedef struct {
	AV **col;
	AV **res_sum;
	AV *res_first, *res_last, *atom_hv, *model_nums;
	HV *elements, *chain_elements, *want_chain;
	//the tally of the chain the current residue belongs to; see chain_tally()
	HV *cur_chain;
	UV n_hydrogen, n_water_atom, bn;
	NV bmin, bmax, bsum;
	NV xmin, ymin, zmin, xmax, ymax, zmax;
	bool have_bbox;
	NV rsx, rsy, rsz, rsb;
	UV rnc, rnb;
	//counted up only; the numbers read out of the file stay signed, as in parse_buf
	UV n_models, n_anisou, n_skipped, n_atom, n_atom_rec, n_het_rec;
	IV want_model, cur_model;
	bool have_res, seen_model, keep_h, keep_water, keep_het, keep_anisou, build_atoms;
	char p_chain[8], p_icode[4], p_resname[8];
	IV p_resseq, p_model;
	bool p_het;
} cif_state;

/*one row of _atom_site.  v[] holds a pointer and a length per known column, or
NULL where the row said `.' or `?' or the file has no such column, so that "the
field was not given" is one test wherever it is asked.

The order of the tests below is the order the PDB reader uses -- model, then
HETATM, then water, then chain, then hydrogen -- because n_skipped counts what
was thrown away and the two readers have to agree on the count as well as on
the structure.*/
static void cif_atom_row(pTHX_ cif_state *CSP_RESTRICT st,
                         const char **CSP_RESTRICT v, STRLEN *CSP_RESTRICT vn)
{
	char resname[8], chain[8], icode[4], elbuf[4], chgbuf[4];
	STRLEN resname_len, chain_len, ellen, chg_n;
	const char *CSP_RESTRICT nm_s, *CSP_RESTRICT alt_s;
	STRLEN nm_n, alt_n;
	res_info ri;
	bool known, het, changed;
	IV rs = 0, serial = 0, model = 1;
	bool have_rs, have_serial, have_occ, have_xyz, have_b;
	NV xv = 0, yv = 0, zv = 0, ov = 0, bv = 0;

	//which model this row belongs to, counted before anything can skip the row
	if (v[A_MODEL] && str2iv(v[A_MODEL], vn[A_MODEL], &model)) {
		if (!st->seen_model || model != st->cur_model) {
			st->cur_model = model;
			st->seen_model = TRUE;
			st->n_models++;
			av_push(st->model_nums, newSViv(model));
		}
	} else if (!st->seen_model) {
		st->cur_model = 1;
		st->seen_model = TRUE;
		st->n_models++;
		av_push(st->model_nums, newSViv(1));
	}

	het = v[A_GROUP] && vn[A_GROUP] >= 6 && cif_iskw(v[A_GROUP], vn[A_GROUP], "hetatm", 6);
	if (het) st->n_het_rec++; else st->n_atom_rec++;

	if (st->want_model >= 0 && st->cur_model != st->want_model) { st->n_skipped++; return; }
	if (het && !st->keep_het) { st->n_skipped++; return; }

	//residue name: auth_comp_id is the name the PDB record carried
	{
		const char *s = v[A_ACOMP] ? v[A_ACOMP] : v[A_LCOMP];
		STRLEN n = v[A_ACOMP] ? vn[A_ACOMP] : (v[A_LCOMP] ? vn[A_LCOMP] : 0);
		if (!s) { s = ""; n = 0; }
		resname_len = n > sizeof(resname) - 1 ? sizeof(resname) - 1 : n;
		memcpy(resname, s, resname_len);
		resname[resname_len] = '\0';
	}
	known = res_lookup(resname, resname_len, &ri);
	if (!st->keep_water && known && ri.type == RT_WATER) { st->n_skipped++; return; }

	{
		const char *s = v[A_AASYM] ? v[A_AASYM] : v[A_LASYM];
		STRLEN n = v[A_AASYM] ? vn[A_AASYM] : (v[A_LASYM] ? vn[A_LASYM] : 0);
		if (!s) { s = ""; n = 0; }
		chain_len = n > sizeof(chain) - 1 ? sizeof(chain) - 1 : n;
		memcpy(chain, s, chain_len);
		chain[chain_len] = '\0';
	}
	if (st->want_chain && !hv_exists(st->want_chain, chain, (I32)chain_len)) {
		st->n_skipped++;
		return;
	}

	nm_s = v[A_AATOM] ? v[A_AATOM] : v[A_LATOM];
	nm_n = v[A_AATOM] ? vn[A_AATOM] : (v[A_LATOM] ? vn[A_LATOM] : 0);
	if (!nm_s) { nm_s = ""; nm_n = 0; }
	alt_s = v[A_ALT] ? v[A_ALT] : "";
	alt_n = v[A_ALT] ? vn[A_ALT] : 0;

	if (v[A_SYMBOL]) {
		STRLEN k;
		ellen = vn[A_SYMBOL] > 2 ? 2 : vn[A_SYMBOL];
		for (k = 0; k < ellen; k++) elbuf[k] = (char)toupper((unsigned char)v[A_SYMBOL][k]);
	} else {
		ellen = cif_guess_element(nm_s, nm_n, elbuf);
	}
	//settled, and spelled as the PDB reader spells it; see elem_case()
	elem_case(elbuf, ellen);
	if (!st->keep_h && ellen == 1 && (elbuf[0] == 'H' || elbuf[0] == 'D')) {
		st->n_skipped++;
		return;
	}

	//kept.  auth_seq_id is the residue number the PDB record carried; label_seq_id
	//is the position in the entity, and is null for everything that is not polymer
	have_rs = FALSE;
	if (v[A_ASEQ])      have_rs = str2iv(v[A_ASEQ], vn[A_ASEQ], &rs);
	if (!have_rs && v[A_LSEQ]) have_rs = str2iv(v[A_LSEQ], vn[A_LSEQ], &rs);
	if (!have_rs) rs = 0;
	icode[0] = (v[A_ICODE] && vn[A_ICODE]) ? v[A_ICODE][0] : '\0';
	icode[1] = '\0';

	changed = !st->have_res
	        || rs != st->p_resseq
	        || st->p_model != st->cur_model
	        || st->p_het != het
	        || strcmp(st->p_chain, chain) != 0
	        || strcmp(st->p_icode, icode) != 0
	        || strcmp(st->p_resname, resname) != 0;
	if (changed) {
		if (st->have_res) {
			av_push(st->res_last, newSVuv(st->n_atom - 1));
			flush_residue(aTHX_ st->res_sum, st->rsx, st->rsy, st->rsz, st->rnc, st->rsb, st->rnb);
		}
		st->rsx = st->rsy = st->rsz = st->rsb = 0;
		st->rnc = st->rnb = 0;
		av_push(st->res_first, newSVuv(st->n_atom));
		st->have_res = TRUE;
		st->p_resseq = rs;
		st->p_model  = st->cur_model;
		st->p_het    = het;
		my_strlcpy(st->p_chain, chain, sizeof(st->p_chain));
		my_strlcpy(st->p_icode, icode, sizeof(st->p_icode));
		my_strlcpy(st->p_resname, resname, sizeof(st->p_resname));
		st->cur_chain = chain_tally(aTHX_ st->chain_elements, st->cur_model,
		                            chain, chain_len);
	}

	have_xyz = v[A_X] && v[A_Y] && v[A_Z]
	         && str2nv(v[A_X], vn[A_X], &xv)
	         && str2nv(v[A_Y], vn[A_Y], &yv)
	         && str2nv(v[A_Z], vn[A_Z], &zv);
	have_occ    = v[A_OCC] && str2nv(v[A_OCC], vn[A_OCC], &ov);
	have_b      = v[A_B]   && str2nv(v[A_B],   vn[A_B],   &bv);
	have_serial = v[A_ID]  && str2iv(v[A_ID],  vn[A_ID],  &serial);
	chg_n = v[A_CHARGE] ? cif_charge(v[A_CHARGE], vn[A_CHARGE], chgbuf) : 0;

	if (have_xyz) {
		st->rsx += xv; st->rsy += yv; st->rsz += zv; st->rnc++;
		if (!st->have_bbox) {
			st->xmin = st->xmax = xv; st->ymin = st->ymax = yv; st->zmin = st->zmax = zv;
			st->have_bbox = TRUE;
		} else {
			if (xv < st->xmin) st->xmin = xv; else if (xv > st->xmax) st->xmax = xv;
			if (yv < st->ymin) st->ymin = yv; else if (yv > st->ymax) st->ymax = yv;
			if (zv < st->zmin) st->zmin = zv; else if (zv > st->zmax) st->zmax = zv;
		}
	}
	if (have_b) {
		st->rsb += bv; st->rnb++;
		if (!st->bn++) { st->bmin = st->bmax = bv; }
		else if (bv < st->bmin) st->bmin = bv;
		else if (bv > st->bmax) st->bmax = bv;
		st->bsum += bv;
	}
	if (ellen) {
		elem_bump(aTHX_ st->elements, elbuf, ellen);
		if (st->cur_chain) elem_bump(aTHX_ st->cur_chain, elbuf, ellen);
		if (ellen == 1 && (elbuf[0] == 'H' || elbuf[0] == 'D')) st->n_hydrogen++;
	}
	if (known && ri.type == RT_WATER) st->n_water_atom++;

	if (st->build_atoms) {
		HV *a = newHV();
		(void)hv_stores(a, "name",   newSVpvn(nm_s, nm_n));
		(void)hv_stores(a, "serial", have_serial ? newSViv(serial) : newSVsv(&PL_sv_undef));
		(void)hv_stores(a, "altloc", newSVpvn(alt_s, alt_n));
		(void)hv_stores(a, "x", have_xyz ? newSVnv(xv) : newSVsv(&PL_sv_undef));
		(void)hv_stores(a, "y", have_xyz ? newSVnv(yv) : newSVsv(&PL_sv_undef));
		(void)hv_stores(a, "z", have_xyz ? newSVnv(zv) : newSVsv(&PL_sv_undef));
		(void)hv_stores(a, "occupancy", have_occ ? newSVnv(ov) : newSVsv(&PL_sv_undef));
		(void)hv_stores(a, "bfactor",   have_b   ? newSVnv(bv) : newSVsv(&PL_sv_undef));
		(void)hv_stores(a, "element", newSVpvn(elbuf, ellen));
		(void)hv_stores(a, "charge",  newSVpvn(chgbuf, chg_n));
		(void)hv_stores(a, "hetero",  newSViv(het));
		av_push(st->atom_hv, newRV_noinc((SV *)a));
	} else {
		av_push(st->col[C_SERIAL], have_serial ? newSViv(serial) : newSVsv(&PL_sv_undef));
		av_push(st->col[C_NAME],   newSVpvn(nm_s, nm_n));
		av_push(st->col[C_ALTLOC], newSVpvn(alt_s, alt_n));
		av_push(st->col[C_X], have_xyz ? newSVnv(xv) : newSVsv(&PL_sv_undef));
		av_push(st->col[C_Y], have_xyz ? newSVnv(yv) : newSVsv(&PL_sv_undef));
		av_push(st->col[C_Z], have_xyz ? newSVnv(zv) : newSVsv(&PL_sv_undef));
		av_push(st->col[C_OCC], have_occ ? newSVnv(ov) : newSVsv(&PL_sv_undef));
		av_push(st->col[C_B],   have_b   ? newSVnv(bv) : newSVsv(&PL_sv_undef));
		av_push(st->col[C_ELEMENT], newSVpvn(elbuf, ellen));
		av_push(st->col[C_CHARGE],  newSVpvn(chgbuf, chg_n));
	}

	av_push(st->col[C_RESNAME], newSVpvn(resname, resname_len));
	av_push(st->col[C_CHAIN],   newSVpvn(chain, chain_len));
	av_push(st->col[C_RESSEQ],  have_rs ? newSViv(rs) : newSVsv(&PL_sv_undef));
	av_push(st->col[C_ICODE],   newSVpvn(icode, icode[0] ? 1 : 0));
	av_push(st->col[C_HET],     newSViv(het));
	av_push(st->col[C_MODEL],   newSViv(st->cur_model));
	st->n_atom++;
}

//a token as the SV a caller reads: a null field is undef, not an empty string
static SV *cif_sv(pTHX_ const cif_tok *CSP_RESTRICT t)
{
	return cif_null(t) ? newSVsv(&PL_sv_undef) : newSVpvn(t->s, t->n);
}

//lowercased, so that a caller looking a tag up never has to guess at its case
static SV *cif_key(pTHX_ const char *CSP_RESTRICT s, STRLEN n)
{
	SV *sv = newSVpvn(s, n);
	char *CSP_RESTRICT p = SvPVX(sv);
	STRLEN i;
	for (i = 0; i < n; i++) p[i] = (char)tolower((unsigned char)p[i]);
	return sv;
}

static void cif_store_lc(pTHX_ HV *CSP_RESTRICT h, const char *CSP_RESTRICT k, STRLEN kn,
                         SV *CSP_RESTRICT val)
{
	SV *key = cif_key(aTHX_ k, kn);
	(void)hv_store_ent(h, key, val, 0);
	SvREFCNT_dec(key);
}

static AV *cif_loop_av(pTHX_ HV *CSP_RESTRICT loops, const char *CSP_RESTRICT cat, STRLEN catn)
{
	SV *key = cif_key(aTHX_ cat, catn);
	HE *he = hv_fetch_ent(loops, key, 0, 0);
	AV *av;
	if (he && HeVAL(he) && SvROK(HeVAL(he))) {
		av = (AV *)SvRV(HeVAL(he));
	} else {
		av = newAV();
		(void)hv_store_ent(loops, key, newRV_noinc((SV *)av), 0);
	}
	SvREFCNT_dec(key);
	return av;
}

/*parse_cif_buf() -- the whole file, once through, into the same hash
parse_buf() fills.  The atom table goes to cif_atom_row(); everything else is
kept as it was written, single tags in `cif' and loops in `cif_loops', for Perl
to make sense of.  Nothing here decides what a category means: that is the same
division of labour the PDB reader keeps, where C reads the coordinates and Perl
reads the header.*/
static HV *parse_cif_buf(pTHX_ const char *CSP_RESTRICT buf, STRLEN len, HV *CSP_RESTRICT opts)
{
	HV *out = newHV();
	HV *meta = newHV(), *cif = newHV(), *loops = newHV();
	AV *col[NCOL], *res_sum[NRSUM];
	cif_state st;
	cif_lex lx;
	cif_tok tk;
	size_t i;
	bool keep_meta;
	UV lineno = 0;
	SV *block = NULL;
	const char *v[A_NFIELD];
	STRLEN vn[A_NFIELD];
	//an _atom_site table written as plain tags rather than as a loop_
	const char *sv_[A_NFIELD];
	STRLEN svn_[A_NFIELD];
	bool have_single = FALSE;

	Zero(&st, 1, cif_state);
	for (i = 0; i < NCOL; i++)  col[i] = newAV();
	for (i = 0; i < NRSUM; i++) res_sum[i] = newAV();
	for (i = 0; i < A_NFIELD; i++) { sv_[i] = NULL; svn_[i] = 0; }

	st.col        = col;
	st.res_sum    = res_sum;
	st.res_first  = newAV();
	st.res_last   = newAV();
	st.atom_hv    = newAV();
	st.model_nums = newAV();
	st.elements   = newHV();
	st.chain_elements = newHV();
	st.cur_model  = 1;
	st.want_model = opt_iv(aTHX_ opts, "model", 1);
	st.keep_h     = opt_bool(aTHX_ opts, "hydrogens", TRUE);
	st.keep_water = opt_bool(aTHX_ opts, "waters", TRUE);
	st.keep_het   = opt_bool(aTHX_ opts, "hetatm", TRUE);
	st.keep_anisou= opt_bool(aTHX_ opts, "anisou", FALSE);
	st.build_atoms= opt_bool(aTHX_ opts, "atom_hashes", FALSE);
	keep_meta     = opt_bool(aTHX_ opts, "meta", TRUE);
	{
		SV *c = opt_get(aTHX_ opts, "chains");
		if (c && SvROK(c) && SvTYPE(SvRV(c)) == SVt_PVHV) st.want_chain = (HV *)SvRV(c);
	}

	{	//lines are not what this format is made of, but callers still count them
		const char *CSP_RESTRICT p = buf, *CSP_RESTRICT e = buf + len;
		while (p < e) {
			const char *CSP_RESTRICT nl = (const char *)memchr(p, '\n', (STRLEN)(e - p));
			lineno++;
			if (!nl) break;
			p = nl + 1;
		}
	}

	lx.buf = buf;
	lx.end = buf + len;
	lx.p   = buf;

	for (;;) {
		cif_next(&lx, &tk);
		if (tk.kind == CT_EOF) break;

		if (tk.kind == CT_DATA) {
			if (!block && tk.n) block = newSVpvn(tk.s, tk.n);
			continue;
		}
		if (tk.kind == CT_STOP || tk.kind == CT_SAVE) continue;

		if (tk.kind == CT_LOOP) {
			const char **CSP_RESTRICT tags = NULL;
			STRLEN *CSP_RESTRICT tagn = NULL;
			short int *CSP_RESTRICT fld_of = NULL; //a column this reader ignores is -1
			size_t ntags = 0, cap = 16; //an _atom_site loop has 21 tags: one Renew
			bool is_atom = FALSE, is_aniso = FALSE, keep_loop = FALSE;
			const char *CSP_RESTRICT cat = NULL;
			STRLEN catn = 0;
			AV *rows = NULL;

			Newx(tags, cap, const char *);
			Newx(tagn, cap, STRLEN);
			for (;;) {
				const char *CSP_RESTRICT save = lx.p;
				cif_next(&lx, &tk);
				if (tk.kind != CT_TAG) { lx.p = save; break; }
				if (ntags == cap) {
					cap *= 2;
					Renew(tags, cap, const char *);
					Renew(tagn, cap, STRLEN);
				}
				tags[ntags] = tk.s;
				tagn[ntags] = tk.n;
				ntags++;
			}
			if (ntags == 0) { Safefree(tags); Safefree(tagn); continue; }

			{
				const char *item;
				STRLEN itemn;
				catn = cif_split_tag(tags[0], tagn[0], &item, &itemn);
				cat  = tags[0];
				is_atom  = (catn == 10 && cif_iskw(cat, catn, "_atom_site", 10));
				is_aniso = (catn == 20 && cif_iskw(cat, catn, "_atom_site_anisotrop", 20));
			}
			/*_atom_site_anisotrop is as long as the atom table and is wanted about
			as often as ANISOU is, which is to say hardly ever; it is counted
			always and kept only when the caller asked for it*/
			keep_loop = keep_meta && !is_atom && (!is_aniso || st.keep_anisou);
			if (keep_loop) rows = cif_loop_av(aTHX_ loops, cat, catn);

			Newx(fld_of, ntags, short int);
			for (i = 0; i < ntags; i++) {
				const char *item;
				STRLEN itemn;
				(void)cif_split_tag(tags[i], tagn[i], &item, &itemn);
				fld_of[i] = is_atom ? cif_atom_field(item, itemn) : -1;
			}

			for (;;) {
				const char *CSP_RESTRICT save = lx.p;
				size_t col_i;
				HV *row = NULL;
				cif_next(&lx, &tk);
				if (tk.kind != CT_VALUE) { lx.p = save; break; }
				if (is_atom) for (i = 0; i < A_NFIELD; i++) { v[i] = NULL; vn[i] = 0; }
				if (keep_loop) row = newHV();

				for (col_i = 0; ; col_i++) {
					if (col_i > 0) {
						save = lx.p;
						cif_next(&lx, &tk);
						if (tk.kind != CT_VALUE) { lx.p = save; break; }
					}
					if (is_atom) {
						short int f = fld_of[col_i];
						if (f >= 0 && !cif_null(&tk)) { v[f] = tk.s; vn[f] = tk.n; }
					} else if (row) {
						const char *item;
						STRLEN itemn;
						(void)cif_split_tag(tags[col_i], tagn[col_i], &item, &itemn);
						cif_store_lc(aTHX_ row, item, itemn, cif_sv(aTHX_ &tk));
					}
					if (col_i + 1 >= ntags) break;
				}
				if (is_atom)      cif_atom_row(aTHX_ &st, v, vn);
				else if (is_aniso) st.n_anisou++;
				if (row) av_push(rows, newRV_noinc((SV *)row));
			}
			Safefree(tags);
			Safefree(tagn);
			Safefree(fld_of);
			continue;
		}

		if (tk.kind == CT_TAG) {
			const char *tag = tk.s, *item;
			STRLEN tagn = tk.n, itemn, catn;
			const char *save = lx.p;
			catn = cif_split_tag(tag, tagn, &item, &itemn);
			cif_next(&lx, &tk);
			if (tk.kind != CT_VALUE) { lx.p = save; continue; }   //a tag with no value
			if (catn == 10 && cif_iskw(tag, catn, "_atom_site", 10)) {
				short int f = cif_atom_field(item, itemn);
				if (f >= 0 && !cif_null(&tk)) { sv_[f] = tk.s; svn_[f] = tk.n; have_single = TRUE; }
				continue;
			}
			if (catn == 20 && cif_iskw(tag, catn, "_atom_site_anisotrop", 20)) {
				st.n_anisou++;
				continue;
			}
			if (keep_meta) cif_store_lc(aTHX_ cif, tag, tagn, cif_sv(aTHX_ &tk));
			continue;
		}
		//a stray value with no tag in front of it: nothing to attach it to
	}

	if (have_single) cif_atom_row(aTHX_ &st, sv_, svn_);

	if (st.have_res) {
		av_push(st.res_last, newSVuv(st.n_atom - 1));
		flush_residue(aTHX_ res_sum, st.rsx, st.rsy, st.rsz, st.rnc, st.rsb, st.rnb);
	}

	(void)hv_stores(out, "atoms", newRV_noinc((SV *)st.atom_hv));
	for (i = 0; i < NRSUM; i++)
		(void)hv_store(out, res_sum_name[i], (I32)strlen(res_sum_name[i]),
		               newRV_noinc((SV *)res_sum[i]), 0);
	(void)hv_stores(out, "elements",      newRV_noinc((SV *)st.elements));
	(void)hv_stores(out, "chain_elements", newRV_noinc((SV *)st.chain_elements));
	(void)hv_stores(out, "n_hydrogens",   newSVuv(st.n_hydrogen));
	(void)hv_stores(out, "n_water_atoms", newSVuv(st.n_water_atom));
	if (st.bn) {
		HV *b = newHV();
		(void)hv_stores(b, "min",  newSVnv(st.bmin));
		(void)hv_stores(b, "max",  newSVnv(st.bmax));
		(void)hv_stores(b, "mean", newSVnv(st.bsum / (NV)st.bn));
		(void)hv_stores(b, "n",    newSVuv(st.bn));
		(void)hv_stores(out, "bfactor_stats", newRV_noinc((SV *)b));
	} else {
		(void)hv_stores(out, "bfactor_stats", newSVsv(&PL_sv_undef));
	}
	if (st.have_bbox) {
		HV *bb = newHV();
		AV *ctr = newAV();
		(void)hv_stores(bb, "xmin", newSVnv(st.xmin));
		(void)hv_stores(bb, "ymin", newSVnv(st.ymin));
		(void)hv_stores(bb, "zmin", newSVnv(st.zmin));
		(void)hv_stores(bb, "xmax", newSVnv(st.xmax));
		(void)hv_stores(bb, "ymax", newSVnv(st.ymax));
		(void)hv_stores(bb, "zmax", newSVnv(st.zmax));
		av_push(ctr, newSVnv((st.xmin + st.xmax) / 2));
		av_push(ctr, newSVnv((st.ymin + st.ymax) / 2));
		av_push(ctr, newSVnv((st.zmin + st.zmax) / 2));
		(void)hv_stores(out, "bbox",   newRV_noinc((SV *)bb));
		(void)hv_stores(out, "center", newRV_noinc((SV *)ctr));
	} else {
		(void)hv_stores(out, "bbox",   newSVsv(&PL_sv_undef));
		(void)hv_stores(out, "center", newSVsv(&PL_sv_undef));
	}

	for (i = 0; i < NCOL; i++) {
		if (i == C_LINENO) { SvREFCNT_dec((SV *)col[i]); continue; }
		(void)hv_store(out, col_name[i], (I32)strlen(col_name[i]),
		               newRV_noinc((SV *)col[i]), 0);
	}
	(void)hv_stores(out, "res_first",     newRV_noinc((SV *)st.res_first));
	(void)hv_stores(out, "res_last",      newRV_noinc((SV *)st.res_last));
	//TER is a PDB record and has no mmCIF counterpart; the key is here so that
	//a caller reading it does not have to know which format the file was
	(void)hv_stores(out, "ter",           newRV_noinc((SV *)newAV()));
	(void)hv_stores(out, "meta",          newRV_noinc((SV *)meta));
	(void)hv_stores(out, "cif",           newRV_noinc((SV *)cif));
	(void)hv_stores(out, "cif_loops",     newRV_noinc((SV *)loops));
	(void)hv_stores(out, "data_block",    block ? block : newSVsv(&PL_sv_undef));
	(void)hv_stores(out, "model_numbers", newRV_noinc((SV *)st.model_nums));
	(void)hv_stores(out, "n_atoms",       newSVuv(st.n_atom));
	(void)hv_stores(out, "n_residues",    newSVuv((UV)(av_len(st.res_first) + 1)));
	(void)hv_stores(out, "n_models",      newSVuv(st.n_models ? st.n_models : 1));
	(void)hv_stores(out, "n_anisou",      newSVuv(st.n_anisou));
	(void)hv_stores(out, "n_skipped",     newSVuv(st.n_skipped));
	(void)hv_stores(out, "n_atom_records",   newSVuv(st.n_atom_rec));
	(void)hv_stores(out, "n_hetatm_records", newSVuv(st.n_het_rec));
	(void)hv_stores(out, "n_lines",       newSVuv(lineno));
	(void)hv_stores(out, "format",        newSVpvs("mmcif"));
	return out;
}

/*slurp() -- read the whole file.  Chunked rather than stat-then-read so that
a named pipe or /dev/stdin works the same as a file on disk.*/
static char *slurp(pTHX_ const char *CSP_RESTRICT path, STRLEN *CSP_RESTRICT lenp)
{
	FILE *CSP_RESTRICT fh = fopen(path, "rb");
	char *CSP_RESTRICT buf;
	/*1 MB to start and doubling from there: an ordinary entry is one
	allocation, and the 33 MB largest in PDBbind v2020 is six.*/
	STRLEN cap = 1 << 20, len = 0;
	if (!fh) croak("Chem::Structure::Parser: cannot read '%s': %s", path, Strerror(errno));
	Newx(buf, cap, char);
	for (;;) {
		size_t got;
		if (len == cap) {
			cap *= 2;
			Renew(buf, cap, char);
		}
		got = fread(buf + len, 1, cap - len, fh);
		len += got;
		if (got == 0) break;
	}
	if (ferror(fh)) {
		int e = errno;
		Safefree(buf);
		fclose(fh);
		croak("Chem::Structure::Parser: error reading '%s': %s", path, Strerror(e));
	}
	fclose(fh);
	*lenp = len;
	return buf;
}

/*asking a parsed chain a question.

Reading the hash of hashes back is Perl's job nine times in ten -- one lookup
is one line of Perl and nothing to be gained by crossing into C for it.  A
question asked once per chain of every chain of every file in a set is the
exception: a structure with a dozen ions has more chains that are one ion than
chains that are a polymer, and the loop that has to tell them apart runs the
whole way down the set.

Nothing here re-does chemistry, and nothing here consults it either: the answer
is in the shape of the chain.*/

/*the chain hash the caller means, from either shape of argument: a chain on
its own, or the structure and a chain id.  Every wrong argument that has an
obvious right one nearby says so, since the three hashes in a parsed structure
-- structure, chain, residue -- all look alike from inside a hash reference and
a false answer would be taken at face value.*/
static HV *chain_arg(pTHX_ SV *CSP_RESTRICT first, SV *CSP_RESTRICT id, const char *CSP_RESTRICT who)
{
	HV *h;
	if (!SvROK(first) || SvTYPE(SvRV(first)) != SVt_PVHV)
		croak("%s: expected the chain hash reference from structure_info()", who);
	h = (HV *)SvRV(first);
	if (SvOK(id)) {
		HE *he;
		SV **p = hv_fetchs(h, "chains", 0);
		if (!p || !*p || !SvROK(*p) || SvTYPE(SvRV(*p)) != SVt_PVHV)
			croak("%s: expected the hash reference from structure_info()", who);
		he = hv_fetch_ent((HV *)SvRV(*p), id, 0, 0);
		if (!he || !SvROK(HeVAL(he)) || SvTYPE(SvRV(HeVAL(he))) != SVt_PVHV)
			croak("%s: no chain '%s'", who, SvPV_nolen(id));
		h = (HV *)SvRV(HeVAL(he));
	}
	//residue_order is what every chain has and nothing else in the structure does
	if (!hv_exists(h, "residue_order", 13)) {
		if (hv_exists(h, "chains", 6))
			croak("%s: that is the whole structure, not a chain: pass one of $info->{chains}, or call %s($info, $chain)", who, who);
		if (hv_exists(h, "atom_order", 10))
			croak("%s: that is a residue, not a chain: a residue answers for itself, in $res->{type}", who);
		croak("%s: expected the chain hash reference from structure_info()", who);
	}
	return h;
}

/*is_single_ion() -- true when a chain holds exactly one residue.

An ion is often numbered into the chain it sits in -- the zinc of a zinc finger
is residue 202 of chain A -- and just as often given a chain of its own, which
is a chain with one residue in it and no sequence to read.  A caller walking
chain_order wants those out of the way before it asks for a sequence or an
alignment, and the count is the whole of the question: single means one residue
in the chain, so a chain of two zincs is not one and a zinc numbered into a
protein chain does not make that chain one.

Single says nothing about how many atoms the residue has.  A sulphate is five
atoms and a perchlorate is five different ones, and both are one residue, so
both answer the same -- which they would not if the count were of atoms, and
would not either if this asked whether the residue types as an ion.  That type
comes off a table of names, and a table of names cannot be complete: SO4 is on
the module's list and BF4 is not, and the difference between them is which ions
someone got round to writing down, not anything about the file.  Counting
residues asks nothing of the table and so cannot inherit its gaps.

The residue is therefore not asked what it is, and the answer is about the shape
of the chain rather than the chemistry in it: a chain that is one sugar, one
buffer molecule, one water or one free amino acid answers true as well.  In a
real file those are rare next to the ions and are the same nuisance to a caller
walking chains -- one residue, no sequence -- but a caller who needs to tell
them apart has it a lookup away, in the residue's own type.*/
static bool chain_is_single_ion(pTHX_ HV *CSP_RESTRICT c)
{
	SV **p = hv_fetchs(c, "residue_order", 0);
	if (!p || !*p || !SvROK(*p) || SvTYPE(SvRV(*p)) != SVt_PVAV) return FALSE;
	/*residue_order rather than n_residues or the residues hash: it is the
	chain's own list of what is in it, in the order the file had them, and
	counting it needs no iterator reset that a caller mid-each() would feel.
	av_len() is the last index, so 0 is one residue and -1 is none.*/
	return av_len((AV *)SvRV(*p)) == 0;
}

/*Physical properties of a structure.

Everything from here to the MODULE line answers questions about a structure
that has already been read: how much of it the solvent can touch, how big it
is, how heavy it is, and which of its aromatic rings are stacked on each
other.  None of it is part of the parse -- it is asked for by name, on an
$info that already exists -- but all of it has to touch every atom, so by the
rule in notes.txt it belongs here rather than in Perl.

The walk is the reason.  A structure is a hash of chains of residues of atoms,
and the calculations below want it as flat arrays of coordinates: doing that
in Perl would mean building a second copy of the coordinate section, which is
the mistake the second round of benchmarking found and took out.  set_build()
walks the hash once into four NV arrays and hands them to whichever of the
calculations were asked for, so a caller who wants all of them pays for one
walk.

Provenance -- every number below came from somebody else's implementation, and
the tests compare against those implementations rather than against this one:

  SASA          Shrake, A; Rupley, JA (1973) J Mol Biol 79(2):351-71, as
                mdtraj 1.11's mdtraj.geometry.shrake_rupley (sasa.py and the
                C kernel behind mdtraj.geometry._geometry._sasa): golden
                section spiral sphere points, a probe rolled over the van der
                Waals surface, the accessible fraction of each atom's points.
  vdW radii     mdtraj's _ATOMIC_RADII (mdtraj/geometry/sasa.py), which is
                Bondi, A (1964) J Phys Chem 68:441 as extended by Mantina, M
                et al. (2009) J Phys Chem A 113:5806, with Shannon, R D (1976)
                Acta Cryst A32:751 ionic radii substituted for the ions that
                are always ionised in a biophysical setting (Li Na K Cs Be Mg
                Ca Ba Cl), and 2.0 A for an element with no measured value.
  masses        mdtraj/core/element.py.  The five elements it names no mass
                for are marked below.
  pi-stacking   mdtraj.geometry.pi_stacking (mdtraj/geometry/pi_stacking.py),
                which is itself ProLIF's FaceToFace and EdgeToFace geometry.
  max ASA       Tien, M Z; Meyer, A G; Sydykova, D K; Spielman, S J; Wilke, C O
                (2013) PLoS ONE 8(11):e80635, Table 1, the "Theoretical"
                column, for the relative accessibility of a residue.

mdtraj works in nanometres and float32; this file works in angstrom and NV,
which is what the PDB and mmCIF formats are written in and what the rest of
the module already returns.  The formulae are the same ones -- an area is
4*pi*r^2 times the accessible fraction whichever unit r is in -- so the answers
agree to the width of a float32, which is what t/features.t measures rather
than assumes.*/

/*elem_prop_of() -- the van der Waals radius and standard atomic weight of an
element symbol, or false for a field that spells no element.

Keyed the way elem_case() is, on the symbol's one or two bytes packed into a
U32 and upper-cased first, so that a symbol arriving as "ZN", "Zn" or "zn" is
one lookup.  The four heaviest are spelled Nh, Mc, Ts and Og -- IUPAC's 2016
names, which is what elem_case() puts in an $info -- where mdtraj's tables
still carry the placeholder Uut, Uup, Uus and Uuo.*/
typedef struct {
	NV vdw;  //van der Waals radius, angstrom
	NV mass; //standard atomic weight, dalton; 0.0 where mdtraj names none
	U32 key; //the symbol packed and upper-cased, so a caller can ask which it was
} elem_prop;

//mdtraj's own fallback: "Where no van der Waals value is known, a default of
//2 angstroms is used"
#define CSP_VDW_DEFAULT 2.0

static bool elem_prop_of(const char *CSP_RESTRICT s, STRLEN len,
                         elem_prop *CSP_RESTRICT out)
{
	//initialised because the compiler cannot see that every case assigns and
	//the default returns; the values are never the ones used
	NV r = CSP_VDW_DEFAULT, m = 0.0;
	U32 key;
	while (len && (*s == ' ' || *s == '\t')) { s++; len--; }
	while (len && (s[len - 1] == ' ' || s[len - 1] == '\t')) len--;
	if (len == 0 || len > 2) return FALSE;
	key = (len == 1)
	    ? K2((unsigned char)toupper((unsigned char)s[0]), 0)
	    : K2((unsigned char)toupper((unsigned char)s[0]),
	         (unsigned char)toupper((unsigned char)s[1]));
	switch (key) {
		case K2('B',0):       r =   1.92; m = 10.8117;     break;
		case K2('C',0):       r =   1.70; m = 12.01078;    break;
		case K2('F',0):       r =   1.47; m = 18.99840325; break;
		case K2('H',0):       r =   1.20; m = 1.007947;    break;
		case K2('I',0):       r =   1.98; m = 126.904473;  break;
		case K2('K',0):       r =   1.38; m = 39.09831;    break;
		case K2('N',0):       r =   1.55; m = 14.00672;    break;
		case K2('O',0):       r =   1.52; m = 15.99943;    break;
		case K2('P',0):       r =   1.80; m = 30.9737622;  break;
		case K2('S',0):       r =   1.80; m = 32.0655;     break;
		case K2('U',0):       r =   1.86; m = 238.028913;  break;
		case K2('V',0):       r =   2.00; m = 50.94151;    break;
		case K2('W',0):       r =   2.00; m = 183.841;     break;
		case K2('Y',0):       r =   2.00; m = 88.905852;   break;
		case K2('A','C'):     r =   2.00; m = 227.0;       break;
		case K2('A','G'):     r =   1.72; m = 107.86822;   break;
		case K2('A','L'):     r =   1.84; m = 26.98153868; break;
		case K2('A','M'):     r =   2.00; m = 243.0;       break;
		case K2('A','R'):     r =   1.88; m = 39.9481;     break;
		case K2('A','S'):     r =   1.85; m = 74.921602;   break;
		case K2('A','T'):     r =   2.02; m = 210.0;       break;
		case K2('A','U'):     r =   1.66; m = 196.9665694; break;
		case K2('B','A'):     r =   1.49; m = 137.3277;    break;
		case K2('B','E'):     r =   0.59; m = 9.0121823;   break;
		case K2('B','H'):     r =   2.00; m = 264.0;       break;
		case K2('B','I'):     r =   2.07; m = 208.980401;  break;
		case K2('B','K'):     r =   2.00; m = 247.0;       break;
		case K2('B','R'):     r =   1.85; m = 79.9041;     break;
		case K2('C','A'):     r =   1.14; m = 40.0784;     break;
		case K2('C','D'):     r =   1.58; m = 112.4118;    break;
		case K2('C','E'):     r =   2.00; m = 140.1161;    break;
		case K2('C','F'):     r =   2.00; m = 251.0;       break;
		case K2('C','L'):     r =   1.81; m = 35.4532;     break;
		case K2('C','M'):     r =   2.00; m = 247.0;       break;
		case K2('C','N'):     r =   2.00; m = 0.0;         break;
		case K2('C','O'):     r =   2.00; m = 58.9331955;  break;
		case K2('C','R'):     r =   2.00; m = 51.99616;    break;
		case K2('C','S'):     r =   1.67; m = 132.90545192; break;
		case K2('C','U'):     r =   1.40; m = 63.5463;     break;
		case K2('D','B'):     r =   2.00; m = 262.0;       break;
		case K2('D','S'):     r =   2.00; m = 281.0;       break;
		case K2('D','Y'):     r =   2.00; m = 162.5001;    break;
		case K2('E','R'):     r =   2.00; m = 167.2593;    break;
		case K2('E','S'):     r =   2.00; m = 252.0;       break;
		case K2('E','U'):     r =   2.00; m = 151.9641;    break;
		case K2('F','E'):     r =   2.00; m = 55.8452;     break;
		case K2('F','L'):     r =   2.00; m = 0.0;         break;
		case K2('F','M'):     r =   2.00; m = 257.0;       break;
		case K2('F','R'):     r =   3.48; m = 223.0;       break;
		case K2('G','A'):     r =   1.87; m = 69.7231;     break;
		case K2('G','D'):     r =   2.00; m = 157.253;     break;
		case K2('G','E'):     r =   2.11; m = 72.641;      break;
		case K2('H','E'):     r =   1.40; m = 4.003;       break;
		case K2('H','F'):     r =   2.00; m = 178.492;     break;
		case K2('H','G'):     r =   1.55; m = 200.592;     break;
		case K2('H','O'):     r =   2.00; m = 164.930322;  break;
		case K2('H','S'):     r =   2.00; m = 269.0;       break;
		case K2('I','N'):     r =   1.93; m = 114.8183;    break;
		case K2('I','R'):     r =   2.00; m = 192.2173;    break;
		case K2('K','R'):     r =   2.02; m = 83.7982;     break;
		case K2('L','A'):     r =   2.00; m = 138.905477;  break;
		case K2('L','I'):     r =   0.76; m = 6.9412;      break;
		case K2('L','R'):     r =   2.00; m = 262.0;       break;
		case K2('L','U'):     r =   2.00; m = 174.9671;    break;
		case K2('L','V'):     r =   2.00; m = 0.0;         break;
		case K2('M','C'):     r =   2.00; m = 288.0;       break;
		case K2('M','D'):     r =   2.00; m = 258.0;       break;
		case K2('M','G'):     r =   0.86; m = 24.30506;    break;
		case K2('M','N'):     r =   2.00; m = 54.9380455;  break;
		case K2('M','O'):     r =   2.00; m = 95.942;      break;
		case K2('M','T'):     r =   2.00; m = 268.0;       break;
		case K2('N','A'):     r =   1.02; m = 22.989769282; break;
		case K2('N','B'):     r =   2.00; m = 92.906382;   break;
		case K2('N','D'):     r =   2.00; m = 144.2423;    break;
		case K2('N','E'):     r =   1.54; m = 20.17976;    break;
		case K2('N','H'):     r =   2.00; m = 284.0;       break;
		case K2('N','I'):     r =   1.63; m = 58.69342;    break;
		case K2('N','O'):     r =   2.00; m = 259.0;       break;
		case K2('N','P'):     r =   2.00; m = 237.0;       break;
		case K2('O','G'):     r =   2.00; m = 0.0;         break;
		case K2('O','S'):     r =   2.00; m = 190.233;     break;
		case K2('P','A'):     r =   2.00; m = 231.035882;  break;
		case K2('P','B'):     r =   2.02; m = 207.21;      break;
		case K2('P','D'):     r =   1.63; m = 106.421;     break;
		case K2('P','M'):     r =   2.00; m = 145.0;       break;
		case K2('P','O'):     r =   1.97; m = 209.0;       break;
		case K2('P','R'):     r =   2.00; m = 140.907652;  break;
		case K2('P','T'):     r =   1.75; m = 195.0849;    break;
		case K2('P','U'):     r =   2.00; m = 244.0;       break;
		case K2('R','A'):     r =   2.83; m = 226.0;       break;
		case K2('R','B'):     r =   3.03; m = 85.46783;    break;
		case K2('R','E'):     r =   2.00; m = 186.2071;    break;
		case K2('R','F'):     r =   2.00; m = 261.0;       break;
		case K2('R','G'):     r =   2.00; m = 272.0;       break;
		case K2('R','H'):     r =   2.00; m = 102.905502;  break;
		case K2('R','N'):     r =   2.20; m = 222.018;     break;
		case K2('R','U'):     r =   2.00; m = 101.072;     break;
		case K2('S','B'):     r =   2.06; m = 121.7601;    break;
		case K2('S','C'):     r =   2.11; m = 44.9559126;  break;
		case K2('S','E'):     r =   1.90; m = 78.963;      break;
		case K2('S','G'):     r =   2.00; m = 266.0;       break;
		case K2('S','I'):     r =   2.10; m = 28.08553;    break;
		case K2('S','M'):     r =   2.00; m = 150.362;     break;
		case K2('S','N'):     r =   2.17; m = 118.7107;    break;
		case K2('S','R'):     r =   2.49; m = 87.621;      break;
		case K2('T','A'):     r =   2.00; m = 180.947882;  break;
		case K2('T','B'):     r =   2.00; m = 158.925352;  break;
		case K2('T','C'):     r =   2.00; m = 98.0;        break;
		case K2('T','E'):     r =   2.06; m = 127.603;     break;
		case K2('T','H'):     r =   2.00; m = 232.038062;  break;
		case K2('T','I'):     r =   2.00; m = 47.8671;     break;
		case K2('T','L'):     r =   1.96; m = 204.38332;   break;
		case K2('T','M'):     r =   2.00; m = 168.934212;  break;
		case K2('T','S'):     r =   2.00; m = 0.0;         break;
		case K2('X','E'):     r =   2.16; m = 131.2936;    break;
		case K2('Y','B'):     r =   2.00; m = 173.043;     break;
		case K2('Z','N'):     r =   1.39; m = 65.4094;     break;
		case K2('Z','R'):     r =   2.00; m = 91.2242;     break;
		default: return FALSE; //not a symbol: the caller decides what to do
	}
	out->vdw  = r;
	out->mass = m;
	out->key  = key;
	return TRUE;
}

/*max_asa[] -- the largest solvent-accessible surface a residue of each kind
can have, in angstrom^2, indexed by its single-letter code less 'A'.

Dividing a residue's SASA by this gives its relative accessibility, which is
the number the buried/exposed question is actually asked of: 130 A^2 is most of
an alanine and a sliver of a tryptophan.  The values are the "Theoretical"
column of Table 1 of Tien et al. (2013) -- a Gly-X-Gly tripeptide extended to
its maximum -- rather than the empirical column of the same table, because the
empirical one is the largest value seen in a particular set of structures and
so goes stale as the archive grows.

0.0 means the letter has no value: the six ambiguity and placeholder codes
(B J O U X Z), and the twenty letters that are not amino acid codes at all.
A residue with no value gets no relative accessibility rather than a wrong
one.*/
static const NV max_asa[26] = {
	129.0,   0.0, 167.0, 193.0, 223.0, 240.0, 104.0, 224.0, 197.0,   0.0, //A-J
	236.0, 201.0, 224.0, 195.0,   0.0, 159.0, 225.0, 274.0, 155.0, 172.0, //K-T
	  0.0, 174.0, 285.0,   0.0, 263.0,   0.0                              //U-Z
};

/*Aromatic rings, by residue name.

A ring is a list of atom names, and the pi-stacking geometry needs it in a
fixed order: mdtraj's compute_ring_normal() takes the plane's normal from the
cross product of the first two atoms' offsets from the centroid, so two readers
that list the same ring in a different order get normals that differ by
whatever the ring departs from planarity in a real structure.  The order below
is the one t/data/features.py hands mdtraj, which is what makes the two
comparable.

Only the residues whose ring is fixed by the format's own atom naming are
here: the four aromatic amino acids and the nucleobases.  A ligand's rings
would need bond perception, which this module does not do -- it never reads a
CONECT record or guesses a bond -- so a ligand contributes no rings and
structure_pi_stacking() says so rather than quietly finding none.

HID/HIE/HIP and HSD/HSE/HSP are histidine under the names AMBER and CHARMM
give its three protonation states; a structure that has been through either
comes back with those spellings and the same five ring atoms.*/
#define RING_MAX 6 //the most atoms any ring in the table below has

typedef struct {
	const char *const *atom; //ring atom names, normal taken from the first two
	unsigned short int n;    //how many
	char label;              //'6' six-membered, '5' five-membered
} ring_def;

static const char *const ring_phe[6] = { "CG", "CD1", "CD2", "CE1", "CE2", "CZ" };
static const char *const ring_his[5] = { "CG", "ND1", "CE1", "NE2", "CD2" };
static const char *const ring_trp5[5] = { "CG", "CD1", "NE1", "CE2", "CD2" };
static const char *const ring_trp6[6] = { "CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3" };
static const char *const ring_pur6[6] = { "N1", "C2", "N3", "C4", "C5", "C6" };
static const char *const ring_pur5[5] = { "C4", "C5", "N7", "C8", "N9" };
static const char *const ring_pyr6[6] = { "N1", "C2", "N3", "C4", "C5", "C6" };

static const ring_def rings_phe[1] = { { ring_phe, 6, '6' } };
static const ring_def rings_his[1] = { { ring_his, 5, '5' } };
static const ring_def rings_trp[2] = { { ring_trp6, 6, '6' }, { ring_trp5, 5, '5' } };
static const ring_def rings_pur[2] = { { ring_pur6, 6, '6' }, { ring_pur5, 5, '5' } };
static const ring_def rings_pyr[1] = { { ring_pyr6, 6, '6' } };

//how many rings a residue name has, and where they are; 0 for everything else
static unsigned short int ring_defs(U32 key, const ring_def *CSP_RESTRICT *CSP_RESTRICT out)
{
	switch (key) {
		case K3('P','H','E'): case K3('T','Y','R'):
			*out = rings_phe; return 1;
		case K3('H','I','S'):
		case K3('H','I','D'): case K3('H','I','E'): case K3('H','I','P'):
		case K3('H','S','D'): case K3('H','S','E'): case K3('H','S','P'):
			*out = rings_his; return 1;
		case K3('T','R','P'):
			*out = rings_trp; return 2;
		//purines: adenine and guanine, DNA and RNA
		case K3(' ',' ','A'): case K3(' ',' ','G'):
		case K3(' ','D','A'): case K3(' ','D','G'):
			*out = rings_pur; return 2;
		//pyrimidines: cytosine, thymine and uracil, DNA and RNA
		case K3(' ',' ','C'): case K3(' ',' ','T'): case K3(' ',' ','U'):
		case K3(' ','D','C'): case K3(' ','D','T'): case K3(' ','D','U'):
			*out = rings_pyr; return 1;
		default: *out = NULL; return 0;
	}
}

//hash field accessors, for walking an $info that Perl built
static SV *hvf_sv(pTHX_ HV *CSP_RESTRICT h, const char *CSP_RESTRICT k, STRLEN klen)
{
	SV **p = hv_fetch(h, k, (I32)klen, 0);
	return (p && *p && SvOK(*p)) ? *p : NULL;
}

static AV *hvf_av(pTHX_ HV *CSP_RESTRICT h, const char *CSP_RESTRICT k, STRLEN klen)
{
	SV *s = hvf_sv(aTHX_ h, k, klen);
	return (s && SvROK(s) && SvTYPE(SvRV(s)) == SVt_PVAV) ? (AV *)SvRV(s) : NULL;
}

static HV *hvf_hv(pTHX_ HV *CSP_RESTRICT h, const char *CSP_RESTRICT k, STRLEN klen)
{
	SV *s = hvf_sv(aTHX_ h, k, klen);
	return (s && SvROK(s) && SvTYPE(SvRV(s)) == SVt_PVHV) ? (HV *)SvRV(s) : NULL;
}

//the hash an array element refers to, or NULL if it is not a hash reference
static HV *hvf_ent_hv(pTHX_ HV *CSP_RESTRICT h, SV *CSP_RESTRICT key)
{
	HE *e = hv_fetch_ent(h, key, 0, 0);
	SV *v = e ? HeVAL(e) : NULL;
	return (v && SvROK(v) && SvTYPE(SvRV(v)) == SVt_PVHV) ? (HV *)SvRV(v) : NULL;
}

/*The coordinate section of an $info, flattened.

Atoms are in walk order, which is chain by chain, residue by residue and, inside
a residue, the order the file wrote them: so a residue's atoms are one
contiguous run and a chain's residues are another, and the roll-ups below are
sums over ranges rather than lookups.  Same shape the parse itself hands back,
for the same reason.*/
typedef struct {
	//per atom
	NV *x, *y, *z;
	NV *rad;      //van der Waals radius plus the probe, angstrom
	NV *area;     //solvent-accessible surface, angstrom^2; NULL until computed
	NV *alone;    //the same with the other chains taken away; NULL unless asked for
	NV *mass;     //dalton; 0.0 for an atom whose element has no mass
	unsigned char *apolar; //1 = carbon or sulphur, 0 = everything else
	HV **atom_hv; //the atom hash the row was read from, for writing back
	UV n_atom;
	//per residue: the atom range [first, last)
	HV **res_hv;
	UV *res_first, *res_last;
	char *res_one; //single-letter code, '\0' when the residue has none
	unsigned char *res_type; //RT_*, for the questions only an amino acid answers
	unsigned char *res_std;  //1 = one of the twenty, as the residue hash says
	U32 *res_key;  //packed residue name, for ring_defs()
	UV n_res;
	/*Kabsch-Sander: each residue's two best acceptors and their energies, laid
	out two per residue.  An entry of n_res means there is no bond there.*/
	UV *ks_acc;
	NV *ks_e;
	//per chain: the residue range [first, last)
	HV **chain_hv;
	UV *chain_first, *chain_last;
	UV n_chain;
	//tallies gathered on the way past
	NV mass_total;
	UV n_no_element; //atoms whose element field spells no element
	/*How many names the residues' atom_order lists held, which is not n_atom:
	that one counts the atoms a position could be read from.  The two differ for
	a file whose coordinate columns are unreadable, and telling them apart is
	what keeps such a file from being mistaken for one read with atoms => 0.*/
	UV n_atom_slots;
} structset;

static void set_free(pTHX_ structset *CSP_RESTRICT s)
{
	Safefree(s->x);        Safefree(s->y);         Safefree(s->z);
	Safefree(s->rad);      Safefree(s->area);      Safefree(s->alone);
	Safefree(s->mass);
	Safefree(s->apolar);   Safefree(s->atom_hv);
	Safefree(s->res_hv);   Safefree(s->res_first); Safefree(s->res_last);
	Safefree(s->res_one);  Safefree(s->res_type);  Safefree(s->res_std);
	Safefree(s->res_key);
	Safefree(s->ks_acc);   Safefree(s->ks_e);
	Safefree(s->chain_hv); Safefree(s->chain_first); Safefree(s->chain_last);
	Zero(s, 1, structset);
}

/*set_build() -- walk an $info into the arrays above.

Sized from what the structure already says about itself rather than by counting
first: a chain's n_atoms counts records, so it is an upper bound on how many
distinct atoms its residues hold (alternate conformers of one atom are one
entry in the residue's atoms hash and several records), and the residue counts
are exact.  One over-allocation, trimmed with Renew afterwards, against a
second walk of the same hashes.

Croaks rather than returning an empty set when there are atom records to be had
and no atom hashes to read them from: that is a structure read with atoms => 0,
and the only thing wrong with it is that nobody said so.*/
static void set_build(pTHX_ HV *CSP_RESTRICT info, structset *CSP_RESTRICT s,
                      NV probe, const char *CSP_RESTRICT who)
{
	AV *order = hvf_av(aTHX_ info, "chain_order", 11);
	HV *chains = hvf_hv(aTHX_ info, "chains", 6);
	UV cap_atom = 0, cap_res = 0, n_chain = 0;
	SSize_t ci, nci;

	Zero(s, 1, structset);
	if (!order || !chains) return;
	nci = av_len(order) + 1;

	//first pass: the chains only, which is a few dozen fetches, for the sizes
	for (ci = 0; ci < nci; ci++) {
		SV **cs = av_fetch(order, ci, 0);
		HV *c;
		AV *ro;
		SV *na;
		if (!cs || !*cs || !SvOK(*cs)) continue;
		c = hvf_ent_hv(aTHX_ chains, *cs);
		if (!c) continue;
		n_chain++;
		ro = hvf_av(aTHX_ c, "residue_order", 13);
		if (ro) cap_res += (UV)(av_len(ro) + 1);
		na = hvf_sv(aTHX_ c, "n_atoms", 7);
		if (na) cap_atom += SvUV(na);
	}
	if (n_chain == 0) return;

	Newx(s->x,   cap_atom ? cap_atom : 1, NV);
	Newx(s->y,   cap_atom ? cap_atom : 1, NV);
	Newx(s->z,   cap_atom ? cap_atom : 1, NV);
	Newx(s->rad, cap_atom ? cap_atom : 1, NV);
	Newx(s->mass, cap_atom ? cap_atom : 1, NV);
	Newx(s->apolar, cap_atom ? cap_atom : 1, unsigned char);
	Newx(s->atom_hv, cap_atom ? cap_atom : 1, HV *);
	Newx(s->res_hv,    cap_res ? cap_res : 1, HV *);
	Newx(s->res_first, cap_res ? cap_res : 1, UV);
	Newx(s->res_last,  cap_res ? cap_res : 1, UV);
	Newx(s->res_one,   cap_res ? cap_res : 1, char);
	Newx(s->res_type,  cap_res ? cap_res : 1, unsigned char);
	Newx(s->res_std,   cap_res ? cap_res : 1, unsigned char);
	Newx(s->res_key,   cap_res ? cap_res : 1, U32);
	Newx(s->chain_hv,    n_chain, HV *);
	Newx(s->chain_first, n_chain, UV);
	Newx(s->chain_last,  n_chain, UV);

	for (ci = 0; ci < nci; ci++) {
		SV **cs = av_fetch(order, ci, 0);
		HV *c;
		AV *ro;
		HV *residues;
		SSize_t ri, nri;
		if (!cs || !*cs || !SvOK(*cs)) continue;
		c = hvf_ent_hv(aTHX_ chains, *cs);
		if (!c) continue;
		s->chain_hv[s->n_chain] = c;
		s->chain_first[s->n_chain] = s->n_res;
		ro = hvf_av(aTHX_ c, "residue_order", 13);
		residues = hvf_hv(aTHX_ c, "residues", 8);
		nri = (ro && residues) ? av_len(ro) + 1 : 0;
		for (ri = 0; ri < nri; ri++) {
			SV **rs = av_fetch(ro, ri, 0);
			HV *r;
			AV *ao;
			HV *atoms;
			SV *one;
			SV *rn;
			SSize_t ai, nai;
			if (!rs || !*rs || !SvOK(*rs)) continue;
			r = hvf_ent_hv(aTHX_ residues, *rs);
			if (!r || s->n_res >= cap_res) continue;
			s->res_hv[s->n_res] = r;
			s->res_first[s->n_res] = s->n_atom;
			one = hvf_sv(aTHX_ r, "one", 3);
			if (one) {
				STRLEN olen;
				const char *op = SvPV_const(one, olen);
				//a residue with no single-letter code has an empty string here,
				//not a missing key, so length is what tells them apart
				s->res_one[s->n_res] = (olen == 1) ? *op : '\0';
			} else {
				s->res_one[s->n_res] = '\0';
			}
			rn = hvf_sv(aTHX_ r, "resname", 7);
			if (rn) {
				STRLEN rl;
				res_info ri;
				const char *rp = SvPV_const(rn, rl);
				s->res_key[s->n_res] = res_key(rp, rl);
				s->res_type[s->n_res] = res_lookup(rp, rl, &ri) ? ri.type : RT_OTHER;
			} else {
				s->res_key[s->n_res] = 0;
				s->res_type[s->n_res] = RT_OTHER;
			}
			{
				SV *st = hvf_sv(aTHX_ r, "standard", 8);
				s->res_std[s->n_res] = (st && SvTRUE(st)) ? 1 : 0;
			}
			ao = hvf_av(aTHX_ r, "atom_order", 10);
			atoms = hvf_hv(aTHX_ r, "atoms", 5);
			nai = (ao && atoms) ? av_len(ao) + 1 : 0;
			s->n_atom_slots += (UV)nai;
			for (ai = 0; ai < nai; ai++) {
				SV **as = av_fetch(ao, ai, 0);
				HV *a;
				SV *xs, *ys, *zs, *es;
				elem_prop ep;
				unsigned char apolar;
				if (!as || !*as || !SvOK(*as)) continue;
				a = hvf_ent_hv(aTHX_ atoms, *as);
				if (!a || s->n_atom >= cap_atom) continue;
				xs = hvf_sv(aTHX_ a, "x", 1);
				ys = hvf_sv(aTHX_ a, "y", 1);
				zs = hvf_sv(aTHX_ a, "z", 1);
				//an atom whose line was truncated before the coordinates has
				//no position, and a position is what every one of these
				//calculations is about; it is left out and counted nowhere
				if (!xs || !ys || !zs) continue;
				es = hvf_sv(aTHX_ a, "element", 7);
				apolar = 0;
				if (es) {
					STRLEN el;
					const char *ep_s = SvPV_const(es, el);
					/*Carbon and sulphur are the apolar surface and everything
					else is the polar one, which is the split Chothia, C (1974)
					Nature 248:338 made when he first added a protein's buried
					surface up.  Only the symbol decides it, so a sulphur that
					is part of a sulphate counts as apolar here; a caller who
					wants a chemistry-aware split has the per-atom areas.

					Asked of the packed key rather than of the field, so that a
					field written " C" or "c" answers the same as "C" -- the
					same normalisation the radius lookup just did.*/
					if (!elem_prop_of(ep_s, el, &ep)) {
						ep.vdw = CSP_VDW_DEFAULT;
						ep.mass = 0.0;
						s->n_no_element++;
					} else if (ep.key == K2('C', 0) || ep.key == K2('S', 0)) {
						apolar = 1;
					}
				} else {
					ep.vdw = CSP_VDW_DEFAULT;
					ep.mass = 0.0;
					s->n_no_element++;
				}
				s->x[s->n_atom] = SvNV(xs);
				s->y[s->n_atom] = SvNV(ys);
				s->z[s->n_atom] = SvNV(zs);
				s->rad[s->n_atom] = ep.vdw + probe;
				s->mass[s->n_atom] = ep.mass;
				s->apolar[s->n_atom] = apolar;
				s->mass_total += ep.mass;
				s->atom_hv[s->n_atom] = a;
				s->n_atom++;
			}
			s->res_last[s->n_res] = s->n_atom;
			s->n_res++;
		}
		s->chain_last[s->n_chain] = s->n_res;
		s->n_chain++;
	}
	/*No atom hashes at all, in a structure whose chains count atom records:
	that is atoms => 0, and the only thing wrong with it is that nobody said so.
	A structure whose atoms are all there but have no readable coordinates is a
	different thing and is not an error -- t/errors.t reads one -- so the test is
	on the names the residues listed, not on the positions that came back.*/
	if (s->n_atom_slots == 0 && cap_atom > 0) {
		set_free(aTHX_ s);
		croak("%s: this structure has no atom hashes to work from; "
		      "read it again without atoms => 0", who);
	}
	//give the over-allocation back before the calculations pile their own
	//working arrays on top of it
	if (s->n_atom < cap_atom && s->n_atom > 0) {
		Renew(s->x, s->n_atom, NV);
		Renew(s->y, s->n_atom, NV);
		Renew(s->z, s->n_atom, NV);
		Renew(s->rad, s->n_atom, NV);
		Renew(s->mass, s->n_atom, NV);
		Renew(s->apolar, s->n_atom, unsigned char);
		Renew(s->atom_hv, s->n_atom, HV *);
	}
}

/*A uniform grid over a set of points, for "which points are within cut of this
one?".

The Shrake-Rupley kernel asks that question once per atom and mdtraj answers it
by comparing every atom against every other, which is fine for the few thousand
atoms an MD frame holds and is not fine here: the largest entry in PDBbind
v2020 is ~400,000 atoms, and 1.6e11 distance comparisons is not a wait anybody
would sit through.  A grid of cells one cutoff wide gives the same neighbour
set from the 27 cells around each atom, in time proportional to the number of
atoms rather than its square.

Same set, not a similar one: an atom is a neighbour when the distance is below
the cutoff, and the cutoff is twice the largest radius, so every neighbour is
inside the 27-cell block whatever the cell size is as long as it is at least
the cutoff.  Which is why the loop that widens the cells to keep their number
down is safe -- it can only make the block bigger.  The order neighbours come
back in does change, and does not matter: the kernel asks whether any of them
covers a point and stops at the first that does.

The cell count is capped at eight per atom so that the grid cannot cost more
memory than the coordinates it indexes -- a thin, extended structure in a large
box would otherwise want more cells than there are atoms to put in them.*/
typedef struct {
	UV *start; //ncell + 1 offsets into idx
	UV *idx;   //point indices, grouped by cell
	NV cell;   //edge length, angstrom
	NV x0, y0, z0;
	UV nx, ny, nz;
} cell_grid;

static void grid_free(pTHX_ cell_grid *CSP_RESTRICT g)
{
	Safefree(g->start);
	Safefree(g->idx);
	Zero(g, 1, cell_grid);
}

//which cell along one axis a coordinate falls in, clamped to the grid.  Written
//so that no non-finite value is ever cast to a UV, which is undefined
//behaviour: a NaN fails `d > 0.0' and an infinity fails `k < n'.
static UV grid_axis(NV d, NV cell, UV n)
{
	NV k;
	if (n <= 1) return 0;
	if (!(d > 0.0)) return 0;
	k = d / cell;
	if (!(k < (NV)n)) return n - 1;
	return (UV)k; //truncation is floor here: k is positive
}

static void grid_build(pTHX_ cell_grid *CSP_RESTRICT g,
                       const NV *CSP_RESTRICT x, const NV *CSP_RESTRICT y,
                       const NV *CSP_RESTRICT z, UV n, NV cut)
{
	NV xmax, ymax, zmax, cells = 0.0;
	UV i, ncell, c;
	Zero(g, 1, cell_grid);
	if (n == 0) return;
	g->x0 = xmax = x[0];
	g->y0 = ymax = y[0];
	g->z0 = zmax = z[0];
	for (i = 1; i < n; i++) {
		if (x[i] < g->x0) g->x0 = x[i]; else if (x[i] > xmax) xmax = x[i];
		if (y[i] < g->y0) g->y0 = y[i]; else if (y[i] > ymax) ymax = y[i];
		if (z[i] < g->z0) g->z0 = z[i]; else if (z[i] > zmax) zmax = z[i];
	}
	g->cell = (cut > 0.0) ? cut : 1.0;
	/*Widen the cells until there are at most eight per point.  Doubling, so a
	box a kilometre across settles in a few dozen turns; the counter is there for
	the case the extent is not a finite number at all, where the test below can
	never come true and the grid falls back to a single cell.*/
	{
		unsigned short int turn;
		NV ex = 0.0, ey = 0.0, ez = 0.0;
		for (turn = 0; turn < 4096; turn++) {
			ex = (xmax - g->x0) / g->cell;
			ey = (ymax - g->y0) / g->cell;
			ez = (zmax - g->z0) / g->cell;
			cells = (ex + 1.0) * (ey + 1.0) * (ez + 1.0);
			if (cells <= (NV)n * 8.0 + 1024.0) break;
			g->cell *= 2.0;
		}
		//the cast is safe only under that test, which bounds each extent by the
		//cell count it just passed
		if (cells <= (NV)n * 8.0 + 1024.0) {
			g->nx = (UV)ex + 1;
			g->ny = (UV)ey + 1;
			g->nz = (UV)ez + 1;
		} else {
			g->nx = g->ny = g->nz = 1;
		}
	}
	ncell = g->nx * g->ny * g->nz;
	Newxz(g->start, ncell + 1, UV);
	Newx(g->idx, n, UV);
	for (i = 0; i < n; i++) {
		c = (grid_axis(x[i] - g->x0, g->cell, g->nx) * g->ny
		   + grid_axis(y[i] - g->y0, g->cell, g->ny)) * g->nz
		   + grid_axis(z[i] - g->z0, g->cell, g->nz);
		g->start[c]++;
	}
	//counting sort: prefix the counts into starts, fill, then shift the ends
	//back into starts, which is one array rather than a second cursor one
	{
		UV run = 0;
		for (c = 0; c < ncell; c++) { UV t = g->start[c]; g->start[c] = run; run += t; }
	}
	for (i = 0; i < n; i++) {
		c = (grid_axis(x[i] - g->x0, g->cell, g->nx) * g->ny
		   + grid_axis(y[i] - g->y0, g->cell, g->ny)) * g->nz
		   + grid_axis(z[i] - g->z0, g->cell, g->nz);
		g->idx[g->start[c]++] = i;
	}
	for (c = ncell; c > 0; c--) g->start[c] = g->start[c - 1];
	g->start[0] = 0;
}

/*sasa_sphere() -- n points spread over the unit sphere by the golden section
spiral, which is what mdtraj's kernel uses.

The spiral is cheap and very nearly even, which is all the algorithm needs: the
area an atom contributes is 4*pi*r^2 times the fraction of its points that no
other atom covers, and that fraction is only as good as the points are evenly
spread.  mdtraj's own note says as much -- points that repelled each other to an
energy minimum would be better and would cost more than the rest of the
calculation.*/
static void sasa_sphere(NV *CSP_RESTRICT pts, UV n)
{
	const NV inc = CSP_PI * (3.0 - nv_sqrt((NV)5.0));
	const NV offset = 2.0 / (NV)n;
	UV i;
	for (i = 0; i < n; i++) {
		NV y = (NV)i * offset - 1.0 + offset / 2.0;
		NV t = 1.0 - y * y;
		NV r = (t > 0.0) ? nv_sqrt(t) : 0.0; //the half-offset keeps |y| < 1; this is in
		NV phi = (NV)i * inc;                //case rounding at n = 1 says otherwise
		pts[3 * i]     = nv_cos(phi) * r;
		pts[3 * i + 1] = y;
		pts[3 * i + 2] = nv_sin(phi) * r;
	}
}

/*sasa_compute() -- Shrake and Rupley, as mdtraj computes it.

Each atom gets a sphere of npts points at its van der Waals radius plus the
probe radius; a point is accessible when no other atom's sphere covers it; the
atom's area is 4*pi*r^2 times the accessible fraction.

Two departures from a literal transcription of mdtraj's kernel, neither of
which changes an answer:

The neighbour search is the grid above rather than a scan of every other atom.

The distance test is squared -- d2 < (ri+rj)^2 rather than sqrt(d2) < ri+rj --
which saves a square root per pair.  The two differ only for a pair whose
separation is within an ulp of the sum of their radii, and a neighbour at
exactly that separation covers none of the other's points: its sphere touches
at one point and covers an open region of measure zero around it.

The rotation of the neighbour list is mdtraj's and is kept: consecutive points
on the spiral are close together, so the atom that covered the last one is the
one most likely to cover this one, and starting the scan there rather than at
the beginning is most of the kernel's speed.*/
static void sasa_kernel(pTHX_ const NV *CSP_RESTRICT x, const NV *CSP_RESTRICT y,
                        const NV *CSP_RESTRICT z, const NV *CSP_RESTRICT rad, UV n,
                        const NV *CSP_RESTRICT pts, UV npts, NV *CSP_RESTRICT area)
{
	UV *CSP_RESTRICT nbr = NULL;
	/*grown by doubling as needed; 64 covers every atom of a protein, where the
	neighbour count runs to the low tens, without a single reallocation*/
	UV nbr_cap = 64;
	cell_grid g;
	NV rmax = 0.0, constant;
	UV i;

	if (n == 0) return;
	for (i = 0; i < n; i++) if (rad[i] > rmax) rmax = rad[i];
	grid_build(aTHX_ &g, x, y, z, n, 2.0 * rmax);
	Newx(nbr, nbr_cap, UV);
	constant = 4.0 * CSP_PI / (NV)npts;

	for (i = 0; i < n; i++) {
		const NV xi = x[i], yi = y[i], zi = z[i], ri = rad[i];
		UV n_nbr = 0, acc = 0, k_closest = 0, j;
		UV cx = grid_axis(xi - g.x0, g.cell, g.nx);
		UV cy = grid_axis(yi - g.y0, g.cell, g.ny);
		UV cz = grid_axis(zi - g.z0, g.cell, g.nz);
		UV ax0 = cx ? cx - 1 : 0, ax1 = (cx + 1 < g.nx) ? cx + 1 : g.nx - 1;
		UV ay0 = cy ? cy - 1 : 0, ay1 = (cy + 1 < g.ny) ? cy + 1 : g.ny - 1;
		UV az0 = cz ? cz - 1 : 0, az1 = (cz + 1 < g.nz) ? cz + 1 : g.nz - 1;
		UV bx, by, bz;
		for (bx = ax0; bx <= ax1; bx++)
		for (by = ay0; by <= ay1; by++)
		for (bz = az0; bz <= az1; bz++) {
			UV c = (bx * g.ny + by) * g.nz + bz;
			UV p;
			for (p = g.start[c]; p < g.start[c + 1]; p++) {
				UV a = g.idx[p];
				NV dx, dy, dz, sum;
				if (a == i) continue;
				dx = x[a] - xi; dy = y[a] - yi; dz = z[a] - zi;
				sum = ri + rad[a];
				if (dx * dx + dy * dy + dz * dz >= sum * sum) continue;
				if (n_nbr == nbr_cap) { nbr_cap *= 2; Renew(nbr, nbr_cap, UV); }
				nbr[n_nbr++] = a;
			}
		}
		for (j = 0; j < npts; j++) {
			const NV px = xi + ri * pts[3 * j];
			const NV py = yi + ri * pts[3 * j + 1];
			const NV pz = zi + ri * pts[3 * j + 2];
			bool open = TRUE;
			UV k;
			for (k = 0; k < n_nbr; k++) {
				UV kp = k_closest + k;
				UV a;
				NV dx, dy, dz, ra;
				if (kp >= n_nbr) kp -= n_nbr;
				a = nbr[kp];
				dx = px - x[a]; dy = py - y[a]; dz = pz - z[a];
				ra = rad[a];
				if (dx * dx + dy * dy + dz * dz < ra * ra) {
					k_closest = kp;
					open = FALSE;
					break;
				}
			}
			if (open) acc++;
		}
		area[i] = (NV)acc * constant * ri * ri;
	}
	Safefree(nbr);
	grid_free(aTHX_ &g);
}

/*sasa_compute() -- the whole structure, and then each chain on its own.

The second half is what an interface costs: a chain's surface with the rest of
the structure taken away, less the surface it has with the rest of the structure
there, is the area the two bury between them.  It is one more kernel run per
chain, and a chain's run touches only that chain's atoms, so all of them
together cost about what the first run cost -- not the number of chains times
it.  The chains' atoms are one contiguous run each, which is what makes the
range the kernel takes enough to say which.*/
static void sasa_compute(pTHX_ structset *CSP_RESTRICT s, UV npts, bool want_iface)
{
	NV *CSP_RESTRICT pts = NULL;
	UV c;

	Newxz(s->area, s->n_atom ? s->n_atom : 1, NV);
	if (s->n_atom == 0) return;
	Newx(pts, npts * 3, NV);
	sasa_sphere(pts, npts);
	sasa_kernel(aTHX_ s->x, s->y, s->z, s->rad, s->n_atom, pts, npts, s->area);

	if (want_iface && s->n_chain > 1) {
		Newxz(s->alone, s->n_atom, NV);
		for (c = 0; c < s->n_chain; c++) {
			UV r0 = s->chain_first[c], r1 = s->chain_last[c], first, last;
			if (r0 >= r1) continue;
			first = s->res_first[r0];
			last  = s->res_last[r1 - 1];
			if (last <= first) continue;
			sasa_kernel(aTHX_ s->x + first, s->y + first, s->z + first,
			            s->rad + first, last - first, pts, npts, s->alone + first);
		}
	} else if (want_iface) {
		/*One chain buries nothing against anything: its isolated surface is the
		surface it already has, and saying so is cheaper than computing it
		again to arrive at zero.*/
		Newx(s->alone, s->n_atom, NV);
		Copy(s->area, s->alone, s->n_atom, NV);
	}
	Safefree(pts);
}

/*Eigenvalues of a symmetric 3x3 matrix, ascending.

The gyration tensor's three eigenvalues are the principal moments, and the shape
descriptors below are all built out of them.  numpy reaches for LAPACK's
eigvalsh here; this is the cyclic Jacobi method, which for a 3x3 is a handful of
plane rotations and needs nothing linked in.

Jacobi rather than the closed form.  A symmetric 3x3's eigenvalues can be
written with an arccosine, and that expression loses most of its digits when two
of the three are close -- which for a globular structure, whose tensor is nearly
isotropic, is the ordinary case rather than the awkward one.  Jacobi is
backward-stable and converges quadratically; the sweep count below is a bound
that is never reached in practice, not a tuning parameter.*/
static void sym3_eigenvalues(const NV m[3][3], NV *CSP_RESTRICT out)
{
	NV a[3][3];
	unsigned short int sweep, p, q, i;
	for (p = 0; p < 3; p++)
		for (q = 0; q < 3; q++)
			a[p][q] = m[p][q];
	/*Twelve sweeps: each annihilates all three off-diagonal entries, and the
	classical bound for a 3x3 is met in five or six.  The loop leaves early on
	any matrix that is already diagonal to working precision.*/
	for (sweep = 0; sweep < 12; sweep++) {
		NV off = nv_fabs(a[0][1]) + nv_fabs(a[0][2]) + nv_fabs(a[1][2]);
		NV scale = nv_fabs(a[0][0]) + nv_fabs(a[1][1]) + nv_fabs(a[2][2]);
		if (!(off > scale * NV_EPSILON)) break;
		for (p = 0; p < 2; p++) {
			for (q = (unsigned short int)(p + 1); q < 3; q++) {
				NV apq = a[p][q], theta, t, c, sn;
				if (!(nv_fabs(apq) > 0.0)) continue;
				theta = (a[q][q] - a[p][p]) / (2.0 * apq);
				/*the smaller root, which is the rotation that does not swap the
				two diagonal entries and so keeps the iteration converging*/
				t = (theta >= 0.0 ? 1.0 : -1.0)
				  / (nv_fabs(theta) + nv_sqrt(1.0 + theta * theta));
				c = 1.0 / nv_sqrt(1.0 + t * t);
				sn = t * c;
				{
					NV app = a[p][p], aqq = a[q][q];
					unsigned short int r = (unsigned short int)(3 - p - q);
					NV apr = a[p][r], aqr = a[q][r];
					a[p][p] = app - t * apq;
					a[q][q] = aqq + t * apq;
					a[p][q] = a[q][p] = 0.0;
					a[p][r] = a[r][p] = c * apr - sn * aqr;
					a[q][r] = a[r][q] = sn * apr + c * aqr;
				}
			}
		}
	}
	out[0] = a[0][0]; out[1] = a[1][1]; out[2] = a[2][2];
	//three values: an insertion sort is the whole of it
	for (i = 1; i < 3; i++) {
		NV v = out[i];
		short int j = (short int)(i - 1);
		while (j >= 0 && out[j] > v) { out[j + 1] = out[j]; j--; }
		out[j + 1] = v;
	}
}

/*Aromatic rings, found and given a plane.

The centroid is the mean of the ring atoms' positions and the normal is the
cross product of the first two atoms' offsets from it, normalised -- which is
mdtraj's compute_centroid() and compute_ring_normal() exactly, including the
part where only two of the ring's atoms decide the plane.  A ring in a real
structure is not quite planar, so a least-squares plane through all six atoms
would be a different vector; it would also be a different answer from the one
the reference implementation gives, and the reference implementation is what
this is tested against.

The direction the normal points in is arbitrary -- it flips if the ring is
listed the other way round -- and nothing below depends on it, because every
angle is folded into 0..90 degrees before it is compared with anything.*/
typedef struct {
	NV cx, cy, cz; //centroid
	NV nx, ny, nz; //unit normal
	UV res;        //which residue it belongs to
	char label;    //'6' or '5': the ring size, which is what names it
} ring_t;

//the coordinates of a named atom of a residue, false when the residue has no
//atom of that name (a side chain modelled only as far as CB, most often).
//Named for what it does rather than for the first thing that wanted it: the
//rings, the disulfides and the dihedrals all reach for atoms by name.
static bool atom_xyz(pTHX_ HV *CSP_RESTRICT atoms, const char *CSP_RESTRICT name,
                      NV *CSP_RESTRICT px, NV *CSP_RESTRICT py, NV *CSP_RESTRICT pz)
{
	STRLEN nlen = strlen(name);
	SV **slot = hv_fetch(atoms, name, (I32)nlen, 0);
	HV *a;
	SV *xs, *ys, *zs;
	if (!slot || !*slot || !SvROK(*slot) || SvTYPE(SvRV(*slot)) != SVt_PVHV) return FALSE;
	a = (HV *)SvRV(*slot);
	xs = hvf_sv(aTHX_ a, "x", 1);
	ys = hvf_sv(aTHX_ a, "y", 1);
	zs = hvf_sv(aTHX_ a, "z", 1);
	if (!xs || !ys || !zs) return FALSE;
	*px = SvNV(xs); *py = SvNV(ys); *pz = SvNV(zs);
	return TRUE;
}

static UV rings_find(pTHX_ structset *CSP_RESTRICT s, ring_t *CSP_RESTRICT *CSP_RESTRICT out)
{
	ring_t *rings = NULL;
	UV cap = 0, n = 0, r;
	*out = NULL;
	for (r = 0; r < s->n_res; r++) {
		const ring_def *defs;
		unsigned short int nd = ring_defs(s->res_key[r], &defs), d;
		HV *atoms;
		if (nd == 0) continue;
		atoms = hvf_hv(aTHX_ s->res_hv[r], "atoms", 5);
		if (!atoms) continue;
		for (d = 0; d < nd; d++) {
			//the largest ring in the table is six-membered; the array is the
			//size of the table's largest entry, not of anything the file says
			NV px[RING_MAX], py[RING_MAX], pz[RING_MAX];
			NV cx = 0.0, cy = 0.0, cz = 0.0, ux, uy, uz, vx, vy, vz, wx, wy, wz, wl;
			unsigned short int i, na = defs[d].n;
			bool whole = TRUE;
			for (i = 0; i < na; i++) {
				if (!atom_xyz(aTHX_ atoms, defs[d].atom[i], &px[i], &py[i], &pz[i])) {
					whole = FALSE;
					break;
				}
				cx += px[i]; cy += py[i]; cz += pz[i];
			}
			if (!whole) continue; //an incomplete ring has no plane to speak of
			cx /= (NV)na; cy /= (NV)na; cz /= (NV)na;
			ux = px[0] - cx; uy = py[0] - cy; uz = pz[0] - cz;
			vx = px[1] - cx; vy = py[1] - cy; vz = pz[1] - cz;
			wx = uy * vz - uz * vy;
			wy = uz * vx - ux * vz;
			wz = ux * vy - uy * vx;
			wl = nv_sqrt(wx * wx + wy * wy + wz * wz);
			//three collinear points define no plane; a structure with two ring
			//atoms deposited at the same position is the way that happens
			if (!(wl > 0.0)) continue;
			if (n == cap) {
				cap = cap ? cap * 2 : 32;
				Renew(rings, cap, ring_t);
			}
			rings[n].cx = cx; rings[n].cy = cy; rings[n].cz = cz;
			rings[n].nx = wx / wl; rings[n].ny = wy / wl; rings[n].nz = wz / wl;
			rings[n].res = r;
			rings[n].label = defs[d].label;
			n++;
		}
	}
	*out = rings;
	return n;
}

//the angle between two vectors, folded into 0..pi/2 -- mdtraj's compute_angles()
//followed by its cap_angle(), which is how a plane's two-sided normal is made to
//mean one thing.  The cosine is clamped before acos() sees it: a dot product of
//two unit vectors can land a fraction of an ulp outside [-1, 1] and acos() of
//that is a NaN, which numpy warns about and carries and this would not.
static NV vec_angle_capped(NV ax, NV ay, NV az, NV bx, NV by, NV bz)
{
	NV la = nv_sqrt(ax * ax + ay * ay + az * az);
	NV lb = nv_sqrt(bx * bx + by * by + bz * bz);
	NV c, a;
	if (!(la > 0.0) || !(lb > 0.0)) return -1.0; //no angle: the caller drops the pair
	c = (ax * bx + ay * by + az * bz) / (la * lb);
	if (c > 1.0) c = 1.0; else if (c < -1.0) c = -1.0;
	a = nv_acos(c);
	return (a > CSP_PI / 2.0) ? CSP_PI - a : a;
}

static NV det3(NV a, NV b, NV c, NV d, NV e, NV f, NV g, NV h, NV i)
{
	return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
}

/*The distance from the nearer of two rings' centroids to the line where their
planes meet -- mdtraj's calculate_intersection_point() followed by the two
norms it is fed to.

An edge-to-face stack is one ring pointing its edge at the other's face, and
what tells it from two rings merely at right angles some way apart is that the
line their planes share passes close to both of them.  Solving [n1; n2; n1xn2]
p = [n1.c1, n2.c2, 0] puts p on that line, projecting c1 onto the line moves it
to the nearest point, and the smaller of the two centroid distances from there
is what gets compared with the radius.

False when the planes are parallel and there is no line.  mdtraj tests the
determinant against zero and lets a NaN fall through its threshold; this is
only ever reached with the planes 50 degrees or more apart, where the
determinant is |n1 x n2|^2 >= sin^2(50) = 0.587, so the test is a formality
either way.*/
static bool ring_intersect(const ring_t *CSP_RESTRICT a, const ring_t *CSP_RESTRICT b,
                           NV *CSP_RESTRICT dmin)
{
	NV dx = a->ny * b->nz - a->nz * b->ny;
	NV dy = a->nz * b->nx - a->nx * b->nz;
	NV dz = a->nx * b->ny - a->ny * b->nx;
	NV det = det3(a->nx, a->ny, a->nz, b->nx, b->ny, b->nz, dx, dy, dz);
	NV o1, o2, px, py, pz, dl, vx, vy, vz, proj, ix, iy, iz, da, db;
	if (nv_fabs(det) <= NV_EPSILON) return FALSE;
	o1 = a->nx * a->cx + a->ny * a->cy + a->nz * a->cz;
	o2 = b->nx * b->cx + b->ny * b->cy + b->nz * b->cz;
	//Cramer's rule on [n_a; n_b; n_a x n_b] p = [o1, o2, 0]
	px = det3(o1, a->ny, a->nz, o2, b->ny, b->nz, 0.0, dy, dz) / det;
	py = det3(a->nx, o1, a->nz, b->nx, o2, b->nz, dx, 0.0, dz) / det;
	pz = det3(a->nx, a->ny, o1, b->nx, b->ny, o2, dx, dy, 0.0) / det;
	dl = nv_sqrt(dx * dx + dy * dy + dz * dz);
	if (!(dl > 0.0)) return FALSE;
	dx /= dl; dy /= dl; dz /= dl;
	vx = a->cx - px; vy = a->cy - py; vz = a->cz - pz;
	proj = dx * vx + dy * vy + dz * vz;
	ix = px + dx * proj; iy = py + dy * proj; iz = pz + dz * proj;
	da = nv_sqrt((a->cx - ix) * (a->cx - ix) + (a->cy - iy) * (a->cy - iy)
	           + (a->cz - iz) * (a->cz - iz));
	db = nv_sqrt((b->cx - ix) * (b->cx - ix) + (b->cy - iy) * (b->cy - iy)
	           + (b->cz - iz) * (b->cz - iz));
	*dmin = (da < db) ? da : db;
	return TRUE;
}

/*Where a pair of rings has to be for the pair to count as stacked.

Angles are radians here and degrees in the option hash, distances are angstrom
here and nanometres in mdtraj.  The defaults Perl passes down are mdtraj's, in
this file's units, with the one exception argued at the head of
structure_pi_stacking() in the Perl: mdtraj's face-to-face centroid distance is
5.5 in a function whose every other distance is nanometres, which is 55 A -- far
enough that any two aromatic rings in a small protein would qualify.  ProLIF,
which mdtraj's geometry is taken from, has 5.5 A, and mdtraj's three other
distances are ProLIF's converted to nanometres.  So this takes 5.5 A as what
was meant, and the option is there for a caller who wants the number mdtraj
actually ships.*/
typedef struct {
	NV face_dist;                    //centroid separation, angstrom
	NV face_plane_lo, face_plane_hi; //angle between the ring planes, radians
	NV face_norm_lo, face_norm_hi;   //normal to centroid-to-centroid, radians
	NV edge_dist;
	NV edge_plane_lo, edge_plane_hi;
	NV edge_norm_lo, edge_norm_hi;
	NV edge_radius;                  //centroid to the planes' shared line, angstrom
} pi_opt;

//one field of the residue a ring sits in, copied into the interaction hash
static void pi_field(pTHX_ HV *CSP_RESTRICT out, const char *CSP_RESTRICT key, STRLEN klen,
                     HV *CSP_RESTRICT res, const char *CSP_RESTRICT from, STRLEN flen)
{
	SV *v = hvf_sv(aTHX_ res, from, flen);
	(void)hv_store(out, key, (I32)klen, v ? newSVsv(v) : newSVpvn("", 0), 0);
}

/*pi_stacking() -- every stacked pair of rings in the structure.

mdtraj asks the question of one named ligand group against one named receptor
group; this asks it of every pair of rings in the structure, which is the same
geometry run over a different set of pairs.  Two of them are dropped: a ring
paired with itself, which mdtraj drops too because the vector between the two
centroids is the zero vector and every angle taken from it is a NaN, and the
two rings of one tryptophan or one purine, which mdtraj would happily call a
face-to-face stack because they are coplanar and 1.2 A apart.  Two rings fused
along a bond are one aromatic system, not two systems stacked.

The pairs are found through the same grid the SASA uses, on the centroids: a
structure with a few thousand aromatic rings has a few million pairs, and all
but a handful of them are nowhere near each other.*/
static AV *pi_stacking(pTHX_ structset *CSP_RESTRICT s, const pi_opt *CSP_RESTRICT po)
{
	AV *out = newAV();
	ring_t *rings = NULL;
	NV *cx = NULL, *cy = NULL, *cz = NULL;
	NV cut = (po->face_dist > po->edge_dist) ? po->face_dist : po->edge_dist;
	const NV to_deg = 180.0 / CSP_PI;
	cell_grid g;
	UV n, i;

	n = rings_find(aTHX_ s, &rings);
	if (n == 0) { Safefree(rings); return out; }
	Newx(cx, n, NV); Newx(cy, n, NV); Newx(cz, n, NV);
	for (i = 0; i < n; i++) { cx[i] = rings[i].cx; cy[i] = rings[i].cy; cz[i] = rings[i].cz; }
	grid_build(aTHX_ &g, cx, cy, cz, n, cut);

	for (i = 0; i < n; i++) {
		UV bx, by, bz;
		UV ci = grid_axis(cx[i] - g.x0, g.cell, g.nx);
		UV cj = grid_axis(cy[i] - g.y0, g.cell, g.ny);
		UV ck = grid_axis(cz[i] - g.z0, g.cell, g.nz);
		UV ax0 = ci ? ci - 1 : 0, ax1 = (ci + 1 < g.nx) ? ci + 1 : g.nx - 1;
		UV ay0 = cj ? cj - 1 : 0, ay1 = (cj + 1 < g.ny) ? cj + 1 : g.ny - 1;
		UV az0 = ck ? ck - 1 : 0, az1 = (ck + 1 < g.nz) ? ck + 1 : g.nz - 1;
		for (bx = ax0; bx <= ax1; bx++)
		for (by = ay0; by <= ay1; by++)
		for (bz = az0; bz <= az1; bz++) {
			UV cell = (bx * g.ny + by) * g.nz + bz, p;
			for (p = g.start[cell]; p < g.start[cell + 1]; p++) {
				UV j = g.idx[p];
				NV vx, vy, vz, d, plane, ang_i, ang_j, inter = 0.0;
				bool face, edge = FALSE;
				const char *kind;
				HV *h;
				if (j <= i) continue; //each pair once
				if (rings[i].res == rings[j].res) continue; //fused, not stacked
				vx = cx[j] - cx[i]; vy = cy[j] - cy[i]; vz = cz[j] - cz[i];
				d = nv_sqrt(vx * vx + vy * vy + vz * vz);
				if (d > cut) continue;
				plane = vec_angle_capped(rings[i].nx, rings[i].ny, rings[i].nz,
				                         rings[j].nx, rings[j].ny, rings[j].nz);
				/*One capped angle serves for both of mdtraj's: it takes the
				angle of each ring's normal against the vector from that ring to
				the other, and the two vectors differ only in sign, which
				cap_angle() folds away.*/
				ang_i = vec_angle_capped(rings[i].nx, rings[i].ny, rings[i].nz, vx, vy, vz);
				ang_j = vec_angle_capped(rings[j].nx, rings[j].ny, rings[j].nz, vx, vy, vz);
				if (plane < 0.0 || ang_i < 0.0 || ang_j < 0.0) continue;
				face = (d <= po->face_dist
				     && plane >= po->face_plane_lo && plane <= po->face_plane_hi
				     && ((ang_i >= po->face_norm_lo && ang_i <= po->face_norm_hi)
				      || (ang_j >= po->face_norm_lo && ang_j <= po->face_norm_hi)));
				if (!face) {
					edge = (d <= po->edge_dist
					     && plane >= po->edge_plane_lo && plane <= po->edge_plane_hi
					     && ((ang_i >= po->edge_norm_lo && ang_i <= po->edge_norm_hi)
					      || (ang_j >= po->edge_norm_lo && ang_j <= po->edge_norm_hi))
					     && ring_intersect(&rings[i], &rings[j], &inter)
					     && inter <= po->edge_radius);
				}
				if (!face && !edge) continue;
				kind = face ? "face" : "edge";
				h = newHV();
				(void)hv_stores(h, "type", newSVpv(kind, 0));
				pi_field(aTHX_ h, "chain1", 6, s->res_hv[rings[i].res], "chain", 5);
				pi_field(aTHX_ h, "residue1", 8, s->res_hv[rings[i].res], "key", 3);
				pi_field(aTHX_ h, "resname1", 8, s->res_hv[rings[i].res], "resname", 7);
				(void)hv_stores(h, "ring1", newSVpvn(&rings[i].label, 1));
				pi_field(aTHX_ h, "chain2", 6, s->res_hv[rings[j].res], "chain", 5);
				pi_field(aTHX_ h, "residue2", 8, s->res_hv[rings[j].res], "key", 3);
				pi_field(aTHX_ h, "resname2", 8, s->res_hv[rings[j].res], "resname", 7);
				(void)hv_stores(h, "ring2", newSVpvn(&rings[j].label, 1));
				(void)hv_stores(h, "distance", newSVnv(d));
				(void)hv_stores(h, "plane_angle", newSVnv(plane * to_deg));
				(void)hv_stores(h, "normal_angle1", newSVnv(ang_i * to_deg));
				(void)hv_stores(h, "normal_angle2", newSVnv(ang_j * to_deg));
				if (edge) (void)hv_stores(h, "intersect_distance", newSVnv(inter));
				av_push(out, newRV_noinc((SV *)h));
			}
		}
	}
	grid_free(aTHX_ &g);
	Safefree(cx); Safefree(cy); Safefree(cz);
	Safefree(rings);
	return out;
}

/*Watson-Crick base pairs, and the G-U wobble that stands in for one.

A base pair is not something the file declares; like the disulfides below it is
something the coordinates show.  What is looked for is the canonical pairing and
only that -- the three geometries a double helix is built from -- because those
are the ones whose hydrogen bonds are fixed by the pair itself rather than by
whatever else the two bases happen to be doing.  The atom pairs are the ones
Watson, J D; Crick, F H C (1953) Nature 171(4356):737-8 drew, and the numbers
beside them are Saenger, W (1984) Principles of Nucleic Acid Structure, chapter
6, where the twenty-eight ways two bases can pair are tabulated and these are
numbers 19, 20 and 28.  That is the same numbering the archive's own annotation
uses, in _ndb_struct_na_base_pair.hbond_type_28.

**Where the two thresholds come from.**  No reader on this machine finds base
pairs -- gemmi, mdtraj and Biopython all stop short of it -- so the rule was
measured against the annotation the wwPDB deposits with the entry, which is
3DNA's.  The forty entries below are the ones of those read that carry an
_ndb_struct_na_base_pair loop, and between them they hold 1372 pairs of type 19,
20 or 28 in model 1 of the asymmetric unit:

  1B94 1B97 1BNA 1D23 1D61 1D8X 1D9R 1EHZ 1FIR 1GID 1I9V 1JJ2 1KD3 1MSY
  1NJW 1NJX 1NK0 1NK4 1NK7 1NK8 1NK9 1NKC 1NKE 1SDR 1U8D 1Y26 1ZBI 2GIS
  2O1I 2R8S 355D 3AU6 3DNB 3SWP 3V9D 456D 4GLX 4OFA 5YTY 5YTZ

Two bases are a pair when both of the following hold.

  - Every hydrogen bond of that pair type is at most CSP_BP_HBOND long.  The
    longest one in the 1354 annotated pairs built from unmodified bases is
    3.4941 A and the shortest one in a candidate the annotation does not call a
    pair is 3.5144 A, so 3.5 -- which is also the conventional heavy-atom
    hydrogen bond distance -- falls between the two, with 0.006 A over the
    population it has to cover.
  - The two bases are coplanar to within CSP_BP_STAGGER, measured as the offset
    between their ring centroids projected on their ring normals, whichever of
    the two projections is larger.  That is what tells a pair from the stacked
    neighbour above or below it, which brings the same atoms within reach but
    sits a rise away rather than beside it: the largest stagger in a real pair
    is 2.5408 A, the smallest in a rejected candidate 2.6983 A, and the stacked
    population proper begins near 3.0 A.

Run that way the rule finds all 1354 and misses none.  The other eighteen of the
1372 are not pairs this can see at all, for two reasons that are both about the
file rather than the geometry.  Eleven have a modified base on one side -- 5MC,
2MG, BRU, DDG -- which res_lookup() does not spell, so they have no
single-letter code to match the table with.  The other seven are in 1D61 and
2O1I, whose asymmetric unit holds one strand of a self-complementary duplex and
whose second strand is a crystallographic symmetry mate: the annotation names
both strands with the same chain and the same residue numbers, and the two are
thirty and forty angstrom apart in the coordinates as deposited, which is what
this reads.

Two pairs are found that the annotation does not list, and both are worth
reading rather than a fault: 1JJ2's C2542-G2617, a G-C with all three bonds
under 2.94 A and half an angstrom of stagger that appears nowhere in that
entry's 1121 annotated rows, and 3SWP's DT4-DA24, where the entry is 4.11 A, the
annotation pairs its DT4 with DA25 instead, and cannot put a Saenger number on
that pair either.

Every base came out of the set with at most one partner, over all 1354 pairs.
That is what the geometry did, not something imposed, so a residue's answer is a
list as its disulfides are: nothing in the rule forbids a second.*/
#define CSP_BP_HBOND   3.5 //angstrom: the longest hydrogen bond a pair may have
#define CSP_BP_STAGGER 2.6 //angstrom: the furthest out of plane it may be
/*How far from its own ring centroid any atom named in the table below can sit,
which is what makes the centroid grid a safe screen: if every bond is at most
CSP_BP_HBOND then the centroids are at most CSP_BP_HBOND + 2 * CSP_BASE_REACH
apart.  The furthest such atom over the 11,525 in the forty entries above is an
adenine N6 at 2.7734 A.*/
#define CSP_BASE_REACH 3.0
#define CSP_BP_BONDS 3 //the most hydrogen bonds any pair in the table has

typedef struct {
	char a, b;                  //the two bases, the purine first
	unsigned short int saenger; //which of Saenger's twenty-eight this is
	unsigned short int n;       //how many of the bonds below define it
	const char *atom[CSP_BP_BONDS][2]; //the a atom and the b atom of each bond
} bp_def;

static const bp_def bp_kinds[] = {
	{ 'G', 'C', 19, 3, { { "O6", "N4" }, { "N1", "N3" }, { "N2", "O2" } } },
	{ 'A', 'U', 20, 2, { { "N6", "O4" }, { "N1", "N3" }, { NULL, NULL } } },
	{ 'A', 'T', 20, 2, { { "N6", "O4" }, { "N1", "N3" }, { NULL, NULL } } },
	//the wobble, and the G-T that is the same pair with a methyl on it: the
	//archive numbers that one 28 as well, in 1NJW, 1NJX, 1NK8 and 1NKC
	{ 'G', 'U', 28, 2, { { "O6", "N3" }, { "N1", "O2" }, { NULL, NULL } } },
	{ 'G', 'T', 28, 2, { { "O6", "N3" }, { "N1", "O2" }, { NULL, NULL } } }
};

typedef struct {
	NV cx, cy, cz; //the six-membered ring's centroid
	NV nx, ny, nz; //and its normal
	HV *atoms;     //the residue's atoms, for reaching the pairing atoms by name
	UV res;        //which residue it belongs to
	char one;      //A, C, G, T or U
} base_t;

/*The bases, out of the rings pi_stacking() already knows how to find.  A purine
has two rings and only the six-membered one is here, because it is the one both
purines and pyrimidines have and the one every atom in the table above hangs
off.  A modified base has no single-letter code and so is not a base here; that
is the eighteen pairs the header comment accounts for.*/
static UV bases_find(pTHX_ structset *CSP_RESTRICT s, base_t *CSP_RESTRICT *CSP_RESTRICT out)
{
	ring_t *rings = NULL;
	base_t *b = NULL;
	UV n, k, m = 0;
	*out = NULL;
	n = rings_find(aTHX_ s, &rings);
	if (n == 0) { Safefree(rings); return 0; }
	Newx(b, n, base_t);
	for (k = 0; k < n; k++) {
		const UV r = rings[k].res;
		const char one = s->res_one[r];
		if (rings[k].label != '6') continue;
		if (s->res_type[r] != RT_NUC) continue;
		if (one != 'A' && one != 'C' && one != 'G' && one != 'T' && one != 'U') continue;
		b[m].cx = rings[k].cx; b[m].cy = rings[k].cy; b[m].cz = rings[k].cz;
		b[m].nx = rings[k].nx; b[m].ny = rings[k].ny; b[m].nz = rings[k].nz;
		b[m].atoms = hvf_hv(aTHX_ s->res_hv[r], "atoms", 5);
		b[m].res = r;
		b[m].one = one;
		if (b[m].atoms) m++;
	}
	Safefree(rings);
	*out = b;
	return m;
}

//which entry of bp_kinds two bases could be, and which way round; NULL for a
//combination that is not one of the canonical pairs
static const bp_def *bp_kind_of(char x, char y, bool *CSP_RESTRICT swapped)
{
	unsigned short int i;
	for (i = 0; i < (unsigned short int)C_ARRAY_LENGTH(bp_kinds); i++) {
		if (bp_kinds[i].a == x && bp_kinds[i].b == y) { *swapped = FALSE; return &bp_kinds[i]; }
		if (bp_kinds[i].a == y && bp_kinds[i].b == x) { *swapped = TRUE;  return &bp_kinds[i]; }
	}
	return NULL;
}

//the partner half of one pair, pushed onto a base's own list
static void bp_note(pTHX_ HV *CSP_RESTRICT res, HV *CSP_RESTRICT partner,
                    const char *CSP_RESTRICT type, unsigned short int saenger)
{
	SV **slot = hv_fetch(res, "base_pair", 9, 0);
	AV *list;
	HV *e;
	if (slot && *slot && SvROK(*slot) && SvTYPE(SvRV(*slot)) == SVt_PVAV) {
		list = (AV *)SvRV(*slot);
	} else {
		list = newAV();
		(void)hv_stores(res, "base_pair", newRV_noinc((SV *)list));
	}
	e = newHV();
	pi_field(aTHX_ e, "chain", 5, partner, "chain", 5);
	pi_field(aTHX_ e, "residue", 7, partner, "key", 3);
	pi_field(aTHX_ e, "resname", 7, partner, "resname", 7);
	(void)hv_stores(e, "type", newSVpv(type, 0));
	(void)hv_stores(e, "saenger", newSVuv((UV)saenger));
	av_push(list, newRV_noinc((SV *)e));
}

static AV *base_pairs(pTHX_ structset *CSP_RESTRICT s, NV hb, NV stagger, bool store)
{
	AV *out = newAV();
	base_t *bases = NULL;
	NV *cx = NULL, *cy = NULL, *cz = NULL;
	const NV cut = hb + 2.0 * CSP_BASE_REACH;
	const NV to_deg = 180.0 / CSP_PI;
	cell_grid g;
	UV n, i;

	n = bases_find(aTHX_ s, &bases);
	if (n < 2) { Safefree(bases); return out; }
	Newx(cx, n, NV); Newx(cy, n, NV); Newx(cz, n, NV);
	for (i = 0; i < n; i++) { cx[i] = bases[i].cx; cy[i] = bases[i].cy; cz[i] = bases[i].cz; }
	/*Clear what an earlier call left on these residues, so that asking twice
	replaces the answer rather than adding to it.  Which residues are candidates
	is decided by the residue's name and its ring atoms alone, never by either
	threshold, so these are exactly the ones a previous call could have written
	to.*/
	if (store)
		for (i = 0; i < n; i++)
			(void)hv_delete(s->res_hv[bases[i].res], "base_pair", 9, G_DISCARD);
	grid_build(aTHX_ &g, cx, cy, cz, n, cut);

	for (i = 0; i < n; i++) {
		UV bx, by, bz;
		UV ci = grid_axis(cx[i] - g.x0, g.cell, g.nx);
		UV cj = grid_axis(cy[i] - g.y0, g.cell, g.ny);
		UV ck = grid_axis(cz[i] - g.z0, g.cell, g.nz);
		UV ax0 = ci ? ci - 1 : 0, ax1 = (ci + 1 < g.nx) ? ci + 1 : g.nx - 1;
		UV ay0 = cj ? cj - 1 : 0, ay1 = (cj + 1 < g.ny) ? cj + 1 : g.ny - 1;
		UV az0 = ck ? ck - 1 : 0, az1 = (ck + 1 < g.nz) ? ck + 1 : g.nz - 1;
		for (bx = ax0; bx <= ax1; bx++)
		for (by = ay0; by <= ay1; by++)
		for (bz = az0; bz <= az1; bz++) {
			UV cell = (bx * g.ny + by) * g.nz + bz, p;
			for (p = g.start[cell]; p < g.start[cell + 1]; p++) {
				UV j = g.idx[p];
				const bp_def *kind;
				bool swapped = FALSE;
				NV vx, vy, vz, d, plane, sag_i, sag_j, sag;
				NV bond[CSP_BP_BONDS];
				unsigned short int k;
				bool whole = TRUE;
				char tname[4];
				HV *h;
				AV *hbs;
				if (j <= i) continue; //each pair once
				//one residue's two rings cannot pair with each other, and
				//neither can two conformers sharing a residue
				if (bases[i].res == bases[j].res) continue;
				kind = bp_kind_of(bases[i].one, bases[j].one, &swapped);
				if (!kind) continue;
				vx = cx[j] - cx[i]; vy = cy[j] - cy[i]; vz = cz[j] - cz[i];
				d = nv_sqrt(vx * vx + vy * vy + vz * vz);
				if (d > cut) continue;
				for (k = 0; k < kind->n; k++) {
					/*the table names the purine's atom first, so a pair whose
					purine is the later residue reads its two names backwards*/
					const char *an = kind->atom[k][swapped ? 1 : 0];
					const char *bn = kind->atom[k][swapped ? 0 : 1];
					NV ax, ay, az, bx2, by2, bz2, dx, dy, dz;
					if (!atom_xyz(aTHX_ bases[i].atoms, an, &ax, &ay, &az)
					 || !atom_xyz(aTHX_ bases[j].atoms, bn, &bx2, &by2, &bz2)) {
						whole = FALSE;
						break;
					}
					dx = bx2 - ax; dy = by2 - ay; dz = bz2 - az;
					bond[k] = nv_sqrt(dx * dx + dy * dy + dz * dz);
					if (bond[k] > hb) { whole = FALSE; break; }
				}
				if (!whole) continue;
				sag_i = nv_fabs(vx * bases[i].nx + vy * bases[i].ny + vz * bases[i].nz);
				sag_j = nv_fabs(vx * bases[j].nx + vy * bases[j].ny + vz * bases[j].nz);
				sag = (sag_i > sag_j) ? sag_i : sag_j;
				if (sag > stagger) continue;
				plane = vec_angle_capped(bases[i].nx, bases[i].ny, bases[i].nz,
				                         bases[j].nx, bases[j].ny, bases[j].nz);
				if (plane < 0.0) continue;
				//the pair is named in the order it is reported, so that the
				//two letters and the two resname fields read the same way round
				tname[0] = bases[i].one; tname[1] = '-';
				tname[2] = bases[j].one; tname[3] = '\0';
				h = newHV();
				(void)hv_stores(h, "type", newSVpvn(tname, 3));
				(void)hv_stores(h, "saenger", newSVuv((UV)kind->saenger));
				pi_field(aTHX_ h, "chain1", 6, s->res_hv[bases[i].res], "chain", 5);
				pi_field(aTHX_ h, "residue1", 8, s->res_hv[bases[i].res], "key", 3);
				pi_field(aTHX_ h, "resname1", 8, s->res_hv[bases[i].res], "resname", 7);
				pi_field(aTHX_ h, "chain2", 6, s->res_hv[bases[j].res], "chain", 5);
				pi_field(aTHX_ h, "residue2", 8, s->res_hv[bases[j].res], "key", 3);
				pi_field(aTHX_ h, "resname2", 8, s->res_hv[bases[j].res], "resname", 7);
				(void)hv_stores(h, "distance", newSVnv(d));
				(void)hv_stores(h, "plane_angle", newSVnv(plane * to_deg));
				(void)hv_stores(h, "stagger", newSVnv(sag));
				hbs = newAV();
				for (k = 0; k < kind->n; k++) {
					HV *e = newHV();
					(void)hv_stores(e, "atom1", newSVpv(kind->atom[k][swapped ? 1 : 0], 0));
					(void)hv_stores(e, "atom2", newSVpv(kind->atom[k][swapped ? 0 : 1], 0));
					(void)hv_stores(e, "distance", newSVnv(bond[k]));
					av_push(hbs, newRV_noinc((SV *)e));
				}
				(void)hv_stores(h, "hbonds", newRV_noinc((SV *)hbs));
				av_push(out, newRV_noinc((SV *)h));
				if (store) {
					bp_note(aTHX_ s->res_hv[bases[i].res], s->res_hv[bases[j].res],
					        tname, kind->saenger);
					tname[0] = bases[j].one; tname[2] = bases[i].one;
					bp_note(aTHX_ s->res_hv[bases[j].res], s->res_hv[bases[i].res],
					        tname, kind->saenger);
				}
			}
		}
	}
	grid_free(aTHX_ &g);
	Safefree(cx); Safefree(cy); Safefree(cz);
	Safefree(bases);
	return out;
}

/*Base stacking, scored the way Condon's tetramer benchmark scores it.

base_pairs() above answers whether two bases lie side by side in a plane; this
answers the other question a nucleic acid poses, which is whether one base lies
on the face of another.  pi_stacking() answers a version of it too, but as a
yes/no over a ring's centroid and normal, with mdtraj's thresholds; this is the
question asked the way the RNA literature asks it -- three geometric variables
and one number between -100% and 100% built out of them -- so that a caller can
see how stacked a pair is rather than only that it passed.

Provenance.  The definition and every constant below are

  Condon, D E; Kennedy, S D; Mort, B C; Kierzek, R; Yildirim, I; Turner, D H
  (2015) "Stacking in RNA: NMR of Four Tetramers Benchmark Molecular
  Dynamics", J Chem Theory Comput 11(6):2729-2742, doi:10.1021/ct501025q,
  section 2.4 "Criteria for Stacking", equations 9 and 10 and Figures 4 and 5,

and the implementation those numbers were produced with, which is the author's
own PDB_stacker (https://github.com/hhg7/PDB_stacker, pdb_stacking.pl).  The
paper is the definition of record and the script is what was run; where the two
disagree the script is what the published percentages mean, and each such place
is marked below with the measurement that settles it.  The measurement is
Figure 4, which is the one worked example either gives: residues 13 (C) and 14
(G) of chain B of PDB entry 157D, for which the caption reports d0 = 4.5 A,
omega = 40.7 degrees and Xi = 17.3 degrees.  Those residues are t/data/aform.pdb
and t/stacking.t is where the three numbers are checked.

**The base's frame.**  Each base gets a centre of mass over the heavy atoms
listed in nucbase_of() and two vectors a and b from that centre to two named
atoms far apart on the ring, chosen so that the pair spans the base and so that
out-of-plane distortion moves the cross product as little as possible.  a x b
and b x a are the base's two normal vectors, one above the plane and one below;
which is which is arbitrary and nothing here depends on it, because every
quantity below is a minimum over the two.

The normals are *not* normalised, and that is load-bearing rather than an
oversight: |a x b| is the d1 of equation 9, so the length of the cross product
enters omega as a lever arm.  It has units of area and is used as a length,
which is what the reference implementation does and what the published angles
were computed with.

**The three variables**, for an ordered pair of bases (the 5' one first):

  d0     the distance between the two centres of mass.
  omega  "oh-mega", for overlap: the angle at the 5' centre of mass in the
         triangle whose sides are d0, d1 = |a x b| of the 5' base, and d2, the
         distance from the tip of whichever of the 5' base's two normals lands
         nearer to the 3' base's centre of mass.  Equation 9, by the law of
         cosines.  It is small when the 3' base sits over the 5' base's face
         and large when it sits beside it, like the angle between the steps of
         a staircase.
  Xi     the angle between the two bases' normal vectors, 0 when the planes are
         parallel.  Equation 10 and Figure 5: Xi near 0 is a parallel stack and
         Xi near 90 a T-shape, which is a different interaction and is scored
         negative rather than dropped.

**Which base is the 5' one.**  omega is not symmetric -- it is measured from
one base's normal against the other's centre of mass -- so the pair has to be
ordered.  The order is the order the file lists the two residues in, which for
a strand written 5' to 3' as the formats and the archive write it is the
chemical order.  Two bases in different chains have no 5'/3' relation at all;
there the order is the chain order, and the reported omega is the one measured
from the base the pair names first.

**Where this departs from the paper's text, and why.**

  - Equation 10 as printed is min(arcsin(|(a5 x b5) x (a3 x b3)| / ...),
    arcsin(|(a5 x b5) x (b3 x a3)| / ...)).  The two arguments differ only in
    the sign of the second cross product, so the two arcsines are equal and the
    minimum is a formality; and taken literally the expression is the angle
    between the two normals with no reference to where the bases are.  The
    script instead compares the 5' base's normal with the 3' base's normal
    *translated so that both are drawn from the 5' base's centre of mass* --
    its comment says "pretend that both vectors are centered on the 5'
    nucleotide's Center of Mass" -- and takes the smaller of the two answers
    the 3' base's two normals give.  On the Figure 4 pair the printed equation
    gives 9.5 degrees and the script gives 17.31; the caption says 17.3, so the
    script is what Xi means and is what is implemented here.
  - Guanine's centre of mass is taken over ten atoms and not eleven: N2, the
    exocyclic amino nitrogen, is not in the script's list, though adenine's N6
    and cytosine's N4 are in theirs.  It reads like an omission, but it is the
    omission the published numbers were computed with -- on the Figure 4 pair
    including N2 gives d0 = 4.77 A and omega = 43.56 degrees against the
    caption's 4.5 and 40.7, and leaving it out gives 4.53 and 40.74.
  - Criterion I in the paper's text puts the distance knee at 3.5 A ("If
    d0 <= 3.5 A, the stacking score is incremented +1"); the script's
    $DISTANCE_MIN is 4.  4 is used here, again because it is the number behind
    the published percentages.  It is not an option: the two cutoffs that
    decide whether a pair is reported at all are options below, and the shape
    of the ramp between them is the definition.
  - The masses are this module's own standard atomic weights, from
    elem_prop_of(), where the script uses the mass numbers 12, 14 and 16.  One
    mass table in the file is worth more than agreement in the last place the
    difference reaches: on the Figure 4 pair the two centres of mass differ by
    enough to move d0 from 4.528 A to 4.529 A, and omega and Xi not at all to
    two decimal places.

**What it does not do.**  A base whose letter is not one of the six below, or
that is missing any atom the letter's entry names, contributes no frame and so
appears in no pair -- the same rule rings_find() applies to an incomplete ring.
That covers a 4-thiouridine, whose O4 is a sulphur, as well as a base whose
density ran out; no geometry is invented for either.*/

/*The distance beyond which two bases are not stacked, from CCSD(T) calculations
on stacked uracil and adenine dimers -- the paper's criterion I, citing its
references 63 and 45.  This and the angle below are defaults: the caller's two
cutoffs are what base_stacks() builds its grid on and scores against.*/
#define CSP_STACK_DIST 5.0
//the angle beyond which they are not stacked either, and Xi is not computed;
//criterion II, chosen there from X-ray statistics
#define CSP_STACK_OMEGA 50.0
//the knees of the two ramps: full marks at or below these, falling to nothing
//at the cutoffs above.  PDB_stacker's $DISTANCE_MIN and $OMEGA_MIN.
#define CSP_STACK_DIST_KNEE  4.0
#define CSP_STACK_OMEGA_KNEE 25.0
//Xi above this is a T-shape rather than a stack, and flips the score's sign.
//The paper's "45.0 < chi < 135.0" over an angle its own equation confines to
//0..90, so only the lower bound can ever be reached.
#define CSP_STACK_T 45.0

typedef struct {
	const char *const *atom;  //the heavy atoms the centre of mass is taken over
	unsigned short int n;     //how many
	const char *a, *b;        //the two atoms the normal's cross product is built from
} nucbase_def;

/*Adenine's ten heavy atoms, guanine's ten (see N2, above), cytosine's eight and
uracil's eight, exactly as PDB_stacker's %atoms lists them.  Every name in every
list begins with the atom's element, which is what lets the masses be looked up
from the first byte alone; that is true of these fourteen names and is not a
general fact about PDB atom names.*/
static const char *const base_ade[10] = {
	"N1", "C2", "N3", "C4", "C5", "C6", "N6", "N7", "C8", "N9"
};
static const char *const base_gua[10] = {
	"N1", "C2", "N3", "C4", "C5", "C6", "O6", "N7", "C8", "N9"
};
static const char *const base_cyt[8] = {
	"N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"
};
static const char *const base_ura[8] = {
	"N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"
};

/*The frame for a base, by the single-letter code res_lookup() gave its residue,
or false for a letter that names no base this can build one for.

Six letters out of the four entries.  Inosine has guanine's ring and guanine's
O6 and simply has no N2, so the ten atoms guanine is measured over are ten
inosine also has; thymine is uracil with a methyl on C5, and takes uracil's
eight, which leaves C7 out of the centre of mass exactly as guanine's N2 is left
out.  N -- the archive's unknown nucleotide -- names no base and gets none.

Every modified base res_lookup() spells reaches its parent's entry through its
letter: a pseudouridine is measured as a uridine, a 7-methylguanosine as a
guanosine, and so on, which is right as long as the modification leaves the
named atoms where they were.  Where it does not, the atom is missing under that
name and the base is dropped rather than mismeasured.*/
static bool nucbase_of(char one, const nucbase_def *CSP_RESTRICT *CSP_RESTRICT out)
{
	//the two vectors are the paper's: (CoM -> C8, CoM -> N6) for adenine,
	//(CoM -> C8, CoM -> O6) for guanine, (CoM -> O2, CoM -> N4) for cytosine
	//and (CoM -> O2, CoM -> O4) for uracil
	static const nucbase_def defs[4] = {
		{ base_ade, 10, "C8", "N6" },
		{ base_gua, 10, "C8", "O6" },
		{ base_cyt,  8, "O2", "N4" },
		{ base_ura,  8, "O2", "O4" }
	};
	switch (one) {
		case 'A': *out = &defs[0]; return TRUE;
		case 'G': case 'I': *out = &defs[1]; return TRUE;
		case 'C': *out = &defs[2]; return TRUE;
		case 'U': case 'T': *out = &defs[3]; return TRUE;
		default:  *out = NULL;     return FALSE;
	}
}

typedef struct {
	NV cx, cy, cz; //centre of mass over the base's heavy atoms
	NV nx, ny, nz; //a x b, the normal above the plane; NOT a unit vector
	NV nlen;       //|a x b|, which is the d1 of equation 9
	UV res;        //which residue it belongs to
	char one;      //A, C, G, I, T or U
} nucbase_t; //nucbase_t and not stack_t: POSIX <signal.h> has that name

/*Every base in the structure that has a frame, in walk order -- which is chain
order and then the order the file listed the residues, so a lower index is the
5' base of any pair drawn from one strand.*/
static UV stacks_find(pTHX_ structset *CSP_RESTRICT s, nucbase_t *CSP_RESTRICT *CSP_RESTRICT out)
{
	nucbase_t *b = NULL;
	UV n = 0, r;
	*out = NULL;
	if (s->n_res == 0) return 0;
	Newx(b, s->n_res, nucbase_t);
	for (r = 0; r < s->n_res; r++) {
		const nucbase_def *def;
		HV *atoms;
		NV cx = 0.0, cy = 0.0, cz = 0.0, tot = 0.0;
		NV ax, ay, az, bx, by, bz, wx, wy, wz, wl;
		unsigned short int i;
		bool whole = TRUE;
		if (s->res_type[r] != RT_NUC) continue;
		if (!nucbase_of(s->res_one[r], &def)) continue;
		atoms = hvf_hv(aTHX_ s->res_hv[r], "atoms", 5);
		if (!atoms) continue;
		//summed as it goes: unlike rings_find(), which needs its first two
		//atoms again to take the plane from, nothing here looks at an atom a
		//second time -- a and b are fetched by name below
		for (i = 0; i < def->n; i++) {
			elem_prop ep;
			NV px, py, pz;
			if (!atom_xyz(aTHX_ atoms, def->atom[i], &px, &py, &pz)
			 || !elem_prop_of(def->atom[i], 1, &ep) || !(ep.mass > 0.0)) {
				whole = FALSE;
				break;
			}
			cx += ep.mass * px; cy += ep.mass * py; cz += ep.mass * pz;
			tot += ep.mass;
		}
		//an incomplete base has no centre of mass to speak of, the same way an
		//incomplete ring has no plane
		if (!whole || !(tot > 0.0)) continue;
		cx /= tot; cy /= tot; cz /= tot;
		if (!atom_xyz(aTHX_ atoms, def->a, &ax, &ay, &az)) continue;
		if (!atom_xyz(aTHX_ atoms, def->b, &bx, &by, &bz)) continue;
		ax -= cx; ay -= cy; az -= cz;
		bx -= cx; by -= cy; bz -= cz;
		wx = ay * bz - az * by;
		wy = az * bx - ax * bz;
		wz = ax * by - ay * bx;
		wl = nv_sqrt(wx * wx + wy * wy + wz * wz);
		//a and b collinear with the centre of mass define no plane, and d1
		//would be zero in the denominator of equation 9
		if (!(wl > 0.0)) continue;
		b[n].cx = cx; b[n].cy = cy; b[n].cz = cz;
		b[n].nx = wx; b[n].ny = wy; b[n].nz = wz;
		b[n].nlen = wl;
		b[n].res = r;
		b[n].one = s->res_one[r];
		n++;
	}
	*out = b;
	return n;
}

/*The stacking score, as a percentage: the paper's criteria I, II and III and
PDB_stacker's stacking_score(), which award a point for the distance and a point
for the overlap, flip the sign for a T-shape, and report the two points as
100%.

  d0    <= knee            1
        knee < d0 <= cut   1 / (d0 - knee + 1)^3, the paper's "decreased as
                           r^-3 from 1 to 0"
  omega <= knee            1
        knee < om <= cut   linear from 1 at the knee to 0 at the cutoff
  Xi    > 45 degrees       the whole score is multiplied by -1

Both cutoffs are the caller's, so that a caller who wants the paper's 3.5 A
criterion or a different overlap limit has them; the knees are not.  A caller
who moves a cutoff below its knee gets a ramp that never runs, which is the
honest answer to the question asked rather than a croak.

Reached only for a pair already inside both cutoffs, which is where the score
is defined: outside them the paper does not compute the angles at all, and the
score is 0 without being calculated.*/
static NV stack_score(NV d0, NV omega, NV xi, NV dist_cut, NV omega_cut)
{
	NV score = 0.0;
	if (d0 <= CSP_STACK_DIST_KNEE) {
		score += 1.0;
	} else if (d0 <= dist_cut) {
		const NV t = d0 - CSP_STACK_DIST_KNEE + 1.0;
		score += 1.0 / (t * t * t);
	}
	if (omega <= CSP_STACK_OMEGA_KNEE) {
		score += 1.0;
	} else if (omega <= omega_cut && omega_cut > CSP_STACK_OMEGA_KNEE) {
		score += (omega_cut - omega) / (omega_cut - CSP_STACK_OMEGA_KNEE);
	}
	if (xi > CSP_STACK_T) score = -score;
	return 100.0 * score / 2.0;
}

/*The partner half of one stack, pushed onto a base's own list, as bp_note()
does for a pair.  The three variables go on with it, because a residue's own
list is the only place a caller walking residues will look for them.

Both residues get the same d0, omega, Xi and score -- they are the pair's, not
either base's.  What differs is `side', which says which end of the pair this
residue is: omega is measured from the 5' base against the 3' base's centre of
mass, so a residue that reads its own list needs to know which one it was.*/
static void stack_note(pTHX_ HV *CSP_RESTRICT res, HV *CSP_RESTRICT partner,
                       const char *CSP_RESTRICT type, NV d0, NV omega,
                       const NV *CSP_RESTRICT xi, NV score, bool five_prime)
{
	SV **slot = hv_fetch(res, "base_stack", 10, 0);
	AV *list;
	HV *e;
	if (slot && *slot && SvROK(*slot) && SvTYPE(SvRV(*slot)) == SVt_PVAV) {
		list = (AV *)SvRV(*slot);
	} else {
		list = newAV();
		(void)hv_stores(res, "base_stack", newRV_noinc((SV *)list));
	}
	e = newHV();
	pi_field(aTHX_ e, "chain", 5, partner, "chain", 5);
	pi_field(aTHX_ e, "residue", 7, partner, "key", 3);
	pi_field(aTHX_ e, "resname", 7, partner, "resname", 7);
	(void)hv_stores(e, "type", newSVpv(type, 0));
	(void)hv_stores(e, "distance", newSVnv(d0));
	(void)hv_stores(e, "omega", newSVnv(omega));
	if (xi) (void)hv_stores(e, "xi", newSVnv(*xi));
	(void)hv_stores(e, "score", newSVnv(score));
	//which end of the pair this residue is, because omega is measured from the
	//5' base and means nothing read the other way round
	(void)hv_stores(e, "side", newSVpv(five_prime ? "5'" : "3'", 0));
	av_push(list, newRV_noinc((SV *)e));
}

/*base_stacks() -- d0, omega, Xi and the score for every pair of bases close
enough for the paper to compute them.

Every pair within the distance cutoff is reported, stacked or not, because the
three variables are the answer and the score is a summary of them: a caller
benchmarking a force field wants the pair that scored 12% as much as the one
that scored 98%.  A pair whose omega is past its cutoff carries no Xi, which is
criterion II -- "the base is not considered stacked and Xi is not computed" --
and scores 0.

The pairs come off the same centroid grid base_pairs() uses, on the centres of
mass.*/
static AV *base_stacks(pTHX_ structset *CSP_RESTRICT s, NV dist_cut, NV omega_cut,
                       bool store)
{
	AV *out = newAV();
	nucbase_t *bases = NULL;
	NV *cx = NULL, *cy = NULL, *cz = NULL;
	const NV to_deg = 180.0 / CSP_PI;
	cell_grid g;
	UV n, i;

	n = stacks_find(aTHX_ s, &bases);
	if (n < 2) { Safefree(bases); return out; }
	Newx(cx, n, NV); Newx(cy, n, NV); Newx(cz, n, NV);
	for (i = 0; i < n; i++) { cx[i] = bases[i].cx; cy[i] = bases[i].cy; cz[i] = bases[i].cz; }
	/*Clear what an earlier call left behind, so that asking twice replaces the
	answer rather than adding to it.  Which residues are candidates is decided
	by the residue's name and its atoms alone, never by either cutoff, so these
	are exactly the ones a previous call could have written to.*/
	if (store)
		for (i = 0; i < n; i++)
			(void)hv_delete(s->res_hv[bases[i].res], "base_stack", 10, G_DISCARD);
	grid_build(aTHX_ &g, cx, cy, cz, n, dist_cut);

	for (i = 0; i < n; i++) {
		UV bx, by, bz;
		UV ci = grid_axis(cx[i] - g.x0, g.cell, g.nx);
		UV cj = grid_axis(cy[i] - g.y0, g.cell, g.ny);
		UV ck = grid_axis(cz[i] - g.z0, g.cell, g.nz);
		UV ax0 = ci ? ci - 1 : 0, ax1 = (ci + 1 < g.nx) ? ci + 1 : g.nx - 1;
		UV ay0 = cj ? cj - 1 : 0, ay1 = (cj + 1 < g.ny) ? cj + 1 : g.ny - 1;
		UV az0 = ck ? ck - 1 : 0, az1 = (ck + 1 < g.nz) ? ck + 1 : g.nz - 1;
		for (bx = ax0; bx <= ax1; bx++)
		for (by = ay0; by <= ay1; by++)
		for (bz = az0; bz <= az1; bz++) {
			UV cell = (bx * g.ny + by) * g.nz + bz, p;
			for (p = g.start[cell]; p < g.start[cell + 1]; p++) {
				UV j = g.idx[p];
				NV vx, vy, vz, d0, d1, d2, cosw, omega, xi = 0.0, score;
				bool have_xi = FALSE;
				char tname[4];
				HV *h;
				unsigned short int k;
				if (j <= i) continue; //each pair once, the earlier residue first
				//two conformers of one residue are not two bases
				if (bases[i].res == bases[j].res) continue;
				vx = cx[j] - cx[i]; vy = cy[j] - cy[i]; vz = cz[j] - cz[i];
				d0 = nv_sqrt(vx * vx + vy * vy + vz * vz);
				if (d0 > dist_cut || !(d0 > 0.0)) continue;
				/*d2: the 5' base's two normals are its centre of mass plus and
				minus a x b, and the one that lands nearer the 3' base's centre
				of mass is the one equation 9 wants -- the base has two faces
				and the 3' base is over one of them.*/
				d1 = bases[i].nlen;
				d2 = 0.0;
				for (k = 0; k < 2; k++) {
					const NV sign = k ? -1.0 : 1.0;
					const NV ex = vx - sign * bases[i].nx;
					const NV ey = vy - sign * bases[i].ny;
					const NV ez = vz - sign * bases[i].nz;
					const NV e = nv_sqrt(ex * ex + ey * ey + ez * ez);
					if (k == 0 || e < d2) d2 = e;
				}
				//equation 9, the law of cosines.  The cosine is clamped before
				//acos() sees it for the reason vec_angle_capped() clamps its
				//own: a value a fraction of an ulp outside [-1, 1] is a NaN
				//that would travel all the way out to the caller.
				cosw = (d0 * d0 + d1 * d1 - d2 * d2) / (2.0 * d0 * d1);
				if (cosw > 1.0) cosw = 1.0; else if (cosw < -1.0) cosw = -1.0;
				omega = nv_acos(cosw) * to_deg;
				if (omega <= omega_cut) {
					/*equation 10 as PDB_stacker computes it: the 3' base's
					normal is redrawn from the 5' base's centre of mass, so the
					vector compared is (CoM3 - CoM5) +/- n3, and the smaller of
					the two answers is Xi.  atan2(|u x v|, |u . v|) rather than
					the arcsine the paper prints: it is the same angle for
					vectors this cannot fold past 90 degrees, and it keeps its
					precision near 0, which is where a good stack sits and where
					an arcsine of a quotient near 0 has the least of it.*/
					for (k = 0; k < 2; k++) {
						const NV sign = k ? -1.0 : 1.0;
						const NV ux = bases[i].nx, uy = bases[i].ny, uz = bases[i].nz;
						const NV wx2 = vx + sign * bases[j].nx;
						const NV wy2 = vy + sign * bases[j].ny;
						const NV wz2 = vz + sign * bases[j].nz;
						const NV qx = uy * wz2 - uz * wy2;
						const NV qy = uz * wx2 - ux * wz2;
						const NV qz = ux * wy2 - uy * wx2;
						const NV cr = nv_sqrt(qx * qx + qy * qy + qz * qz);
						const NV dp = ux * wx2 + uy * wy2 + uz * wz2;
						NV t;
						if (!(nv_fabs(dp) > 0.0) && !(cr > 0.0)) continue;
						t = nv_atan2(cr, nv_fabs(dp)) * to_deg;
						if (!have_xi || t < xi) { xi = t; have_xi = TRUE; }
					}
				}
				score = have_xi ? stack_score(d0, omega, xi, dist_cut, omega_cut) : 0.0;
				//named in the order it is reported, so the two letters and the
				//two resname fields read the same way round
				tname[0] = bases[i].one; tname[1] = '-';
				tname[2] = bases[j].one; tname[3] = '\0';
				h = newHV();
				(void)hv_stores(h, "type", newSVpvn(tname, 3));
				pi_field(aTHX_ h, "chain1", 6, s->res_hv[bases[i].res], "chain", 5);
				pi_field(aTHX_ h, "residue1", 8, s->res_hv[bases[i].res], "key", 3);
				pi_field(aTHX_ h, "resname1", 8, s->res_hv[bases[i].res], "resname", 7);
				pi_field(aTHX_ h, "chain2", 6, s->res_hv[bases[j].res], "chain", 5);
				pi_field(aTHX_ h, "residue2", 8, s->res_hv[bases[j].res], "key", 3);
				pi_field(aTHX_ h, "resname2", 8, s->res_hv[bases[j].res], "resname", 7);
				(void)hv_stores(h, "distance", newSVnv(d0));
				(void)hv_stores(h, "omega", newSVnv(omega));
				if (have_xi) (void)hv_stores(h, "xi", newSVnv(xi));
				(void)hv_stores(h, "score", newSVnv(score));
				(void)hv_stores(h, "stacked", newSViv(score > 50.0 ? 1 : 0));
				av_push(out, newRV_noinc((SV *)h));
				if (store) {
					stack_note(aTHX_ s->res_hv[bases[i].res], s->res_hv[bases[j].res],
					           tname, d0, omega, have_xi ? &xi : NULL, score, TRUE);
					tname[0] = bases[j].one; tname[2] = bases[i].one;
					stack_note(aTHX_ s->res_hv[bases[j].res], s->res_hv[bases[i].res],
					           tname, d0, omega, have_xi ? &xi : NULL, score, FALSE);
				}
			}
		}
	}
	grid_free(aTHX_ &g);
	Safefree(cx); Safefree(cy); Safefree(cz);
	Safefree(bases);
	return out;
}

/*Disulfide bonds, found by geometry.

$info->{ssbond} is what the file *says* -- an SSBOND record, or an mmCIF
_struct_conn row of type disulf -- and a deposited file is not obliged to be
right or to say anything at all.  This is the other half of that: the bonds the
coordinates show, which a caller can compare with the declared ones.  Neither is
the authority; a disagreement is a fact about the entry, and 1AHW's four bonds
are declared and found alike while 1A4K declares four of the eight it has,
because the depositor listed one copy of a Fab that is in the file twice.

The rule is mdtraj's, from Topology.create_disulfide_bonds() in
mdtraj/core/topology.py: a residue named CYS that has an SG and has no HG,
paired with another when the two SG atoms are less than 0.3 nm apart.  The HG
test is what separates a cysteine whose thiol hydrogen was modelled -- so it is
reduced, and holds no bond -- from one that was not; a crystal structure with no
hydrogens in it has no HG anywhere and every cysteine is a candidate, which is
right.

**mdtraj's own implementation of that rule never fires**, and t/features.t
cannot call it.  PDBTrajectoryFile hands create_disulfide_bonds() the positions
it read out of the file, which are angstrom, and the function compares them with
0.3 -- "this is supposed to be nm. I think we're good", says the comment beside
it.  So it tests for an SG-SG separation under 0.3 A, which no pair of atoms has,
and finds nothing on any file: measured on 1A4K, 8 pairs under 0.3 nm and 0
under 0.3 A.  The SG-SG bonds that do turn up in mdtraj's topology.bonds come
from the file's own CONECT records instead, which is the declared answer wearing
a different hat.  So what t/data/features.txt freezes is mdtraj's rule
transcribed into numpy with its units made consistent -- the same thing the
float64 surface column is, and for the same reason.

CYS and not CYX: AMBER and CHARMM rename a bonded cysteine, and matching those
spellings would be matching something the reference implementation does not.  A
structure that has been through a force field needs its residues named the way
the archive names them.*/
typedef struct {
	NV x, y, z; //where the SG is
	UV res;     //which residue it belongs to
} sulfur_t;

static UV sulfur_find(pTHX_ structset *CSP_RESTRICT s, sulfur_t *CSP_RESTRICT *CSP_RESTRICT out)
{
	sulfur_t *sg = NULL;
	UV cap = 0, n = 0, r;
	*out = NULL;
	for (r = 0; r < s->n_res; r++) {
		HV *atoms;
		NV x, y, z;
		if (s->res_key[r] != K3('C','Y','S')) continue;
		atoms = hvf_hv(aTHX_ s->res_hv[r], "atoms", 5);
		if (!atoms) continue;
		//a modelled thiol hydrogen means a reduced cysteine, which holds no bond
		if (hv_exists(atoms, "HG", 2)) continue;
		if (!atom_xyz(aTHX_ atoms, "SG", &x, &y, &z)) continue;
		if (n == cap) {
			cap = cap ? cap * 2 : 16;
			Renew(sg, cap, sulfur_t);
		}
		sg[n].x = x; sg[n].y = y; sg[n].z = z; sg[n].res = r;
		n++;
	}
	*out = sg;
	return n;
}

/*disulfides() -- every pair of those within the cutoff.

Through the same grid the surface uses, on the sulphurs alone.  A structure has
a few dozen cysteines and comparing every pair would cost nothing, which is what
mdtraj does; the grid is here because it is already written and because a
ribosome-sized entry has thousands.

A cysteine can hold only one disulfide, so a second partner means something is
wrong with the entry -- alternate conformers refined into each other, most often
-- and both pairs are reported rather than one of them silently dropped.*/
//the partner half of one bond, pushed onto a residue's own list
static void ss_note(pTHX_ HV *CSP_RESTRICT res, HV *CSP_RESTRICT partner, NV d)
{
	SV **slot = hv_fetch(res, "disulfide", 9, 0);
	AV *list;
	HV *e;
	if (slot && *slot && SvROK(*slot) && SvTYPE(SvRV(*slot)) == SVt_PVAV) {
		list = (AV *)SvRV(*slot);
	} else {
		list = newAV();
		(void)hv_stores(res, "disulfide", newRV_noinc((SV *)list));
	}
	e = newHV();
	pi_field(aTHX_ e, "chain", 5, partner, "chain", 5);
	pi_field(aTHX_ e, "residue", 7, partner, "key", 3);
	(void)hv_stores(e, "distance", newSVnv(d));
	av_push(list, newRV_noinc((SV *)e));
}

static AV *disulfides(pTHX_ structset *CSP_RESTRICT s, NV cut, bool store)
{
	AV *out = newAV();
	sulfur_t *sg = NULL;
	NV *sx = NULL, *sy = NULL, *sz = NULL;
	cell_grid g;
	UV n, i;

	n = sulfur_find(aTHX_ s, &sg);
	if (n < 2) { Safefree(sg); return out; }
	Newx(sx, n, NV); Newx(sy, n, NV); Newx(sz, n, NV);
	for (i = 0; i < n; i++) { sx[i] = sg[i].x; sy[i] = sg[i].y; sz[i] = sg[i].z; }
	/*Clear what an earlier call left on these residues, so that asking twice
	replaces the answer rather than adding to it.  The candidate set is decided
	by residue name and atom names alone, never by the cutoff, so the residues
	cleared here are exactly the ones a previous call could have written to.*/
	if (store)
		for (i = 0; i < n; i++)
			(void)hv_delete(s->res_hv[sg[i].res], "disulfide", 9, G_DISCARD);
	grid_build(aTHX_ &g, sx, sy, sz, n, cut);

	for (i = 0; i < n; i++) {
		UV bx, by, bz;
		UV ci = grid_axis(sx[i] - g.x0, g.cell, g.nx);
		UV cj = grid_axis(sy[i] - g.y0, g.cell, g.ny);
		UV ck = grid_axis(sz[i] - g.z0, g.cell, g.nz);
		UV ax0 = ci ? ci - 1 : 0, ax1 = (ci + 1 < g.nx) ? ci + 1 : g.nx - 1;
		UV ay0 = cj ? cj - 1 : 0, ay1 = (cj + 1 < g.ny) ? cj + 1 : g.ny - 1;
		UV az0 = ck ? ck - 1 : 0, az1 = (ck + 1 < g.nz) ? ck + 1 : g.nz - 1;
		for (bx = ax0; bx <= ax1; bx++)
		for (by = ay0; by <= ay1; by++)
		for (bz = az0; bz <= az1; bz++) {
			UV cell = (bx * g.ny + by) * g.nz + bz, p;
			for (p = g.start[cell]; p < g.start[cell + 1]; p++) {
				UV j = g.idx[p];
				NV dx, dy, dz, d;
				HV *h;
				if (j <= i) continue; //each pair once
				//two SG atoms of one residue are alternate conformers of the
				//same sulphur, not two sulphurs bonded to each other
				if (sg[i].res == sg[j].res) continue;
				dx = sx[j] - sx[i]; dy = sy[j] - sy[i]; dz = sz[j] - sz[i];
				d = nv_sqrt(dx * dx + dy * dy + dz * dz);
				if (!(d < cut)) continue; //`<', as the reference has it
				h = newHV();
				pi_field(aTHX_ h, "chain1", 6, s->res_hv[sg[i].res], "chain", 5);
				pi_field(aTHX_ h, "residue1", 8, s->res_hv[sg[i].res], "key", 3);
				pi_field(aTHX_ h, "chain2", 6, s->res_hv[sg[j].res], "chain", 5);
				pi_field(aTHX_ h, "residue2", 8, s->res_hv[sg[j].res], "key", 3);
				(void)hv_stores(h, "distance", newSVnv(d));
				av_push(out, newRV_noinc((SV *)h));
				if (store) {
					ss_note(aTHX_ s->res_hv[sg[i].res], s->res_hv[sg[j].res], d);
					ss_note(aTHX_ s->res_hv[sg[j].res], s->res_hv[sg[i].res], d);
				}
			}
		}
	}
	grid_free(aTHX_ &g);
	Safefree(sx); Safefree(sy); Safefree(sz);
	Safefree(sg);
	return out;
}

/*Backbone and side chain torsion angles.

phi, psi and omega say how the backbone is folded and chi1..chi5 how each side
chain is rotated, and between them they are what a Ramachandran plot and a
rotamer library are drawn from.  The atoms of each are mdtraj's, from
PHI_ATOMS, PSI_ATOMS, OMEGA_ATOMS and CHI1_ATOMS..CHI5_ATOMS in
mdtraj/geometry/dihedral.py:

  phi    C of the residue before, then N, CA, C
  psi    N, CA, C, then N of the residue after
  omega  CA, C, then N and CA of the residue after

A chi has alternatives, because the fourth atom is whatever that side chain
happens to be built from: chi1's is CG for most residues, CG1 for valine and
isoleucine, SG for cysteine, OG for serine and OG1 for threonine.  The first
alternative all of whose atoms are present is the one used, which is how mdtraj
picks too, and each chi is looked for on its own -- a residue that has chi2 and
not chi1 gets chi2 and a hole where chi1 would be, rather than neither.

**Where this departs from mdtraj, and why.**  mdtraj takes the residue before
and after to be the ones next to this one in the file, within the same chain,
and asks nothing else of them.  A chain with a gap in it -- residues nobody
could see, which is most crystal structures -- therefore gets a phi and a psi
measured across the gap, between two residues that are not bonded and may be
half the protein apart.  That is a number rather than an answer.  So the link is
tested here: the C of the earlier residue must be within 1.8 A of the N of the
later one, which is the distance Biopython's Bio.PDB.Polypeptide.PPBuilder uses
to decide the same question (its `radius' argument, default 1.8).  A chi needs
no link, being inside one residue, and is unaffected.

t/features.t holds both halves of that: mdtraj's angle for every residue whose
neighbours really are bonded, and nothing at all for the residues where it is
measuring across a gap.*/
#define CHI_MAX 5 //chi1 to chi5, which is as many as any side chain has

typedef struct {
	const char *const *atom; //the four atom names
	unsigned short int n;    //how many alternatives follow this one
} chi_alt;

static const char *const chi1_a[5][4] = {
	{ "N", "CA", "CB", "CG"  }, { "N", "CA", "CB", "CG1" }, { "N", "CA", "CB", "SG" },
	{ "N", "CA", "CB", "OG"  }, { "N", "CA", "CB", "OG1" }
};
static const char *const chi2_a[6][4] = {
	{ "CA", "CB", "CG", "CD"  }, { "CA", "CB", "CG", "CD1" }, { "CA", "CB", "CG1", "CD1" },
	{ "CA", "CB", "CG", "OD1" }, { "CA", "CB", "CG", "ND1" }, { "CA", "CB", "CG", "SD"  }
};
static const char *const chi3_a[4][4] = {
	{ "CB", "CG", "CD", "NE" }, { "CB", "CG", "CD", "CE" },
	{ "CB", "CG", "CD", "OE1" }, { "CB", "CG", "SD", "CE" }
};
static const char *const chi4_a[2][4] = {
	{ "CG", "CD", "NE", "CZ" }, { "CG", "CD", "CE", "NZ" }
};
static const char *const chi5_a[1][4] = { { "CD", "NE", "CZ", "NH1" } };

/*The torsion angle about the p1-p2 bond, in radians, by the usual construction:
project the two outer bonds onto the plane across the middle one and take the
angle between them, signed by which side of that plane the second falls on.

False when there is no such angle, which is when either outer bond lies along
the middle one -- four collinear atoms have no torsion, and every value is as
good as every other.  A structure with real coordinates never gets close, and
one built on a straight line (t/data/mini.pdb, whose atoms are laid out along a
line on purpose) is entirely made of them.  The alternative is to return the
zero that atan2(0, 0) gives, which is a number where there is no answer.

The test is on the shorter projected length against the bond it was projected
from, which is the sine of the angle between them, and the two cases are nowhere
near each other.  Measured: over phi, psi, omega, chi1 and chi2 of 1A42 the
smallest that ratio gets is 0.736, because a real bond angle is nowhere near
straight; over the deliberately collinear t/data/mini.pdb it runs from 5.4e-7 to
1.1e-6, which is float32 rounding on coordinates that are collinear exactly.
1e-3 sits a thousand times above the second and seven hundred times below the
first.  t/data/features.py applies the same rule with the same constant, so the
frozen answer and this agree about which angles exist.*/
#define CSP_COLLINEAR 1e-3

static bool dihedral4(const NV *CSP_RESTRICT p0, const NV *CSP_RESTRICT p1,
                      const NV *CSP_RESTRICT p2, const NV *CSP_RESTRICT p3,
                      NV *CSP_RESTRICT out)
{
	NV b0[3], b1[3], b2[3], v[3], w[3], cr[3];
	NV len, d0, d2, x, y;
	unsigned short int k;
	for (k = 0; k < 3; k++) {
		b0[k] = p0[k] - p1[k];
		b1[k] = p2[k] - p1[k];
		b2[k] = p3[k] - p2[k];
	}
	len = nv_sqrt(b1[0] * b1[0] + b1[1] * b1[1] + b1[2] * b1[2]);
	if (!(len > 0.0)) return FALSE; //two atoms at one position: no bond, no angle
	for (k = 0; k < 3; k++) b1[k] /= len;
	d0 = b0[0] * b1[0] + b0[1] * b1[1] + b0[2] * b1[2];
	d2 = b2[0] * b1[0] + b2[1] * b1[1] + b2[2] * b1[2];
	for (k = 0; k < 3; k++) {
		v[k] = b0[k] - d0 * b1[k];
		w[k] = b2[k] - d2 * b1[k];
	}
	{
		NV l0 = nv_sqrt(b0[0] * b0[0] + b0[1] * b0[1] + b0[2] * b0[2]);
		NV l2 = nv_sqrt(b2[0] * b2[0] + b2[1] * b2[1] + b2[2] * b2[2]);
		NV lv = nv_sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
		NV lw = nv_sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
		if (!(l0 > 0.0) || !(l2 > 0.0)) return FALSE;
		if (lv <= l0 * CSP_COLLINEAR || lw <= l2 * CSP_COLLINEAR) return FALSE;
	}
	cr[0] = b1[1] * v[2] - b1[2] * v[1];
	cr[1] = b1[2] * v[0] - b1[0] * v[2];
	cr[2] = b1[0] * v[1] - b1[1] * v[0];
	x = v[0] * w[0] + v[1] * w[1] + v[2] * w[2];
	y = cr[0] * w[0] + cr[1] * w[1] + cr[2] * w[2];
	*out = nv_atan2(y, x);
	return TRUE;
}

//the four named atoms of one residue, or of a pair of them: `which' says where
//each comes from, 0 for the first residue's hash and 1 for the second's
static bool four_atoms(pTHX_ HV *CSP_RESTRICT a, HV *CSP_RESTRICT b,
                       const char *const *CSP_RESTRICT name,
                       const unsigned short int *CSP_RESTRICT which,
                       NV p[4][3])
{
	unsigned short int k;
	for (k = 0; k < 4; k++) {
		HV *from = which[k] ? b : a;
		if (!from) return FALSE;
		if (!atom_xyz(aTHX_ from, name[k], &p[k][0], &p[k][1], &p[k][2])) return FALSE;
	}
	return TRUE;
}

//Biopython's PPBuilder radius: two residues are peptide-bonded when the C of
//the earlier is this close to the N of the later
#define CSP_PEPTIDE_BOND 1.8

static bool peptide_linked(pTHX_ HV *CSP_RESTRICT prev, HV *CSP_RESTRICT next, NV cut)
{
	HV *pa, *na;
	NV cx, cy, cz, nx, ny, nz, dx, dy, dz;
	if (!prev || !next) return FALSE;
	pa = hvf_hv(aTHX_ prev, "atoms", 5);
	na = hvf_hv(aTHX_ next, "atoms", 5);
	if (!pa || !na) return FALSE;
	if (!atom_xyz(aTHX_ pa, "C", &cx, &cy, &cz)) return FALSE;
	if (!atom_xyz(aTHX_ na, "N", &nx, &ny, &nz)) return FALSE;
	dx = nx - cx; dy = ny - cy; dz = nz - cz;
	return (dx * dx + dy * dy + dz * dz) < cut * cut;
}

//one chi, from whichever of its alternatives the residue is built from
static bool chi_of(pTHX_ HV *CSP_RESTRICT atoms, const char *const alt[][4],
                   unsigned short int n_alt, NV *CSP_RESTRICT out)
{
	static const unsigned short int here[4] = { 0, 0, 0, 0 };
	unsigned short int k;
	for (k = 0; k < n_alt; k++) {
		NV p[4][3];
		if (!four_atoms(aTHX_ atoms, NULL, alt[k], here, p)) continue;
		return dihedral4(p[0], p[1], p[2], p[3], out);
	}
	return FALSE;
}

/*Nucleic acid torsion angles, and the shape of the sugar ring.

A nucleotide's conformation is six backbone torsions, one glycosidic torsion and
the pucker of its ribose, and between them they are what tells an A-form helix
from a B-form one and either from what a loop does.  The names and the atoms of
each are the IUPAC-IUB Joint Commission on Biochemical Nomenclature's (1983)
"Abbreviations and symbols for the description of conformations of
polynucleotide chains", Eur J Biochem 131:9-15:

  alpha    O3' of the residue before, then P, O5', C5'
  beta     P, O5', C5', C4'
  gamma    O5', C5', C4', C3'
  delta    C5', C4', C3', O3'
  epsilon  C4', C3', O3', then P of the residue after
  zeta     C3', O3', then P and O5' of the residue after
  chi      O4', C1', then N9 and C4 of a purine or N1 and C2 of a pyrimidine
  nu0..nu4 the ring itself, C4'-O4'-C1'-C2' round to C3'-C4'-O4'-C1'

chi is stored under the same key an amino acid's side chain torsions are, and no
residue has both -- an amino acid has no C1' and a nucleotide has no CB -- but an
amino acid's chi is a list of up to five and a nucleotide's is one number.

Which four atoms it is turns on N9 rather than on trying one quadruple and then
the other: only a purine has an N9, while a purine has N1 and C2 as well as N9
and C4, so a purine missing its C4 would fall through to the pyrimidine atoms
and answer a different torsion under the same name.  Asking for the atom that
separates the two cases costs one lookup and cannot do that.

Nothing here asks whether the residue is a nucleotide: the atoms decide, as they
do for phi and psi.  A protein residue has no O3' and gets none of these, and a
nucleotide sitting in a binding site as a ligand rather than in a chain gets the
ones its own atoms define, which is the right answer rather than an accident.

alpha needs the residue before and epsilon and zeta the one after, and both
links are tested rather than assumed, for the reason phi and psi test theirs: a
gap in the model puts two residues next to each other that are not bonded, and a
torsion across one is a number rather than an answer.  The test is gemmi's, from
are_connected() in gemmi/polyheur.hpp -- an O3'-to-P separation under 1.5 times
the 1.6 A ideal bond, so 2.4 A -- which is the polynucleotide half of the same
function whose polypeptide half supplies the 1.8 A used above.

**The pucker.**  Altona, C; Sundaralingam, M (1972) J Am Chem Soc 94(23):8205-12
describes the ring with two numbers instead of five, on the observation that the
five nu are one sinusoid sampled at five points: a phase angle P saying which of
the ring's ten envelope shapes it is nearest, and an amplitude nu_max saying how
far from flat it is.  Their equations 1 and 2:

    tan P = ((nu4 + nu1) - (nu3 + nu0)) / (2 * nu2 * (sin 36 + sin 72))
    nu_max = nu2 / cos P

Written as an atan2 of the two halves rather than a division, which puts P in
the right half of the circle without a separate test: the denominator carries
the sign of nu2, so cos P and nu2 always agree in sign and nu_max comes out
positive, which is the paper's convention.  Both sides of the fraction are
linear in the nu, so the units cancel and the radians this file works in give
the same P as the degrees the paper is written in.

P is reported in 0 to 360 rather than the -180 to 180 every other angle here is
in, because the cycle is conventionally read that way and the pucker names are
ranges of it.  They are the ten envelope forms, each the 36 degrees of P centred
on it, as tabulated by Altona and Sundaralingam and by Saenger, W (1984)
Principles of Nucleic Acid Structure, chapter 2.  C3'-endo (P near 18) is what
an A-form helix and almost every ribose are in and C2'-endo (P near 162) what a
B-form helix is, so the name is often the whole answer.

**glycosidic** is which side of the sugar the base is turned to, read off chi:
syn within 90 degrees of 0, anti within 90 of 180, which is Saenger's division
in the same chapter and the one DSSR reports.  Nearly everything is anti, which
is what makes a syn guanosine worth being told about.  It is a bisection and not
a classification with a gap, so the band around -90 that the literature calls
high-anti comes back here as syn: chi itself is beside it for a caller who wants
the third name.*/

//gemmi's are_connected(): 1.5 times the 1.6 A ideal O3'-P bond
#define CSP_PHOSPHO_BOND 2.4

/*The ten envelope forms, each named for the ring atom furthest out of the plane
of the other four and which face of the ring it is on.  Entry k covers the phase
angle from 36k to 36(k+1) degrees.*/
static const char *const nuc_pucker[10] = {
	"C3'-endo", "C4'-exo",  "O4'-endo", "C1'-exo",  "C2'-endo",
	"C3'-exo",  "C4'-endo", "O4'-exo",  "C1'-endo", "C2'-exo"
};

//the keys below, so that asking twice replaces the answer rather than leaving
//half of an older one beside it
static const char *const nuc_keys[] = {
	"alpha", "beta", "gamma", "delta", "epsilon", "zeta", "glycosidic",
	"nu", "pucker", "pucker_phase", "pucker_amplitude"
};

//gemmi's nucleotide bond: the O3' of the earlier residue this close to the P of
//the later one
static bool nucleotide_linked(pTHX_ HV *CSP_RESTRICT prev, HV *CSP_RESTRICT next, NV cut)
{
	HV *pa, *na;
	NV ox, oy, oz, px, py, pz, dx, dy, dz;
	if (!prev || !next) return FALSE;
	pa = hvf_hv(aTHX_ prev, "atoms", 5);
	na = hvf_hv(aTHX_ next, "atoms", 5);
	if (!pa || !na) return FALSE;
	if (!atom_xyz(aTHX_ pa, "O3'", &ox, &oy, &oz)) return FALSE;
	if (!atom_xyz(aTHX_ na, "P",   &px, &py, &pz)) return FALSE;
	dx = px - ox; dy = py - oy; dz = pz - oz;
	return (dx * dx + dy * dy + dz * dz) < cut * cut;
}

//one nucleotide's torsions, onto its residue hash.  Degrees, and a torsion
//whose atoms are not all present simply has no key.
static void nuc_torsions(pTHX_ HV *CSP_RESTRICT me, HV *CSP_RESTRICT my_at,
                         HV *CSP_RESTRICT pv_at, HV *CSP_RESTRICT nx_at,
                         bool prev_linked, bool next_linked)
{
	static const unsigned short int here[4]   = { 0, 0, 0, 0 };
	static const unsigned short int prev_w[4] = { 0, 1, 1, 1 };
	static const unsigned short int eps_w[4]  = { 0, 0, 0, 1 };
	static const unsigned short int zeta_w[4] = { 0, 0, 1, 1 };
	static const char *const alpha_n[4] = { "O3'", "P", "O5'", "C5'" };
	static const char *const beta_n[4]  = { "P", "O5'", "C5'", "C4'" };
	static const char *const gamma_n[4] = { "O5'", "C5'", "C4'", "C3'" };
	static const char *const delta_n[4] = { "C5'", "C4'", "C3'", "O3'" };
	static const char *const eps_n[4]   = { "C4'", "C3'", "O3'", "P" };
	static const char *const zeta_n[4]  = { "C3'", "O3'", "P", "O5'" };
	static const char *const chi_pur[4] = { "O4'", "C1'", "N9", "C4" };
	static const char *const chi_pyr[4] = { "O4'", "C1'", "N1", "C2" };
	static const char *const nu_n[5][4] = {
		{ "C4'", "O4'", "C1'", "C2'" }, { "O4'", "C1'", "C2'", "C3'" },
		{ "C1'", "C2'", "C3'", "C4'" }, { "C2'", "C3'", "C4'", "O4'" },
		{ "C3'", "C4'", "O4'", "C1'" }
	};
	const NV to_deg = 180.0 / CSP_PI;
	NV p[4][3], v, nu[5];
	bool have_nu = TRUE;
	unsigned short int k;

	for (k = 0; k < (unsigned short int)C_ARRAY_LENGTH(nuc_keys); k++)
		(void)hv_delete(me, nuc_keys[k], (I32)strlen(nuc_keys[k]), G_DISCARD);

	if (prev_linked && four_atoms(aTHX_ pv_at, my_at, alpha_n, prev_w, p)
	    && dihedral4(p[0], p[1], p[2], p[3], &v))
		(void)hv_stores(me, "alpha", newSVnv(v * to_deg));
	if (four_atoms(aTHX_ my_at, NULL, beta_n, here, p)
	    && dihedral4(p[0], p[1], p[2], p[3], &v))
		(void)hv_stores(me, "beta", newSVnv(v * to_deg));
	if (four_atoms(aTHX_ my_at, NULL, gamma_n, here, p)
	    && dihedral4(p[0], p[1], p[2], p[3], &v))
		(void)hv_stores(me, "gamma", newSVnv(v * to_deg));
	if (four_atoms(aTHX_ my_at, NULL, delta_n, here, p)
	    && dihedral4(p[0], p[1], p[2], p[3], &v))
		(void)hv_stores(me, "delta", newSVnv(v * to_deg));
	if (next_linked) {
		if (four_atoms(aTHX_ my_at, nx_at, eps_n, eps_w, p)
		    && dihedral4(p[0], p[1], p[2], p[3], &v))
			(void)hv_stores(me, "epsilon", newSVnv(v * to_deg));
		if (four_atoms(aTHX_ my_at, nx_at, zeta_n, zeta_w, p)
		    && dihedral4(p[0], p[1], p[2], p[3], &v))
			(void)hv_stores(me, "zeta", newSVnv(v * to_deg));
	}
	{
		NV n9[3];
		const char *const *CSP_RESTRICT chi_n =
			atom_xyz(aTHX_ my_at, "N9", &n9[0], &n9[1], &n9[2]) ? chi_pur : chi_pyr;
		if (four_atoms(aTHX_ my_at, NULL, chi_n, here, p)
		    && dihedral4(p[0], p[1], p[2], p[3], &v)) {
			const NV deg = v * to_deg;
			(void)hv_stores(me, "chi", newSVnv(deg));
			(void)hv_stores(me, "glycosidic",
			                newSVpv(nv_fabs(deg) <= 90.0 ? "syn" : "anti", 0));
		}
	}

	for (k = 0; k < 5; k++) {
		if (four_atoms(aTHX_ my_at, NULL, nu_n[k], here, p)
		    && dihedral4(p[0], p[1], p[2], p[3], &nu[k])) continue;
		have_nu = FALSE;
		break;
	}
	if (have_nu) {
		/*sin 36 + sin 72, the constant of Altona and Sundaralingam's equation 1.
		Computed rather than written out, so that it is exact to whatever width
		this perl's NV is rather than to however many digits were typed.*/
		const NV k1 = nv_sin(CSP_PI / 5.0) + nv_sin(2.0 * CSP_PI / 5.0);
		const NV num = (nu[4] + nu[1]) - (nu[3] + nu[0]);
		const NV den = 2.0 * nu[2] * k1;
		AV *list = newAV();
		for (k = 0; k < 5; k++) av_push(list, newSVnv(nu[k] * to_deg));
		(void)hv_stores(me, "nu", newRV_noinc((SV *)list));
		/*A ring flat to the last bit has no phase -- every P describes it as
		well as every other, and atan2(0, 0) would answer 0 as though one of
		them did.  Real geometry never gets there; a fixture built on a straight
		line has no nu at all, dihedral4() having refused them.*/
		if (num != 0.0 || den != 0.0) {
			const NV pr = nv_atan2(num, den);
			const NV cp = nv_cos(pr);
			NV phase = pr * to_deg;
			if (phase < 0.0) phase += 360.0;
			//an angle a hair below zero adds the whole turn and rounds to it,
			//and 360 is 0: the cycle closes, and the range is documented as
			//0 to 360 with the upper end not in it
			if (phase >= 360.0) phase = 0.0;
			(void)hv_stores(me, "pucker_phase", newSVnv(phase));
			(void)hv_stores(me, "pucker",
			                newSVpv(nuc_pucker[(unsigned short int)(phase / 36.0)], 0));
			//cos P is zero only when nu2 is exactly zero, where the amplitude
			//is the other half of the fraction; there is no ring like that in
			//a real file and none is invented here
			if (cp != 0.0)
				(void)hv_stores(me, "pucker_amplitude", newSVnv(nu[2] / cp * to_deg));
		}
	}
}

/*chain_torsions() -- the same angles again, as one array per torsion.

The block above writes each angle onto the residue it belongs to, which is where
a question about one residue is asked.  A question about the chain -- a
Ramachandran plot, a rotamer census, the mean of an angle along a helix -- wants
them as columns instead, so they are gathered onto the chain as well, under
`torsions': $info->{chains}{A}{torsions}{phi} and its relatives.

Each array runs parallel to the chain's `residue_order', one element per residue
in the order the walk saw them, and a residue with no such torsion holds an
undef rather than being left out -- the position is what says which residue a
value came from, and a compacted array would lose that.  A key no residue in the
chain has at all is absent instead, so a protein chain carries phi, psi, omega
and chi and a nucleic acid one alpha to zeta, and neither carries eleven arrays
of nothing.

The elements are copies of the residues' scalars, so writing to one does not
reach back; `chi' and `nu' copy the reference, which means the chain and the
residue name one list between them, as two references to a thing always do.*/

/*every key the torsions above write onto a residue: the amino acid ones, then
the nucleotide ones, which are nuc_keys[] with the chi the two kinds share
lifted out of it*/
static const char *const tors_name[] = {
	"phi", "psi", "omega", "chi",
	"alpha", "beta", "gamma", "delta", "epsilon", "zeta",
	"nu", "pucker", "pucker_phase", "pucker_amplitude", "glycosidic"
};
#define NTORS ((unsigned short int)(sizeof tors_name / sizeof *tors_name))

static void chain_torsions(pTHX_ HV *CSP_RESTRICT chain,
                           HV *CSP_RESTRICT *CSP_RESTRICT res, UV n)
{
	HV *out = newHV();
	AV *col[NTORS];
	STRLEN len[NTORS];
	bool any[NTORS];   //FALSE until some residue turns out to have the key
	unsigned short int k;
	UV i;

	for (k = 0; k < NTORS; k++) {
		col[k] = newAV();
		len[k] = strlen(tors_name[k]);
		any[k] = FALSE;
		//n undefs, of which the loop below fills in the ones there are
		if (n) av_fill(col[k], (SSize_t)n - 1);
	}
	for (i = 0; i < n; i++)
		for (k = 0; k < NTORS; k++) {
			SV *v = hvf_sv(aTHX_ res[i], tors_name[k], len[k]);
			if (!v) continue;
			(void)av_store(col[k], (SSize_t)i, newSVsv(v));
			any[k] = TRUE;
		}
	for (k = 0; k < NTORS; k++) {
		if (any[k])
			(void)hv_store(out, tors_name[k], (I32)len[k],
			               newRV_noinc((SV *)col[k]), 0);
		else
			SvREFCNT_dec((SV *)col[k]);
	}
	//replaced whole rather than added to, so that asking twice gives one answer
	(void)hv_stores(chain, "torsions", newRV_noinc((SV *)out));
}

/*dihedrals() -- every torsion angle of every residue, written onto the residue.

Degrees, in -180 to 180, which is what the rest of the module reports an angle
in.  A residue that has no phi -- the first of a chain, or the one after a gap
-- simply has no phi key, rather than a zero that would plot.  The one exception
is the pucker phase, which nuc_torsions() reports in 0 to 360 because that is
how the pseudorotation cycle is read.

One loop for both kinds of residue: the protein torsions and the nucleic ones
want the same three residue hashes and the same two link tests, and a residue is
one or the other, so asking twice would be walking twice.  `bond' is the peptide
cutoff and `phospho' the phosphodiester one.

Each chain is then handed to chain_torsions(), which collects what its residues
were just given into one array per angle; that is the same walk again, on one
chain's worth of residue hashes, rather than a second pass over the structure.*/
static void dihedrals(pTHX_ structset *CSP_RESTRICT s, NV bond, NV phospho)
{
	static const char *const phi_n[4]   = { "C", "N", "CA", "C" };
	static const unsigned short int phi_w[4] = { 0, 1, 1, 1 };
	static const char *const psi_n[4]   = { "N", "CA", "C", "N" };
	static const unsigned short int psi_w[4] = { 0, 0, 0, 1 };
	static const char *const omg_n[4]   = { "CA", "C", "N", "CA" };
	static const unsigned short int omg_w[4] = { 0, 0, 1, 1 };
	const NV to_deg = 180.0 / CSP_PI;
	UV c, r;

	for (c = 0; c < s->n_chain; c++) {
		for (r = s->chain_first[c]; r < s->chain_last[c]; r++) {
			HV *me = s->res_hv[r];
			HV *prev = (r > s->chain_first[c]) ? s->res_hv[r - 1] : NULL;
			HV *next = (r + 1 < s->chain_last[c]) ? s->res_hv[r + 1] : NULL;
			HV *my_at = hvf_hv(aTHX_ me, "atoms", 5);
			HV *pv_at = prev ? hvf_hv(aTHX_ prev, "atoms", 5) : NULL;
			HV *nx_at = next ? hvf_hv(aTHX_ next, "atoms", 5) : NULL;
			NV p[4][3], v;
			unsigned short int k, last = 0;
			NV chi[CHI_MAX];
			bool have[CHI_MAX];
			if (!my_at) continue;
			//everything written here is replaced rather than added to, so that
			//asking twice gives one answer
			(void)hv_delete(me, "phi",   3, G_DISCARD);
			(void)hv_delete(me, "psi",   3, G_DISCARD);
			(void)hv_delete(me, "omega", 5, G_DISCARD);
			(void)hv_delete(me, "chi",   3, G_DISCARD);

			if (peptide_linked(aTHX_ prev, me, bond)
			    && four_atoms(aTHX_ pv_at, my_at, phi_n, phi_w, p)
			    && dihedral4(p[0], p[1], p[2], p[3], &v))
				(void)hv_stores(me, "phi", newSVnv(v * to_deg));
			if (peptide_linked(aTHX_ me, next, bond)) {
				if (four_atoms(aTHX_ my_at, nx_at, psi_n, psi_w, p)
				    && dihedral4(p[0], p[1], p[2], p[3], &v))
					(void)hv_stores(me, "psi", newSVnv(v * to_deg));
				if (four_atoms(aTHX_ my_at, nx_at, omg_n, omg_w, p)
				    && dihedral4(p[0], p[1], p[2], p[3], &v))
					(void)hv_stores(me, "omega", newSVnv(v * to_deg));
			}

			have[0] = chi_of(aTHX_ my_at, chi1_a, 5, &chi[0]);
			have[1] = chi_of(aTHX_ my_at, chi2_a, 6, &chi[1]);
			have[2] = chi_of(aTHX_ my_at, chi3_a, 4, &chi[2]);
			have[3] = chi_of(aTHX_ my_at, chi4_a, 2, &chi[3]);
			have[4] = chi_of(aTHX_ my_at, chi5_a, 1, &chi[4]);
			for (k = 0; k < CHI_MAX; k++) if (have[k]) last = (unsigned short int)(k + 1);
			if (last) {
				AV *list = newAV();
				for (k = 0; k < last; k++)
					av_push(list, have[k] ? newSVnv(chi[k] * to_deg) : newSV(0));
				(void)hv_stores(me, "chi", newRV_noinc((SV *)list));
			}

			nuc_torsions(aTHX_ me, my_at, pv_at, nx_at,
			             nucleotide_linked(aTHX_ prev, me, phospho),
			             nucleotide_linked(aTHX_ me, next, phospho));
		}
		//and the chain's own view of what the residues above were given
		chain_torsions(aTHX_ s->chain_hv[c], s->res_hv + s->chain_first[c],
		               s->chain_last[c] - s->chain_first[c]);
	}
}

/*Residue contacts, and how much of a residue's own half of the world is full.

Two things that both answer "what is near this residue", from two different
implementations, and both through the grid the surface already builds.

**Contacts** are mdtraj's compute_contacts() with its default `closest-heavy'
scheme: the distance between two residues is the shortest distance between any
pair of their heavy atoms, hydrogens left out.  mdtraj's `all' pairs up
residues in the same chain that are three or more apart in it, which is the
question a contact map of one chain asks; this reports those and the pairs that
cross a chain as well, because a caller looking at a complex wants the interface
and dropping it would be a strange thing to do silently.  t/features.t compares
the subset mdtraj has an opinion about.

**Half-sphere exposure** is Biopython's Bio.PDB.HSExposure.HSExposureCB: draw
the vector from a residue's CA to its CB, count the CA atoms of other amino acid
residues within 12 A, and split that count by which side of the plane through CA
they fall -- hse_up towards the side chain, hse_down away from it.  It says
something the accessible surface does not: a residue can be buried and still
have its side chain pointing into a cavity.  Glycine has no CB, and gets the
virtual one Biopython builds, the N position rotated 120 degrees the other way
about the CA-C axis.

The radius, the offset of 0 (no flanking residues skipped) and the CB variant
rather than the CA one are Biopython's defaults.*/
#define CSP_HSE_RADIUS 12.0 //HSExposureCB's `radius'

//is this atom name a hydrogen or deuterium?  Contacts are between heavy atoms,
//and the element is what says so -- a name beginning H is not enough, HG being
//mercury as often as a hydrogen
static bool atom_is_h(pTHX_ HV *CSP_RESTRICT a)
{
	SV *e = hvf_sv(aTHX_ a, "element", 7);
	STRLEN n;
	const char *p;
	if (!e) return FALSE;
	p = SvPV_const(e, n);
	return n == 1 && (*p == 'H' || *p == 'D');
}

//glycine's CB, where there is none: the N position turned -120 degrees about
//the CA-to-C axis, which is Biopython's _get_gly_cb_vector().  Returns the
//CA-to-CB vector, as the real one does.
static bool gly_cb(pTHX_ HV *CSP_RESTRICT atoms, NV *CSP_RESTRICT out)
{
	NV nx, ny, nz, cx, cy, cz, ax, ay, az;
	NV k[3], v[3], len, c, sn, dot;
	const NV ang = -CSP_PI * 120.0 / 180.0;
	unsigned short int i;
	if (!atom_xyz(aTHX_ atoms, "N",  &nx, &ny, &nz)) return FALSE;
	if (!atom_xyz(aTHX_ atoms, "C",  &cx, &cy, &cz)) return FALSE;
	if (!atom_xyz(aTHX_ atoms, "CA", &ax, &ay, &az)) return FALSE;
	v[0] = nx - ax; v[1] = ny - ay; v[2] = nz - az;
	k[0] = cx - ax; k[1] = cy - ay; k[2] = cz - az;
	len = nv_sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);
	if (!(len > 0.0)) return FALSE;
	for (i = 0; i < 3; i++) k[i] /= len;
	c = nv_cos(ang);
	sn = nv_sin(ang);
	dot = k[0] * v[0] + k[1] * v[1] + k[2] * v[2];
	//Rodrigues: v cos + (k x v) sin + k (k.v)(1 - cos)
	out[0] = v[0] * c + (k[1] * v[2] - k[2] * v[1]) * sn + k[0] * dot * (1.0 - c);
	out[1] = v[1] * c + (k[2] * v[0] - k[0] * v[2]) * sn + k[1] * dot * (1.0 - c);
	out[2] = v[2] * c + (k[0] * v[1] - k[1] * v[0]) * sn + k[2] * dot * (1.0 - c);
	return TRUE;
}

//one residue's running list of the residues it touches and how closely
typedef struct {
	UV res;
	NV d;
} touch;

static AV *contacts_find(pTHX_ structset *CSP_RESTRICT s, NV cut, bool store)
{
	AV *out = newAV();
	unsigned char *CSP_RESTRICT heavy = NULL;
	NV *CSP_RESTRICT hx = NULL, *CSP_RESTRICT hy = NULL, *CSP_RESTRICT hz = NULL;
	UV *CSP_RESTRICT hres = NULL, *CSP_RESTRICT hidx = NULL;
	touch *CSP_RESTRICT near = NULL;
	UV near_cap = 32;
	cell_grid g;
	UV n = 0, i, r;

	if (s->n_atom == 0) return out;
	Newx(heavy, s->n_atom, unsigned char);
	for (i = 0; i < s->n_atom; i++)
		heavy[i] = atom_is_h(aTHX_ s->atom_hv[i]) ? 0 : 1;
	Newx(hx, s->n_atom, NV); Newx(hy, s->n_atom, NV); Newx(hz, s->n_atom, NV);
	Newx(hres, s->n_atom, UV); Newx(hidx, s->n_atom, UV);
	for (r = 0; r < s->n_res; r++) {
		for (i = s->res_first[r]; i < s->res_last[r]; i++) {
			if (!heavy[i]) continue;
			hx[n] = s->x[i]; hy[n] = s->y[i]; hz[n] = s->z[i];
			hres[n] = r;
			hidx[n] = i;
			n++;
		}
	}
	Safefree(heavy);
	if (n == 0) {
		Safefree(hx); Safefree(hy); Safefree(hz); Safefree(hres); Safefree(hidx);
		return out;
	}
	grid_build(aTHX_ &g, hx, hy, hz, n, cut);
	Newx(near, near_cap, touch);

	if (store)
		for (r = 0; r < s->n_res; r++)
			(void)hv_stores(s->res_hv[r], "n_contacts", newSVuv(0));

	for (r = 0; r < s->n_res; r++) {
		UV n_near = 0, k;
		for (i = 0; i < n; i++) {
			UV bx, by, bz, ci, cj, ck;
			if (hres[i] != r) continue;
			ci = grid_axis(hx[i] - g.x0, g.cell, g.nx);
			cj = grid_axis(hy[i] - g.y0, g.cell, g.ny);
			ck = grid_axis(hz[i] - g.z0, g.cell, g.nz);
			for (bx = ci ? ci - 1 : 0; bx <= (ci + 1 < g.nx ? ci + 1 : g.nx - 1); bx++)
			for (by = cj ? cj - 1 : 0; by <= (cj + 1 < g.ny ? cj + 1 : g.ny - 1); by++)
			for (bz = ck ? ck - 1 : 0; bz <= (ck + 1 < g.nz ? ck + 1 : g.nz - 1); bz++) {
				UV cell = (bx * g.ny + by) * g.nz + bz, p;
				for (p = g.start[cell]; p < g.start[cell + 1]; p++) {
					UV j = g.idx[p];
					NV dx, dy, dz, d2, d;
					if (hres[j] <= r) continue; //each residue pair once
					dx = hx[j] - hx[i]; dy = hy[j] - hy[i]; dz = hz[j] - hz[i];
					d2 = dx * dx + dy * dy + dz * dz;
					if (d2 >= cut * cut) continue;
					d = nv_sqrt(d2);
					for (k = 0; k < n_near; k++)
						if (near[k].res == hres[j]) break;
					if (k < n_near) {
						if (d < near[k].d) near[k].d = d;
					} else {
						if (n_near == near_cap) { near_cap *= 2; Renew(near, near_cap, touch); }
						near[n_near].res = hres[j];
						near[n_near].d = d;
						n_near++;
					}
				}
			}
		}
		for (k = 0; k < n_near; k++) {
			HV *h = newHV();
			pi_field(aTHX_ h, "chain1", 6, s->res_hv[r], "chain", 5);
			pi_field(aTHX_ h, "residue1", 8, s->res_hv[r], "key", 3);
			pi_field(aTHX_ h, "chain2", 6, s->res_hv[near[k].res], "chain", 5);
			pi_field(aTHX_ h, "residue2", 8, s->res_hv[near[k].res], "key", 3);
			(void)hv_stores(h, "distance", newSVnv(near[k].d));
			av_push(out, newRV_noinc((SV *)h));
			if (!store) continue;
			{
				SV **a = hv_fetch(s->res_hv[r], "n_contacts", 10, 0);
				SV **b = hv_fetch(s->res_hv[near[k].res], "n_contacts", 10, 0);
				if (a && *a) sv_setuv(*a, SvUV(*a) + 1);
				if (b && *b) sv_setuv(*b, SvUV(*b) + 1);
			}
		}
	}
	Safefree(near);
	grid_free(aTHX_ &g);
	Safefree(hx); Safefree(hy); Safefree(hz); Safefree(hres); Safefree(hidx);
	return out;
}

//half-sphere exposure, onto each amino acid residue that has a CA and a CB
static void hse_compute(pTHX_ structset *CSP_RESTRICT s, NV radius)
{
	NV *CSP_RESTRICT cax = NULL, *CSP_RESTRICT cay = NULL, *CSP_RESTRICT caz = NULL;
	UV *CSP_RESTRICT cares = NULL;
	cell_grid g;
	UV n = 0, r;

	if (s->n_res == 0) return;
	Newx(cax, s->n_res, NV); Newx(cay, s->n_res, NV); Newx(caz, s->n_res, NV);
	Newx(cares, s->n_res, UV);
	for (r = 0; r < s->n_res; r++) {
		HV *at;
		/*One of the twenty, not merely an amino acid.  Biopython's
		CaPPBuilder.build_peptides() runs with aa_only => 1, which keeps only
		the standard names, so a selenomethionine is invisible to HSExposureCB
		-- it neither gets a figure nor counts towards anybody else's.  Matching
		that is what lets t/features.t hold this to it; the cost is that a
		structure full of modified residues gets a count that ignores them.*/
		if (s->res_type[r] != RT_AA || !s->res_std[r]) continue;
		at = hvf_hv(aTHX_ s->res_hv[r], "atoms", 5);
		if (!at) continue;
		if (!atom_xyz(aTHX_ at, "CA", &cax[n], &cay[n], &caz[n])) continue;
		cares[n] = r;
		n++;
	}
	if (n == 0) {
		Safefree(cax); Safefree(cay); Safefree(caz); Safefree(cares);
		return;
	}
	grid_build(aTHX_ &g, cax, cay, caz, n, radius);

	for (r = 0; r < n; r++) {
		HV *me = s->res_hv[cares[r]];
		HV *at = hvf_hv(aTHX_ me, "atoms", 5);
		NV pcb[3], bx, by, bz;
		UV up = 0, down = 0, ax, ay, az;
		if (!at) continue;
		if (atom_xyz(aTHX_ at, "CB", &bx, &by, &bz)) {
			pcb[0] = bx - cax[r]; pcb[1] = by - cay[r]; pcb[2] = bz - caz[r];
		} else if (!gly_cb(aTHX_ at, pcb)) {
			continue; //no side chain direction: no half to be on
		}
		{
			UV ci = grid_axis(cax[r] - g.x0, g.cell, g.nx);
			UV cj = grid_axis(cay[r] - g.y0, g.cell, g.ny);
			UV ck = grid_axis(caz[r] - g.z0, g.cell, g.nz);
			for (ax = ci ? ci - 1 : 0; ax <= (ci + 1 < g.nx ? ci + 1 : g.nx - 1); ax++)
			for (ay = cj ? cj - 1 : 0; ay <= (cj + 1 < g.ny ? cj + 1 : g.ny - 1); ay++)
			for (az = ck ? ck - 1 : 0; az <= (ck + 1 < g.nz ? ck + 1 : g.nz - 1); az++) {
				UV cell = (ax * g.ny + ay) * g.nz + az, p;
				for (p = g.start[cell]; p < g.start[cell + 1]; p++) {
					UV j = g.idx[p];
					NV dx, dy, dz;
					if (j == r) continue;
					dx = cax[j] - cax[r]; dy = cay[j] - cay[r]; dz = caz[j] - caz[r];
					if (dx * dx + dy * dy + dz * dz >= radius * radius) continue;
					if (dx * pcb[0] + dy * pcb[1] + dz * pcb[2] > 0.0) up++;
					else down++;
				}
			}
		}
		(void)hv_stores(me, "hse_up", newSVuv(up));
		(void)hv_stores(me, "hse_down", newSVuv(down));
	}
	grid_free(aTHX_ &g);
	Safefree(cax); Safefree(cay); Safefree(caz); Safefree(cares);
}

/*Backbone hydrogen bonds, by Kabsch and Sander's electrostatic definition.

Kabsch, W; Sander, C (1983) Biopolymers 22(12):2577-637 put a charge of -0.42 e
on the carbonyl oxygen and +0.20 e on the amide hydrogen, their opposites on the
carbon and the nitrogen, and called the four-way Coulomb sum the bond energy:

    E = 0.42 * 0.20 * 332 * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN)   kcal/mol

with the distances in angstrom.  Anything below -0.5 kcal/mol is a bond.  The
implementation is mdtraj's kabsch_sander(), and every rule below was checked
against it pair by pair on four real entries -- 1A42, 1A22, 1AHW and 3AU6 --
where the two agree on which bonds exist and on their energies to 4e-5 kcal/mol,
which is float32 rounding on mdtraj's side:

  * the amide hydrogen is not read from the file, it is placed: one angstrom
    from N along the direction of the previous residue's C=O.  A crystal
    structure has no hydrogens, and this is the construction that lets the
    definition be used on one anyway.
  * proline has no amide hydrogen at all and can never donate.
  * a donor's own carbonyl, and the one its hydrogen was built from, are not
    acceptors: acceptor == donor and acceptor == donor - 1 are both skipped.
    The second is the important one -- an "H-bond" to the very carbonyl that
    placed the hydrogen is an artefact of the construction.
  * only pairs whose CA atoms are within 9 A are considered, which is DSSP's
    own screen and is what the grid here is built at.
  * a nitrogen has two lone pairs' worth of hydrogen to give and no more, so
    each donor keeps its two best acceptors and no others.
  * the energy is clipped at -9.9 kcal/mol, for the pathological case of two
    atoms refined on top of each other.

**Where this departs from mdtraj**, in the same way and for the same reason as
the torsion angles: mdtraj takes "the previous residue" to be whichever one
comes before this one in the file and builds the hydrogen from its carbonyl,
even across a chain break, which puts the hydrogen somewhere the chemistry never
put it.  Here the two must be peptide-bonded first.  t/features.t freezes
mdtraj's answer together with whether the link is real, and holds this to the
first where the second is true.*/
#define CSP_KS_Q      27.888 //0.42 * 0.20 * 332, kcal A / mol
#define CSP_KS_CUTOFF (-0.5) //kcal/mol: below this it is a bond
#define CSP_KS_FLOOR  (-9.9) //kcal/mol: no bond is reported as stronger
#define CSP_KS_CA     9.0    //angstrom: DSSP's CA-to-CA screen
#define CSP_KS_KEEP   2      //acceptors per donor

//the backbone atoms one residue offers, and whether it has them all
typedef struct {
	NV n[3], ca[3], c[3], o[3];
	bool whole;   //all four present, so it can donate and accept
	bool proline; //no amide hydrogen: it can accept but never donate
	bool linked;  //peptide-bonded to the residue before it in the walk
} backbone;

static void backbone_read(pTHX_ structset *CSP_RESTRICT s, backbone *CSP_RESTRICT bb, NV bond)
{
	UV r;
	for (r = 0; r < s->n_res; r++) {
		HV *at = hvf_hv(aTHX_ s->res_hv[r], "atoms", 5);
		backbone *b = &bb[r];
		b->whole = FALSE;
		b->proline = (s->res_key[r] == K3('P','R','O'));
		b->linked = FALSE;
		if (!at) continue;
		if (!atom_xyz(aTHX_ at, "N",  &b->n[0],  &b->n[1],  &b->n[2]))  continue;
		if (!atom_xyz(aTHX_ at, "CA", &b->ca[0], &b->ca[1], &b->ca[2])) continue;
		if (!atom_xyz(aTHX_ at, "C",  &b->c[0],  &b->c[1],  &b->c[2]))  continue;
		if (!atom_xyz(aTHX_ at, "O",  &b->o[0],  &b->o[1],  &b->o[2]))  continue;
		b->whole = TRUE;
	}
	for (r = 1; r < s->n_res; r++)
		bb[r].linked = peptide_linked(aTHX_ s->res_hv[r - 1], s->res_hv[r], bond);
}

//the energy of the bond from acceptor a's C=O to donor d's N-H, with the
//hydrogen already placed
static NV ks_energy(const backbone *CSP_RESTRICT a, const NV *CSP_RESTRICT h,
                    const NV *CSP_RESTRICT nd)
{
	NV r_on = 0.0, r_ch = 0.0, r_oh = 0.0, r_cn = 0.0, e;
	unsigned short int k;
	for (k = 0; k < 3; k++) {
		NV d;
		d = a->o[k] - nd[k];   r_on += d * d;
		d = a->c[k] - h[k];    r_ch += d * d;
		d = a->o[k] - h[k];    r_oh += d * d;
		d = a->c[k] - nd[k];   r_cn += d * d;
	}
	r_on = nv_sqrt(r_on); r_ch = nv_sqrt(r_ch);
	r_oh = nv_sqrt(r_oh); r_cn = nv_sqrt(r_cn);
	if (!(r_on > 0.0) || !(r_ch > 0.0) || !(r_oh > 0.0) || !(r_cn > 0.0))
		return CSP_KS_FLOOR;
	e = CSP_KS_Q * (1.0 / r_on + 1.0 / r_ch - 1.0 / r_oh - 1.0 / r_cn);
	return (e < CSP_KS_FLOOR) ? CSP_KS_FLOOR : e;
}

/*ks_compute() -- each residue's two best acceptors, and their energies.

Kept on the structset rather than returned, because the secondary structure
below reads the same table: DSSP is these bonds and nothing else, read for the
patterns they fall into.  acc[r][k] is the residue whose carbonyl residue r's
amide hydrogen is bonded to, or the count of residues when there is none.*/
static void ks_compute(pTHX_ structset *CSP_RESTRICT s, const backbone *CSP_RESTRICT bb)
{
	NV *CSP_RESTRICT cx = NULL, *CSP_RESTRICT cy = NULL, *CSP_RESTRICT cz = NULL;
	UV *CSP_RESTRICT which = NULL;
	cell_grid g;
	UV n = 0, r, k;

	Newx(s->ks_acc, s->n_res ? s->n_res * CSP_KS_KEEP : 1, UV);
	Newx(s->ks_e,   s->n_res ? s->n_res * CSP_KS_KEEP : 1, NV);
	for (r = 0; r < s->n_res * CSP_KS_KEEP; r++) {
		s->ks_acc[r] = s->n_res;  //"none": one past the last residue
		s->ks_e[r] = 0.0;
	}
	if (s->n_res == 0) return;

	Newx(cx, s->n_res, NV); Newx(cy, s->n_res, NV); Newx(cz, s->n_res, NV);
	Newx(which, s->n_res, UV);
	for (r = 0; r < s->n_res; r++) {
		if (!bb[r].whole) continue;
		cx[n] = bb[r].ca[0]; cy[n] = bb[r].ca[1]; cz[n] = bb[r].ca[2];
		which[n] = r;
		n++;
	}
	if (n == 0) {
		Safefree(cx); Safefree(cy); Safefree(cz); Safefree(which);
		return;
	}
	grid_build(aTHX_ &g, cx, cy, cz, n, CSP_KS_CA);

	for (k = 0; k < n; k++) {
		UV d = which[k];
		NV h[3], len;
		UV bx, by, bz, ci, cj, ck;
		unsigned short int j;
		if (bb[d].proline || !bb[d].linked) continue;
		if (d == 0 || !bb[d - 1].whole) continue;
		len = 0.0;
		for (j = 0; j < 3; j++) {
			h[j] = bb[d - 1].c[j] - bb[d - 1].o[j];
			len += h[j] * h[j];
		}
		len = nv_sqrt(len);
		if (!(len > 0.0)) continue;
		for (j = 0; j < 3; j++) h[j] = bb[d].n[j] + h[j] / len;

		ci = grid_axis(cx[k] - g.x0, g.cell, g.nx);
		cj = grid_axis(cy[k] - g.y0, g.cell, g.ny);
		ck = grid_axis(cz[k] - g.z0, g.cell, g.nz);
		for (bx = ci ? ci - 1 : 0; bx <= (ci + 1 < g.nx ? ci + 1 : g.nx - 1); bx++)
		for (by = cj ? cj - 1 : 0; by <= (cj + 1 < g.ny ? cj + 1 : g.ny - 1); by++)
		for (bz = ck ? ck - 1 : 0; bz <= (ck + 1 < g.nz ? ck + 1 : g.nz - 1); bz++) {
			UV cell = (bx * g.ny + by) * g.nz + bz, p;
			for (p = g.start[cell]; p < g.start[cell + 1]; p++) {
				UV q = g.idx[p], a = which[q];
				NV dx, dy, dz, e;
				if (a == d || (d > 0 && a == d - 1)) continue;
				dx = cx[q] - cx[k]; dy = cy[q] - cy[k]; dz = cz[q] - cz[k];
				if (dx * dx + dy * dy + dz * dz >= CSP_KS_CA * CSP_KS_CA) continue;
				e = ks_energy(&bb[a], h, bb[d].n);
				if (!(e < CSP_KS_CUTOFF)) continue;
				//keep the two best, best first
				if (e < s->ks_e[d * CSP_KS_KEEP] || s->ks_acc[d * CSP_KS_KEEP] == s->n_res) {
					s->ks_acc[d * CSP_KS_KEEP + 1] = s->ks_acc[d * CSP_KS_KEEP];
					s->ks_e[d * CSP_KS_KEEP + 1]   = s->ks_e[d * CSP_KS_KEEP];
					s->ks_acc[d * CSP_KS_KEEP] = a;
					s->ks_e[d * CSP_KS_KEEP]   = e;
				} else if (e < s->ks_e[d * CSP_KS_KEEP + 1]
				           || s->ks_acc[d * CSP_KS_KEEP + 1] == s->n_res) {
					s->ks_acc[d * CSP_KS_KEEP + 1] = a;
					s->ks_e[d * CSP_KS_KEEP + 1]   = e;
				}
			}
		}
	}
	grid_free(aTHX_ &g);
	Safefree(cx); Safefree(cy); Safefree(cz); Safefree(which);
}

//the bonds ks_compute() found, as a list a caller can read
static AV *hbond_list(pTHX_ structset *CSP_RESTRICT s)
{
	AV *out = newAV();
	UV d;
	unsigned short int k;
	for (d = 0; d < s->n_res; d++) {
		for (k = 0; k < CSP_KS_KEEP; k++) {
			UV a = s->ks_acc[d * CSP_KS_KEEP + k];
			HV *h;
			if (a >= s->n_res) continue;
			h = newHV();
			//the donor gives the hydrogen, the acceptor the carbonyl
			pi_field(aTHX_ h, "donor_chain",    11, s->res_hv[d], "chain", 5);
			pi_field(aTHX_ h, "donor_residue",  13, s->res_hv[d], "key", 3);
			pi_field(aTHX_ h, "acceptor_chain",   14, s->res_hv[a], "chain", 5);
			pi_field(aTHX_ h, "acceptor_residue", 16, s->res_hv[a], "key", 3);
			(void)hv_stores(h, "energy", newSVnv(s->ks_e[d * CSP_KS_KEEP + k]));
			av_push(out, newRV_noinc((SV *)h));
		}
	}
	return out;
}

/*Secondary structure, by the Kabsch-Sander dictionary.

DSSP is the hydrogen bonds above and nothing else: it reads the patterns they
fall into and gives every residue one of eight letters.  Kabsch, W; Sander, C
(1983) Biopolymers 22(12):2577-637 is the definition and mdtraj's compute_dssp()
is the implementation checked against.

  H  alpha helix: a 4-turn beginning at this residue and the one before it
  G  3-10 helix, from 3-turns; I  pi helix, from 5-turns
  E  extended strand: a residue in a bridge that has a neighbour in one
  B  isolated beta bridge: a bridge with no neighbouring bridge
  T  hydrogen-bonded turn: in an n-turn but in none of the above
  S  bend: the chain turns through more than 70 degrees here
     coil: none of them, written as a space

An n-turn at i is a bond from the C=O of i to the N-H of i+n.  A bridge between
i and j is the four-bond pattern that makes two strands run beside each other,
parallel or antiparallel.  The order above is the priority: a residue that
qualifies for two letters gets the earlier one.

Bridges are searched from the bonds rather than over every pair of residues.
Each of the four patterns names a bond, so every bridge has at least one bond in
it and the candidates are the handful of pairs each bond can belong to -- which
turns a search that is quadratic in the residues into one that is linear in the
bonds.  A ribosome has half a million residues and a few hundred thousand bonds.

Chain boundaries are respected, as mdtraj's DSSP does (it is handed the chain of
every residue): no turn, bridge or bend reaches across one.*/
#define SS_COIL  ' '
#define SS_BEND  'S'
#define SS_TURN  'T'
#define SS_PI    'I'
#define SS_G310  'G'
#define SS_BRIDGE 'B'
#define SS_STRAND 'E'
#define SS_HELIX  'H'

//is the C=O of `a' bonded to the N-H of `d'?  The table keeps two acceptors per
//donor, so this is two comparisons rather than a search
static bool ks_bonded(const structset *CSP_RESTRICT s, UV a, UV d)
{
	if (a >= s->n_res || d >= s->n_res) return FALSE;
	return s->ks_acc[d * CSP_KS_KEEP] == a || s->ks_acc[d * CSP_KS_KEEP + 1] == a;
}

//two residues are in the same chain when the walk put them in one
static bool same_chain(const structset *CSP_RESTRICT s, UV a, UV b)
{
	UV c;
	if (a >= s->n_res || b >= s->n_res) return FALSE;
	for (c = 0; c < s->n_chain; c++)
		if (a >= s->chain_first[c] && a < s->chain_last[c])
			return b >= s->chain_first[c] && b < s->chain_last[c];
	return FALSE;
}

//an n-turn at i: the C=O of i reaches the N-H of i+n, both in one chain
static bool n_turn(const structset *CSP_RESTRICT s, UV i, UV n)
{
	return same_chain(s, i, i + n) && ks_bonded(s, i, i + n);
}

#define BR_NONE  0
#define BR_PARA  1 //the two strands run the same way
#define BR_ANTI  2 //they run opposite ways

static unsigned char bridge_kind(const structset *CSP_RESTRICT s, UV i, UV j)
{
	if (i == 0 || j == 0) return BR_NONE;
	if (i + 1 >= s->n_res || j + 1 >= s->n_res) return BR_NONE;
	//DSSP asks for three residues between them, and for one chain
	if (!same_chain(s, i, j)) return BR_NONE;
	if (!(same_chain(s, i - 1, i + 1) && same_chain(s, j - 1, j + 1))) return BR_NONE;
	if ((i > j ? i - j : j - i) < 3) return BR_NONE;
	if ((ks_bonded(s, i - 1, j) && ks_bonded(s, j, i + 1))
	 || (ks_bonded(s, j - 1, i) && ks_bonded(s, i, j + 1)))
		return BR_PARA;
	if ((ks_bonded(s, i, j) && ks_bonded(s, j, i))
	 || (ks_bonded(s, i - 1, j + 1) && ks_bonded(s, j - 1, i + 1)))
		return BR_ANTI;
	return BR_NONE;
}

/*dssp_compute() -- the eight-letter assignment, onto each residue.

`ss' is the letter; `ss_simple' is mdtraj's three-state reduction of it, H for
the three helices, E for the two sheet letters and C for everything else.  A
residue with no backbone to speak of gets neither key rather than a letter
saying it is coil, because it is not coil, it is not protein.*/
static void dssp_compute(pTHX_ structset *CSP_RESTRICT s,
                         const backbone *CSP_RESTRICT bb, bool store)
{
	char *CSP_RESTRICT ss = NULL;
	unsigned char *CSP_RESTRICT bridge = NULL;
	UV i, d;
	unsigned short int k;

	if (s->n_res == 0) return;
	Newx(ss, s->n_res, char);
	Newxz(bridge, s->n_res, unsigned char);
	for (i = 0; i < s->n_res; i++) ss[i] = SS_COIL;

	/*The order matters, and it is not the order of the priority list.  A 3-10 or
	a pi helix is laid down only where the whole run of it is free -- if any one
	residue of the run is already an alpha helix or a sheet, the whole run is
	refused rather than trimmed, which is what the reference implementations do
	and is why 1A42's residues 13 and 14 come out T and not G: a 3-turn pair
	starts at 12, and 12 is already the last residue of an alpha helix.

	So: alpha helices, then sheets, then the two weaker helices where nothing
	stronger sits, then turns, then bends.*/
	for (i = 1; i < s->n_res; i++) {
		UV j;
		if (n_turn(s, i - 1, 4) && n_turn(s, i, 4))
			for (j = i; j < i + 4 && j < s->n_res; j++) ss[j] = SS_HELIX;
	}

	/*Bridges, from the bonds rather than from every pair.  Each of the four
	patterns is a pair of bonds, so reading one bond's two residues back through
	every place it can appear gives all eight pairs a bridge could be built on --
	and every bridge has a bond in it, so no bridge is missed.  That turns a
	search quadratic in the residues into one linear in the bonds.*/
	for (d = 0; d < s->n_res; d++) {
		for (k = 0; k < CSP_KS_KEEP; k++) {
			UV a = s->ks_acc[d * CSP_KS_KEEP + k];
			UV cand[8][2];
			unsigned short int c;
			UV am1 = (a > 0) ? a - 1 : s->n_res;
			UV dm1 = (d > 0) ? d - 1 : s->n_res;
			if (a >= s->n_res) continue;
			(void)am1;
			//parallel, first clause:  Hbond(i-1,j) and Hbond(j,i+1)
			cand[0][0] = a + 1; cand[0][1] = d;
			cand[1][0] = dm1;   cand[1][1] = a;
			//parallel, second clause: Hbond(j-1,i) and Hbond(i,j+1)
			cand[2][0] = d;     cand[2][1] = a + 1;
			cand[3][0] = a;     cand[3][1] = dm1;
			//antiparallel, first clause:  Hbond(i,j) and Hbond(j,i)
			cand[4][0] = a;     cand[4][1] = d;
			cand[5][0] = d;     cand[5][1] = a;
			//antiparallel, second clause: Hbond(i-1,j+1) and Hbond(j-1,i+1)
			cand[6][0] = a + 1; cand[6][1] = dm1;
			cand[7][0] = dm1;   cand[7][1] = a + 1;
			for (c = 0; c < 8; c++) {
				UV bi = cand[c][0], bj = cand[c][1];
				if (bi >= s->n_res || bj >= s->n_res) continue;
				if (bridge_kind(s, bi, bj) == BR_NONE) continue;
				bridge[bi] = 1;
				bridge[bj] = 1;
			}
		}
	}
	/*A bridge next to another bridge is a strand; one on its own is an isolated
	bridge.  Both outrank every letter but the alpha helix.*/
	for (i = 0; i < s->n_res; i++) {
		bool run;
		if (!bridge[i]) continue;
		if (ss[i] == SS_HELIX) continue;
		run = (i > 0 && bridge[i - 1] && same_chain(s, i - 1, i))
		   || (i + 1 < s->n_res && bridge[i + 1] && same_chain(s, i, i + 1));
		ss[i] = run ? SS_STRAND : SS_BRIDGE;
	}

	//the 3-10 helix, and then the pi helix, each only where its whole run is free
	for (i = 1; i < s->n_res; i++) {
		UV j;
		bool free_run = TRUE;
		if (!(n_turn(s, i - 1, 3) && n_turn(s, i, 3))) continue;
		for (j = i; j < i + 3 && j < s->n_res; j++)
			if (ss[j] != SS_COIL) free_run = FALSE;
		if (!free_run) continue;
		for (j = i; j < i + 3 && j < s->n_res; j++) ss[j] = SS_G310;
	}
	for (i = 1; i < s->n_res; i++) {
		UV j;
		bool free_run = TRUE;
		if (!(n_turn(s, i - 1, 5) && n_turn(s, i, 5))) continue;
		for (j = i; j < i + 5 && j < s->n_res; j++)
			if (ss[j] != SS_COIL) free_run = FALSE;
		if (!free_run) continue;
		for (j = i; j < i + 5 && j < s->n_res; j++) ss[j] = SS_PI;
	}

	//turns: any residue spanned by an n-turn that nothing stronger claimed
	for (i = 0; i < s->n_res; i++) {
		for (k = 3; k <= 5; k++) {
			UV j;
			if (!n_turn(s, i, k)) continue;
			for (j = i + 1; j < i + k && j < s->n_res; j++)
				if (ss[j] == SS_COIL) ss[j] = SS_TURN;
		}
	}

	//bends, last: the chain turns through more than 70 degrees between the two
	//CA atoms two either side of this one
	for (i = 2; i + 2 < s->n_res; i++) {
		NV u[3], v[3], lu = 0.0, lv = 0.0, dot = 0.0;
		if (ss[i] != SS_COIL) continue;
		if (!bb[i].whole || !bb[i - 2].whole || !bb[i + 2].whole) continue;
		if (!same_chain(s, i - 2, i + 2)) continue;
		for (k = 0; k < 3; k++) {
			u[k] = bb[i].ca[k]     - bb[i - 2].ca[k];
			v[k] = bb[i + 2].ca[k] - bb[i].ca[k];
			lu += u[k] * u[k];
			lv += v[k] * v[k];
			dot += u[k] * v[k];
		}
		lu = nv_sqrt(lu); lv = nv_sqrt(lv);
		if (!(lu > 0.0) || !(lv > 0.0)) continue;
		dot /= lu * lv;
		if (dot > 1.0) dot = 1.0; else if (dot < -1.0) dot = -1.0;
		if (nv_acos(dot) * 180.0 / CSP_PI > 70.0) ss[i] = SS_BEND;
	}

	if (store) {
		for (i = 0; i < s->n_res; i++) {
			char simple;
			(void)hv_delete(s->res_hv[i], "ss", 2, G_DISCARD);
			(void)hv_delete(s->res_hv[i], "ss_simple", 9, G_DISCARD);
			if (!bb[i].whole) continue; //not protein: it has no secondary structure
			simple = (ss[i] == SS_HELIX || ss[i] == SS_G310 || ss[i] == SS_PI) ? 'H'
			       : (ss[i] == SS_STRAND || ss[i] == SS_BRIDGE) ? 'E' : 'C';
			(void)hv_stores(s->res_hv[i], "ss", newSVpvn(&ss[i], 1));
			(void)hv_stores(s->res_hv[i], "ss_simple", newSVpvn(&simple, 1));
		}
	}
	Safefree(ss);
	Safefree(bridge);
}

//an option that is a number, with the default the caller wrote down
static NV opt_nv(pTHX_ HV *CSP_RESTRICT o, const char *CSP_RESTRICT k, NV dflt)
{
	SV *v = opt_get(aTHX_ o, k);
	return v ? SvNV(v) : dflt;
}

/*set_free() as a scope destructor, so that a croak between set_build() and the
end of features_do() -- which at this point can only be an allocation failing --
does not walk off with the coordinate arrays.

No CSP_RESTRICT on the parameter: the prototype is SAVEDESTRUCTOR_X's, not this
file's, and a qualifier there would say something about a pointer perl chose.*/
static void set_free_cb(pTHX_ void *p)
{
	set_free(aTHX_ (structset *)p);
}

//degrees in the option hash, radians in the comparison
static NV opt_rad(pTHX_ HV *CSP_RESTRICT o, const char *CSP_RESTRICT k, NV dflt_deg)
{
	return opt_nv(aTHX_ o, k, dflt_deg) * CSP_PI / 180.0;
}

/*features_do() -- one walk of the structure, every property that was asked for.

The roll-ups are sums over the contiguous ranges set_build() left behind: a
residue's atoms are one run and a chain's residues are another, so the
per-residue and per-chain areas are the same loop the per-atom areas came out
of rather than a second pass keyed on anything.

What is written back into $info -- an atom's sasa, a residue's sasa and rsa, a
chain's sasa -- is written here rather than returned and grafted on in Perl,
for the reason the whole file exists: a 400,000-atom structure would otherwise
mean 400,000 more SVs and a Perl loop to put them in place.*/
static HV *features_do(pTHX_ HV *CSP_RESTRICT info, HV *CSP_RESTRICT o,
                       const char *CSP_RESTRICT who)
{
	structset s;
	HV *out;
	HV *sasa_hv = NULL;
	AV *pi = NULL;
	AV *ss = NULL;
	AV *bp = NULL;
	AV *bs = NULL;
	AV *cont = NULL;
	AV *hb = NULL;
	const bool want_sasa  = opt_bool(aTHX_ o, "sasa", TRUE);
	const bool want_iface = opt_bool(aTHX_ o, "interface", TRUE);
	const bool want_shape = opt_bool(aTHX_ o, "shape", TRUE);
	const bool want_pi   = opt_bool(aTHX_ o, "pi_stacking", TRUE);
	const bool want_ss   = opt_bool(aTHX_ o, "disulfides", TRUE);
	const bool want_bp   = opt_bool(aTHX_ o, "base_pairs", TRUE);
	const bool want_bs   = opt_bool(aTHX_ o, "base_stacks", TRUE);
	const bool want_tors = opt_bool(aTHX_ o, "dihedrals", TRUE);
	const bool want_cont = opt_bool(aTHX_ o, "contacts", TRUE);
	const bool want_hse  = opt_bool(aTHX_ o, "exposure", TRUE);
	const bool want_hb   = opt_bool(aTHX_ o, "hbonds", TRUE);
	const bool want_ssq  = opt_bool(aTHX_ o, "secondary", TRUE);
	const bool store     = opt_bool(aTHX_ o, "store", TRUE);
	const NV probe = opt_nv(aTHX_ o, "probe", 1.4);
	const IV points = opt_iv(aTHX_ o, "points", 960);
	NV total = 0.0, apolar = 0.0, buried = 0.0, mass_total = 0.0;
	NV gyr[3][3], moment[3];
	NV cx = 0.0, cy = 0.0, cz = 0.0, mx = 0.0, my = 0.0, mz = 0.0;
	NV rg = 0.0, rg_mass = 0.0;
	bool have_rg_mass = FALSE;
	UV i, r, c;

	if (probe < 0.0) croak("%s: probe must not be negative", who);
	//The upper bound is memory, not accuracy: the sphere points are three NVs
	//each, so ten million of them is 240 MB on a double perl and four times
	//that on a quadmath one, for an atom whose area is one number.
	if (points < 1 || points > 10000000)
		croak("%s: points must be between 1 and 10000000", who);

	ENTER;
	Zero(&s, 1, structset);
	SAVEDESTRUCTOR_X(set_free_cb, &s);
	set_build(aTHX_ info, &s, probe, who);

	if (want_sasa) {
		sasa_compute(aTHX_ &s, (UV)points, want_iface);
		for (i = 0; i < s.n_atom; i++) {
			total += s.area[i];
			if (s.apolar[i]) apolar += s.area[i];
			if (s.alone) buried += s.alone[i] - s.area[i];
			if (store) (void)hv_stores(s.atom_hv[i], "sasa", newSVnv(s.area[i]));
		}
		for (c = 0; c < s.n_chain; c++) {
			NV c_area = 0.0, c_alone = 0.0;
			for (r = s.chain_first[c]; r < s.chain_last[c]; r++) {
				NV r_area = 0.0;
				for (i = s.res_first[r]; i < s.res_last[r]; i++) {
					r_area += s.area[i];
					if (s.alone) c_alone += s.alone[i];
				}
				c_area += r_area;
				if (!store) continue;
				(void)hv_stores(s.res_hv[r], "sasa", newSVnv(r_area));
				/*Relative accessibility, and only for an amino acid: the
				single-letter codes of the nucleotides are amino acid codes
				too, and dividing a guanine's area by glycine's maximum would
				be a number rather than an answer.*/
				if (s.res_type[r] == RT_AA && s.res_one[r] >= 'A' && s.res_one[r] <= 'Z') {
					NV maxa = max_asa[s.res_one[r] - 'A'];
					if (maxa > 0.0)
						(void)hv_stores(s.res_hv[r], "rsa", newSVnv(r_area / maxa));
				}
			}
			if (!store) continue;
			(void)hv_stores(s.chain_hv[c], "sasa", newSVnv(c_area));
			if (s.alone) {
				(void)hv_stores(s.chain_hv[c], "sasa_alone", newSVnv(c_alone));
				(void)hv_stores(s.chain_hv[c], "buried", newSVnv(c_alone - c_area));
			}
		}
	}

	/*Size and weight.  The unweighted radius of gyration is mdtraj's
	compute_rg() with its default masses -- the root mean square distance from
	the centroid -- and is comparable with it directly.  The mass-weighted one
	is taken about the centre of mass, which is what the quantity means and is
	not what mdtraj does: compute_rg(traj, masses=m) weights the distances by
	mass but still measures them from the geometric centroid, so the two answers
	differ by however far the two centres are apart.  t/features.t pins this one
	against the corrected expression rather than against that call.*/
	for (i = 0; i < s.n_atom; i++) {
		cx += s.x[i]; cy += s.y[i]; cz += s.z[i];
		mass_total += s.mass[i];
		mx += s.mass[i] * s.x[i]; my += s.mass[i] * s.y[i]; mz += s.mass[i] * s.z[i];
	}
	Zero(gyr, 1, NV[3][3]);
	Zero(moment, 3, NV);
	if (s.n_atom) {
		NV sum = 0.0;
		cx /= (NV)s.n_atom; cy /= (NV)s.n_atom; cz /= (NV)s.n_atom;
		for (i = 0; i < s.n_atom; i++) {
			NV dx = s.x[i] - cx, dy = s.y[i] - cy, dz = s.z[i] - cz;
			sum += dx * dx + dy * dy + dz * dz;
			/*The gyration tensor, in the same pass: it is the outer product of
			the same offsets whose squares the radius of gyration sums, so its
			trace is that sum and rg*rg is the trace over n.  Six entries rather
			than nine, because it is symmetric by construction.*/
			gyr[0][0] += dx * dx; gyr[1][1] += dy * dy; gyr[2][2] += dz * dz;
			gyr[0][1] += dx * dy; gyr[0][2] += dx * dz; gyr[1][2] += dy * dz;
		}
		rg = nv_sqrt(sum / (NV)s.n_atom);
		gyr[0][0] /= (NV)s.n_atom; gyr[1][1] /= (NV)s.n_atom; gyr[2][2] /= (NV)s.n_atom;
		gyr[0][1] /= (NV)s.n_atom; gyr[0][2] /= (NV)s.n_atom; gyr[1][2] /= (NV)s.n_atom;
		gyr[1][0] = gyr[0][1]; gyr[2][0] = gyr[0][2]; gyr[2][1] = gyr[1][2];
		sym3_eigenvalues(gyr, moment);
	}
	if (mass_total > 0.0) {
		NV sum = 0.0;
		mx /= mass_total; my /= mass_total; mz /= mass_total;
		for (i = 0; i < s.n_atom; i++) {
			NV dx = s.x[i] - mx, dy = s.y[i] - my, dz = s.z[i] - mz;
			sum += s.mass[i] * (dx * dx + dy * dy + dz * dz);
		}
		rg_mass = nv_sqrt(sum / mass_total);
		have_rg_mass = TRUE;
	}

	if (want_pi) {
		pi_opt po;
		po.face_dist     = opt_nv(aTHX_ o, "face_distance", 5.5);
		po.face_plane_lo = opt_rad(aTHX_ o, "face_plane_min", 0.0);
		po.face_plane_hi = opt_rad(aTHX_ o, "face_plane_max", 35.0);
		po.face_norm_lo  = opt_rad(aTHX_ o, "face_normal_min", 0.0);
		po.face_norm_hi  = opt_rad(aTHX_ o, "face_normal_max", 33.0);
		po.edge_dist     = opt_nv(aTHX_ o, "edge_distance", 6.5);
		po.edge_plane_lo = opt_rad(aTHX_ o, "edge_plane_min", 50.0);
		po.edge_plane_hi = opt_rad(aTHX_ o, "edge_plane_max", 90.0);
		po.edge_norm_lo  = opt_rad(aTHX_ o, "edge_normal_min", 0.0);
		po.edge_norm_hi  = opt_rad(aTHX_ o, "edge_normal_max", 30.0);
		po.edge_radius   = opt_nv(aTHX_ o, "edge_radius", 1.5);
		pi = pi_stacking(aTHX_ &s, &po);
		sv_2mortal((SV *)pi);
	}

	if (want_tors && store) {
		const NV bond = opt_nv(aTHX_ o, "peptide_bond", CSP_PEPTIDE_BOND);
		const NV phospho = opt_nv(aTHX_ o, "phosphodiester_bond", CSP_PHOSPHO_BOND);
		if (bond <= 0.0) croak("%s: peptide_bond must be a positive number", who);
		if (phospho <= 0.0)
			croak("%s: phosphodiester_bond must be a positive number", who);
		dihedrals(aTHX_ &s, bond, phospho);
	}

	if (want_hb || want_ssq) {
		backbone *CSP_RESTRICT bb = NULL;
		const NV bond = opt_nv(aTHX_ o, "peptide_bond", CSP_PEPTIDE_BOND);
		Newxz(bb, s.n_res ? s.n_res : 1, backbone);
		backbone_read(aTHX_ &s, bb, bond);
		ks_compute(aTHX_ &s, bb);
		if (want_ssq) dssp_compute(aTHX_ &s, bb, store);
		Safefree(bb);
		if (want_hb) {
			hb = hbond_list(aTHX_ &s);
			sv_2mortal((SV *)hb);
		}
	}

	if (want_cont) {
		const NV cut = opt_nv(aTHX_ o, "contact_distance", 4.5);
		if (cut <= 0.0) croak("%s: contact_distance must be a positive number", who);
		cont = contacts_find(aTHX_ &s, cut, store);
		sv_2mortal((SV *)cont);
	}
	if (want_hse && store) hse_compute(aTHX_ &s, CSP_HSE_RADIUS);

	if (want_ss) {
		const NV cut = opt_nv(aTHX_ o, "disulfide_distance", 3.0);
		if (cut < 0.0) croak("%s: disulfide_distance must not be negative", who);
		ss = disulfides(aTHX_ &s, cut, store);
		sv_2mortal((SV *)ss);
	}

	if (want_bp) {
		const NV hb  = opt_nv(aTHX_ o, "base_pair_hbond", CSP_BP_HBOND);
		const NV sag = opt_nv(aTHX_ o, "base_pair_stagger", CSP_BP_STAGGER);
		if (hb <= 0.0) croak("%s: base_pair_hbond must be a positive number", who);
		if (sag < 0.0) croak("%s: base_pair_stagger must not be negative", who);
		bp = base_pairs(aTHX_ &s, hb, sag, store);
		sv_2mortal((SV *)bp);
	}

	if (want_bs) {
		const NV dc = opt_nv(aTHX_ o, "base_stack_distance", CSP_STACK_DIST);
		const NV oc = opt_nv(aTHX_ o, "base_stack_omega", CSP_STACK_OMEGA);
		if (dc <= 0.0) croak("%s: base_stack_distance must be a positive number", who);
		if (oc < 0.0 || oc > 180.0)
			croak("%s: base_stack_omega must be between 0 and 180", who);
		bs = base_stacks(aTHX_ &s, dc, oc, store);
		sv_2mortal((SV *)bs);
	}

	//the result hash is built last, so that nothing between here and the return
	//can croak with it half-filled and unreferenced
	out = newHV();
	(void)hv_stores(out, "n_atoms",      newSVuv(s.n_atom));
	(void)hv_stores(out, "n_residues",   newSVuv(s.n_res));
	(void)hv_stores(out, "n_chains",     newSVuv(s.n_chain));
	(void)hv_stores(out, "n_no_element", newSVuv(s.n_no_element));
	(void)hv_stores(out, "mass",         newSVnv(mass_total));
	if (s.n_atom) {
		AV *ctr = newAV();
		av_push(ctr, newSVnv(cx)); av_push(ctr, newSVnv(cy)); av_push(ctr, newSVnv(cz));
		(void)hv_stores(out, "center", newRV_noinc((SV *)ctr));
		(void)hv_stores(out, "rg", newSVnv(rg));
	}
	if (have_rg_mass) {
		AV *com = newAV();
		av_push(com, newSVnv(mx)); av_push(com, newSVnv(my)); av_push(com, newSVnv(mz));
		(void)hv_stores(out, "center_of_mass", newRV_noinc((SV *)com));
		(void)hv_stores(out, "rg_mass", newSVnv(rg_mass));
	}
	if (want_sasa) {
		sasa_hv = newHV();
		(void)hv_stores(sasa_hv, "total",  newSVnv(total));
		(void)hv_stores(sasa_hv, "apolar", newSVnv(apolar));
		(void)hv_stores(sasa_hv, "polar",  newSVnv(total - apolar));
		(void)hv_stores(sasa_hv, "probe",  newSVnv(probe));
		(void)hv_stores(sasa_hv, "points", newSVuv((UV)points));
		/*What the chains bury against each other: the surface they have apart
		less the surface they have together.  A two-body interface is usually
		quoted as half of this, because the number is counted once on each side
		of it, and a caller with more than two chains has the per-chain figures
		to add up however the question was asked.*/
		if (s.alone) (void)hv_stores(sasa_hv, "buried", newSVnv(buried));
		(void)hv_stores(out, "sasa", newRV_noinc((SV *)sasa_hv));
	}
	if (want_shape && s.n_atom) {
		/*The gyration tensor's eigenvalues, and the three numbers mdtraj builds
		out of them in mdtraj/geometry/shape.py: asphericity is how far the
		largest moment stands above the mean of the other two, acylindricity how
		far the middle one stands above the smallest, and the relative shape
		anisotropy is 0 for a sphere and 1 for a straight line.  All in
		angstrom^2, where mdtraj's are nm^2; the anisotropy has no units.*/
		HV *sh = newHV();
		AV *pm = newAV();
		AV *tensor = newAV();
		NV tr = moment[0] + moment[1] + moment[2];
		unsigned short int a, b;
		for (a = 0; a < 3; a++) av_push(pm, newSVnv(moment[a]));
		for (a = 0; a < 3; a++) {
			AV *row = newAV();
			for (b = 0; b < 3; b++) av_push(row, newSVnv(gyr[a][b]));
			av_push(tensor, newRV_noinc((SV *)row));
		}
		(void)hv_stores(sh, "gyration_tensor",  newRV_noinc((SV *)tensor));
		(void)hv_stores(sh, "principal_moments", newRV_noinc((SV *)pm));
		(void)hv_stores(sh, "asphericity",   newSVnv(moment[2] - (moment[0] + moment[1]) / 2.0));
		(void)hv_stores(sh, "acylindricity", newSVnv(moment[1] - moment[0]));
		if (tr > 0.0)
			(void)hv_stores(sh, "anisotropy", newSVnv(
				1.5 * (moment[0] * moment[0] + moment[1] * moment[1] + moment[2] * moment[2])
				/ (tr * tr) - 0.5));
		(void)hv_stores(out, "shape", newRV_noinc((SV *)sh));
	}
	if (want_pi) (void)hv_stores(out, "pi_stacking", newRV_inc((SV *)pi));
	if (want_ss) (void)hv_stores(out, "disulfides", newRV_inc((SV *)ss));
	if (want_bp) (void)hv_stores(out, "base_pairs", newRV_inc((SV *)bp));
	if (want_bs) (void)hv_stores(out, "base_stacks", newRV_inc((SV *)bs));
	if (want_cont) (void)hv_stores(out, "contacts", newRV_inc((SV *)cont));
	if (want_hb) (void)hv_stores(out, "hbonds", newRV_inc((SV *)hb));
	LEAVE;
	return out;
}

MODULE = Chem::Structure::Parser		PACKAGE = Chem::Structure::Parser

PROTOTYPES: DISABLE

SV * _parse_file(path, opts = &PL_sv_undef)
	SV *path
	SV *opts
	PREINIT:
		char *CSP_RESTRICT buf;
		STRLEN len;
		HV *o = NULL, *res;
		const char *CSP_RESTRICT p;
	CODE:
		if (!SvOK(path)) croak("Chem::Structure::Parser: file name is undefined");
		if (SvOK(opts)) {
			if (!SvROK(opts) || SvTYPE(SvRV(opts)) != SVt_PVHV)
				croak("Chem::Structure::Parser: options must be a hash reference");
			o = (HV *)SvRV(opts);
		}
		p = SvPV_nolen(path);
		buf = slurp(aTHX_ p, &len);
		res = parse_buf(aTHX_ buf, len, o);
		Safefree(buf);
		(void)hv_stores(res, "file", newSVsv(path));
		RETVAL = newRV_noinc((SV *)res);
	OUTPUT:
		RETVAL

SV * _parse_string(text, opts = &PL_sv_undef)
	SV *text
	SV *opts
	PREINIT:
		STRLEN len;
		const char *CSP_RESTRICT buf;
		HV *o = NULL, *res;
	CODE:
		if (!SvOK(text)) croak("Chem::Structure::Parser: PDB text is undefined");
		if (SvOK(opts)) {
			if (!SvROK(opts) || SvTYPE(SvRV(opts)) != SVt_PVHV)
				croak("Chem::Structure::Parser: options must be a hash reference");
			o = (HV *)SvRV(opts);
		}
		buf = SvPV_const(text, len);
		res = parse_buf(aTHX_ buf, len, o);
		RETVAL = newRV_noinc((SV *)res);
	OUTPUT:
		RETVAL

SV *
_parse_cif_file(path, opts = &PL_sv_undef)
	SV *path
	SV *opts
	PREINIT:
		char *CSP_RESTRICT buf;
		STRLEN len;
		HV *o = NULL, *res;
		const char *CSP_RESTRICT p;
	CODE:
		if (!SvOK(path)) croak("Chem::Structure::Parser: file name is undefined");
		if (SvOK(opts)) {
			if (!SvROK(opts) || SvTYPE(SvRV(opts)) != SVt_PVHV)
				croak("Chem::Structure::Parser: options must be a hash reference");
			o = (HV *)SvRV(opts);
		}
		p = SvPV_nolen(path);
		buf = slurp(aTHX_ p, &len);
		res = parse_cif_buf(aTHX_ buf, len, o);
		Safefree(buf);
		(void)hv_stores(res, "file", newSVsv(path));
		RETVAL = newRV_noinc((SV *)res);
	OUTPUT:
		RETVAL

SV *
_parse_cif_string(text, opts = &PL_sv_undef)
	SV *text
	SV *opts
	PREINIT:
		STRLEN len;
		const char *CSP_RESTRICT buf;
		HV *o = NULL, *res;
	CODE:
		if (!SvOK(text)) croak("Chem::Structure::Parser: mmCIF text is undefined");
		if (SvOK(opts)) {
			if (!SvROK(opts) || SvTYPE(SvRV(opts)) != SVt_PVHV)
				croak("Chem::Structure::Parser: options must be a hash reference");
			o = (HV *)SvRV(opts);
		}
		buf = SvPV_const(text, len);
		res = parse_cif_buf(aTHX_ buf, len, o);
		RETVAL = newRV_noinc((SV *)res);
	OUTPUT:
		RETVAL

void
_str2nv_paths(text)
	SV *text
	PREINIT:
		STRLEN n;
		const char *CSP_RESTRICT s;
		NV a = 0, b = 0;
		bool ok_a, ok_b;
	PPCODE:
		if (!SvOK(text)) croak("_str2nv_paths: string is undefined");
		s = SvPV_const(text, n);
		ok_a = str2nv_fixed(s, n, &a);
		ok_b = str2nv_slow(s, n, &b);
		EXTEND(SP, 4);
		PUSHs(sv_2mortal(newSViv(ok_a)));
		PUSHs(ok_a ? sv_2mortal(newSVnv(a)) : &PL_sv_undef);
		PUSHs(sv_2mortal(newSViv(ok_b)));
		PUSHs(ok_b ? sv_2mortal(newSVnv(b)) : &PL_sv_undef);

SV *
aa3to1(name)
	SV *name
	PREINIT:
		STRLEN n;
		const char *s;
		res_info ri;
	CODE:
		if (!SvOK(name)) croak("aa3to1: residue name is undefined");
		s = SvPV_const(name, n);
		RETVAL = (res_lookup(s, n, &ri) && ri.type == RT_AA)
		       ? newSVpvn(&ri.one, 1)
		       : newSVpvn("", 0);
	OUTPUT:
		RETVAL

SV *
aa1to3(one)
	SV *one
	PREINIT:
		STRLEN n;
		const char *s;
		const char *name;
	CODE:
		if (!SvOK(one)) croak("aa1to3: single-letter code is undefined");
		s = SvPV_const(one, n);
		name = aa1to3_lookup(s, n);
		RETVAL = name ? newSVpvn(name, 3) : newSVpvn("", 0);
	OUTPUT:
		RETVAL

SV *
res1(name)
	SV *name
	PREINIT:
		STRLEN n;
		const char *s;
		res_info ri;
	CODE:
		if (!SvOK(name)) croak("res1: residue name is undefined");
		s = SvPV_const(name, n);
		RETVAL = (res_lookup(s, n, &ri) && ri.one)
		       ? newSVpvn(&ri.one, 1)
		       : newSVpvn("", 0);
	OUTPUT:
		RETVAL

SV *
res_type(name)
	SV *name
	PREINIT:
		STRLEN n;
		const char *CSP_RESTRICT s;
		res_info ri;
		const char *CSP_RESTRICT t;
	CODE:
		if (!SvOK(name)) croak("res_type: residue name is undefined");
		s = SvPV_const(name, n);
		if (!res_lookup(s, n, &ri)) t = "other";
		else if (ri.type == RT_AA)    t = "amino_acid";
		else if (ri.type == RT_NUC)   t = "nucleotide";
		else if (ri.type == RT_WATER) t = "water";
		else                          t = "other";
		RETVAL = newSVpv(t, 0);
	OUTPUT:
		RETVAL

bool
is_single_ion(chain, id = &PL_sv_undef)
	SV *chain
	SV *id
	CODE:
		//two arguments mean the structure and a chain: an undefined id there is
		//a chain that was not given, not a chain hash passed on its own
		if (items > 1 && !SvOK(id)) croak("is_single_ion: no chain given");
		RETVAL = chain_is_single_ion(aTHX_ chain_arg(aTHX_ chain, id, "is_single_ion"));
	OUTPUT:
		RETVAL

SV *
_features(info, opts, who)
	SV *info
	SV *opts
	SV *who
	PREINIT:
		HV *o = NULL;
	CODE:
		if (!SvROK(info) || SvTYPE(SvRV(info)) != SVt_PVHV)
			croak("Chem::Structure::Parser: structure must be a hash reference");
		if (SvOK(opts)) {
			if (!SvROK(opts) || SvTYPE(SvRV(opts)) != SVt_PVHV)
				croak("Chem::Structure::Parser: options must be a hash reference");
			o = (HV *)SvRV(opts);
		}
		RETVAL = newRV_noinc((SV *)features_do(aTHX_ (HV *)SvRV(info), o,
		                                       SvOK(who) ? SvPV_nolen(who) : "structure_features"));
	OUTPUT:
		RETVAL
