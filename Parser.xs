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

This layer is the same one Stats::LikeR carries, cut down to the five functions
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
static unsigned short int ring_defs(U32 key, const ring_def *CSP_RESTRICT *out)
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
	NV *mass;     //dalton; 0.0 for an atom whose element has no mass
	unsigned char *apolar; //1 = carbon or sulphur, 0 = everything else
	HV **atom_hv; //the atom hash the row was read from, for writing back
	UV n_atom;
	//per residue: the atom range [first, last)
	HV **res_hv;
	UV *res_first, *res_last;
	char *res_one; //single-letter code, '\0' when the residue has none
	unsigned char *res_type; //RT_*, for the questions only an amino acid answers
	U32 *res_key;  //packed residue name, for ring_defs()
	UV n_res;
	//per chain: the residue range [first, last)
	HV **chain_hv;
	UV *chain_first, *chain_last;
	UV n_chain;
	//tallies gathered on the way past
	NV mass_total;
	UV n_no_element; //atoms whose element field spells no element
} structset;

static void set_free(pTHX_ structset *CSP_RESTRICT s)
{
	Safefree(s->x);        Safefree(s->y);         Safefree(s->z);
	Safefree(s->rad);      Safefree(s->area);      Safefree(s->mass);
	Safefree(s->apolar);   Safefree(s->atom_hv);
	Safefree(s->res_hv);   Safefree(s->res_first); Safefree(s->res_last);
	Safefree(s->res_one);  Safefree(s->res_type);  Safefree(s->res_key);
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
			ao = hvf_av(aTHX_ r, "atom_order", 10);
			atoms = hvf_hv(aTHX_ r, "atoms", 5);
			nai = (ao && atoms) ? av_len(ao) + 1 : 0;
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
	if (s->n_atom == 0 && cap_atom > 0) {
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
static void sasa_compute(pTHX_ structset *CSP_RESTRICT s, UV npts)
{
	const UV n = s->n_atom;
	NV *CSP_RESTRICT pts = NULL;
	UV *CSP_RESTRICT nbr = NULL;
	//grown by doubling as needed; 64 covers every atom of a protein, where the
	//neighbour count runs to the low tens, without a single reallocation
	UV nbr_cap = 64;
	cell_grid g;
	NV rmax = 0.0, constant;
	UV i;

	Newxz(s->area, n ? n : 1, NV);
	if (n == 0) return;
	for (i = 0; i < n; i++) if (s->rad[i] > rmax) rmax = s->rad[i];
	Newx(pts, npts * 3, NV);
	sasa_sphere(pts, npts);
	grid_build(aTHX_ &g, s->x, s->y, s->z, n, 2.0 * rmax);
	Newx(nbr, nbr_cap, UV);
	constant = 4.0 * CSP_PI / (NV)npts;

	for (i = 0; i < n; i++) {
		const NV xi = s->x[i], yi = s->y[i], zi = s->z[i], ri = s->rad[i];
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
				dx = s->x[a] - xi; dy = s->y[a] - yi; dz = s->z[a] - zi;
				sum = ri + s->rad[a];
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
				dx = px - s->x[a]; dy = py - s->y[a]; dz = pz - s->z[a];
				ra = s->rad[a];
				if (dx * dx + dy * dy + dz * dz < ra * ra) {
					k_closest = kp;
					open = FALSE;
					break;
				}
			}
			if (open) acc++;
		}
		s->area[i] = (NV)acc * constant * ri * ri;
	}
	Safefree(nbr);
	Safefree(pts);
	grid_free(aTHX_ &g);
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
//atom of that name (a side chain modelled only as far as CB, most often)
static bool ring_atom(pTHX_ HV *CSP_RESTRICT atoms, const char *CSP_RESTRICT name,
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

static UV rings_find(pTHX_ structset *CSP_RESTRICT s, ring_t *CSP_RESTRICT *out)
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
				if (!ring_atom(aTHX_ atoms, defs[d].atom[i], &px[i], &py[i], &pz[i])) {
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

//an option that is a number, with the default the caller wrote down
static NV opt_nv(pTHX_ HV *CSP_RESTRICT o, const char *CSP_RESTRICT k, NV dflt)
{
	SV *v = opt_get(aTHX_ o, k);
	return v ? SvNV(v) : dflt;
}

//set_free() as a scope destructor, so that a croak between set_build() and the
//end of features_do() -- which at this point can only be an allocation failing
//-- does not walk off with the coordinate arrays
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
	const bool want_sasa = opt_bool(aTHX_ o, "sasa", TRUE);
	const bool want_pi   = opt_bool(aTHX_ o, "pi_stacking", TRUE);
	const bool store     = opt_bool(aTHX_ o, "store", TRUE);
	const NV probe = opt_nv(aTHX_ o, "probe", 1.4);
	const IV points = opt_iv(aTHX_ o, "points", 960);
	NV total = 0.0, apolar = 0.0, mass_total = 0.0;
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
		sasa_compute(aTHX_ &s, (UV)points);
		for (i = 0; i < s.n_atom; i++) {
			total += s.area[i];
			if (s.apolar[i]) apolar += s.area[i];
			if (store) (void)hv_stores(s.atom_hv[i], "sasa", newSVnv(s.area[i]));
		}
		for (c = 0; c < s.n_chain; c++) {
			NV c_area = 0.0;
			for (r = s.chain_first[c]; r < s.chain_last[c]; r++) {
				NV r_area = 0.0;
				for (i = s.res_first[r]; i < s.res_last[r]; i++) r_area += s.area[i];
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
			if (store) (void)hv_stores(s.chain_hv[c], "sasa", newSVnv(c_area));
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
	if (s.n_atom) {
		NV sum = 0.0;
		cx /= (NV)s.n_atom; cy /= (NV)s.n_atom; cz /= (NV)s.n_atom;
		for (i = 0; i < s.n_atom; i++) {
			NV dx = s.x[i] - cx, dy = s.y[i] - cy, dz = s.z[i] - cz;
			sum += dx * dx + dy * dy + dz * dz;
		}
		rg = nv_sqrt(sum / (NV)s.n_atom);
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
		(void)hv_stores(out, "sasa", newRV_noinc((SV *)sasa_hv));
	}
	if (want_pi) (void)hv_stores(out, "pi_stacking", newRV_inc((SV *)pi));
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
