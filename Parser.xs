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

/*residue table.  One table serves aa3to1(), res1() and res_type(), so the
three can never disagree about what a residue is.  Keyed on the three-byte
name packed into an int, left-padded with blanks, so that "ALA", " DA" and
"  A" are all one switch on an integer.*/
#define RT_OTHER 0
#define RT_AA    1
#define RT_NUC   2
#define RT_WATER 3

#define K3(a, b, c) (((unsigned)(a) << 16) | ((unsigned)(b) << 8) | (unsigned)(c))

typedef struct {
	char one;           //single-letter code; '\0' when there isn't one
	unsigned char type; //RT_*
} res_info;

//pack a (possibly blank-padded, possibly short) residue name into a key
static unsigned res_key(const char *CSP_RESTRICT s, STRLEN len)
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

/*res_lookup() -- classify a residue name.  1 when the name is known, 0 when
it is not (a ligand, an ion, a sugar: anything the caller has to work out for
itself from the atoms).*/
static int res_lookup(const char *CSP_RESTRICT s, STRLEN len, res_info *CSP_RESTRICT out)
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

		default: return 0;
	}
	out->one  = one;
	out->type = type;
	return 1;
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
static int str2iv(const char *CSP_RESTRICT s, STRLEN n, IV *CSP_RESTRICT out)
{
	int neg = 0, seen = 0;
	IV v = 0;
	STRLEN i = 0;
	if (i < n && (s[i] == '+' || s[i] == '-')) { neg = (s[i] == '-'); i++; }
	for (; i < n; i++) {
		if (!isdigit((unsigned char)s[i])) return 0;
		v = v * 10 + (s[i] - '0');
		seen = 1;
	}
	if (!seen) return 0;
	*out = neg ? -v : v;
	return 1;
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
static int str2nv_slow(const char *CSP_RESTRICT s, STRLEN n, NV *CSP_RESTRICT out)
{
	char buf[64], *end;
	NV v;
	if (n == 0 || n >= sizeof(buf)) return 0;
	memcpy(buf, s, n);
	buf[n] = '\0';
	errno = 0;
	v = STR2NV_STRTOD(buf, &end);
	if (end == buf) return 0;
	*out = v;
	return 1;
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

static int str2nv_fixed(const char *CSP_RESTRICT s, STRLEN n, NV *CSP_RESTRICT out)
{
	NV m = 0;
	STRLEN i = 0, nd = 0, frac = 0;
	int neg = 0, dot = 0;

	if (n == 0) return 0;
	if (s[0] == '+' || s[0] == '-') { neg = (s[0] == '-'); i = 1; }
	for (; i < n; i++) {
		const char c = s[i];
		if (c >= '0' && c <= '9') {
			if (++nd > STR2NV_MAX_DIGITS) return 0;
			m = m * 10 + (NV)(c - '0');
			if (dot) frac++;
		} else if (c == '.' && !dot) {
			dot = 1;
		} else {
			return 0; //an exponent, a stray letter, the '*****' of an overflowed field
		}
	}
	if (nd == 0) return 0; //"", "-", ".", "+."
	if (frac) m /= str2nv_pow10[frac];
	*out = neg ? -m : m;
	return 1;
}

//str2nv() -- coordinates and the like
static int str2nv(const char *CSP_RESTRICT s, STRLEN n, NV *CSP_RESTRICT out)
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

static int fld_iv(const char *CSP_RESTRICT line, STRLEN llen, STRLEN from, STRLEN to,
                  IV *CSP_RESTRICT out)
{
	const char *s;
	STRLEN n = fld(line, llen, from, to, &s);
	return n ? str2iv(s, n, out) : 0;
}

static int fld_nv(const char *CSP_RESTRICT line, STRLEN llen, STRLEN from, STRLEN to,
                  NV *CSP_RESTRICT out)
{
	const char *s;
	STRLEN n = fld(line, llen, from, to, &s);
	return n ? str2nv(s, n, out) : 0;
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

static int opt_bool(pTHX_ HV *CSP_RESTRICT o, const char *CSP_RESTRICT k, int dflt)
{
	SV *v = opt_get(aTHX_ o, k);
	return v ? (SvTRUE(v) ? 1 : 0) : dflt;
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

static void flush_residue(pTHX_ AV **CSP_RESTRICT rs, NV sx, NV sy, NV sz, IV nc, NV sb, IV nb)
{
	av_push(rs[R_SX],   nc ? newSVnv(sx) : newSVsv(&PL_sv_undef));
	av_push(rs[R_SY],   nc ? newSVnv(sy) : newSVsv(&PL_sv_undef));
	av_push(rs[R_SZ],   nc ? newSVnv(sz) : newSVsv(&PL_sv_undef));
	av_push(rs[R_NXYZ], newSViv(nc));
	av_push(rs[R_SB],   nb ? newSVnv(sb) : newSVsv(&PL_sv_undef));
	av_push(rs[R_NB],   newSViv(nb));
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
	AV *atom_hv = newAV();
	AV *col[NCOL];
	AV *res_sum[NRSUM];
	AV *res_first = newAV(), *res_last = newAV();
	AV *ter = newAV(), *model_nums = newAV();
	HV *want_chain = NULL;
	//whole-structure statistics over the atoms that were kept
	IV n_hydrogen = 0, n_water_atom = 0, bn = 0;
	NV bmin = 0, bmax = 0, bsum = 0;
	NV xmin = 0, ymin = 0, zmin = 0, xmax = 0, ymax = 0, zmax = 0;
	int have_bbox = 0;
	//the residue being accumulated
	NV rsx = 0, rsy = 0, rsz = 0, rsb = 0;
	IV rnc = 0, rnb = 0;
	IV want_model, n_models = 0, n_anisou = 0, n_skipped = 0, n_atom = 0;
	IV n_atom_rec = 0, n_het_rec = 0, cur_model = 1, lineno = 0;
	int keep_h, keep_water, keep_het, keep_meta, keep_anisou, keep_lineno, build_atoms;
	int i, have_res = 0;
	//previous kept atom's residue identity, for boundary detection
	char p_chain[8], p_icode[4], p_resname[8];
	IV p_resseq = 0, p_model = 0;
	int p_het = 0;
	STRLEN pos = 0;

	/*a negative model number means every model.  Zero cannot be the sentinel:
	an ensemble whose models are numbered from 0 is unusual but legal, and
	"model 0" would then quietly mean "all of them".*/
	want_model  = opt_iv(aTHX_ opts, "model", 1);
	keep_h      = opt_bool(aTHX_ opts, "hydrogens", 1);
	keep_water  = opt_bool(aTHX_ opts, "waters", 1);
	keep_het    = opt_bool(aTHX_ opts, "hetatm", 1);
	keep_meta   = opt_bool(aTHX_ opts, "meta", 1);
	keep_anisou = opt_bool(aTHX_ opts, "anisou", 0);
	keep_lineno = opt_bool(aTHX_ opts, "lineno", 0);
	build_atoms = opt_bool(aTHX_ opts, "atom_hashes", 0);
	{
		SV *c = opt_get(aTHX_ opts, "chains");
		if (c && SvROK(c) && SvTYPE(SvRV(c)) == SVt_PVHV) want_chain = (HV *)SvRV(c);
	}
	p_chain[0] = p_icode[0] = p_resname[0] = '\0';

	for (i = 0; i < NCOL; i++)  col[i] = newAV();
	for (i = 0; i < NRSUM; i++) res_sum[i] = newAV();

	while (pos < len) {
		const char *line = buf + pos;
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
			int het = (line[0] == 'H');
			res_info ri;
			int known;
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
			} else {
				ellen = guess_element(nm_raw, nm_rawlen, elbuf);
			}
			if (!keep_h && ellen == 1 && (elbuf[0] == 'H' || elbuf[0] == 'D')) { n_skipped++; continue; }

			/*kept.  The residue boundary is settled first, before this atom has
			been added to anything, so that closing the previous residue is a
			matter of flushing what is already there.  Everything the test needs
			is read from the record's own columns.*/
			{
				IV rs = 0;
				int changed;
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
						av_push(res_last, newSViv(n_atom - 1));
						flush_residue(aTHX_ res_sum, rsx, rsy, rsz, rnc, rsb, rnb);
					}
					rsx = rsy = rsz = rsb = 0; rnc = rnb = 0;
					av_push(res_first, newSViv(n_atom));
					have_res = 1;
					p_resseq = rs;
					p_model  = cur_model;
					p_het    = het;
					my_strlcpy(p_chain, chain, sizeof(p_chain));
					my_strlcpy(p_icode, icode, sizeof(p_icode));
					my_strlcpy(p_resname, resname, sizeof(p_resname));
				}
			}

			//the atom's own fields
			{
				const char *nm_s, *alt_s, *chg_s;
				STRLEN nm_n, alt_n, chg_n;
				IV serial = 0;
				int have_serial, have_occ;
				NV xv, yv, zv, bv, ov;
				int have_xyz = fld_nv(line, llen, 30, 37, &xv)
				             & fld_nv(line, llen, 38, 45, &yv)
				             & fld_nv(line, llen, 46, 53, &zv);
				int have_b = fld_nv(line, llen, 60, 65, &bv);

				have_occ    = fld_nv(line, llen, 54, 59, &ov);
				have_serial = fld_iv(line, llen, 6, 10, &serial);
				nm_n  = fld(line, llen, 12, 15, &nm_s);
				alt_n = fld(line, llen, 16, 16, &alt_s);
				chg_n = fld(line, llen, 78, 79, &chg_s);

				/*The sums and extremes are gathered here rather than in Perl
				because they have to touch every atom whether or not the caller
				wanted per-atom hashes.  Doing them in the loop that is already
				reading the numbers costs nothing; doing them again in Perl
				costs more than the whole parse.*/
				if (have_xyz) {
					rsx += xv; rsy += yv; rsz += zv; rnc++;
					if (!have_bbox) {
						xmin = xmax = xv; ymin = ymax = yv; zmin = zmax = zv;
						have_bbox = 1;
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
					SV **cnt = hv_fetch(elements, elbuf, (I32)ellen, 1);
					if (cnt && *cnt) sv_inc(*cnt);
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
			if (keep_lineno) av_push(col[C_LINENO], newSViv(lineno));

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
			cur_model = fld_iv(line, llen, 10, 13, &m) ? m : n_models;
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
		av_push(res_last, newSViv(n_atom - 1));
		flush_residue(aTHX_ res_sum, rsx, rsy, rsz, rnc, rsb, rnb);
	}

	(void)hv_stores(out, "atoms", newRV_noinc((SV *)atom_hv));
	for (i = 0; i < NRSUM; i++)
		(void)hv_store(out, res_sum_name[i], (I32)strlen(res_sum_name[i]),
		               newRV_noinc((SV *)res_sum[i]), 0);
	(void)hv_stores(out, "elements",     newRV_noinc((SV *)elements));
	(void)hv_stores(out, "n_hydrogens",  newSViv(n_hydrogen));
	(void)hv_stores(out, "n_water_atoms", newSViv(n_water_atom));
	if (bn) {
		HV *b = newHV();
		(void)hv_stores(b, "min",  newSVnv(bmin));
		(void)hv_stores(b, "max",  newSVnv(bmax));
		(void)hv_stores(b, "mean", newSVnv(bsum / (NV)bn));
		(void)hv_stores(b, "n",    newSViv(bn));
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

	for (i = 0; i < NCOL; i++) {
		if (i == C_LINENO && !keep_lineno) { SvREFCNT_dec((SV *)col[i]); continue; }
		(void)hv_store(out, col_name[i], (I32)strlen(col_name[i]), newRV_noinc((SV *)col[i]), 0);
	}
	(void)hv_stores(out, "res_first",    newRV_noinc((SV *)res_first));
	(void)hv_stores(out, "res_last",     newRV_noinc((SV *)res_last));
	(void)hv_stores(out, "ter",          newRV_noinc((SV *)ter));
	(void)hv_stores(out, "meta",         newRV_noinc((SV *)meta));
	(void)hv_stores(out, "model_numbers", newRV_noinc((SV *)model_nums));
	(void)hv_stores(out, "n_atoms",      newSViv(n_atom));
	(void)hv_stores(out, "n_residues",   newSViv(av_len(res_first) + 1));
	(void)hv_stores(out, "n_models",     newSViv(n_models ? n_models : 1));
	(void)hv_stores(out, "n_anisou",     newSViv(n_anisou));
	(void)hv_stores(out, "n_skipped",    newSViv(n_skipped));
	(void)hv_stores(out, "n_atom_records",   newSViv(n_atom_rec));
	(void)hv_stores(out, "n_hetatm_records", newSViv(n_het_rec));
	(void)hv_stores(out, "n_lines",      newSViv(lineno));
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

#define CT_EOF   0
#define CT_VALUE 1
#define CT_TAG   2
#define CT_LOOP  3
#define CT_DATA  4
#define CT_SAVE  5
#define CT_STOP  6

typedef struct {
	const char *buf;
	const char *end;
	const char *p;
} cif_lex;

typedef struct {
	const char *s;
	STRLEN n;
	int kind;
	/*a quoted value is text and nothing else: '.' in quotes is a full stop, and
	only a bare '.' is the null the format means*/
	int quoted;
} cif_tok;

static int cif_iskw(const char *CSP_RESTRICT s, STRLEN n, const char *CSP_RESTRICT kw,
                    STRLEN kwn)
{
	STRLEN i;
	if (n < kwn) return 0;
	for (i = 0; i < kwn; i++)
		if (tolower((unsigned char)s[i]) != kw[i]) return 0;
	return 1;
}

static int cif_space(char c)
{
	return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}

//a bare '.' (not applicable) or '?' (unknown): the field was not filled in
static int cif_null(const cif_tok *CSP_RESTRICT t)
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
	t->quoted = 0;
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
			const char *nl = (const char *)memchr(q, '\n', (STRLEN)(end - q));
			if (!nl) { q = end; lx->p = end; break; }
			if (nl + 1 < end && nl[1] == ';') { q = nl; lx->p = nl + 2; break; }
			if (nl + 1 >= end) { q = nl; lx->p = end; break; }
			q = nl + 1;
		}
		t->s = s;
		t->n = (STRLEN)(q - s);
		t->kind = CT_VALUE;
		t->quoted = 1;
		return;
	}

	if (*p == '\'' || *p == '"') {
		/*the closing quote is one followed by whitespace or the end of the file.
		Anything else is a quote inside the value, which is how a CIF writes
		O5' without escaping it.*/
		char qc = *p;
		const char *s = p + 1, *q = s;
		for (;;) {
			while (q < end && *q != qc) q++;
			if (q >= end) break;
			if (q + 1 >= end || cif_space(q[1])) break;
			q++;
		}
		t->s = s;
		t->n = (STRLEN)(q - s);
		t->kind = CT_VALUE;
		t->quoted = 1;
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

static int cif_atom_field(const char *CSP_RESTRICT item, STRLEN n)
{
	int i;
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
	HV *elements, *want_chain;
	IV n_hydrogen, n_water_atom, bn;
	NV bmin, bmax, bsum;
	NV xmin, ymin, zmin, xmax, ymax, zmax;
	int have_bbox;
	NV rsx, rsy, rsz, rsb;
	IV rnc, rnb;
	IV want_model, n_models, n_anisou, n_skipped, n_atom, n_atom_rec, n_het_rec;
	IV cur_model;
	int keep_h, keep_water, keep_het, keep_anisou, build_atoms;
	int have_res, seen_model;
	char p_chain[8], p_icode[4], p_resname[8];
	IV p_resseq, p_model;
	int p_het;
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
	const char *nm_s, *alt_s;
	STRLEN nm_n, alt_n;
	res_info ri;
	int known, het, changed;
	IV rs = 0, serial = 0, model = 1;
	int have_rs, have_serial, have_occ, have_xyz, have_b;
	NV xv = 0, yv = 0, zv = 0, ov = 0, bv = 0;

	//which model this row belongs to, counted before anything can skip the row
	if (v[A_MODEL] && str2iv(v[A_MODEL], vn[A_MODEL], &model)) {
		if (!st->seen_model || model != st->cur_model) {
			st->cur_model = model;
			st->seen_model = 1;
			st->n_models++;
			av_push(st->model_nums, newSViv(model));
		}
	} else if (!st->seen_model) {
		st->cur_model = 1;
		st->seen_model = 1;
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
	if (!st->keep_h && ellen == 1 && (elbuf[0] == 'H' || elbuf[0] == 'D')) {
		st->n_skipped++;
		return;
	}

	//kept.  auth_seq_id is the residue number the PDB record carried; label_seq_id
	//is the position in the entity, and is null for everything that is not polymer
	have_rs = 0;
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
			av_push(st->res_last, newSViv(st->n_atom - 1));
			flush_residue(aTHX_ st->res_sum, st->rsx, st->rsy, st->rsz, st->rnc, st->rsb, st->rnb);
		}
		st->rsx = st->rsy = st->rsz = st->rsb = 0;
		st->rnc = st->rnb = 0;
		av_push(st->res_first, newSViv(st->n_atom));
		st->have_res = 1;
		st->p_resseq = rs;
		st->p_model  = st->cur_model;
		st->p_het    = het;
		my_strlcpy(st->p_chain, chain, sizeof(st->p_chain));
		my_strlcpy(st->p_icode, icode, sizeof(st->p_icode));
		my_strlcpy(st->p_resname, resname, sizeof(st->p_resname));
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
			st->have_bbox = 1;
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
		SV **cnt = hv_fetch(st->elements, elbuf, (I32)ellen, 1);
		if (cnt && *cnt) sv_inc(*cnt);
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
	char *p = SvPVX(sv);
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
	int i, keep_meta;
	IV lineno = 0;
	SV *block = NULL;
	const char *v[A_NFIELD];
	STRLEN vn[A_NFIELD];
	//an _atom_site table written as plain tags rather than as a loop_
	const char *sv_[A_NFIELD];
	STRLEN svn_[A_NFIELD];
	int have_single = 0;

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
	st.cur_model  = 1;
	st.want_model = opt_iv(aTHX_ opts, "model", 1);
	st.keep_h     = opt_bool(aTHX_ opts, "hydrogens", 1);
	st.keep_water = opt_bool(aTHX_ opts, "waters", 1);
	st.keep_het   = opt_bool(aTHX_ opts, "hetatm", 1);
	st.keep_anisou= opt_bool(aTHX_ opts, "anisou", 0);
	st.build_atoms= opt_bool(aTHX_ opts, "atom_hashes", 0);
	keep_meta     = opt_bool(aTHX_ opts, "meta", 1);
	{
		SV *c = opt_get(aTHX_ opts, "chains");
		if (c && SvROK(c) && SvTYPE(SvRV(c)) == SVt_PVHV) st.want_chain = (HV *)SvRV(c);
	}

	{	//lines are not what this format is made of, but callers still count them
		const char *p = buf, *e = buf + len;
		while (p < e) {
			const char *nl = (const char *)memchr(p, '\n', (STRLEN)(e - p));
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
			const char **tags = NULL;
			STRLEN *tagn = NULL;
			int *fld_of = NULL;
			int ntags = 0, cap = 16, is_atom = 0, is_aniso = 0, keep_loop = 0;
			const char *cat = NULL;
			STRLEN catn = 0;
			AV *rows = NULL;

			Newx(tags, cap, const char *);
			Newx(tagn, cap, STRLEN);
			for (;;) {
				const char *save = lx.p;
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

			Newx(fld_of, ntags, int);
			for (i = 0; i < ntags; i++) {
				const char *item;
				STRLEN itemn;
				(void)cif_split_tag(tags[i], tagn[i], &item, &itemn);
				fld_of[i] = is_atom ? cif_atom_field(item, itemn) : -1;
			}

			for (;;) {
				const char *save = lx.p;
				int col_i;
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
						int f = fld_of[col_i];
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
				int f = cif_atom_field(item, itemn);
				if (f >= 0 && !cif_null(&tk)) { sv_[f] = tk.s; svn_[f] = tk.n; have_single = 1; }
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
		av_push(st.res_last, newSViv(st.n_atom - 1));
		flush_residue(aTHX_ res_sum, st.rsx, st.rsy, st.rsz, st.rnc, st.rsb, st.rnb);
	}

	(void)hv_stores(out, "atoms", newRV_noinc((SV *)st.atom_hv));
	for (i = 0; i < NRSUM; i++)
		(void)hv_store(out, res_sum_name[i], (I32)strlen(res_sum_name[i]),
		               newRV_noinc((SV *)res_sum[i]), 0);
	(void)hv_stores(out, "elements",      newRV_noinc((SV *)st.elements));
	(void)hv_stores(out, "n_hydrogens",   newSViv(st.n_hydrogen));
	(void)hv_stores(out, "n_water_atoms", newSViv(st.n_water_atom));
	if (st.bn) {
		HV *b = newHV();
		(void)hv_stores(b, "min",  newSVnv(st.bmin));
		(void)hv_stores(b, "max",  newSVnv(st.bmax));
		(void)hv_stores(b, "mean", newSVnv(st.bsum / (NV)st.bn));
		(void)hv_stores(b, "n",    newSViv(st.bn));
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
	(void)hv_stores(out, "n_atoms",       newSViv(st.n_atom));
	(void)hv_stores(out, "n_residues",    newSViv(av_len(st.res_first) + 1));
	(void)hv_stores(out, "n_models",      newSViv(st.n_models ? st.n_models : 1));
	(void)hv_stores(out, "n_anisou",      newSViv(st.n_anisou));
	(void)hv_stores(out, "n_skipped",     newSViv(st.n_skipped));
	(void)hv_stores(out, "n_atom_records",   newSViv(st.n_atom_rec));
	(void)hv_stores(out, "n_hetatm_records", newSViv(st.n_het_rec));
	(void)hv_stores(out, "n_lines",       newSViv(lineno));
	(void)hv_stores(out, "format",        newSVpvs("mmcif"));
	return out;
}

/*slurp() -- read the whole file.  Chunked rather than stat-then-read so that
a named pipe or /dev/stdin works the same as a file on disk.*/
static char *slurp(pTHX_ const char *restrict CSP_RESTRICT path, STRLEN *restrict CSP_RESTRICT lenp)
{
	FILE *restrict fh = fopen(path, "rb");
	char *restrict buf;
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

MODULE = Chem::Structure::Parser		PACKAGE = Chem::Structure::Parser

PROTOTYPES: DISABLE

SV * _parse_file(path, opts = &PL_sv_undef)
	SV *path
	SV *opts
	PREINIT:
		char *restrict buf;
		STRLEN len;
		HV *o = NULL, *res;
		const char *restrict p;
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
		const char *restrict buf;
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
		char *restrict buf;
		STRLEN len;
		HV *o = NULL, *res;
		const char *restrict p;
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
		const char *restrict buf;
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
		const char *restrict s;
		NV a = 0, b = 0;
		int ok_a, ok_b;
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
		const char *s;
		res_info ri;
		const char *t;
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
