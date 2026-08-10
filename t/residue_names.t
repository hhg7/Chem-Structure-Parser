#!/usr/bin/env perl
# aa3to1(), res1() and res_type() -- the residue name table in the XS.
require 5.010;
use strict;
use warnings FATAL => 'all';
use Structure::Info;
use Test::Exception;
use Test::More;

#--------
# the twenty
#--------
my %STANDARD = (
	ALA => 'A', ARG => 'R', ASN => 'N', ASP => 'D', CYS => 'C',
	GLN => 'Q', GLU => 'E', GLY => 'G', HIS => 'H', ILE => 'I',
	LEU => 'L', LYS => 'K', MET => 'M', PHE => 'F', PRO => 'P',
	SER => 'S', THR => 'T', TRP => 'W', TYR => 'Y', VAL => 'V',
);
is(aa3to1($_), $STANDARD{$_}, "aa3to1: $_ is $STANDARD{$_}") for sort keys %STANDARD;

# every one of the twenty gets a different letter, which is the one property
# the table cannot get away with breaking
my %seen;
$seen{ aa3to1($_) }++ for keys %STANDARD;
is(scalar keys %seen, 20, 'aa3to1: the twenty map to twenty distinct letters');

#--------
# the codes the format allows for an ambiguous or unusual residue
#--------
is(aa3to1('ASX'), 'B', 'aa3to1: ASX is B (ASP or ASN)');
is(aa3to1('GLX'), 'Z', 'aa3to1: GLX is Z (GLU or GLN)');
is(aa3to1('XLE'), 'J', 'aa3to1: XLE is J (LEU or ILE)');
is(aa3to1('SEC'), 'U', 'aa3to1: SEC is U (selenocysteine)');
is(aa3to1('PYL'), 'O', 'aa3to1: PYL is O (pyrrolysine)');
is(aa3to1('UNK'), 'X', 'aa3to1: UNK is X');

#--------
# modified residues map to the residue they were made from.  A structure
# solved with selenomethionine has the same sequence as one without it, and
# an X every seventh position is no use to anyone.
#--------
my %MODIFIED = (
	MSE => 'M', FME => 'M',            # selenomethionine, formylmethionine
	SEP => 'S', TPO => 'T', PTR => 'Y',# phosphoserine/threonine/tyrosine
	CSO => 'C', CME => 'C', OCS => 'C',
	MLY => 'K', KCX => 'K', LLP => 'K',
	HYP => 'P', PCA => 'E', CGU => 'E',
	HSD => 'H', HSE => 'H', HIP => 'H',# the names force fields use
	SAR => 'G', ABA => 'A', NLE => 'L',
	DAL => 'A', DPN => 'F', DTY => 'Y',# D-amino acids
	CIR => 'R',                        # citrulline
);
is(aa3to1($_), $MODIFIED{$_}, "aa3to1: $_ maps to its parent $MODIFIED{$_}") for sort keys %MODIFIED;

#--------
# things that are not amino acids
#--------
is(aa3to1($_), '', "aa3to1: $_ is not an amino acid") for qw(HOH WAT NAG ZN SO4 DA A ATP);

#--------
# whitespace and case: names arrive straight out of columns 18-20, which are
# blank padded, and files are not consistent about case
#--------
is(aa3to1('  ALA'), 'A', 'aa3to1: leading blanks are ignored');
is(aa3to1('ALA  '), 'A', 'aa3to1: trailing blanks are ignored');
is(aa3to1(' ALA '), 'A', 'aa3to1: blanks on both sides are ignored');
is(aa3to1('ala'),   'A', 'aa3to1: lower case is accepted');
is(aa3to1('Ala'),   'A', 'aa3to1: mixed case is accepted');
is(aa3to1(''),      '',  'aa3to1: the empty string is not an amino acid');
is(aa3to1('    '),  '',  'aa3to1: blanks alone are not an amino acid');
is(aa3to1('TOOLONG'), '', 'aa3to1: a name longer than three characters is unknown');

#--------
# res1() widens the same table to nucleotides
#--------
is(res1('ALA'), 'A', 'res1: amino acids answer as aa3to1 does');
is(res1('MSE'), 'M', 'res1: modified amino acids too');
is(res1(' DA'), 'A', 'res1: DA is deoxyadenosine');
is(res1('DA'),  'A', 'res1: DA without its padding blank');
is(res1('DC'),  'C', 'res1: DC');
is(res1('DG'),  'G', 'res1: DG');
is(res1('DT'),  'T', 'res1: DT');
is(res1('  A'), 'A', 'res1: A is adenosine');
is(res1('U'),   'U', 'res1: U is uridine');
is(res1('PSU'), 'U', 'res1: pseudouridine is a U');
is(res1('5MC'), 'C', 'res1: 5-methylcytidine is a C');
is(res1('HOH'), '',  'res1: water has no single-letter code');
is(res1('NAG'), '',  'res1: a sugar has no single-letter code');

# CYS is an amino acid and CYT is cytosine: three-letter names that differ in
# one character must not collide in the packed key
is(aa3to1('CYS'),   'C', 'CYS is cysteine');
is(res_type('CYT'), 'nucleotide', 'CYT is cytosine, not cysteine');

#--------
# res_type()
#--------
is(res_type('ALA'), 'amino_acid', 'res_type: ALA');
is(res_type('MSE'), 'amino_acid', 'res_type: MSE is still an amino acid');
is(res_type('DA'),  'nucleotide', 'res_type: DA');
is(res_type('U'),   'nucleotide', 'res_type: U');
is(res_type($_), 'water', "res_type: $_ is water") for qw(HOH WAT DOD H2O SOL);
is(res_type($_), 'other', "res_type: $_ is neither polymer nor water") for qw(NAG ZN SO4 ATP HEM);

#--------
# an undefined name is a mistake worth hearing about: it means a column was
# read out of a record that did not have one
#--------
throws_ok { aa3to1(undef)   } qr/undefined/, 'aa3to1: undef dies';
throws_ok { res1(undef)     } qr/undefined/, 'res1: undef dies';
throws_ok { res_type(undef) } qr/undefined/, 'res_type: undef dies';

done_testing();
