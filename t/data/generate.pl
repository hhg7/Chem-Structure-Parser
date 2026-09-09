#!/usr/bin/env perl
# Writes the fixtures in this directory.
#
# The fixtures are generated rather than typed because a PDB record is a
# fixed-column format: an atom name in the wrong column is a different
# element, and a residue number shifted by one is a different residue.  Hand
# editing gets that wrong silently.  Run this from t/data/ after changing a
# fixture, and commit both the script and what it wrote.
require 5.010;
use strict;
use warnings FATAL => 'all';
use autodie ':default';

# atom_line() -- one ATOM/HETATM record, in the columns the format wants.
#
# The atom name is the fiddly part: a one-letter element is right-justified
# from column 14 (" CA " is a carbon alpha), a two-letter element starts in
# column 13 ("CA  " is a calcium).  That is the rule the parser's element
# guess relies on, so the fixtures have to follow it exactly.
sub atom_line {
	my (%a) = @_;
	my $name = length($a{element}) == 2 || length($a{name}) == 4
	         ? sprintf('%-4s', $a{name})
	         : sprintf(' %-3s', $a{name});
	return sprintf(
		'%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%-2s',
		$a{record}, $a{serial}, $name, ($a{altloc} // ''), $a{resname},
		$a{chain}, $a{resseq}, ($a{icode} // ''),
		$a{x}, $a{y}, $a{z}, ($a{occ} // 1), ($a{b} // 20),
		$a{element}, ($a{charge} // '')
	);
}

my $serial = 0;
sub atom {
	my ($rec, $chain, $resname, $resseq, $icode, $atoms, $base) = @_;
	my @out;
	my $i = 0;
	for my $a (@$atoms) {
		my ($name, $element, $altloc, $occ) = @$a;
		$serial++;
		push @out, atom_line(
			record => $rec, serial => $serial, name => $name, element => $element,
			altloc => $altloc, resname => $resname, chain => $chain,
			resseq => $resseq, icode => $icode,
			x => $base->[0] + $i * 1.5, y => $base->[1] + $i, z => $base->[2] + $i * 0.5,
			occ => $occ, b => 15 + $i,
		);
		$i++;
	}
	return @out;
}

# backbone of a residue, plus whatever side-chain atoms are named
sub bb { my @s = @_; return ([ 'N', 'N' ], [ 'CA', 'C' ], [ 'C', 'C' ], [ 'O', 'O' ], @s) }

# cols() -- build a record by column number rather than by counting spaces.
# Takes (start, width, value) triples, with start being the 1-based column the
# format specification gives, optionally followed by 'R' for a right-justified
# field.  Records below can then be checked against the spec by reading the
# numbers rather than counting anything.
sub cols {
	my @spec = @_;
	my $line = ' ' x 80;
	while (@spec) {
		my ($start, $width, $value) = splice @spec, 0, 3;
		# a start column is always a number, so an 'R' here is this field's
		# justification and not the beginning of the next one
		my $just = (@spec && defined $spec[0] && $spec[0] =~ /\A[LR]\z/) ? shift(@spec) : 'L';
		$value = '' unless defined $value;
		$value = $just eq 'R' ? sprintf('%*s', $width, $value)
		                      : sprintf('%-*s', $width, $value);
		substr($line, $start - 1, $width) = substr($value, 0, $width);
	}
	$line =~ s/\s+\z//;
	return $line;
}

# --- mini.pdb -- one of everything the reader knows how to look at ---------
my @mini = (
'HEADER    HYDROLASE/PEPTIDE INHIBITOR             01-JAN-20   9XYZ              ',
'TITLE     A SMALL TEST STRUCTURE WITH A GAP, AN INSERTION CODE, AN              ',
'TITLE    2 ALTERNATE CONFORMER AND A LIGAND                                     ',
'COMPND    MOL_ID: 1;                                                            ',
'COMPND   2 MOLECULE: TEST PROTEIN;                                              ',
'COMPND   3 CHAIN: A;                                                            ',
'COMPND   4 EC: 3.4.21.5;                                                        ',
'COMPND   5 ENGINEERED: YES;                                                     ',
'COMPND   6 MOL_ID: 2;                                                           ',
'COMPND   7 MOLECULE: TEST DNA;                                                  ',
'COMPND   8 CHAIN: B                                                             ',
'SOURCE    MOL_ID: 1;                                                            ',
'SOURCE   2 ORGANISM_SCIENTIFIC: HOMO SAPIENS;                                   ',
'SOURCE   3 ORGANISM_COMMON: HUMAN;                                              ',
'SOURCE   4 ORGANISM_TAXID: 9606;                                                ',
'SOURCE   5 EXPRESSION_SYSTEM: ESCHERICHIA COLI;                                 ',
'SOURCE   6 MOL_ID: 2;                                                           ',
'SOURCE   7 SYNTHETIC: YES                                                       ',
'KEYWDS    HYDROLASE, TEST STRUCTURE, COMPLEX (HYDROLASE-                        ',
'KEYWDS   2 PEPTIDE)                                                             ',
'EXPDTA    X-RAY DIFFRACTION                                                     ',
'NUMMDL    1                                                                     ',
'AUTHOR    D.E.CONDON,A.N.OTHER                                                  ',
'REVDAT   1   01-JAN-20 9XYZ    0                                                ',
'JRNL        AUTH   D.E.CONDON,A.N.OTHER                                         ',
'JRNL        TITL   A STRUCTURE MADE UP FOR A TEST SUITE, AND WHAT IT            ',
'JRNL        TITL 2 CONTAINS                                                     ',
'JRNL        REF    J.INVENTED.RES.               V.  10    42 2020              ',
'JRNL        PMID   12345678                                                     ',
'JRNL        DOI    10.1000/INVENTED.2020.42                                     ',
'REMARK   2                                                                      ',
'REMARK   2 RESOLUTION.    1.85 ANGSTROMS.                                       ',
'REMARK   3                                                                      ',
'REMARK   3   R VALUE            (WORKING SET) : 0.174                           ',
'REMARK   3   FREE R VALUE                     : 0.219                           ',
'REMARK   3   BIN FREE R VALUE                    : 0.999                        ',
'REMARK 200   TEMPERATURE           (KELVIN) : 100.0                             ',
'REMARK 200   PH                             : 7.5                               ',
'REMARK 465   MISSING RESIDUES                                                   ',
);

# The fixed-field annotation records, placed by the column numbers in the PDB
# format specification (v3.3).  Written this way so that a record can be
# checked against the spec by reading the numbers, not by counting spaces --
# which is how the LINK record in an earlier draft of this file ended up two
# columns to the left and silently parsed its chain ids as blanks.
push @mini,
	# DBREF: idCode 8-11, chain 13, seqBegin 15-18, seqEnd 21-24,
	#        database 27-32, dbAccession 34-41, dbIdCode 43-54,
	#        dbseqBegin 56-60, dbseqEnd 63-67
	cols(1,6,'DBREF', 8,4,'9XYZ', 13,1,'A', 15,4,1,'R', 21,4,11,'R',
	     27,6,'UNP', 34,8,'P12345', 43,12,'TEST_HUMAN', 56,5,1,'R', 63,5,11,'R'),
	# SEQADV: resName 13-15, chain 17, seqNum 19-22, database 25-28,
	#         dbAccession 30-38, dbRes 40-42, dbSeq 44-48, conflict 50-70
	cols(1,6,'SEQADV', 8,4,'9XYZ', 13,3,'MSE', 17,1,'A', 19,4,7,'R',
	     25,4,'UNP', 30,9,'P12345', 40,3,'MET', 44,5,7,'R', 50,21,'MODIFIED RESIDUE'),
	# SEQRES: serNum 8-10, chain 12, numRes 14-17, residues from 20 in 4s
	cols(1,6,'SEQRES', 8,3,1,'R', 12,1,'A', 14,4,11,'R',
	     20,51,join('', map { sprintf('%3s ', $_) } qw(MET ALA GLY LEU LYS CYS MSE HIS HIS SER CYS))),
	cols(1,6,'SEQRES', 8,3,1,'R', 12,1,'B', 14,4,4,'R',
	     20,51,join('', map { sprintf('%3s ', $_) } qw(DA DC DG DT))),
	# MODRES: idCode 8-11, resName 13-15, chain 17, seqNum 19-22,
	#         stdRes 25-27, comment 30-70
	cols(1,6,'MODRES', 8,4,'9XYZ', 13,3,'MSE', 17,1,'A', 19,4,7,'R',
	     25,3,'MET', 30,41,'SELENOMETHIONINE'),
	# HET: hetID 8-10, chain 13, seqNum 14-17, iCode 18, numHetAtoms 21-25
	cols(1,6,'HET', 8,3,'NAG','R', 13,1,'A', 14,4,201,'R', 21,5,14,'R'),
	cols(1,6,'HET', 8,3,'ZN', 'R', 13,1,'A', 14,4,202,'R', 21,5,1, 'R'),
	# HETNAM: continuation 9-10, hetID 12-14, text 16-70
	cols(1,6,'HETNAM', 12,3,'NAG','R', 16,55,'2-ACETAMIDO-2-DEOXY-BETA-D-GLUCOPYRANOSE'),
	cols(1,6,'HETNAM', 12,3,'ZN', 'R', 16,55,'ZINC ION'),
	# FORMUL: compNum 9-10, hetID 13-15, continuation 17-18, asterisk 19, text 20-70
	cols(1,6,'FORMUL', 9,2,3,'R', 13,3,'NAG','R', 20,51,'C8 H15 N O6'),
	cols(1,6,'FORMUL', 9,2,4,'R', 13,3,'ZN', 'R', 20,51,'ZN 2+'),
	cols(1,6,'FORMUL', 9,2,5,'R', 13,3,'HOH','R', 19,1,'*', 20,51,'2(H2 O)'),
	# HELIX: serNum 8-10, helixID 12-14, initResName 16-18, initChain 20,
	#        initSeqNum 22-25, endResName 28-30, endChain 32, endSeqNum 34-37,
	#        helixClass 39-40, length 72-76
	cols(1,6,'HELIX', 8,3,1,'R', 12,3,'AA1', 16,3,'MET', 20,1,'A', 22,4,1,'R',
	     28,3,'GLY', 32,1,'A', 34,4,3,'R', 39,2,1,'R', 72,5,3,'R'),
	# SHEET: strand 8-10, sheetID 12-14, numStrands 15-16, initResName 18-20,
	#        initChain 22, initSeqNum 23-26, endResName 29-31, endChain 33,
	#        endSeqNum 34-37, sense 39-40
	cols(1,6,'SHEET', 8,3,1,'R', 12,3,'AA1', 15,2,2,'R', 18,3,'CYS', 22,1,'A',
	     23,4,6,'R', 29,3,'HIS', 33,1,'A', 34,4,8,'R', 39,2,0,'R'),
	# SSBOND: serNum 8-10, CYS 12-14, chain1 16, seqNum1 18-21,
	#         CYS 26-28, chain2 30, seqNum2 32-35, sym1 60-65, sym2 67-72, length 74-78
	cols(1,6,'SSBOND', 8,3,1,'R', 12,3,'CYS', 16,1,'A', 18,4,6,'R',
	     26,3,'CYS', 30,1,'A', 32,4,10,'R', 60,6,'1555', 67,6,'1555', 74,5,'2.03','R'),
	# LINK: name1 13-16, resName1 18-20, chain1 22, resSeq1 23-26,
	#       name2 43-46, resName2 48-50, chain2 52, resSeq2 53-56, length 74-78
	cols(1,6,'LINK', 13,4,'ZN', 18,3,'ZN','R', 22,1,'A', 23,4,202,'R',
	     43,4,' SG', 48,3,'CYS', 52,1,'A', 53,4,6,'R', 74,5,'2.31','R'),
	# CISPEP: serNum 8-10, pep1 12-14, chain1 16, seqNum1 18-21,
	#         pep2 26-28, chain2 30, seqNum2 32-35, modNum 44-46, measure 54-59
	cols(1,6,'CISPEP', 8,3,1,'R', 12,3,'GLY', 16,1,'A', 18,4,3,'R',
	     26,3,'CYS', 30,1,'A', 32,4,6,'R', 44,3,0,'R', 54,6,'-0.42','R'),
	# CRYST1: a 7-15, b 16-24, c 25-33, alpha 34-40, beta 41-47, gamma 48-54,
	#         sGroup 56-66, z 67-70
	cols(1,6,'CRYST1', 7,9,'40.100','R', 16,9,'50.200','R', 25,9,'60.300','R',
	     34,7,'90.00','R', 41,7,'95.50','R', 48,7,'90.00','R',
	     56,11,'P 1 21 1', 67,4,4,'R');

push @mini, atom('ATOM  ', 'A', 'MET',  1, '', [ bb([ 'CB', 'C' ]) ], [ 10, 10, 10 ]);
# an alternate conformer: two CB records, altlocs A and B, unequal occupancy
push @mini, atom('ATOM  ', 'A', 'ALA',  2, '', [ bb() ], [ 13, 11, 11 ]);
push @mini, atom_line(record => 'ATOM  ', serial => ++$serial, name => 'CB', element => 'C',
	altloc => 'A', resname => 'ALA', chain => 'A', resseq => 2, icode => '',
	x => 19.0, y => 15.0, z => 13.0, occ => 0.40, b => 22);
push @mini, atom_line(record => 'ATOM  ', serial => ++$serial, name => 'CB', element => 'C',
	altloc => 'B', resname => 'ALA', chain => 'A', resseq => 2, icode => '',
	x => 19.5, y => 15.5, z => 13.5, occ => 0.60, b => 25);
push @mini, atom('ATOM  ', 'A', 'GLY',  3, '', [ bb() ], [ 16, 12, 12 ]);
# residues 4 and 5 are in SEQRES but were never modelled: a gap
push @mini, atom('ATOM  ', 'A', 'CYS',  6, '', [ bb([ 'SG', 'S' ]) ], [ 19, 13, 13 ]);
# a modified residue, written as HETATM, that is still an M in the sequence
push @mini, atom('HETATM', 'A', 'MSE',  7, '', [ bb([ 'CB', 'C' ], [ 'SE', 'SE' ]) ], [ 22, 14, 14 ]);
push @mini, atom('ATOM  ', 'A', 'HIS',  8, '',  [ bb([ 'CB', 'C' ]) ], [ 25, 15, 15 ]);
push @mini, atom('ATOM  ', 'A', 'HIS',  8, 'A', [ bb([ 'CB', 'C' ]) ], [ 28, 16, 16 ]);
# a hydrogen, so that hydrogens => 0 has something to remove
push @mini, atom('ATOM  ', 'A', 'SER',  9, '', [ bb([ 'CB', 'C' ], [ 'HB2', 'H' ]) ], [ 31, 17, 17 ]);
push @mini, atom('ATOM  ', 'A', 'CYS', 10, '', [ bb([ 'SG', 'S' ]) ], [ 34, 18, 18 ]);
push @mini, sprintf('%-6s%5d      %3s %1s%4d%1s', 'TER', ++$serial, 'CYS', 'A', 10, '');
push @mini, atom('HETATM', 'A', 'NAG', 201, '', [ [ 'C1', 'C' ], [ 'C2', 'C' ], [ 'O5', 'O' ], [ 'N2', 'N' ] ], [ 40, 20, 20 ]);
push @mini, atom('HETATM', 'A', 'ZN',  202, '', [ [ 'ZN', 'ZN' ] ], [ 45, 22, 22 ]);
push @mini, atom('HETATM', 'A', 'HOH', 301, '', [ [ 'O', 'O' ] ], [ 50, 24, 24 ]);
push @mini, atom('HETATM', 'A', 'HOH', 302, '', [ [ 'O', 'O' ] ], [ 52, 25, 25 ]);

# a DNA chain, to have a chain that is not a protein
my $z = 0;
for my $r ([ 'DA', 1 ], [ 'DC', 2 ], [ 'DG', 3 ], [ 'DT', 4 ]) {
	push @mini, atom('ATOM  ', 'B', $r->[0], $r->[1], '',
		[ [ 'P', 'P' ], [ 'OP1', 'O' ], [ "C1'", 'C' ] ], [ 60 + $z, 30, 30 ]);
	$z += 3;
}
push @mini, sprintf('%-6s%5d      %3s %1s%4d%1s', 'TER', ++$serial, 'DT', 'B', 4, '');
push @mini, 'CONECT   57   58   59';
push @mini, 'MASTER      000    0    0    1    1    0    0    6   66    2    0    2          ';
push @mini, 'END';

# --- nmr.pdb -- three models of the same tripeptide -----------------------
my @nmr = (
'HEADER    TEST                                    01-JAN-20   9NMR              ',
'TITLE     A THREE MODEL ENSEMBLE                                                ',
'EXPDTA    SOLUTION NMR                                                          ',
'NUMMDL    3                                                                     ',
'SEQRES   1 A    3  GLY SER TRP                                                  ',
);
for my $m (1 .. 3) {
	$serial = 0;
	push @nmr, sprintf('MODEL     %4d', $m);
	my $i = 0;
	for my $r ([ 'GLY', 1 ], [ 'SER', 2 ], [ 'TRP', 3 ]) {
		push @nmr, atom('ATOM  ', 'A', $r->[0], $r->[1], '', [ bb() ], [ 5 + $m, 5 + $i, 5 ]);
		$i++;
	}
	push @nmr, 'ENDMDL';
}
push @nmr, 'END';

# --- bare.pdb -- coordinates and nothing else -----------------------------
# No header at all, and no element columns, so the element has to be worked
# out from the atom name.  Files this old, and files written by simulation
# programs, both look like this.
my @bare;
$serial = 0;
for my $r ([ 'VAL', 1 ], [ 'LYS', 2 ]) {
	for my $a ([ 'N', 'N' ], [ 'CA', 'C' ], [ 'C', 'C' ], [ 'O', 'O' ], [ 'CB', 'C' ]) {
		$serial++;
		my $l = atom_line(
			record => 'ATOM  ', serial => $serial, name => $a->[0], element => $a->[1],
			resname => $r->[0], chain => 'A', resseq => $r->[1], altloc => '', icode => '',
			x => $serial, y => $serial + 1, z => $serial + 2,
		);
		push @bare, substr($l, 0, 66);    # truncate before the element columns
	}
}
push @bare, 'END';

# --- stack.pdb -- aromatic residues that are actually stacked --------------
#
# Seven residues of 1A42 (human carbonic anhydrase II, 1.05 A), lifted whole
# from /home/con/ui/pepPriML/PPB/PDB/PDBbind.v2020/1a42.ent.pdb chain A with
# their coordinates, occupancies and B factors exactly as deposited.  Nothing
# here is invented: the geometry is what makes them stack, so a made-up ring
# would test the arithmetic against itself.
#
# They were chosen by asking the module for every stacked pair in 1A42 and
# keeping the residues of the four that are in one chain and cover both kinds:
#
#   face  TRP 5 six-ring  / HIS 64          planes 4.6 degrees apart, 3.70 A
#   edge  PHE 66          / PHE 226         planes 78.4 apart, 5.36 A
#   edge  TRP 97 five-ring / PHE 226        planes 65.4 apart, 5.53 A
#   edge  HIS 119         / TRP 209 six-ring planes 67.6 apart, 4.86 A
#
# which between them exercise the six-membered ring of a phenylalanine and a
# tyrosine, both rings of a tryptophan, and the five-membered ring of a
# histidine.  The residues are not consecutive and there is no backbone
# between them; that is deliberate, and is why the file has no TER and no
# SEQRES -- it is a set of residues, not a chain that folded.
my @stack_res = (
	[ 'TRP', 5, [
		[ 'N',    'N',     8.257,   -0.414,   10.635,  1.00,  27.01 ],
		[ 'CA',   'C',     7.481,   -1.387,   11.384,  1.00,  19.87 ],
		[ 'C',    'C',     6.510,   -2.204,   10.563,  1.00,  19.91 ],
		[ 'O',    'O',     6.066,   -1.770,    9.486,  1.00,  19.19 ],
		[ 'CB',   'C',     6.719,   -0.667,   12.505,  1.00,  14.02 ],
		[ 'CG',   'C',     5.438,    0.076,   12.142,  1.00,  12.29 ],
		[ 'CD1',  'C',     5.444,    1.400,   11.760,  1.00,   6.38 ],
		[ 'CD2',  'C',     4.172,   -0.466,   12.199,  1.00,   6.22 ],
		[ 'NE1',  'N',     4.172,    1.704,   11.581,  1.00,   8.90 ],
		[ 'CE2',  'C',     3.380,    0.633,   11.826,  1.00,   7.53 ],
		[ 'CE3',  'C',     3.588,   -1.693,   12.512,  1.00,   4.86 ],
		[ 'CZ2',  'C',     1.995,    0.506,   11.767,  1.00,   3.51 ],
		[ 'CZ3',  'C',     2.199,   -1.821,   12.451,  1.00,   2.00 ],
		[ 'CH2',  'C',     1.408,   -0.726,   12.081,  1.00,   4.15 ],
	] ],
	[ 'HIS', 64, [
		[ 'N',    'N',    -0.871,   -2.013,    6.309,  1.00,  10.84 ],
		[ 'CA',   'C',    -1.016,   -2.605,    7.638,  1.00,  14.16 ],
		[ 'C',    'C',    -2.362,   -3.292,    7.812,  1.00,  13.62 ],
		[ 'O',    'O',    -2.465,   -4.335,    8.449,  1.00,  16.31 ],
		[ 'CB',   'C',    -0.853,   -1.521,    8.690,  1.00,  11.37 ],
		[ 'CG',   'C',     0.551,   -0.960,    8.667,  1.00,  12.92 ],
		[ 'ND1',  'N',     1.686,   -1.623,    8.866,  1.00,  16.68 ],
		[ 'CD2',  'C',     0.873,    0.353,    8.418,  1.00,  13.77 ],
		[ 'CE1',  'C',     2.684,   -0.777,    8.748,  1.00,  13.76 ],
		[ 'NE2',  'N',     2.181,    0.404,    8.481,  1.00,  15.68 ],
	] ],
	[ 'PHE', 66, [
		[ 'N',    'N',    -6.779,   -2.621,    6.033,  1.00,   8.55 ],
		[ 'CA',   'C',    -7.541,   -1.754,    5.181,  1.00,   9.88 ],
		[ 'C',    'C',    -8.423,   -0.958,    6.121,  1.00,  10.30 ],
		[ 'O',    'O',    -8.898,   -1.492,    7.127,  1.00,  12.04 ],
		[ 'CB',   'C',    -8.378,   -2.552,    4.184,  1.00,   2.93 ],
		[ 'CG',   'C',    -9.381,   -3.566,    4.720,  1.00,   5.04 ],
		[ 'CD1',  'C',    -8.974,   -4.859,    5.004,  1.00,   5.57 ],
		[ 'CD2',  'C',   -10.717,   -3.200,    4.911,  1.00,   6.77 ],
		[ 'CE1',  'C',    -9.895,   -5.788,    5.477,  1.00,   9.61 ],
		[ 'CE2',  'C',   -11.629,   -4.135,    5.386,  1.00,   6.91 ],
		[ 'CZ',   'C',   -11.221,   -5.428,    5.669,  1.00,   2.26 ],
	] ],
	[ 'TRP', 97, [
		[ 'N',    'N',    -6.638,   -9.426,   11.785,  1.00,   4.96 ],
		[ 'CA',   'C',    -7.086,  -10.799,   11.614,  1.00,   6.57 ],
		[ 'C',    'C',    -5.877,  -11.699,   11.551,  1.00,   6.50 ],
		[ 'O',    'O',    -4.726,  -11.233,   11.412,  1.00,   6.06 ],
		[ 'CB',   'C',    -7.941,  -10.972,   10.316,  1.00,   9.03 ],
		[ 'CG',   'C',    -7.254,  -10.572,    9.015,  1.00,   7.42 ],
		[ 'CD1',  'C',    -6.549,  -11.475,    8.255,  1.00,   3.72 ],
		[ 'CD2',  'C',    -7.228,   -9.295,    8.510,  1.00,   8.60 ],
		[ 'NE1',  'N',    -6.059,  -10.759,    7.276,  1.00,   4.52 ],
		[ 'CE2',  'C',    -6.434,   -9.467,    7.376,  1.00,   9.06 ],
		[ 'CE3',  'C',    -7.753,   -8.040,    8.824,  1.00,   6.27 ],
		[ 'CZ2',  'C',    -6.157,   -8.381,    6.543,  1.00,   8.08 ],
		[ 'CZ3',  'C',    -7.475,   -6.965,    7.991,  1.00,   5.35 ],
		[ 'CH2',  'C',    -6.683,   -7.135,    6.859,  1.00,   7.95 ],
	] ],
	[ 'HIS', 119, [
		[ 'N',    'N',   -12.075,   -2.497,   14.976,  1.00,   4.93 ],
		[ 'CA',   'C',   -11.839,   -1.316,   15.736,  1.00,   5.12 ],
		[ 'C',    'C',   -12.522,   -0.129,   15.084,  1.00,   4.80 ],
		[ 'O',    'O',   -12.185,    0.229,   13.955,  1.00,   2.03 ],
		[ 'CB',   'C',   -10.335,   -1.119,   15.812,  1.00,   7.51 ],
		[ 'CG',   'C',    -9.582,   -2.139,   16.636,  1.00,   3.45 ],
		[ 'ND1',  'N',    -8.293,   -2.403,   16.550,  1.00,   2.88 ],
		[ 'CD2',  'C',   -10.107,   -2.864,   17.672,  1.00,   2.71 ],
		[ 'CE1',  'C',    -8.001,   -3.249,   17.510,  1.00,   2.00 ],
		[ 'NE2',  'N',    -9.095,   -3.514,   18.162,  1.00,   5.49 ],
	] ],
	[ 'TRP', 209, [
		[ 'N',    'N',   -10.711,    2.259,   25.863,  1.00,   6.60 ],
		[ 'CA',   'C',   -11.092,    1.095,   25.119,  1.00,   5.32 ],
		[ 'C',    'C',   -12.423,    0.570,   25.579,  1.00,   6.87 ],
		[ 'O',    'O',   -12.554,    0.309,   26.771,  1.00,   9.07 ],
		[ 'CB',   'C',   -10.045,    0.056,   25.307,  1.00,   3.91 ],
		[ 'CG',   'C',    -8.962,    0.144,   24.258,  1.00,   6.41 ],
		[ 'CD1',  'C',    -7.752,    0.759,   24.456,  1.00,   5.46 ],
		[ 'CD2',  'C',    -9.110,   -0.395,   23.010,  1.00,   8.08 ],
		[ 'NE1',  'N',    -7.128,    0.612,   23.309,  1.00,   5.55 ],
		[ 'CE2',  'C',    -7.885,   -0.068,   22.421,  1.00,   8.97 ],
		[ 'CE3',  'C',   -10.107,   -1.098,   22.314,  1.00,   8.95 ],
		[ 'CZ2',  'C',    -7.649,   -0.456,   21.101,  1.00,   4.36 ],
		[ 'CZ3',  'C',    -9.868,   -1.477,   20.996,  1.00,   6.86 ],
		[ 'CH2',  'C',    -8.644,   -1.152,   20.408,  1.00,   8.25 ],
	] ],
	[ 'PHE', 226, [
		[ 'N',    'N',   -10.267,  -11.466,    3.537,  1.00,  16.86 ],
		[ 'CA',   'C',    -9.378,  -10.545,    4.252,  1.00,  11.14 ],
		[ 'C',    'C',    -7.959,  -10.982,    4.058,  1.00,   8.34 ],
		[ 'O',    'O',    -7.092,  -10.183,    3.737,  1.00,   9.01 ],
		[ 'CB',   'C',    -9.594,  -10.512,    5.778,  1.00,  11.87 ],
		[ 'CG',   'C',   -10.858,   -9.859,    6.337,  1.00,  15.35 ],
		[ 'CD1',  'C',   -11.949,   -9.518,    5.524,  1.00,  13.24 ],
		[ 'CD2',  'C',   -10.920,   -9.623,    7.715,  1.00,  16.54 ],
		[ 'CE1',  'C',   -13.095,   -8.949,    6.086,  1.00,  19.53 ],
		[ 'CE2',  'C',   -12.064,   -9.056,    8.276,  1.00,  19.19 ],
		[ 'CZ',   'C',   -13.153,   -8.718,    7.466,  1.00,  22.23 ],
	] ],
);
my @stack;
$serial = 0;
for my $r (@stack_res) {
	my ($resname, $resseq, $atoms) = @$r;
	for my $a (@$atoms) {
		$serial++;
		push @stack, atom_line(
			record => 'ATOM  ', serial => $serial, name => $a->[0], element => $a->[1],
			altloc => '', resname => $resname, chain => 'A', resseq => $resseq, icode => '',
			x => $a->[2], y => $a->[3], z => $a->[4], occ => $a->[5], b => $a->[6],
		);
	}
}
push @stack, 'END';

# --- bases.pdb -- a stack of nucleobases -----------------------------------
#
# Chain P of 3AU6 (an archaeal PolX bound to a gapped DNA), lifted whole from
# /home/con/ui/pepPriML/PPB/PDB/PDBbind.v2020/3au6.ent.pdb, again exactly as
# deposited.  Seven nucleotides, all four DNA bases, and eight face-to-face
# stacks between consecutive ones -- which is the other half of the ring table:
# the six- and five-membered rings of a purine and the six-membered ring of a
# pyrimidine.
#
# A separate file from stack.pdb rather than a second chain in it, because the
# two came out of different entries and so are in different frames of
# reference: putting them in one file would invite whatever contacts the two
# coordinate systems happened to produce, which would mean nothing.
my @bases_res = (
	[ 'DC', 1, [
		[ "O5'",  'O',    -9.582,  -36.543,   33.130,  1.00,  61.80 ],
		[ "C5'",  'C',    -8.825,  -36.688,   34.320,  1.00,  61.52 ],
		[ "C4'",  'C',    -7.494,  -37.368,   34.079,  1.00,  61.00 ],
		[ "O4'",  'O',    -7.732,  -38.690,   33.607,  1.00,  60.58 ],
		[ "C3'",  'C',    -6.614,  -36.774,   33.007,  1.00,  61.45 ],
		[ "O3'",  'O',    -5.808,  -35.855,   33.628,  1.00,  63.73 ],
		[ "C2'",  'C',    -5.796,  -37.931,   32.458,  1.00,  60.15 ],
		[ "C1'",  'C',    -6.715,  -39.096,   32.717,  1.00,  60.16 ],
		[ 'N1',   'N',    -7.426,  -39.650,   31.567,  1.00,  58.92 ],
		[ 'C2',   'C',    -6.743,  -40.331,   30.565,  1.00,  59.75 ],
		[ 'O2',   'O',    -5.498,  -40.451,   30.600,  1.00,  59.27 ],
		[ 'N3',   'N',    -7.501,  -40.835,   29.553,  1.00,  59.87 ],
		[ 'C4',   'C',    -8.824,  -40.676,   29.534,  1.00,  56.89 ],
		[ 'N4',   'N',    -9.539,  -41.173,   28.541,  1.00,  56.58 ],
		[ 'C5',   'C',    -9.501,  -40.000,   30.549,  1.00,  57.52 ],
		[ 'C6',   'C',    -8.771,  -39.520,   31.535,  1.00,  58.07 ],
	] ],
	[ 'DA', 2, [
		[ 'P',    'P',    -5.110,  -34.750,   32.754,  1.00,  65.30 ],
		[ 'OP1',  'O',    -5.557,  -33.500,   33.406,  1.00,  66.32 ],
		[ 'OP2',  'O',    -5.407,  -34.972,   31.325,  1.00,  66.30 ],
		[ "O5'",  'O',    -3.554,  -35.016,   32.906,  1.00,  65.14 ],
		[ "C5'",  'C',    -3.079,  -35.975,   33.792,  1.00,  67.61 ],
		[ "C4'",  'C',    -2.072,  -36.873,   33.104,  1.00,  68.92 ],
		[ "O4'",  'O',    -2.713,  -37.691,   32.098,  1.00,  69.05 ],
		[ "C3'",  'C',    -0.965,  -36.135,   32.370,  1.00,  69.18 ],
		[ "O3'",  'O',     0.245,  -36.826,   32.590,  1.00,  70.41 ],
		[ "C2'",  'C',    -1.480,  -36.164,   30.936,  1.00,  69.14 ],
		[ "C1'",  'C',    -1.985,  -37.576,   30.888,  1.00,  69.32 ],
		[ 'N9',   'N',    -2.954,  -37.923,   29.873,  1.00,  70.65 ],
		[ 'C8',   'C',    -4.179,  -37.347,   29.699,  1.00,  71.08 ],
		[ 'N7',   'N',    -4.915,  -37.872,   28.760,  1.00,  69.99 ],
		[ 'C5',   'C',    -4.101,  -38.871,   28.301,  1.00,  70.40 ],
		[ 'C6',   'C',    -4.319,  -39.788,   27.285,  1.00,  71.09 ],
		[ 'N6',   'N',    -5.467,  -39.769,   26.600,  1.00,  69.40 ],
		[ 'N1',   'N',    -3.310,  -40.676,   27.030,  1.00,  73.29 ],
		[ 'C2',   'C',    -2.171,  -40.645,   27.772,  1.00,  72.89 ],
		[ 'N3',   'N',    -1.859,  -39.807,   28.771,  1.00,  72.02 ],
		[ 'C4',   'C',    -2.884,  -38.941,   28.969,  1.00,  71.12 ],
	] ],
	[ 'DG', 3, [
		[ 'P',    'P',     1.505,  -36.494,   31.680,  1.00,  74.08 ],
		[ 'OP1',  'O',     2.769,  -36.955,   32.313,  1.00,  73.14 ],
		[ 'OP2',  'O',     1.367,  -35.053,   31.382,  1.00,  75.77 ],
		[ "O5'",  'O',     1.182,  -37.269,   30.306,  1.00,  74.84 ],
		[ "C5'",  'C',     1.414,  -38.659,   30.106,  1.00,  75.73 ],
		[ "C4'",  'C',     2.231,  -38.953,   28.850,  1.00,  76.52 ],
		[ "O4'",  'O',     1.378,  -39.089,   27.681,  1.00,  75.82 ],
		[ "C3'",  'C',     3.272,  -37.926,   28.404,  1.00,  77.07 ],
		[ "O3'",  'O',     4.294,  -38.652,   27.649,  1.00,  78.57 ],
		[ "C2'",  'C',     2.415,  -36.935,   27.613,  1.00,  75.93 ],
		[ "C1'",  'C',     1.450,  -37.898,   26.931,  1.00,  74.90 ],
		[ 'N9',   'N',     0.080,  -37.454,   26.831,  1.00,  73.34 ],
		[ 'C8',   'C',    -0.609,  -36.487,   27.531,  1.00,  72.20 ],
		[ 'N7',   'N',    -1.852,  -36.382,   27.163,  1.00,  71.06 ],
		[ 'C5',   'C',    -1.967,  -37.343,   26.164,  1.00,  71.80 ],
		[ 'C6',   'C',    -3.063,  -37.713,   25.362,  1.00,  71.73 ],
		[ 'O6',   'O',    -4.189,  -37.249,   25.380,  1.00,  72.03 ],
		[ 'N1',   'N',    -2.777,  -38.738,   24.476,  1.00,  72.15 ],
		[ 'C2',   'C',    -1.552,  -39.332,   24.366,  1.00,  73.18 ],
		[ 'N2',   'N',    -1.487,  -40.288,   23.443,  1.00,  75.14 ],
		[ 'N3',   'N',    -0.489,  -39.011,   25.094,  1.00,  72.21 ],
		[ 'C4',   'C',    -0.790,  -38.014,   25.959,  1.00,  72.22 ],
	] ],
	[ 'DT', 4, [
		[ 'P',    'P',     5.307,  -37.953,   26.623,  1.00,  79.20 ],
		[ 'OP1',  'O',     6.568,  -38.728,   26.549,  1.00,  78.67 ],
		[ 'OP2',  'O',     5.338,  -36.498,   26.943,  1.00,  77.91 ],
		[ "O5'",  'O',     4.515,  -38.220,   25.261,  1.00,  79.72 ],
		[ "C5'",  'C',     4.240,  -39.570,   24.895,  1.00,  80.54 ],
		[ "C4'",  'C',     3.606,  -39.640,   23.514,  1.00,  81.11 ],
		[ "O4'",  'O',     2.225,  -39.166,   23.549,  1.00,  80.87 ],
		[ "C3'",  'C',     4.308,  -38.796,   22.454,  1.00,  80.86 ],
		[ "O3'",  'O',     4.141,  -39.371,   21.161,  1.00,  82.80 ],
		[ "C2'",  'C',     3.539,  -37.504,   22.571,  1.00,  79.81 ],
		[ "C1'",  'C',     2.135,  -38.100,   22.630,  1.00,  79.11 ],
		[ 'N1',   'N',     1.081,  -37.120,   23.011,  1.00,  77.00 ],
		[ 'C2',   'C',    -0.126,  -37.213,   22.369,  1.00,  74.84 ],
		[ 'O2',   'O',    -0.356,  -38.076,   21.552,  1.00,  72.68 ],
		[ 'N3',   'N',    -1.046,  -36.262,   22.747,  1.00,  74.84 ],
		[ 'C4',   'C',    -0.865,  -35.248,   23.676,  1.00,  74.32 ],
		[ 'O4',   'O',    -1.714,  -34.432,   23.980,  1.00,  72.54 ],
		[ 'C5',   'C',     0.419,  -35.203,   24.296,  1.00,  75.30 ],
		[ 'C7',   'C',     0.684,  -34.141,   25.313,  1.00,  75.98 ],
		[ 'C6',   'C',     1.329,  -36.122,   23.946,  1.00,  76.20 ],
	] ],
	[ 'DA', 5, [
		[ 'P',    'P',     5.053,  -38.899,   19.929,  1.00,  84.39 ],
		[ 'OP1',  'O',     6.458,  -39.144,   20.384,  1.00,  82.68 ],
		[ 'OP2',  'O',     4.597,  -37.529,   19.558,  1.00,  82.76 ],
		[ "O5'",  'O',     4.692,  -39.866,   18.681,  1.00,  82.05 ],
		[ "C5'",  'C',     3.391,  -40.489,   18.513,  1.00,  82.76 ],
		[ "C4'",  'C',     2.224,  -39.665,   17.908,  1.00,  82.63 ],
		[ "O4'",  'O',     1.803,  -38.578,   18.774,  1.00,  81.80 ],
		[ "C3'",  'C',     2.428,  -39.010,   16.531,  1.00,  83.07 ],
		[ "O3'",  'O',     1.394,  -39.351,   15.581,  1.00,  83.36 ],
		[ "C2'",  'C',     2.378,  -37.518,   16.834,  1.00,  83.13 ],
		[ "C1'",  'C',     1.322,  -37.578,   17.908,  1.00,  82.01 ],
		[ 'N9',   'N',     1.157,  -36.411,   18.728,  1.00,  81.07 ],
		[ 'C8',   'C',     2.102,  -35.766,   19.464,  1.00,  80.77 ],
		[ 'N7',   'N',     1.626,  -34.752,   20.133,  1.00,  80.74 ],
		[ 'C5',   'C',     0.282,  -34.761,   19.807,  1.00,  80.02 ],
		[ 'C6',   'C',    -0.775,  -33.935,   20.193,  1.00,  79.32 ],
		[ 'N6',   'N',    -0.607,  -32.906,   21.031,  1.00,  79.33 ],
		[ 'N1',   'N',    -1.992,  -34.207,   19.677,  1.00,  78.90 ],
		[ 'C2',   'C',    -2.134,  -35.235,   18.832,  1.00,  79.06 ],
		[ 'N3',   'N',    -1.217,  -36.088,   18.398,  1.00,  79.23 ],
		[ 'C4',   'C',    -0.027,  -35.778,   18.935,  1.00,  80.08 ],
	] ],
	[ 'DT', 6, [
		[ 'P',    'P',     1.616,  -39.002,   14.023,  1.00,  83.82 ],
		[ 'OP1',  'O',     1.063,  -40.106,   13.197,  1.00,  82.71 ],
		[ 'OP2',  'O',     3.008,  -38.483,   13.848,  1.00,  82.37 ],
		[ "O5'",  'O',     0.632,  -37.785,   13.842,  1.00,  83.12 ],
		[ "C5'",  'C',    -0.738,  -37.999,   14.025,  1.00,  84.23 ],
		[ "C4'",  'C',    -1.444,  -36.760,   13.535,  1.00,  84.61 ],
		[ "O4'",  'O',    -1.550,  -35.798,   14.638,  1.00,  83.76 ],
		[ "C3'",  'C',    -0.661,  -36.122,   12.384,  1.00,  84.60 ],
		[ "O3'",  'O',    -1.527,  -35.743,   11.285,  1.00,  85.81 ],
		[ "C2'",  'C',     0.042,  -34.965,   13.103,  1.00,  83.83 ],
		[ "C1'",  'C',    -0.941,  -34.593,   14.229,  1.00,  82.44 ],
		[ 'N1',   'N',    -0.358,  -33.819,   15.394,  1.00,  80.40 ],
		[ 'C2',   'C',    -1.218,  -33.092,   16.181,  1.00,  78.33 ],
		[ 'O2',   'O',    -2.410,  -33.061,   16.020,  1.00,  77.70 ],
		[ 'N3',   'N',    -0.631,  -32.391,   17.189,  1.00,  77.77 ],
		[ 'C4',   'C',     0.709,  -32.338,   17.497,  1.00,  78.12 ],
		[ 'O4',   'O',     1.149,  -31.679,   18.414,  1.00,  78.20 ],
		[ 'C5',   'C',     1.573,  -33.110,   16.668,  1.00,  79.26 ],
		[ 'C7',   'C',     3.047,  -33.095,   16.950,  1.00,  79.32 ],
		[ 'C6',   'C',     1.012,  -33.799,   15.659,  1.00,  80.45 ],
	] ],
);
my @bases;
$serial = 0;
for my $r (@bases_res) {
	my ($resname, $resseq, $atoms) = @$r;
	for my $a (@$atoms) {
		$serial++;
		push @bases, atom_line(
			record => 'ATOM  ', serial => $serial, name => $a->[0], element => $a->[1],
			altloc => '', resname => $resname, chain => 'P', resseq => $resseq, icode => '',
			x => $a->[2], y => $a->[3], z => $a->[4], occ => $a->[5], b => $a->[6],
		);
	}
}
push @bases, 'END';

# --- ss.pdb -- cysteines, two of them bonded and one not ------------------
#
# Five residues of chain A of 1AHW (a Fab bound to tissue factor), lifted whole
# from /home/con/ui/pepPriML/PPB/PDB/PDBbind.v2020/1ahw.ent.pdb with their
# coordinates as deposited.  Chain A was chosen by asking the module for a chain
# holding both cases at once:
#
#   CYS 23  - CYS 88   a disulfide, SG to SG 2.05 A
#   CYS 134 - CYS 194  a disulfide, SG to SG 2.03 A
#   CYS 214            free: its SG is nowhere near another
#
# so the file exercises the rule in both directions without inventing a
# geometry.  The entry declares all four bonds in SSBOND records, which is the
# other half of the comparison: t/features.t checks the computed answer and
# t/real.t checks it against what files say about themselves.
my @ss_res = (
	[ 'CYS', 23, [
		[ 'N',    'N',     7.801,    3.335,   46.174,  1.00,  14.71 ],
		[ 'CA',   'C',     6.431,    3.744,   46.278,  1.00,  14.71 ],
		[ 'C',    'C',     5.941,    3.493,   47.682,  1.00,  14.71 ],
		[ 'O',    'O',     6.073,    2.342,   48.178,  1.00,  13.80 ],
		[ 'CB',   'C',     5.573,    2.937,   45.308,  1.00,  13.80 ],
		[ 'SG',   'S',     5.378,    3.900,   43.821,  1.00,  13.80 ],
	] ],
	[ 'CYS', 88, [
		[ 'N',    'N',     4.046,    2.898,   39.606,  1.00,   5.78 ],
		[ 'CA',   'C',     3.522,    1.876,   40.496,  1.00,   5.78 ],
		[ 'C',    'C',     2.068,    2.241,   40.462,  1.00,   5.78 ],
		[ 'O',    'O',     1.714,    3.408,   40.321,  1.00,   2.48 ],
		[ 'CB',   'C',     4.088,    1.955,   41.921,  1.00,   2.48 ],
		[ 'SG',   'S',     3.641,    3.373,   42.956,  1.00,   2.48 ],
	] ],
	[ 'CYS', 134, [
		[ 'N',    'N',    19.022,   20.078,    9.689,  1.00,  38.85 ],
		[ 'CA',   'C',    20.346,   19.673,   10.165,  1.00,  38.85 ],
		[ 'C',    'C',    19.961,   18.522,   11.113,  1.00,  38.85 ],
		[ 'O',    'O',    19.020,   18.665,   11.878,  1.00,  44.13 ],
		[ 'CB',   'C',    21.047,   20.855,   10.909,  1.00,  44.13 ],
		[ 'SG',   'S',    22.782,   20.725,   11.552,  1.00,  44.13 ],
	] ],
	[ 'CYS', 194, [
		[ 'N',    'N',    24.887,   24.953,    9.823,  1.00,  45.40 ],
		[ 'CA',   'C',    24.569,   24.862,   11.247,  1.00,  45.40 ],
		[ 'C',    'C',    25.827,   24.465,   12.010,  1.00,  45.40 ],
		[ 'O',    'O',    26.671,   23.731,   11.474,  1.00,  38.27 ],
		[ 'CB',   'C',    23.461,   23.822,   11.501,  1.00,  38.27 ],
		[ 'SG',   'S',    23.832,   22.248,   10.689,  1.00,  38.27 ],
	] ],
	[ 'CYS', 214, [
		[ 'N',    'N',    22.758,   26.326,  -10.486,  1.00,  88.86 ],
		[ 'CA',   'C',    23.222,   27.025,  -11.669,  1.00,  88.86 ],
		[ 'C',    'C',    24.024,   25.987,  -12.439,  1.00,  88.86 ],
		[ 'O',    'O',    25.290,   26.053,  -12.471,  1.00,  91.66 ],
		[ 'CB',   'C',    24.081,   28.225,  -11.227,  1.00,  91.66 ],
		[ 'SG',   'S',    23.724,   28.735,   -9.458,  1.00,  91.66 ],
		[ 'OXT',  'O',    23.323,   25.062,  -12.929,  1.00,  91.66 ],
	] ],
);
my @ss;
$serial = 0;
for my $r (@ss_res) {
	my ($resname, $resseq, $atoms) = @$r;
	for my $a (@$atoms) {
		$serial++;
		push @ss, atom_line(
			record => 'ATOM  ', serial => $serial, name => $a->[0], element => $a->[1],
			altloc => '', resname => $resname, chain => 'A', resseq => $resseq, icode => '',
			x => $a->[2], y => $a->[3], z => $a->[4], occ => $a->[5], b => $a->[6],
		);
	}
}
push @ss, 'END';

# --- fold.pdb -- a stretch of protein that is actually folded ---------------
#
# Residues 3 to 62 of chain A of 1A42 (human carbonic anhydrase II, 1.05 A),
# lifted whole from /home/con/ui/pepPriML/PPB/PDB/PDBbind.v2020/1a42.ent.pdb
# with their coordinates as deposited.
#
# The other fixtures are a few residues each, chosen for one question apiece,
# and between them they contain almost no secondary structure and almost no
# backbone hydrogen bonds -- so the parts of the module that read those would go
# untested against mdtraj.  This is a contiguous run long enough to hold an
# alpha helix, a 3-10 helix, two strands and the turns between them, which is
# what t/features.t compares residue by residue.
#
# Sixty consecutive residues, so the phi and psi angles have real neighbours and
# the peptide-bond test has something to say yes to.
my @fold_res = (
	[ 'HIS', 4, [
		[ 'N',    'N',    10.662,   -0.401,    7.977,  1.00,  38.25 ],
		[ 'CA',   'C',    10.127,    0.268,    9.147,  1.00,  36.18 ],
		[ 'C',    'C',     9.137,   -0.769,    9.692,  1.00,  33.60 ],
		[ 'O',    'O',     9.181,   -1.906,    9.193,  1.00,  32.41 ],
		[ 'CB',   'C',     9.415,    1.580,    8.759,  1.00,  40.50 ],
		[ 'CG',   'C',     9.808,    2.795,    9.631,  1.00,  46.39 ],
		[ 'ND1',  'N',     9.962,    2.952,   10.959,  1.00,  46.24 ],
		[ 'CD2',  'C',    10.076,    4.036,    9.071,  1.00,  48.38 ],
		[ 'CE1',  'C',    10.304,    4.208,   11.202,  1.00,  48.85 ],
		[ 'NE2',  'N',    10.369,    4.845,   10.060,  1.00,  49.86 ],
	] ],
	[ 'TRP', 5, [
		[ 'N',    'N',     8.257,   -0.414,   10.635,  1.00,  27.01 ],
		[ 'CA',   'C',     7.481,   -1.387,   11.384,  1.00,  19.87 ],
		[ 'C',    'C',     6.510,   -2.204,   10.563,  1.00,  19.91 ],
		[ 'O',    'O',     6.066,   -1.770,    9.486,  1.00,  19.19 ],
		[ 'CB',   'C',     6.719,   -0.667,   12.505,  1.00,  14.02 ],
		[ 'CG',   'C',     5.438,    0.076,   12.142,  1.00,  12.29 ],
		[ 'CD1',  'C',     5.444,    1.400,   11.760,  1.00,   6.38 ],
		[ 'CD2',  'C',     4.172,   -0.466,   12.199,  1.00,   6.22 ],
		[ 'NE1',  'N',     4.172,    1.704,   11.581,  1.00,   8.90 ],
		[ 'CE2',  'C',     3.380,    0.633,   11.826,  1.00,   7.53 ],
		[ 'CE3',  'C',     3.588,   -1.693,   12.512,  1.00,   4.86 ],
		[ 'CZ2',  'C',     1.995,    0.506,   11.767,  1.00,   3.51 ],
		[ 'CZ3',  'C',     2.199,   -1.821,   12.451,  1.00,   2.00 ],
		[ 'CH2',  'C',     1.408,   -0.726,   12.081,  1.00,   4.15 ],
	] ],
	[ 'GLY', 6, [
		[ 'N',    'N',     6.158,   -3.341,   11.144,  1.00,  13.85 ],
		[ 'CA',   'C',     5.208,   -4.220,   10.536,  1.00,  12.88 ],
		[ 'C',    'C',     4.962,   -5.333,   11.521,  1.00,  14.39 ],
		[ 'O',    'O',     5.048,   -5.127,   12.727,  1.00,  14.80 ],
	] ],
	[ 'TYR', 7, [
		[ 'N',    'N',     4.655,   -6.518,   10.998,  1.00,  13.15 ],
		[ 'CA',   'C',     4.314,   -7.670,   11.805,  1.00,  14.72 ],
		[ 'C',    'C',     5.170,   -8.832,   11.402,  1.00,  14.17 ],
		[ 'O',    'O',     4.839,   -9.977,   11.668,  1.00,  18.90 ],
		[ 'CB',   'C',     2.852,   -8.046,   11.617,  1.00,  10.37 ],
		[ 'CG',   'C',     1.936,   -6.922,   12.029,  1.00,  13.25 ],
		[ 'CD1',  'C',     1.594,   -6.749,   13.379,  1.00,  11.77 ],
		[ 'CD2',  'C',     1.474,   -6.048,   11.048,  1.00,   9.56 ],
		[ 'CE1',  'C',     0.778,   -5.681,   13.754,  1.00,  13.85 ],
		[ 'CE2',  'C',     0.662,   -4.981,   11.411,  1.00,  14.02 ],
		[ 'CZ',   'C',     0.323,   -4.809,   12.759,  1.00,  15.07 ],
		[ 'OH',   'O',    -0.505,   -3.754,   13.077,  1.00,  19.57 ],
	] ],
	[ 'GLY', 8, [
		[ 'N',    'N',     6.301,   -8.588,   10.760,  1.00,  18.22 ],
		[ 'CA',   'C',     7.129,   -9.688,   10.305,  1.00,  20.53 ],
		[ 'C',    'C',     8.175,   -9.968,   11.360,  1.00,  23.45 ],
		[ 'O',    'O',     8.411,   -9.142,   12.250,  1.00,  25.40 ],
	] ],
	[ 'LYS', 9, [
		[ 'N',    'N',     8.886,  -11.088,   11.219,  1.00,  23.74 ],
		[ 'CA',   'C',     9.933,  -11.465,   12.148,  1.00,  25.06 ],
		[ 'C',    'C',    10.881,  -10.309,   12.438,  1.00,  25.45 ],
		[ 'O',    'O',    11.251,  -10.004,   13.579,  1.00,  23.26 ],
		[ 'CB',   'C',    10.716,  -12.653,   11.567,  1.00,  29.42 ],
		[ 'CG',   'C',    11.882,  -13.121,   12.465,  1.00,  34.54 ],
		[ 'CD',   'C',    12.685,  -14.290,   11.892,  1.00,  38.37 ],
		[ 'CE',   'C',    11.984,  -15.633,   12.082,  1.00,  43.35 ],
		[ 'NZ',   'N',    12.002,  -16.063,   13.475,  1.00,  42.96 ],
	] ],
	[ 'HIS', 10, [
		[ 'N',    'N',    11.187,   -9.576,   11.366,  1.00,  22.56 ],
		[ 'CA',   'C',    12.216,   -8.572,   11.501,  1.00,  22.56 ],
		[ 'C',    'C',    11.741,   -7.221,   11.979,  1.00,  20.04 ],
		[ 'O',    'O',    12.529,   -6.480,   12.563,  1.00,  20.12 ],
		[ 'CB',   'C',    12.954,   -8.470,   10.152,  1.00,  25.27 ],
		[ 'CG',   'C',    13.652,   -9.804,    9.857,  1.00,  28.90 ],
		[ 'ND1',  'N',    13.292,  -10.778,    9.010,  1.00,  33.34 ],
		[ 'CD2',  'C',    14.755,  -10.260,   10.542,  1.00,  29.87 ],
		[ 'CE1',  'C',    14.108,  -11.793,    9.170,  1.00,  32.13 ],
		[ 'NE2',  'N',    14.980,  -11.463,   10.092,  1.00,  33.23 ],
	] ],
	[ 'ASN', 11, [
		[ 'N',    'N',    10.459,   -6.922,   11.899,  1.00,  16.88 ],
		[ 'CA',   'C',     9.987,   -5.588,   12.216,  1.00,  15.73 ],
		[ 'C',    'C',     8.678,   -5.658,   12.978,  1.00,  18.45 ],
		[ 'O',    'O',     7.879,   -4.709,   12.945,  1.00,  19.79 ],
		[ 'CB',   'C',     9.790,   -4.824,   10.923,  1.00,  14.73 ],
		[ 'CG',   'C',     8.772,   -5.490,    9.992,  1.00,  16.61 ],
		[ 'OD1',  'O',     8.222,   -6.561,   10.271,  1.00,  13.34 ],
		[ 'ND2',  'N',     8.443,   -4.895,    8.853,  1.00,  14.29 ],
	] ],
	[ 'GLY', 12, [
		[ 'N',    'N',     8.440,   -6.813,   13.616,  1.00,  17.60 ],
		[ 'CA',   'C',     7.215,   -7.040,   14.359,  1.00,  16.90 ],
		[ 'C',    'C',     7.196,   -6.398,   15.736,  1.00,  17.50 ],
		[ 'O',    'O',     8.172,   -5.745,   16.140,  1.00,  15.75 ],
	] ],
	[ 'PRO', 13, [
		[ 'N',    'N',     6.141,   -6.656,   16.529,  1.00,  17.75 ],
		[ 'CA',   'C',     5.915,   -6.089,   17.867,  1.00,  17.55 ],
		[ 'C',    'C',     7.137,   -6.046,   18.796,  1.00,  15.98 ],
		[ 'O',    'O',     7.536,   -5.024,   19.357,  1.00,  14.52 ],
		[ 'CB',   'C',     4.789,   -6.940,   18.400,  1.00,  19.63 ],
		[ 'CG',   'C',     3.965,   -7.296,   17.183,  1.00,  17.16 ],
		[ 'CD',   'C',     5.070,   -7.596,   16.183,  1.00,  17.36 ],
	] ],
	[ 'GLU', 14, [
		[ 'N',    'N',     7.845,   -7.157,   18.850,  1.00,  12.59 ],
		[ 'CA',   'C',     9.054,   -7.327,   19.661,  1.00,  14.57 ],
		[ 'C',    'C',    10.158,   -6.309,   19.385,  1.00,  12.09 ],
		[ 'O',    'O',    11.141,   -6.198,   20.118,  1.00,  15.74 ],
		[ 'CB',   'C',     9.676,   -8.718,   19.446,  1.00,  19.20 ],
		[ 'CG',   'C',     8.851,   -9.870,   18.837,  1.00,  30.27 ],
		[ 'CD',   'C',     8.431,   -9.706,   17.360,  1.00,  36.96 ],
		[ 'OE1',  'O',     9.268,   -9.330,   16.531,  1.00,  40.80 ],
		[ 'OE2',  'O',     7.255,   -9.933,   17.049,  1.00,  36.61 ],
	] ],
	[ 'HIS', 15, [
		[ 'N',    'N',    10.048,   -5.570,   18.287,  1.00,  10.66 ],
		[ 'CA',   'C',    11.036,   -4.612,   17.855,  1.00,   9.85 ],
		[ 'C',    'C',    10.496,   -3.220,   17.995,  1.00,   8.80 ],
		[ 'O',    'O',    11.277,   -2.274,   18.088,  1.00,  12.25 ],
		[ 'CB',   'C',    11.397,   -4.771,   16.384,  1.00,   9.92 ],
		[ 'CG',   'C',    12.137,   -6.042,   16.030,  1.00,  14.18 ],
		[ 'ND1',  'N',    13.443,   -6.287,   16.030,  1.00,  15.65 ],
		[ 'CD2',  'C',    11.495,   -7.186,   15.628,  1.00,  15.97 ],
		[ 'CE1',  'C',    13.612,   -7.530,   15.647,  1.00,  18.19 ],
		[ 'NE2',  'N',    12.433,   -8.059,   15.409,  1.00,  17.69 ],
	] ],
	[ 'TRP', 16, [
		[ 'N',    'N',     9.192,   -3.019,   18.090,  1.00,  12.01 ],
		[ 'CA',   'C',     8.643,   -1.669,   17.993,  1.00,  13.38 ],
		[ 'C',    'C',     9.204,   -0.743,   19.068,  1.00,  15.29 ],
		[ 'O',    'O',     9.472,    0.434,   18.809,  1.00,  19.52 ],
		[ 'CB',   'C',     7.113,   -1.689,   18.107,  1.00,   7.18 ],
		[ 'CG',   'C',     6.326,   -2.385,   16.999,  1.00,  11.80 ],
		[ 'CD1',  'C',     6.891,   -2.781,   15.810,  1.00,  12.22 ],
		[ 'CD2',  'C',     4.991,   -2.715,   17.073,  1.00,   8.99 ],
		[ 'NE1',  'N',     5.924,   -3.374,   15.144,  1.00,  11.17 ],
		[ 'CE2',  'C',     4.781,   -3.361,   15.847,  1.00,   9.47 ],
		[ 'CE3',  'C',     3.947,   -2.571,   17.984,  1.00,   5.18 ],
		[ 'CZ2',  'C',     3.524,   -3.874,   15.521,  1.00,   4.07 ],
		[ 'CZ3',  'C',     2.695,   -3.083,   17.666,  1.00,   4.69 ],
		[ 'CH2',  'C',     2.484,   -3.732,   16.443,  1.00,   8.68 ],
	] ],
	[ 'HIS', 17, [
		[ 'N',    'N',     9.521,   -1.296,   20.248,  1.00,  14.80 ],
		[ 'CA',   'C',    10.020,   -0.491,   21.357,  1.00,  13.57 ],
		[ 'C',    'C',    11.293,    0.264,   20.979,  1.00,  13.18 ],
		[ 'O',    'O',    11.600,    1.269,   21.608,  1.00,  15.46 ],
		[ 'CB',   'C',    10.307,   -1.375,   22.598,  1.00,   9.49 ],
		[ 'CG',   'C',    11.572,   -2.231,   22.491,  1.00,   9.14 ],
		[ 'ND1',  'N',    12.799,   -1.928,   22.917,  1.00,   6.36 ],
		[ 'CD2',  'C',    11.640,   -3.481,   21.898,  1.00,  11.19 ],
		[ 'CE1',  'C',    13.596,   -2.935,   22.613,  1.00,   6.57 ],
		[ 'NE2',  'N',    12.892,   -3.859,   22.001,  1.00,  10.04 ],
	] ],
	[ 'LYS', 18, [
		[ 'N',    'N',    12.061,   -0.204,   19.981,  1.00,  12.69 ],
		[ 'CA',   'C',    13.297,    0.459,   19.611,  1.00,  12.08 ],
		[ 'C',    'C',    13.027,    1.802,   18.954,  1.00,  13.26 ],
		[ 'O',    'O',    13.869,    2.688,   18.988,  1.00,  15.52 ],
		[ 'CB',   'C',    14.078,   -0.439,   18.676,  1.00,  10.72 ],
		[ 'CG',   'C',    14.772,   -1.608,   19.367,  1.00,  10.61 ],
		[ 'CD',   'C',    15.545,   -2.491,   18.406,  1.00,  12.27 ],
		[ 'CE',   'C',    14.603,   -3.305,   17.519,  1.00,  12.87 ],
		[ 'NZ',   'N',    15.340,   -4.109,   16.552,  1.00,   7.96 ],
	] ],
	[ 'ASP', 19, [
		[ 'N',    'N',    11.864,    2.001,   18.350,  1.00,  14.87 ],
		[ 'CA',   'C',    11.558,    3.285,   17.742,  1.00,  17.22 ],
		[ 'C',    'C',    10.411,    3.977,   18.447,  1.00,  15.76 ],
		[ 'O',    'O',    10.266,    5.196,   18.361,  1.00,  12.15 ],
		[ 'CB',   'C',    11.183,    3.135,   16.271,  1.00,  23.23 ],
		[ 'CG',   'C',    12.386,    2.907,   15.369,  1.00,  25.80 ],
		[ 'OD1',  'O',    13.246,    3.791,   15.239,  1.00,  27.96 ],
		[ 'OD2',  'O',    12.447,    1.817,   14.807,  1.00,  32.22 ],
	] ],
	[ 'PHE', 20, [
		[ 'N',    'N',     9.558,    3.231,   19.147,  1.00,  13.15 ],
		[ 'CA',   'C',     8.472,    3.829,   19.884,  1.00,  11.47 ],
		[ 'C',    'C',     8.672,    3.273,   21.283,  1.00,  13.90 ],
		[ 'O',    'O',     8.096,    2.240,   21.663,  1.00,  12.69 ],
		[ 'CB',   'C',     7.126,    3.403,   19.286,  1.00,  13.66 ],
		[ 'CG',   'C',     7.063,    3.667,   17.790,  1.00,  12.24 ],
		[ 'CD1',  'C',     6.793,    4.955,   17.320,  1.00,  16.22 ],
		[ 'CD2',  'C',     7.331,    2.625,   16.895,  1.00,  15.45 ],
		[ 'CE1',  'C',     6.804,    5.202,   15.935,  1.00,  18.11 ],
		[ 'CE2',  'C',     7.338,    2.878,   15.518,  1.00,  19.04 ],
		[ 'CZ',   'C',     7.077,    4.165,   15.033,  1.00,  16.21 ],
	] ],
	[ 'PRO', 21, [
		[ 'N',    'N',     9.514,    3.927,   22.101,  1.00,  12.65 ],
		[ 'CA',   'C',     9.817,    3.501,   23.455,  1.00,  10.98 ],
		[ 'C',    'C',     8.570,    3.332,   24.320,  1.00,  13.87 ],
		[ 'O',    'O',     8.535,    2.486,   25.221,  1.00,  17.59 ],
		[ 'CB',   'C',    10.761,    4.568,   23.930,  1.00,  13.45 ],
		[ 'CG',   'C',    11.475,    5.022,   22.685,  1.00,  12.57 ],
		[ 'CD',   'C',    10.302,    5.098,   21.736,  1.00,   9.93 ],
	] ],
	[ 'ILE', 22, [
		[ 'N',    'N',     7.475,    4.063,   24.038,  1.00,  13.17 ],
		[ 'CA',   'C',     6.229,    3.953,   24.792,  1.00,  11.68 ],
		[ 'C',    'C',     5.631,    2.545,   24.740,  1.00,  13.35 ],
		[ 'O',    'O',     4.730,    2.222,   25.513,  1.00,  15.23 ],
		[ 'CB',   'C',     5.239,    5.022,   24.238,  1.00,  10.83 ],
		[ 'CG1',  'C',     3.986,    5.066,   25.085,  1.00,  11.88 ],
		[ 'CG2',  'C',     4.806,    4.697,   22.806,  1.00,   8.53 ],
		[ 'CD1',  'C',     4.241,    5.674,   26.448,  1.00,  16.93 ],
	] ],
	[ 'ALA', 23, [
		[ 'N',    'N',     6.137,    1.680,   23.852,  1.00,  14.54 ],
		[ 'CA',   'C',     5.714,    0.293,   23.762,  1.00,  14.78 ],
		[ 'C',    'C',     5.860,   -0.437,   25.083,  1.00,  14.86 ],
		[ 'O',    'O',     5.239,   -1.475,   25.277,  1.00,  15.22 ],
		[ 'CB',   'C',     6.545,   -0.472,   22.742,  1.00,  14.19 ],
	] ],
	[ 'LYS', 24, [
		[ 'N',    'N',     6.717,    0.092,   25.971,  1.00,  15.11 ],
		[ 'CA',   'C',     7.018,   -0.488,   27.274,  1.00,  17.68 ],
		[ 'C',    'C',     6.320,    0.277,   28.407,  1.00,  18.95 ],
		[ 'O',    'O',     6.667,    0.157,   29.580,  1.00,  21.32 ],
		[ 'CB',   'C',     8.511,   -0.417,   27.504,  1.00,  16.86 ],
		[ 'CG',   'C',     9.430,   -1.013,   26.486,  1.00,  16.77 ],
		[ 'CD',   'C',     9.456,   -2.501,   26.649,  1.00,  15.76 ],
		[ 'CE',   'C',    10.920,   -2.860,   26.767,  1.00,  20.15 ],
		[ 'NZ',   'N',    11.088,   -4.254,   26.397,  1.00,  27.00 ],
	] ],
	[ 'GLY', 25, [
		[ 'N',    'N',     5.369,    1.141,   28.089,  1.00,  22.15 ],
		[ 'CA',   'C',     4.742,    2.003,   29.070,  1.00,  18.64 ],
		[ 'C',    'C',     3.697,    1.302,   29.917,  1.00,  16.40 ],
		[ 'O',    'O',     3.393,    0.120,   29.799,  1.00,  16.03 ],
	] ],
	[ 'GLU', 26, [
		[ 'N',    'N',     3.068,    2.189,   30.672,  1.00,  17.49 ],
		[ 'CA',   'C',     2.087,    1.853,   31.680,  1.00,  20.28 ],
		[ 'C',    'C',     0.665,    1.561,   31.249,  1.00,  18.33 ],
		[ 'O',    'O',    -0.078,    1.000,   32.038,  1.00,  17.31 ],
		[ 'CB',   'C',     2.064,    2.974,   32.711,  1.00,  25.93 ],
		[ 'CG',   'C',     3.311,    3.025,   33.601,  1.00,  38.03 ],
		[ 'CD',   'C',     3.589,    1.727,   34.363,  1.00,  44.34 ],
		[ 'OE1',  'O',     2.744,    1.281,   35.150,  1.00,  50.25 ],
		[ 'OE2',  'O',     4.658,    1.152,   34.150,  1.00,  49.34 ],
	] ],
	[ 'ARG', 27, [
		[ 'N',    'N',     0.174,    1.998,   30.095,  1.00,  14.96 ],
		[ 'CA',   'C',    -1.168,    1.602,   29.716,  1.00,  12.83 ],
		[ 'C',    'C',    -1.229,    1.244,   28.232,  1.00,  12.50 ],
		[ 'O',    'O',    -1.909,    1.825,   27.394,  1.00,  10.33 ],
		[ 'CB',   'C',    -2.168,    2.727,   30.087,  1.00,  10.19 ],
		[ 'CG',   'C',    -1.782,    4.124,   29.684,  1.00,  12.21 ],
		[ 'CD',   'C',    -2.913,    5.117,   29.883,  1.00,  14.90 ],
		[ 'NE',   'N',    -3.940,    5.079,   28.852,  1.00,   9.72 ],
		[ 'CZ',   'C',    -4.669,    6.166,   28.583,  1.00,  11.81 ],
		[ 'NH1',  'N',    -4.498,    7.316,   29.225,  1.00,  11.52 ],
		[ 'NH2',  'N',    -5.625,    6.103,   27.671,  1.00,  12.77 ],
	] ],
	[ 'GLN', 28, [
		[ 'N',    'N',    -0.433,    0.233,   27.888,  1.00,  12.99 ],
		[ 'CA',   'C',    -0.358,   -0.290,   26.534,  1.00,   9.65 ],
		[ 'C',    'C',    -1.472,   -1.290,   26.313,  1.00,  10.91 ],
		[ 'O',    'O',    -1.915,   -1.978,   27.240,  1.00,  12.30 ],
		[ 'CB',   'C',     1.002,   -0.959,   26.312,  1.00,   7.47 ],
		[ 'CG',   'C',     2.125,    0.088,   26.175,  1.00,   4.90 ],
		[ 'CD',   'C',     2.000,    0.860,   24.878,  1.00,  11.03 ],
		[ 'OE1',  'O',     2.178,    0.296,   23.798,  1.00,  12.48 ],
		[ 'NE2',  'N',     1.615,    2.131,   24.878,  1.00,  11.49 ],
	] ],
	[ 'SER', 29, [
		[ 'N',    'N',    -1.953,   -1.306,   25.071,  1.00,  11.47 ],
		[ 'CA',   'C',    -2.997,   -2.193,   24.590,  1.00,   8.96 ],
		[ 'C',    'C',    -2.462,   -3.139,   23.520,  1.00,   6.80 ],
		[ 'O',    'O',    -1.450,   -2.819,   22.897,  1.00,   5.78 ],
		[ 'CB',   'C',    -4.133,   -1.356,   24.022,  1.00,   8.64 ],
		[ 'OG',   'O',    -4.770,   -0.698,   25.101,  1.00,   4.94 ],
	] ],
	[ 'PRO', 30, [
		[ 'N',    'N',    -3.076,   -4.286,   23.225,  1.00,   5.57 ],
		[ 'CA',   'C',    -4.273,   -4.820,   23.881,  1.00,   7.55 ],
		[ 'C',    'C',    -3.953,   -5.501,   25.237,  1.00,   7.66 ],
		[ 'O',    'O',    -2.782,   -5.540,   25.607,  1.00,   6.22 ],
		[ 'CB',   'C',    -4.809,   -5.731,   22.783,  1.00,   3.72 ],
		[ 'CG',   'C',    -3.544,   -6.410,   22.336,  1.00,   3.03 ],
		[ 'CD',   'C',    -2.606,   -5.219,   22.202,  1.00,   2.00 ],
	] ],
	[ 'VAL', 31, [
		[ 'N',    'N',    -4.919,   -6.046,   26.005,  1.00,   8.25 ],
		[ 'CA',   'C',    -4.700,   -6.743,   27.269,  1.00,   3.84 ],
		[ 'C',    'C',    -5.666,   -7.933,   27.243,  1.00,   4.86 ],
		[ 'O',    'O',    -6.559,   -8.055,   26.389,  1.00,   4.29 ],
		[ 'CB',   'C',    -5.035,   -5.876,   28.561,  1.00,   5.08 ],
		[ 'CG1',  'C',    -4.212,   -4.589,   28.654,  1.00,   2.00 ],
		[ 'CG2',  'C',    -6.500,   -5.516,   28.531,  1.00,   2.21 ],
	] ],
	[ 'ASP', 32, [
		[ 'N',    'N',    -5.522,   -8.816,   28.221,  1.00,   9.02 ],
		[ 'CA',   'C',    -6.410,   -9.950,   28.378,  1.00,  10.39 ],
		[ 'C',    'C',    -7.555,   -9.471,   29.248,  1.00,   9.34 ],
		[ 'O',    'O',    -7.301,   -8.825,   30.260,  1.00,  11.79 ],
		[ 'CB',   'C',    -5.653,  -11.099,   29.058,  1.00,  10.29 ],
		[ 'CG',   'C',    -6.494,  -12.327,   29.398,  1.00,  12.11 ],
		[ 'OD1',  'O',    -7.357,  -12.705,   28.614,  1.00,  13.91 ],
		[ 'OD2',  'O',    -6.283,  -12.923,   30.447,  1.00,  15.55 ],
	] ],
	[ 'ILE', 33, [
		[ 'N',    'N',    -8.810,   -9.670,   28.894,  1.00,  10.18 ],
		[ 'CA',   'C',    -9.889,   -9.355,   29.808,  1.00,  12.61 ],
		[ 'C',    'C',   -10.097,  -10.631,   30.646,  1.00,  14.49 ],
		[ 'O',    'O',   -10.582,  -11.660,   30.141,  1.00,  12.79 ],
		[ 'CB',   'C',   -11.218,   -8.988,   29.061,  1.00,  10.37 ],
		[ 'CG1',  'C',   -11.207,   -7.569,   28.471,  1.00,  13.21 ],
		[ 'CG2',  'C',   -12.344,   -8.956,   30.070,  1.00,   8.33 ],
		[ 'CD1',  'C',   -10.361,   -7.329,   27.201,  1.00,  11.78 ],
	] ],
	[ 'ASP', 34, [
		[ 'N',    'N',    -9.672,  -10.646,   31.919,  1.00,  17.43 ],
		[ 'CA',   'C',    -9.996,  -11.764,   32.806,  1.00,  18.67 ],
		[ 'C',    'C',   -11.372,  -11.480,   33.416,  1.00,  19.66 ],
		[ 'O',    'O',   -11.575,  -10.584,   34.249,  1.00,  19.01 ],
		[ 'CB',   'C',    -8.971,  -11.915,   33.930,  1.00,  23.16 ],
		[ 'CG',   'C',    -9.308,  -13.049,   34.909,  1.00,  25.79 ],
		[ 'OD1',  'O',    -9.715,  -14.160,   34.528,  1.00,  28.60 ],
		[ 'OD2',  'O',    -9.160,  -12.794,   36.087,  1.00,  30.58 ],
	] ],
	[ 'THR', 35, [
		[ 'N',    'N',   -12.333,  -12.298,   33.028,  1.00,  18.43 ],
		[ 'CA',   'C',   -13.701,  -12.046,   33.381,  1.00,  22.62 ],
		[ 'C',    'C',   -13.956,  -12.109,   34.861,  1.00,  24.44 ],
		[ 'O',    'O',   -14.749,  -11.348,   35.396,  1.00,  32.28 ],
		[ 'CB',   'C',   -14.580,  -13.044,   32.614,  1.00,  23.04 ],
		[ 'OG1',  'O',   -14.178,  -14.378,   32.953,  1.00,  21.47 ],
		[ 'CG2',  'C',   -14.485,  -12.765,   31.108,  1.00,  22.29 ],
	] ],
	[ 'HIS', 36, [
		[ 'N',    'N',   -13.239,  -12.983,   35.551,  1.00,  28.79 ],
		[ 'CA',   'C',   -13.424,  -13.164,   36.981,  1.00,  26.47 ],
		[ 'C',    'C',   -12.837,  -12.042,   37.801,  1.00,  26.90 ],
		[ 'O',    'O',   -13.249,  -11.905,   38.942,  1.00,  31.19 ],
		[ 'CB',   'C',   -12.807,  -14.483,   37.393,  1.00,  24.14 ],
		[ 'CG',   'C',   -13.589,  -15.640,   36.797,  1.00,  26.35 ],
		[ 'ND1',  'N',   -14.872,  -15.927,   37.023,  1.00,  26.71 ],
		[ 'CD2',  'C',   -13.101,  -16.545,   35.874,  1.00,  26.61 ],
		[ 'CE1',  'C',   -15.188,  -16.965,   36.275,  1.00,  28.08 ],
		[ 'NE2',  'N',   -14.119,  -17.322,   35.592,  1.00,  28.35 ],
	] ],
	[ 'THR', 37, [
		[ 'N',    'N',   -11.904,  -11.220,   37.347,  1.00,  26.10 ],
		[ 'CA',   'C',   -11.443,  -10.180,   38.231,  1.00,  25.86 ],
		[ 'C',    'C',   -11.847,   -8.806,   37.742,  1.00,  23.68 ],
		[ 'O',    'O',   -11.373,   -7.792,   38.271,  1.00,  26.34 ],
		[ 'CB',   'C',    -9.899,  -10.285,   38.393,  1.00,  27.18 ],
		[ 'OG1',  'O',    -9.278,  -10.207,   37.114,  1.00,  26.39 ],
		[ 'CG2',  'C',    -9.526,  -11.588,   39.096,  1.00,  26.41 ],
	] ],
	[ 'ALA', 38, [
		[ 'N',    'N',   -12.692,   -8.716,   36.724,  1.00,  23.18 ],
		[ 'CA',   'C',   -13.076,   -7.415,   36.203,  1.00,  20.38 ],
		[ 'C',    'C',   -14.146,   -6.896,   37.130,  1.00,  18.52 ],
		[ 'O',    'O',   -15.052,   -7.609,   37.546,  1.00,  17.34 ],
		[ 'CB',   'C',   -13.655,   -7.531,   34.806,  1.00,  21.62 ],
	] ],
	[ 'LYS', 39, [
		[ 'N',    'N',   -14.099,   -5.633,   37.455,  1.00,  18.58 ],
		[ 'CA',   'C',   -15.009,   -5.108,   38.430,  1.00,  21.35 ],
		[ 'C',    'C',   -16.302,   -4.669,   37.788,  1.00,  23.40 ],
		[ 'O',    'O',   -16.278,   -3.785,   36.928,  1.00,  22.41 ],
		[ 'CB',   'C',   -14.303,   -3.961,   39.104,  1.00,  23.23 ],
		[ 'CG',   'C',   -15.089,   -3.231,   40.152,  1.00,  29.50 ],
		[ 'CD',   'C',   -15.262,   -4.062,   41.408,  1.00,  35.46 ],
		[ 'CE',   'C',   -15.902,   -3.182,   42.477,  1.00,  37.51 ],
		[ 'NZ',   'N',   -15.032,   -2.069,   42.833,  1.00,  43.00 ],
	] ],
	[ 'TYR', 40, [
		[ 'N',    'N',   -17.445,   -5.240,   38.157,  1.00,  24.76 ],
		[ 'CA',   'C',   -18.711,   -4.713,   37.677,  1.00,  28.33 ],
		[ 'C',    'C',   -18.904,   -3.319,   38.273,  1.00,  30.75 ],
		[ 'O',    'O',   -18.877,   -3.156,   39.501,  1.00,  34.63 ],
		[ 'CB',   'C',   -19.866,   -5.648,   38.095,  1.00,  31.12 ],
		[ 'CG',   'C',   -21.229,   -4.994,   37.936,  1.00,  34.15 ],
		[ 'CD1',  'C',   -21.694,   -4.585,   36.683,  1.00,  34.89 ],
		[ 'CD2',  'C',   -21.967,   -4.708,   39.082,  1.00,  36.94 ],
		[ 'CE1',  'C',   -22.884,   -3.877,   36.577,  1.00,  35.53 ],
		[ 'CE2',  'C',   -23.159,   -3.994,   38.979,  1.00,  40.03 ],
		[ 'CZ',   'C',   -23.605,   -3.580,   37.728,  1.00,  37.87 ],
		[ 'OH',   'O',   -24.756,   -2.824,   37.661,  1.00,  42.12 ],
	] ],
	[ 'ASP', 41, [
		[ 'N',    'N',   -19.108,   -2.304,   37.441,  1.00,  31.55 ],
		[ 'CA',   'C',   -19.278,   -0.964,   37.931,  1.00,  33.21 ],
		[ 'C',    'C',   -20.718,   -0.595,   37.607,  1.00,  35.68 ],
		[ 'O',    'O',   -21.095,   -0.621,   36.434,  1.00,  39.64 ],
		[ 'CB',   'C',   -18.270,   -0.054,   37.229,  1.00,  31.82 ],
		[ 'CG',   'C',   -18.123,    1.355,   37.812,  1.00,  36.42 ],
		[ 'OD1',  'O',   -19.049,    1.858,   38.443,  1.00,  37.61 ],
		[ 'OD2',  'O',   -17.076,    1.980,   37.630,  1.00,  38.27 ],
	] ],
	[ 'PRO', 42, [
		[ 'N',    'N',   -21.589,   -0.270,   38.583,  1.00,  37.72 ],
		[ 'CA',   'C',   -22.974,    0.152,   38.351,  1.00,  37.04 ],
		[ 'C',    'C',   -23.114,    1.654,   38.115,  1.00,  34.79 ],
		[ 'O',    'O',   -24.213,    2.171,   37.982,  1.00,  37.03 ],
		[ 'CB',   'C',   -23.708,   -0.336,   39.581,  1.00,  35.91 ],
		[ 'CG',   'C',   -22.689,   -0.047,   40.648,  1.00,  37.72 ],
		[ 'CD',   'C',   -21.370,   -0.484,   40.014,  1.00,  37.06 ],
	] ],
	[ 'SER', 43, [
		[ 'N',    'N',   -22.040,    2.427,   38.152,  1.00,  33.64 ],
		[ 'CA',   'C',   -22.114,    3.839,   37.817,  1.00,  32.16 ],
		[ 'C',    'C',   -21.893,    4.037,   36.307,  1.00,  29.84 ],
		[ 'O',    'O',   -21.892,    5.160,   35.795,  1.00,  28.48 ],
		[ 'CB',   'C',   -21.047,    4.568,   38.612,  1.00,  31.61 ],
		[ 'OG',   'O',   -20.807,    3.929,   39.870,  1.00,  40.10 ],
	] ],
	[ 'LEU', 44, [
		[ 'N',    'N',   -21.665,    2.944,   35.573,  1.00,  26.70 ],
		[ 'CA',   'C',   -21.368,    3.004,   34.152,  1.00,  25.83 ],
		[ 'C',    'C',   -22.688,    3.210,   33.411,  1.00,  25.12 ],
		[ 'O',    'O',   -23.648,    2.492,   33.712,  1.00,  26.66 ],
		[ 'CB',   'C',   -20.654,    1.674,   33.775,  1.00,  20.33 ],
		[ 'CG',   'C',   -19.157,    1.641,   33.370,  1.00,  15.58 ],
		[ 'CD1',  'C',   -18.357,    2.583,   34.190,  1.00,  16.59 ],
		[ 'CD2',  'C',   -18.595,    0.242,   33.578,  1.00,  12.12 ],
	] ],
	[ 'LYS', 45, [
		[ 'N',    'N',   -22.788,    4.220,   32.534,  1.00,  23.96 ],
		[ 'CA',   'C',   -23.978,    4.436,   31.712,  1.00,  21.86 ],
		[ 'C',    'C',   -23.968,    3.508,   30.483,  1.00,  23.26 ],
		[ 'O',    'O',   -22.926,    2.942,   30.134,  1.00,  20.36 ],
		[ 'CB',   'C',   -24.040,    5.866,   31.227,  1.00,  22.96 ],
		[ 'CG',   'C',   -24.273,    6.936,   32.287,  1.00,  27.73 ],
		[ 'CD',   'C',   -22.963,    7.290,   32.987,  1.00,  32.33 ],
		[ 'CE',   'C',   -23.171,    8.268,   34.139,  1.00,  34.82 ],
		[ 'NZ',   'N',   -21.909,    8.809,   34.632,  1.00,  31.45 ],
	] ],
	[ 'PRO', 46, [
		[ 'N',    'N',   -25.073,    3.249,   29.775,  1.00,  24.49 ],
		[ 'CA',   'C',   -25.061,    2.582,   28.471,  1.00,  25.08 ],
		[ 'C',    'C',   -24.337,    3.394,   27.388,  1.00,  23.19 ],
		[ 'O',    'O',   -24.157,    4.615,   27.452,  1.00,  21.52 ],
		[ 'CB',   'C',   -26.534,    2.354,   28.129,  1.00,  24.89 ],
		[ 'CG',   'C',   -27.250,    2.396,   29.461,  1.00,  25.49 ],
		[ 'CD',   'C',   -26.454,    3.477,   30.209,  1.00,  28.39 ],
	] ],
	[ 'LEU', 47, [
		[ 'N',    'N',   -23.869,    2.669,   26.392,  1.00,  20.76 ],
		[ 'CA',   'C',   -23.270,    3.273,   25.230,  1.00,  19.28 ],
		[ 'C',    'C',   -24.389,    3.866,   24.434,  1.00,  17.04 ],
		[ 'O',    'O',   -25.485,    3.297,   24.467,  1.00,  17.44 ],
		[ 'CB',   'C',   -22.607,    2.267,   24.319,  1.00,  18.70 ],
		[ 'CG',   'C',   -21.298,    1.715,   24.758,  1.00,  18.86 ],
		[ 'CD1',  'C',   -20.855,    0.649,   23.761,  1.00,  20.33 ],
		[ 'CD2',  'C',   -20.316,    2.864,   24.916,  1.00,  18.28 ],
	] ],
	[ 'SER', 48, [
		[ 'N',    'N',   -24.080,    4.947,   23.729,  1.00,  14.18 ],
		[ 'CA',   'C',   -24.993,    5.512,   22.772,  1.00,  18.24 ],
		[ 'C',    'C',   -24.147,    5.677,   21.533,  1.00,  17.73 ],
		[ 'O',    'O',   -23.122,    6.357,   21.452,  1.00,  17.23 ],
		[ 'CB',   'C',   -25.551,    6.860,   23.214,  1.00,  17.09 ],
		[ 'OG',   'O',   -24.626,    7.579,   24.006,  1.00,  30.98 ],
	] ],
	[ 'VAL', 49, [
		[ 'N',    'N',   -24.569,    4.869,   20.586,  1.00,  20.08 ],
		[ 'CA',   'C',   -23.953,    4.797,   19.288,  1.00,  21.54 ],
		[ 'C',    'C',   -24.977,    5.470,   18.372,  1.00,  21.85 ],
		[ 'O',    'O',   -26.121,    5.020,   18.303,  1.00,  20.74 ],
		[ 'CB',   'C',   -23.741,    3.311,   18.957,  1.00,  21.94 ],
		[ 'CG1',  'C',   -23.044,    3.219,   17.635,  1.00,  22.35 ],
		[ 'CG2',  'C',   -22.872,    2.616,   19.994,  1.00,  22.90 ],
	] ],
	[ 'SER', 50, [
		[ 'N',    'N',   -24.678,    6.590,   17.720,  1.00,  23.21 ],
		[ 'CA',   'C',   -25.623,    7.175,   16.781,  1.00,  21.76 ],
		[ 'C',    'C',   -24.956,    7.062,   15.432,  1.00,  20.16 ],
		[ 'O',    'O',   -24.195,    7.929,   15.018,  1.00,  19.01 ],
		[ 'CB',   'C',   -25.869,    8.625,   17.127,  1.00,  23.34 ],
		[ 'OG',   'O',   -26.300,    8.673,   18.473,  1.00,  28.98 ],
	] ],
	[ 'TYR', 51, [
		[ 'N',    'N',   -25.142,    5.941,   14.756,  1.00,  22.43 ],
		[ 'CA',   'C',   -24.520,    5.758,   13.451,  1.00,  24.77 ],
		[ 'C',    'C',   -25.458,    5.922,   12.260,  1.00,  25.70 ],
		[ 'O',    'O',   -25.093,    5.625,   11.122,  1.00,  27.42 ],
		[ 'CB',   'C',   -23.862,    4.369,   13.420,  1.00,  21.91 ],
		[ 'CG',   'C',   -22.552,    4.201,   14.184,  1.00,  21.53 ],
		[ 'CD1',  'C',   -21.952,    5.279,   14.856,  1.00,  22.86 ],
		[ 'CD2',  'C',   -21.952,    2.938,   14.213,  1.00,  18.44 ],
		[ 'CE1',  'C',   -20.759,    5.110,   15.553,  1.00,  21.20 ],
		[ 'CE2',  'C',   -20.762,    2.757,   14.916,  1.00,  21.12 ],
		[ 'CZ',   'C',   -20.169,    3.847,   15.583,  1.00,  22.50 ],
		[ 'OH',   'O',   -18.968,    3.678,   16.281,  1.00,  23.88 ],
	] ],
	[ 'ASP', 52, [
		[ 'N',    'N',   -26.680,    6.415,   12.485,  1.00,  28.94 ],
		[ 'CA',   'C',   -27.665,    6.564,   11.413,  1.00,  30.42 ],
		[ 'C',    'C',   -27.239,    7.505,   10.291,  1.00,  29.31 ],
		[ 'O',    'O',   -27.458,    7.190,    9.126,  1.00,  28.97 ],
		[ 'CB',   'C',   -29.008,    7.047,   11.987,  1.00,  34.57 ],
		[ 'CG',   'C',   -29.032,    8.482,   12.506,  1.00,  38.93 ],
		[ 'OD1',  'O',   -28.454,    8.764,   13.559,  1.00,  45.06 ],
		[ 'OD2',  'O',   -29.632,    9.320,   11.838,  1.00,  44.23 ],
	] ],
	[ 'GLN', 53, [
		[ 'N',    'N',   -26.596,    8.626,   10.614,  1.00,  25.10 ],
		[ 'CA',   'C',   -26.190,    9.569,    9.605,  1.00,  24.77 ],
		[ 'C',    'C',   -24.707,    9.513,    9.283,  1.00,  21.20 ],
		[ 'O',    'O',   -24.114,   10.506,    8.862,  1.00,  20.23 ],
		[ 'CB',   'C',   -26.623,   10.963,   10.065,  1.00,  31.21 ],
		[ 'CG',   'C',   -26.431,   11.213,   11.551,  1.00,  44.84 ],
		[ 'CD',   'C',   -26.867,   12.577,   12.060,  1.00,  52.21 ],
		[ 'OE1',  'O',   -27.953,   13.064,   11.732,  1.00,  54.57 ],
		[ 'NE2',  'N',   -26.052,   13.222,   12.893,  1.00,  53.94 ],
	] ],
	[ 'ALA', 54, [
		[ 'N',    'N',   -24.052,    8.362,    9.423,  1.00,  19.59 ],
		[ 'CA',   'C',   -22.627,    8.262,    9.100,  1.00,  22.21 ],
		[ 'C',    'C',   -22.337,    8.307,    7.593,  1.00,  21.45 ],
		[ 'O',    'O',   -23.064,    7.743,    6.774,  1.00,  21.70 ],
		[ 'CB',   'C',   -22.031,    6.969,    9.630,  1.00,  18.96 ],
	] ],
	[ 'THR', 55, [
		[ 'N',    'N',   -21.286,    9.022,    7.192,  1.00,  21.32 ],
		[ 'CA',   'C',   -20.930,    9.143,    5.795,  1.00,  17.74 ],
		[ 'C',    'C',   -19.585,    8.501,    5.480,  1.00,  16.62 ],
		[ 'O',    'O',   -18.546,    9.098,    5.795,  1.00,  11.41 ],
		[ 'CB',   'C',   -20.952,   10.640,    5.470,  1.00,  20.84 ],
		[ 'OG1',  'O',   -22.279,   11.088,    5.752,  1.00,  24.44 ],
		[ 'CG2',  'C',   -20.575,   10.942,    4.021,  1.00,  22.48 ],
	] ],
	[ 'SER', 56, [
		[ 'N',    'N',   -19.559,    7.276,    4.932,  1.00,  14.57 ],
		[ 'CA',   'C',   -18.298,    6.680,    4.540,  1.00,  18.61 ],
		[ 'C',    'C',   -17.924,    7.280,    3.202,  1.00,  20.27 ],
		[ 'O',    'O',   -18.770,    7.703,    2.389,  1.00,  22.10 ],
		[ 'CB',   'C',   -18.410,    5.171,    4.412,  1.00,  13.25 ],
		[ 'OG',   'O',   -19.431,    4.736,    3.532,  1.00,  19.48 ],
	] ],
	[ 'LEU', 57, [
		[ 'N',    'N',   -16.618,    7.348,    3.035,  1.00,  19.49 ],
		[ 'CA',   'C',   -16.081,    7.921,    1.838,  1.00,  15.54 ],
		[ 'C',    'C',   -15.251,    6.938,    1.069,  1.00,  15.08 ],
		[ 'O',    'O',   -15.485,    6.802,   -0.126,  1.00,  16.53 ],
		[ 'CB',   'C',   -15.204,    9.096,    2.161,  1.00,  14.41 ],
		[ 'CG',   'C',   -15.857,   10.256,    2.874,  1.00,  18.29 ],
		[ 'CD1',  'C',   -14.766,   11.140,    3.422,  1.00,  18.41 ],
		[ 'CD2',  'C',   -16.831,   10.948,    1.941,  1.00,  13.90 ],
	] ],
	[ 'ARG', 58, [
		[ 'N',    'N',   -14.353,    6.189,    1.718,  1.00,  12.47 ],
		[ 'CA',   'C',   -13.280,    5.519,    1.012,  1.00,   9.99 ],
		[ 'C',    'C',   -12.722,    4.304,    1.725,  1.00,  10.83 ],
		[ 'O',    'O',   -12.715,    4.304,    2.954,  1.00,   8.62 ],
		[ 'CB',   'C',   -12.244,    6.588,    0.827,  1.00,  10.49 ],
		[ 'CG',   'C',   -10.844,    6.185,    0.503,  1.00,  14.59 ],
		[ 'CD',   'C',    -9.975,    7.134,    1.295,  1.00,  12.65 ],
		[ 'NE',   'N',    -8.953,    7.605,    0.405,  1.00,  13.99 ],
		[ 'CZ',   'C',    -7.741,    7.982,    0.768,  1.00,  12.60 ],
		[ 'NH1',  'N',    -7.289,    8.004,    2.025,  1.00,  15.79 ],
		[ 'NH2',  'N',    -6.939,    8.217,   -0.248,  1.00,  17.59 ],
	] ],
	[ 'ILE', 59, [
		[ 'N',    'N',   -12.291,    3.260,    1.021,  1.00,   8.27 ],
		[ 'CA',   'C',   -11.602,    2.165,    1.671,  1.00,   9.17 ],
		[ 'C',    'C',   -10.176,    2.234,    1.140,  1.00,  11.49 ],
		[ 'O',    'O',    -9.943,    2.657,    0.001,  1.00,  13.74 ],
		[ 'CB',   'C',   -12.279,    0.818,    1.348,  1.00,   8.19 ],
		[ 'CG1',  'C',   -11.553,   -0.254,    2.168,  1.00,   2.00 ],
		[ 'CG2',  'C',   -12.328,    0.557,   -0.154,  1.00,   2.00 ],
		[ 'CD1',  'C',   -12.372,   -1.549,    2.216,  1.00,   4.60 ],
	] ],
	[ 'LEU', 60, [
		[ 'N',    'N',    -9.197,    1.858,    1.948,  1.00,  13.54 ],
		[ 'CA',   'C',    -7.790,    2.046,    1.617,  1.00,  15.37 ],
		[ 'C',    'C',    -6.959,    0.903,    2.186,  1.00,  13.46 ],
		[ 'O',    'O',    -7.179,    0.489,    3.325,  1.00,  12.87 ],
		[ 'CB',   'C',    -7.304,    3.376,    2.210,  1.00,  15.32 ],
		[ 'CG',   'C',    -5.820,    3.694,    2.140,  1.00,  22.87 ],
		[ 'CD1',  'C',    -5.465,    4.279,    0.783,  1.00,  20.45 ],
		[ 'CD2',  'C',    -5.481,    4.663,    3.251,  1.00,  17.67 ],
	] ],
	[ 'ASN', 61, [
		[ 'N',    'N',    -6.034,    0.359,    1.412,  1.00,  11.19 ],
		[ 'CA',   'C',    -5.094,   -0.587,    1.948,  1.00,  10.07 ],
		[ 'C',    'C',    -3.908,    0.291,    2.286,  1.00,   9.63 ],
		[ 'O',    'O',    -3.312,    0.925,    1.404,  1.00,   8.87 ],
		[ 'CB',   'C',    -4.688,   -1.605,    0.921,  1.00,   7.26 ],
		[ 'CG',   'C',    -3.685,   -2.605,    1.445,  1.00,   9.46 ],
		[ 'OD1',  'O',    -2.859,   -2.341,    2.330,  1.00,  10.96 ],
		[ 'ND2',  'N',    -3.717,   -3.804,    0.900,  1.00,   6.76 ],
	] ],
	[ 'ASN', 62, [
		[ 'N',    'N',    -3.531,    0.353,    3.563,  1.00,   8.75 ],
		[ 'CA',   'C',    -2.429,    1.216,    3.976,  1.00,   8.07 ],
		[ 'C',    'C',    -1.116,    0.497,    4.270,  1.00,   8.74 ],
		[ 'O',    'O',    -0.204,    1.096,    4.836,  1.00,   6.81 ],
		[ 'CB',   'C',    -2.845,    2.015,    5.204,  1.00,   8.11 ],
		[ 'CG',   'C',    -3.112,    1.194,    6.460,  1.00,   7.99 ],
		[ 'OD1',  'O',    -2.760,    0.024,    6.607,  1.00,   5.57 ],
		[ 'ND2',  'N',    -3.784,    1.771,    7.428,  1.00,  13.33 ],
	] ],
);
my @fold;
$serial = 0;
for my $r (@fold_res) {
	my ($resname, $resseq, $atoms) = @$r;
	for my $a (@$atoms) {
		$serial++;
		push @fold, atom_line(
			record => 'ATOM  ', serial => $serial, name => $a->[0], element => $a->[1],
			altloc => '', resname => $resname, chain => 'A', resseq => $resseq, icode => '',
			x => $a->[2], y => $a->[3], z => $a->[4], occ => $a->[5], b => $a->[6],
		);
	}
}
push @fold, 'END';

# --- the mmCIF twins ------------------------------------------------------
#
# The same structures, written the other way.  The coordinates are converted
# from the records above rather than typed again, because the point of the
# .cif fixtures is that reading either file gives the same answer, and a
# fixture pair that was typed twice tests the typing.
#
# What is deliberately not converted is the naming.  A real mmCIF file carries
# two sets of identifiers -- label_* assigned by the archive, auth_* as the
# depositor numbered them -- and only the auth_* ones match the PDB record.
# So the label_* columns written below are the other ones on purpose: chains
# lettered straight through including the waters, residues numbered from 1
# with no gap and no insertion code.  A reader that reached for label_asym_id
# would produce a structure with six chains in it, and t/cif.t would say so.

# cifq() -- one value, quoted the way the format needs it
sub cifq {
	my ($v) = @_;
	return '?' unless defined $v && length $v;
	return $v unless $v =~ /[\s'"]/ || $v =~ /\A[_\#\$\[\]]/ || $v =~ /\A(?:data|loop|save|stop|global)_/i;
	return "'$v'" if $v !~ /'/;
	return "\"$v\"" if $v !~ /"/;
	return "\n;$v\n;";      # a value holding both quotes has only one way left
}

my @ATOM_ITEM = qw(
	group_PDB id type_symbol label_atom_id label_alt_id label_comp_id
	label_asym_id label_entity_id label_seq_id pdbx_PDB_ins_code
	Cartn_x Cartn_y Cartn_z occupancy B_iso_or_equiv pdbx_formal_charge
	auth_seq_id auth_comp_id auth_asym_id auth_atom_id pdbx_PDB_model_num
);

# atom_site_loop() -- the ATOM/HETATM records of a PDB file as an mmCIF loop.
# Read by column, because that is where a PDB record keeps its fields, and the
# whole point is to carry every one of them across unchanged.
sub atom_site_loop {
	my ($lines, %opt) = @_;
	my (@rows, %asym, $model, $seq);
	$model = 1;
	my %seq_of;      # label_seq_id: the polymer position, counted per label asym
	for my $l (@$lines) {
		if ($l =~ /\AMODEL\s+(\d+)/) { $model = $1; next }
		next unless $l =~ /\A(ATOM  |HETATM)/;
		my %a = (
			group   => ($1 eq 'ATOM  ' ? 'ATOM' : 'HETATM'),
			serial  => _t(substr($l, 6, 5)),
			name    => _t(substr($l, 12, 4)),
			altloc  => _t(substr($l, 16, 1)),
			resname => _t(substr($l, 17, 3)),
			chain   => _t(substr($l, 21, 1)),
			resseq  => _t(substr($l, 22, 4)),
			icode   => _t(substr($l, 26, 1)),
			x       => _t(substr($l, 30, 8)),
			y       => _t(substr($l, 38, 8)),
			z       => _t(substr($l, 46, 8)),
			occ     => _t(substr($l, 54, 6)),
			b       => _t(substr($l, 60, 6)),
			element => (length($l) > 76 ? _t(substr($l, 76, 2)) : ''),
			charge  => (length($l) > 78 ? _t(substr($l, 78, 2)) : ''),
		);
		# label_asym_id: a fresh letter per chain and per kind of thing in it,
		# which is how the archive assigns them and is not the PDB chain id
		my $kind = $a{group} eq 'HETATM' ? ($a{resname} eq 'HOH' ? 'w' : "h$a{resname}") : 'p';
		# lettered from B rather than from A, so that no label_asym_id can
		# coincide with the auth_asym_id of the chain it belongs to and a
		# reader that took the wrong one cannot pass by luck
		my $ak = "$a{chain}/$kind";
		$asym{$ak} = chr(ord('B') + $asym{n}++) unless exists $asym{$ak};
		my $lasym = $asym{$ak};
		my $lseq  = '.';
		if ($kind eq 'p') {
			my $rk = "$ak/$a{resseq}$a{icode}";
			$seq_of{$lasym}{$rk} ||= ++$seq_of{$lasym}{n};
			$lseq = $seq_of{$lasym}{$rk};
		}
		# a PDB charge is "2+", an mmCIF one is 2; converted here so that the
		# reader has the conversion to undo
		my $chg = '?';
		if ($a{charge} =~ /\A(\d)([-+])\z/) { $chg = ($2 eq '-' ? "-$1" : $1) }
		# the two spellings of nothing: '.' where the item does not apply to
		# this row, '?' where it does and the file does not know it.  Both are
		# written, because both have to read back as an empty field.
		push @rows, [
			$a{group}, $a{serial},
			(length $a{element} ? $a{element} : '?'),
			$a{name},
			(length $a{altloc} ? $a{altloc} : '.'),
			$a{resname}, $lasym, ($kind eq 'p' ? 1 : 2), $lseq,
			(length $a{icode} ? $a{icode} : '?'),
			$a{x}, $a{y}, $a{z}, $a{occ}, $a{b}, $chg,
			$a{resseq}, $a{resname}, $a{chain}, $a{name}, $model,
		];
	}
	my @items = @ATOM_ITEM;
	my @keep  = 0 .. $#items;
	if ($opt{no_element}) {      # a file with no type_symbol, as bare.pdb has no element
		@keep  = grep { $items[$_] ne 'type_symbol' } @keep;
	}
	my @out = ('loop_', map { "_atom_site.$items[$_]" } @keep);
	for my $r (@rows) {
		push @out, join ' ', map { ($_ eq '.' || $_ eq '?') ? $_ : cifq($_) } @{$r}[@keep];
	}
	return @out;
}

sub _t { my $s = shift; return '' unless defined $s; $s =~ s/\A\s+//; $s =~ s/\s+\z//; return $s }

# a category with one row, written as plain tags
sub cif_pairs {
	my ($cat, @kv) = @_;
	my @out;
	while (@kv) {
		my ($k, $v) = splice @kv, 0, 2;
		push @out, sprintf('%-52s %s', "_$cat.$k", cifq($v));
	}
	return @out;
}

# a category with several rows, written as a loop_
sub cif_loop {
	my ($cat, $items, @rows) = @_;
	my @out = ('loop_', map { "_$cat.$_" } @$items);
	push @out, join ' ', map { cifq($_) } @$_ for @rows;
	return @out;
}

my @minicif = ('data_9XYZ', '#');
push @minicif,
	cif_pairs('entry', id => '9XYZ'), '#',
	# a semicolon text field, which is the only way the format has of writing
	# a value too long for a line -- and the only token that spans lines
	'_struct.entry_id   9XYZ',
	'_struct.title',
	';A SMALL TEST STRUCTURE WITH A GAP, AN INSERTION CODE, AN ALTERNATE CONFORMER AND A LIGAND',
	';', '#',
	cif_pairs('struct_keywords',
		entry_id      => '9XYZ',
		pdbx_keywords => 'HYDROLASE/PEPTIDE INHIBITOR',
		text          => 'HYDROLASE, TEST STRUCTURE, COMPLEX (HYDROLASE-PEPTIDE)'), '#',
	cif_pairs('pdbx_database_status', entry_id => '9XYZ',
		recvd_initial_deposition_date => '2020-01-01'), '#',
	cif_pairs('exptl', entry_id => '9XYZ', method => 'X-RAY DIFFRACTION'), '#',
	cif_pairs('refine',
		entry_id            => '9XYZ',
		'ls_d_res_high'     => '1.85',
		'ls_R_factor_R_work'=> '0.174',
		'ls_R_factor_R_free'=> '0.219'), '#',
	cif_pairs('diffrn', id => 1, ambient_temp => '100.0'), '#',
	cif_pairs('exptl_crystal_grow', crystal_id => 1, pH => '7.5'), '#',
	cif_loop('audit_author', [qw(name pdbx_ordinal)],
		[ 'Condon, D.E.', 1 ], [ 'Other, A.N.', 2 ]), '#',
	cif_loop('citation',
		[qw(id title journal_abbrev journal_volume page_first year
		    pdbx_database_id_PubMed pdbx_database_id_DOI)],
		[ 'primary', 'A STRUCTURE MADE UP FOR A TEST SUITE, AND WHAT IT CONTAINS',
		  'J.Invented.Res.', 10, 42, 2020, 12345678, '10.1000/INVENTED.2020.42' ]), '#',
	cif_loop('citation_author', [qw(citation_id name ordinal)],
		[ 'primary', 'Condon, D.E.', 1 ], [ 'primary', 'Other, A.N.', 2 ]), '#',
	cif_pairs('cell', entry_id => '9XYZ',
		length_a => '40.100', length_b => '50.200', length_c => '60.300',
		angle_alpha => '90.00', angle_beta => '95.50', angle_gamma => '90.00',
		'Z_PDB' => 4), '#',
	cif_pairs('symmetry', entry_id => '9XYZ', 'space_group_name_H-M' => 'P 1 21 1'), '#',
	cif_loop('entity', [qw(id type src_method pdbx_description pdbx_ec)],
		[ 1, 'polymer',     'man', 'TEST PROTEIN', '3.4.21.5' ],
		[ 2, 'polymer',     'syn', 'TEST DNA',     undef ],
		[ 3, 'non-polymer', 'syn', '2-ACETAMIDO-2-DEOXY-BETA-D-GLUCOPYRANOSE', undef ],
		[ 4, 'non-polymer', 'syn', 'ZINC ION',     undef ],
		[ 5, 'water',       'nat', 'water',        undef ]), '#',
	cif_loop('entity_poly', [qw(entity_id type pdbx_seq_one_letter_code_can pdbx_strand_id)],
		[ 1, 'polypeptide(L)',       'MAGLKCMHHSC', 'A' ],
		[ 2, 'polydeoxyribonucleotide', 'ACGT',      'B' ]), '#',
	cif_loop('entity_poly_seq', [qw(entity_id num mon_id hetero)],
		(map { [ 1, $_->[0], $_->[1], 'n' ] }
		 map { [ $_ + 1, (qw(MET ALA GLY LEU LYS CYS MSE HIS HIS SER CYS))[$_] ] } 0 .. 10),
		(map { [ 2, $_->[0], $_->[1], 'n' ] }
		 map { [ $_ + 1, (qw(DA DC DG DT))[$_] ] } 0 .. 3)), '#',
	cif_pairs('entity_src_gen',
		entity_id                        => 1,
		pdbx_gene_src_scientific_name    => 'HOMO SAPIENS',
		pdbx_gene_src_ncbi_taxonomy_id   => 9606,
		pdbx_host_org_scientific_name    => 'ESCHERICHIA COLI'), '#',
	cif_loop('chem_comp', [qw(id name formula type)],
		[ 'MSE', 'SELENOMETHIONINE', 'C5 H11 N O2 Se', 'L-peptide linking' ],
		[ 'NAG', '2-ACETAMIDO-2-DEOXY-BETA-D-GLUCOPYRANOSE', 'C8 H15 N O6', 'D-saccharide' ],
		[ 'ZN',  'ZINC ION',  'ZN 2+', 'non-polymer' ],
		[ 'HOH', 'WATER',     'H2 O',  'water' ]), '#',
	cif_loop('pdbx_nonpoly_scheme',
		[qw(asym_id entity_id mon_id pdb_strand_id pdb_seq_num pdb_ins_code)],
		[ 'C', 3, 'NAG', 'A', 201, '.' ],
		[ 'D', 4, 'ZN',  'A', 202, '.' ],
		[ 'E', 5, 'HOH', 'A', 301, '.' ],
		[ 'E', 5, 'HOH', 'A', 302, '.' ]), '#',
	cif_loop('pdbx_struct_mod_residue',
		[qw(id label_comp_id auth_comp_id auth_asym_id auth_seq_id parent_comp_id details)],
		[ 1, 'MSE', 'MSE', 'A', 7, 'MET', 'SELENOMETHIONINE' ]), '#',
	cif_loop('struct_conf',
		[qw(conf_type_id id beg_auth_comp_id beg_auth_asym_id beg_auth_seq_id
		    end_auth_comp_id end_auth_asym_id end_auth_seq_id
		    pdbx_PDB_helix_class pdbx_PDB_helix_length)],
		[ 'HELX_P', 'AA1', 'MET', 'A', 1, 'GLY', 'A', 3, 1, 3 ]), '#',
	cif_loop('struct_sheet_range',
		[qw(sheet_id id beg_auth_comp_id beg_auth_asym_id beg_auth_seq_id
		    end_auth_comp_id end_auth_asym_id end_auth_seq_id)],
		[ 'AA1', 1, 'CYS', 'A', 6, 'HIS', 'A', 8 ]), '#',
	cif_loop('struct_conn',
		[qw(id conn_type_id ptnr1_label_atom_id ptnr1_auth_comp_id ptnr1_auth_asym_id
		    ptnr1_auth_seq_id ptnr2_label_atom_id ptnr2_auth_comp_id ptnr2_auth_asym_id
		    ptnr2_auth_seq_id pdbx_dist_value)],
		[ 'disulf1', 'disulf', 'SG', 'CYS', 'A', 6,  'SG', 'CYS', 'A', 10, '2.03' ],
		[ 'covale1', 'covale', 'ZN', 'ZN',  'A', 202, 'SG', 'CYS', 'A', 6,  '2.31' ]), '#',
	cif_loop('struct_mon_prot_cis',
		[qw(pdbx_id auth_comp_id auth_asym_id auth_seq_id pdbx_auth_comp_id_2
		    pdbx_auth_asym_id_2 pdbx_auth_seq_id_2 pdbx_omega_angle)],
		[ 1, 'GLY', 'A', 3, 'CYS', 'A', 6, '-0.42' ]), '#',
	cif_loop('struct_ref', [qw(id db_name db_code pdbx_db_accession entity_id)],
		[ 1, 'UNP', 'TEST_HUMAN', 'P12345', 1 ]), '#',
	cif_loop('struct_ref_seq',
		[qw(align_id ref_id pdbx_strand_id pdbx_auth_seq_align_beg
		    pdbx_auth_seq_align_end pdbx_db_accession db_align_beg db_align_end)],
		[ 1, 1, 'A', 1, 11, 'P12345', 1, 11 ]), '#';
push @minicif, atom_site_loop(\@mini), '#';

my @nmrcif = ('data_9NMR', '#');
push @nmrcif,
	cif_pairs('entry', id => '9NMR'), '#',
	'_struct.title    "A THREE MODEL ENSEMBLE"', '#',
	cif_pairs('exptl', entry_id => '9NMR', method => 'SOLUTION NMR'), '#',
	cif_pairs('pdbx_nmr_ensemble', entry_id => '9NMR',
		conformers_submitted_total_number => 3), '#',
	cif_loop('entity', [qw(id type pdbx_description)], [ 1, 'polymer', 'TEST PEPTIDE' ]), '#',
	cif_loop('entity_poly', [qw(entity_id pdbx_strand_id)], [ 1, 'A' ]), '#',
	cif_loop('entity_poly_seq', [qw(entity_id num mon_id)],
		[ 1, 1, 'GLY' ], [ 1, 2, 'SER' ], [ 1, 3, 'TRP' ]), '#';
push @nmrcif, atom_site_loop(\@nmr), '#';

# bare.cif -- coordinates and nothing else, and no type_symbol, so the element
# has to come out of the atom name here exactly as it does from a PDB record
# with no element columns
my @barecif = ('data_bare', '#', atom_site_loop(\@bare, no_element => 1), '#');

# The two feature fixtures, the other way round.  Both are plain ATOM records
# with nothing else in them, so the conversion is the whole file.
my @stackcif = ('data_9ARO', '#', atom_site_loop(\@stack), '#');
my @basescif = ('data_9NUC', '#', atom_site_loop(\@bases), '#');
my @sscif    = ('data_9SSB', '#', atom_site_loop(\@ss), '#');
my @foldcif  = ('data_9FLD', '#', atom_site_loop(\@fold), '#');

# quirks.cif -- everything about the way the format is written down that a
# reader has to get right, in one file: comments in every position, both kinds
# of quote, a quote inside a value (O5', which is an atom name and not a
# string that someone forgot to close), a semicolon text field, '.' and '?' for
# the two kinds of nothing, a quoted '.' that is a full stop and not a null, a
# formal charge in the mmCIF spelling, a category written as plain tags where
# it is usually a loop, and columns in an order no writer uses.
my @quirks = (
'# a comment before anything at all',
'data_QRK    # and one after the block name',
'#',
'_entry.id   QRK',
'_struct.title     "A file that leans on the syntax"',
"_struct_keywords.text    'one, two, three'",
'_exptl.method',
';SOLUTION NMR',
';',
'#',
'loop_',
'_atom_site.auth_atom_id',          # the columns in a deliberately odd order
'_atom_site.pdbx_formal_charge',
'_atom_site.auth_comp_id',
'_atom_site.group_PDB',
'_atom_site.auth_asym_id',
'_atom_site.Cartn_x',
'_atom_site.Cartn_y',
'_atom_site.Cartn_z',
'_atom_site.auth_seq_id',
'_atom_site.type_symbol',
'_atom_site.id',
'_atom_site.label_alt_id',
'_atom_site.occupancy',
'_atom_site.B_iso_or_equiv',
'_atom_site.pdbx_PDB_ins_code',
"P     ?   G  ATOM   B  1.000 2.000 3.000 1 P  1 . 1.00 10.00 ?",
"OP1   -1  G  ATOM   B  2.000 3.000 4.000 1 O  2 . 1.00 11.00 ?",
"\"O5'\" 0   G  ATOM   B  3.000 4.000 5.000 1 O  3 . 1.00 12.00 ?",
"\"C1'\" 3   G  ATOM   B  4.000 5.000 6.000 1 C  4 . 1.00 13.00 ?",
'# the ion is a HETATM, is charged, and has an insertion code',
"ZN    2   ZN HETATM B  9.000 9.000 9.000 40 ZN 5 . 1.00 14.00 A",
'#',
'# a category that is usually a loop, written as plain tags because it has',
'# one row -- which the format allows and a reader has to accept',
"_chem_comp.id        ZN",
"_chem_comp.name      'ZINC ION'",
"_chem_comp.formula   'ZN 2+'",
'#',
'loop_',
'_citation.id',
'_citation.title',
'_citation.year',
"primary  'A paper with a full stop.  And two sentences.'  2021",
'#',
'_cell.length_a    .',                  # not applicable
'_cell.length_b    ?',                  # unknown
"_pdbx_database_status.recvd_initial_deposition_date  '.'",  # quoted: a value
'#',
);

for my $f ([ 'mini.pdb', \@mini ], [ 'nmr.pdb', \@nmr ], [ 'bare.pdb', \@bare ],
           [ 'mini.cif', \@minicif ], [ 'nmr.cif', \@nmrcif ], [ 'bare.cif', \@barecif ],
           [ 'stack.pdb', \@stack ], [ 'stack.cif', \@stackcif ],
           [ 'bases.pdb', \@bases ], [ 'bases.cif', \@basescif ],
           [ 'ss.pdb', \@ss ], [ 'ss.cif', \@sscif ],
           [ 'fold.pdb', \@fold ], [ 'fold.cif', \@foldcif ],
           [ 'quirks.cif', \@quirks ]) {
	open my $fh, '>', $f->[0];
	print {$fh} "$_\n" for @{ $f->[1] };
	close $fh;
	print "wrote $f->[0] (" . scalar(@{ $f->[1] }) . " lines)\n";
}

# an empty file is a legitimate thing to be handed, and must not die
open my $e, '>', 'empty.pdb';
close $e;
print "wrote empty.pdb\n";
open my $ec, '>', 'empty.cif';
close $ec;
print "wrote empty.cif\n";
