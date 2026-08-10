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

for my $f ([ 'mini.pdb', \@mini ], [ 'nmr.pdb', \@nmr ], [ 'bare.pdb', \@bare ]) {
	open my $fh, '>', $f->[0];
	print {$fh} "$_\n" for @{ $f->[1] };
	close $fh;
	print "wrote $f->[0] (" . scalar(@{ $f->[1] }) . " lines)\n";
}

# an empty file is a legitimate thing to be handed, and must not die
open my $e, '>', 'empty.pdb';
close $e;
print "wrote empty.pdb\n";
