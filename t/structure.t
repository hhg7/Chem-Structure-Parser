#!/usr/bin/env perl
# The assembled hash of hashes: chains, residues, atoms, and the numbers that
# go with them, read from t/data/mini.pdb.
require 5.010;
use strict;
use warnings FATAL => 'all';
use Cwd 'abs_path';
use File::Basename 'dirname';
use Structure::Info;
use Test::Exception;
use Test::More;

my $data = dirname(abs_path(__FILE__)) . '/data';
my $info = structure_info("$data/mini.pdb");

#--------
# the top of the structure
#--------
is($info->{id},       '9XYZ', 'id comes from the HEADER record');
is($info->{format},   'pdb',  'the format is recorded');
is($info->{n_models}, 1,      'a file with no MODEL records has one model');
is($info->{model},    1,      'and the chains were built from it');
like($info->{file}, qr/mini\.pdb\z/, 'the file it was read from is kept');

#--------
# chains
#--------
is_deeply($info->{chain_order}, [ 'A', 'B' ], 'chains are in the order the file has them');
is($info->{chains}{A}{type}, 'protein', 'chain A is a protein');
is($info->{chains}{B}{type}, 'dna',     'chain B is DNA, not RNA: its residues are DA, DC, DG, DT');
is($info->{chains}{A}{molecule}, 'TEST PROTEIN', 'the chain knows what molecule it is, from COMPND');
is($info->{chains}{A}{organism}, 'HOMO SAPIENS', 'and what it came from, from SOURCE');
is($info->{chains}{A}{ec},       '3.4.21.5',     'and its EC number');
is($info->{chains}{A}{mol_id},   '1',            'and which COMPND entity it belongs to');
is($info->{chains}{B}{molecule}, 'TEST DNA',     'the second entity is mapped to its own chain');

#--------
# the sequence, which is the thing this was written for
#--------
is($info->{chains}{A}{sequence}, 'MAGCMHHSC',
	'chain A single-letter sequence: MSE counts as M, and the ligand, ion and waters are not in it');
is($info->{chains}{B}{sequence}, 'ACGT', 'chain B reads as DNA');
is($info->{chains}{A}{seqres}, 'MAGLKCMHHSC', 'SEQRES is the full sequence, including what was not modelled');
is($info->{chains}{A}{seqres_length}, 11, 'and its length');
is($info->{chains}{A}{n_missing},      2, 'two residues in SEQRES have no coordinates');
is($info->{chains}{B}{seqres}, 'ACGT', 'SEQRES for a DNA chain');

is_deeply(structure_sequences($info), { A => 'MAGCMHHSC', B => 'ACGT' },
	'structure_sequences: every chain that has a sequence');
is_deeply(structure_sequences("$data/mini.pdb"), { A => 'MAGCMHHSC', B => 'ACGT' },
	'structure_sequences: a file name is read on the spot');
is_deeply(structure_sequences("$data/mini.pdb", atoms => 0, meta => 0),
	{ A => 'MAGCMHHSC', B => 'ACGT' },
	'structure_sequences: and reads it with the options it was given');
throws_ok { structure_sequences("$data/no-such-file.pdb") } qr/does not exist/,
	'structure_sequences: a file name that is not a file says so';
throws_ok { structure_sequences($info, atoms => 0) } qr/options apply to reading a file/,
	'structure_sequences: options with an already-parsed structure would do nothing, so they are refused';
throws_ok { structure_sequences(undef) } qr/expected a file name/,
	'structure_sequences: and nothing at all is still an error';
is(chain_sequence($info, 'A'), 'MAGCMHHSC', 'chain_sequence: observed by default');
is(chain_sequence($info, 'A', 'seqres'), 'MAGLKCMHHSC', 'chain_sequence: seqres on request');

#--------
# gaps -- where the observed sequence and SEQRES part company
#--------
is($info->{chains}{A}{n_gaps}, 1, 'chain A has one gap');
is_deeply($info->{chains}{A}{gaps}, [ { after => '3', before => '6', missing => 2 } ],
	'the gap names the residues either side of it and how many are missing');
is($info->{chains}{B}{n_gaps}, 0, 'chain B is continuous');

#--------
# counts
#--------
is($info->{chains}{A}{n_residues}, 13, 'chain A: 9 polymer residues, a ligand, an ion and two waters');
is($info->{chains}{A}{n_polymer},   9, 'nine of them are polymer');
is($info->{chains}{A}{n_water},     2, 'two are water');
is($info->{chains}{A}{n_ligand},    2, 'two are heterogens (the sugar and the zinc)');
is($info->{chains}{B}{n_residues},  4, 'chain B has four nucleotides');
is($info->{stats}{n_atoms}, 66, 'every coordinate record was counted');
is($info->{stats}{n_atoms},
	$info->{chains}{A}{n_atoms} + $info->{chains}{B}{n_atoms},
	'and the chains account for all of them');

#--------
# residues
#--------
my $r = $info->{chains}{A}{residues}{6};
is($r->{resname}, 'CYS', 'residue 6 is a cysteine');
is($r->{one},     'C',   'its single-letter code');
is($r->{type},    'amino_acid', 'its type');
is($r->{number},  6,   'its number');
is($r->{icode},   '',  'and no insertion code');
is($r->{hetero},  0,   'it was written as ATOM');
is($r->{standard}, 1,  'it is one of the twenty');
is($r->{modified}, 0,  'and not modified');
is($r->{n_atoms},  5,  'five atoms');
is_deeply($r->{atom_order}, [ qw(N CA C O SG) ], 'in the order the file has them');

# a modified residue is still part of the sequence
my $mse = $info->{chains}{A}{residues}{7};
is($mse->{resname},  'MSE',        'residue 7 is selenomethionine');
is($mse->{one},      'M',          'which reads as M');
is($mse->{type},     'amino_acid', 'and counts as an amino acid');
is($mse->{hetero},   1,            'even though it was written as HETATM');
is($mse->{modified}, 1,            'it is flagged as modified');
is($mse->{standard}, 0,            'and not standard');

# insertion codes: 8 and 8A are two residues, not one
ok(exists $info->{chains}{A}{residues}{'8'},  'residue 8 exists');
ok(exists $info->{chains}{A}{residues}{'8A'}, 'and residue 8A beside it');
is($info->{chains}{A}{residues}{'8A'}{icode},  'A', 'the insertion code is kept');
is($info->{chains}{A}{residues}{'8A'}{number},   8, 'and the number is still 8');
isnt($info->{chains}{A}{residues}{'8'}{atoms}{CA}{serial},
     $info->{chains}{A}{residues}{'8A'}{atoms}{CA}{serial},
     'the two residues hold different atoms');

# a ligand, an ion and water are told apart
is($info->{chains}{A}{residues}{201}{type}, 'ligand', 'NAG is a ligand');
is($info->{chains}{A}{residues}{202}{type}, 'ion',    'a lone ZN is an ion, not a ligand');
is($info->{chains}{A}{residues}{301}{type}, 'water',  'HOH is water');
is($info->{chains}{A}{residues}{201}{one},  '',       'a ligand has no single-letter code');

#--------
# atoms
#--------
my $ca = $info->{chains}{A}{residues}{6}{atoms}{CA};
is($ca->{name},    'CA', 'the atom knows its name');
is($ca->{element}, 'C',  'and its element');
is($ca->{hetero},  0,    'and which record it came from');
ok(!exists $ca->{altlocs}, 'an atom with one conformer carries no altloc list');
cmp_ok(abs($ca->{x} - 20.5), '<', 1e-9, 'x');
cmp_ok(abs($ca->{y} - 14.0), '<', 1e-9, 'y');
cmp_ok(abs($ca->{z} - 13.5), '<', 1e-9, 'z');

is($info->{chains}{A}{residues}{202}{atoms}{ZN}{element}, 'ZN',
	'a two-letter element is read as two letters');

# a residue's centre and mean B-factor
my $c = $info->{chains}{A}{residues}{6}{center};
is(scalar @$c, 3, 'a residue has a centre');
ok(defined $info->{chains}{A}{residues}{6}{b_mean}, 'and a mean B-factor');

#--------
# alternate conformers.  Both are kept on the atom; which one supplies the
# coordinates is the altloc option's business, tested in options.t.
#--------
my $cb = $info->{chains}{A}{residues}{2}{atoms}{CB};
is($cb->{altloc}, 'A', 'the first conformer supplies the coordinates by default');
is(scalar @{ $cb->{altlocs} }, 2, 'but both conformers are recorded');
is_deeply([ map { $_->{altloc} } @{ $cb->{altlocs} } ], [ 'A', 'B' ], 'in file order');
is($cb->{altlocs}[1]{occupancy}, 0.6, 'with the occupancy of each');
is($info->{chains}{A}{residues}{2}{n_atoms}, 6,
	'the residue counts both conformer records');
is(scalar @{ $info->{chains}{A}{residues}{2}{atom_order} }, 5,
	'but has five distinct atom names');

#--------
# whole-structure statistics
#--------
is($info->{stats}{n_hetatm}, 13, 'HETATM records are counted: MSE, NAG, ZN and two waters');
is($info->{stats}{n_water_atoms}, 2, 'water atoms are counted');
is($info->{stats}{n_hydrogens},   1, 'so are hydrogens');
is($info->{stats}{elements}{S},   2, 'elements are tallied');
is($info->{stats}{elements}{SE},  1, 'including two-letter ones');
ok($info->{stats}{bfactor}{min} <= $info->{stats}{bfactor}{mean}, 'B-factor min <= mean');
ok($info->{stats}{bfactor}{mean} <= $info->{stats}{bfactor}{max}, 'B-factor mean <= max');
ok($info->{stats}{bbox}{xmin} < $info->{stats}{bbox}{xmax}, 'the bounding box has a width');
is(scalar @{ $info->{stats}{center} }, 3, 'and a centre');

#--------
# a file with nothing but coordinates in it
#--------
{
	my $bare = structure_info("$data/bare.pdb");
	is($bare->{id}, 'BARE', 'with no HEADER, the id falls back to the file name');
	is($bare->{chains}{A}{sequence}, 'VK', 'and the sequence is still read');
	is($bare->{chains}{A}{seqres}, undef, 'a file with no SEQRES has no seqres sequence');
	is($bare->{chains}{A}{residues}{1}{atoms}{CA}{element}, 'C',
		'and the elements are worked out from the atom names');
	is($bare->{title}, undef, 'a missing TITLE is undef rather than missing');
	is_deeply($bare->{keywords}, [], 'and a missing KEYWDS is an empty list');
}

#--------
# an empty file
#--------
{
	my $empty = structure_info("$data/empty.pdb");
	is_deeply($empty->{chain_order}, [], 'an empty file has no chains');
	is($empty->{stats}{n_atoms}, 0, 'and no atoms');
	lives_ok { structure_summary($empty) } 'and can still be summarised';
	is_deeply(structure_sequences($empty), {}, 'and has no sequences');
}

#--------
# the same file through a string
#--------
{
	open my $fh, '<', "$data/mini.pdb" or die $!;
	my $text = do { local $/; <$fh> };
	close $fh;
	my $s = structure_info_string($text);
	is($s->{id}, '9XYZ', 'structure_info_string: reads a structure from a string');
	is($s->{chains}{A}{sequence}, $info->{chains}{A}{sequence},
		'and gets the same sequence as reading the file');
	is($s->{file}, undef, 'with no file name to record');
}

#--------
# pdb_info() is structure_info() with the format settled
#--------
is_deeply(pdb_info("$data/mini.pdb"), $info, 'pdb_info: the same answer as structure_info');

done_testing();
