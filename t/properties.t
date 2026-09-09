#!/usr/bin/env perl
# The surface of structure_features(), structure_sasa() and
# structure_pi_stacking(): every option, every call form, every croak, and what
# each of them leaves behind in $info.
#
# What the numbers *are* is t/features.t's job -- it compares them against
# mdtraj and Biopython.  This file is about the parts that are this module's own
# and that no other implementation has an opinion on.
require 5.010;
use strict;
use warnings FATAL => 'all';
use Cwd 'abs_path';
use File::Basename 'dirname';
use File::Temp 'tempdir';
use Test::More;
use Test::Exception;
use Chem::Structure::Parser;

my $data = dirname(abs_path(__FILE__)) . '/data';

# ---- what comes back -----------------------------------------------------
{
	my $i = structure_info("$data/stack.pdb");
	my $f = structure_features($i);

	is($f->{n_atoms},   84, 'n_atoms is the atoms the walk found');
	is($f->{n_residues}, 7, 'n_residues likewise');
	is($f->{n_chains},   1, 'and n_chains');
	is($f->{n_no_element}, 0, 'every atom of a modern file has an element');
	ok(exists $f->{sasa} && exists $f->{pi_stacking}, 'both calculations ran by default');
	is(scalar @{ $f->{center} },         3, 'the centroid is a triple');
	is(scalar @{ $f->{center_of_mass} }, 3, 'so is the centre of mass');

	# the surface adds up three ways: over atoms, over residues, and the split
	# by element class
	my ($atoms, $residues) = (0, 0);
	my $c = $i->{chains}{A};
	for my $rk (@{ $c->{residue_order} }) {
		my $r = $c->{residues}{$rk};
		$residues += $r->{sasa};
		$atoms += $r->{atoms}{$_}{sasa} for @{ $r->{atom_order} };
	}
	cmp_ok(abs($atoms - $f->{sasa}{total}), '<', 1e-9,
		'the per-atom surfaces add up to the total');
	cmp_ok(abs($residues - $f->{sasa}{total}), '<', 1e-9,
		'and so do the per-residue ones');
	cmp_ok(abs($c->{sasa} - $f->{sasa}{total}), '<', 1e-9,
		'and the chain carries the same figure');
	cmp_ok(abs($f->{sasa}{apolar} + $f->{sasa}{polar} - $f->{sasa}{total}), '<', 1e-9,
		'apolar and polar account for all of it');
	cmp_ok($f->{sasa}{apolar}, '>', 0, 'some of it is apolar');
	cmp_ok($f->{sasa}{polar},  '>', 0, 'and some of it is polar');
	is($f->{sasa}{probe},  1.4, 'the probe radius used is reported back');
	is($f->{sasa}{points}, 960, 'and the number of sphere points');
}

# ---- relative accessibility ----------------------------------------------
{
	my $i = structure_info("$data/stack.pdb");
	structure_features($i);
	my $c = $i->{chains}{A};
	for my $rk (@{ $c->{residue_order} }) {
		my $r = $c->{residues}{$rk};
		ok(defined $r->{rsa}, "$r->{resname} $rk has a relative accessibility");
		# these seven residues are cut out of their protein and so are close to
		# fully exposed; Tien's maximum is a Gly-X-Gly tripeptide, which a bare
		# residue can exceed a little
		cmp_ok($r->{rsa}, '>', 0.5, "and it is high for a residue with no neighbours");
		cmp_ok($r->{rsa}, '<', 1.5, "and not absurd");
	}
}

# a nucleotide gets no relative accessibility: the single-letter codes of the
# nucleotides are amino acid codes too, and dividing a guanine's surface by
# glycine's maximum would be a number rather than an answer
{
	my $i = structure_info("$data/bases.pdb");
	structure_features($i);
	my $c = $i->{chains}{P};
	my @rsa = grep { defined } map { $c->{residues}{$_}{sasa} } @{ $c->{residue_order} };
	is(scalar @rsa, scalar @{ $c->{residue_order} }, 'every nucleotide has a surface');
	my @has = grep { defined $c->{residues}{$_}{rsa} } @{ $c->{residue_order} };
	is_deeply(\@has, [], 'and none of them has a relative accessibility');
}

# ---- store => 0 ----------------------------------------------------------
{
	my $i = structure_info("$data/stack.pdb");
	my $f = structure_features($i, store => 0);
	cmp_ok($f->{sasa}{total}, '>', 0, 'the total still comes back');
	my $c = $i->{chains}{A};
	ok(!exists $c->{sasa}, 'the chain was left alone');
	ok(!exists $c->{hydropathy}, 'and so was its hydropathy');
	my $r = $c->{residues}{ $c->{residue_order}[0] };
	ok(!exists $r->{sasa} && !exists $r->{rsa}, 'and the residue');
	ok(!exists $r->{atoms}{ $r->{atom_order}[0] }{sasa}, 'and the atom');
}

# ---- structure_sasa and structure_pi_stacking ----------------------------
{
	my $i = structure_info("$data/stack.pdb");
	my $s = structure_features($i)->{sasa};
	my $j = structure_info("$data/stack.pdb");
	my $t = structure_sasa($j);
	is_deeply($t, $s, 'structure_sasa returns what structure_features puts under sasa');
	ok(defined $j->{chains}{A}{sasa}, 'and stores the per-chain figure too');

	my $k = structure_info("$data/stack.pdb");
	my $p = structure_pi_stacking($k);
	is(scalar @$p, 4, 'structure_pi_stacking returns the pairs on their own');
	ok(!exists $k->{chains}{A}{sasa},
		'and computes no surface, because it was not asked for one');
	my %type;
	$type{ $_->{type} }++ for @$p;
	is_deeply(\%type, { face => 1, edge => 3 }, 'one face stack and three edge stacks');
	for my $e (@$p) {
		like($e->{type}, qr/\A(?:face|edge)\z/, 'each is a face or an edge stack');
		like($e->{ring1}, qr/\A[56]\z/, 'each names the size of the first ring');
		like($e->{ring2}, qr/\A[56]\z/, 'and of the second');
		cmp_ok($e->{distance}, '>', 0, 'each has a centroid distance');
		cmp_ok($e->{plane_angle}, '>=', 0, 'and a plane angle in degrees');
		cmp_ok($e->{plane_angle}, '<=', 90, 'folded into the first quadrant');
		cmp_ok($e->{normal_angle1}, '<=', 90, 'as are the two normal angles');
		cmp_ok($e->{normal_angle2}, '<=', 90, '...');
		ok(defined $e->{resname1} && defined $e->{residue1} && defined $e->{chain1},
			'and names the residue each ring belongs to');
		if ($e->{type} eq 'edge') {
			ok(defined $e->{intersect_distance},
				'an edge stack says how far the shared line passes from a centroid');
		} else {
			ok(!exists $e->{intersect_distance},
				'a face stack has no intersection to report');
		}
	}
}

# ---- the options do something -------------------------------------------
{
	my $i = structure_info("$data/stack.pdb");
	my $full  = structure_sasa($i)->{total};
	my $naked = structure_sasa($i, probe => 0)->{total};
	cmp_ok($naked, '<', $full,
		'probe => 0 gives the van der Waals surface, which is smaller than the accessible one');

	my $coarse = structure_sasa($i, points => 100)->{total};
	cmp_ok(abs($coarse - $full) / $full, '<', 0.05,
		'a hundred sphere points is within five percent of nine hundred and sixty');
	isnt($coarse, $full, 'but not the same number');
	is(structure_sasa($i, points => 100)->{points}, 100, 'and says which it used');

	# a bigger probe reaches less far in, so the surface shrinks
	my $fat = structure_sasa($i, probe => 3.0)->{total};
	cmp_ok($fat, '>', $full, 'a larger probe rolls over a larger sphere');
}

{
	my $i = structure_info("$data/stack.pdb");
	is(scalar @{ structure_pi_stacking($i, face_plane_max => 0) }, 3,
		'no ring pair is exactly coplanar, so only the edge stacks survive');
	is(scalar @{ structure_pi_stacking($i, edge_plane_min => 90) }, 1,
		'and requiring exactly perpendicular planes leaves only the face stack');
	cmp_ok(scalar @{ structure_pi_stacking($i, face_distance => 55) }, '>', 4,
		'mdtraj-as-shipped 55 A face cutoff finds pairs that are not stacked at all');
	is(scalar @{ structure_pi_stacking($i, edge_radius => 0) }, 1,
		'an edge stack needs the planes to meet near a ring');
}

# ---- both formats, same answer ------------------------------------------
for my $pair ([ 'stack.pdb', 'stack.cif' ], [ 'bases.pdb', 'bases.cif' ],
              [ 'mini.pdb', 'mini.cif' ]) {
	my ($a, $b) = @$pair;
	my $fa = structure_features(structure_info("$data/$a"));
	my $fb = structure_features(structure_info("$data/$b"));
	# the residue keys travel, so the pi-stacking lists compare directly
	is_deeply($fb, $fa, "$a and $b give the same properties");
}

# ---- structures with nothing in them ------------------------------------
{
	my $i = structure_info("$data/empty.pdb");
	my $f = structure_features($i);
	is($f->{n_atoms}, 0, 'an empty file has no atoms');
	is($f->{sasa}{total}, 0, 'and no surface');
	is_deeply($f->{pi_stacking}, [], 'and no stacked rings');
	ok(!exists $f->{rg}, 'and no radius of gyration, because there is nothing to gyrate');
	ok(!exists $f->{center}, 'and no centroid');
	is($f->{mass}, 0, 'its mass is zero');
	ok(!exists $f->{hydropathy}, 'and it has no sequence to be hydropathic');
}

# an element the table does not know is counted rather than passed over in
# silence: it gets mdtraj's 2.0 A default radius and no mass
{
	my $i = structure_info_string(
		"ATOM      1  XX  UNK A   1       0.000   0.000   0.000  1.00  0.00          Xx\nEND\n");
	my $f = structure_features($i);
	is($f->{n_atoms}, 1, 'the atom is there');
	is($f->{n_no_element}, 1, 'and is counted as having no element this table knows');
	is($f->{mass}, 0, 'so it weighs nothing');
	cmp_ok($f->{sasa}{total}, '>', 0, 'but it still has a surface, at the default radius');
}

# ---- what it refuses to do ----------------------------------------------
{
	my $i = structure_info("$data/stack.pdb", atoms => 0);
	throws_ok { structure_features($i) } qr/atoms => 0/,
		'a structure read without atom hashes says so rather than reporting no surface';
	throws_ok { structure_sasa($i) } qr/atoms => 0/, 'and structure_sasa too';
}

{
	my $i = structure_info("$data/stack.pdb");
	throws_ok { structure_features($i, hydrogen => 1) } qr/unknown option 'hydrogen'/,
		'a misspelt option is a typo, not something to ignore';
	throws_ok { structure_sasa($i, face_distance => 5) } qr/unknown option 'face_distance'/,
		'a pi-stacking threshold means nothing to structure_sasa';
	throws_ok { structure_pi_stacking($i, probe => 1.4) } qr/unknown option 'probe'/,
		'and a probe radius means nothing to structure_pi_stacking';
	throws_ok { structure_features($i, probe => -1) } qr/probe must be a number/,
		'a negative probe is refused';
	throws_ok { structure_features($i, probe => 'wide') } qr/probe must be a number/,
		'and so is one that is not a number';
	throws_ok { structure_features($i, points => 0) } qr/points must be an integer/,
		'zero sphere points is refused';
	throws_ok { structure_features($i, points => 'many') } qr/points must be an integer/,
		'and so is a word';
	throws_ok { structure_features($i, points => 10_000_001) } qr/points must be an integer/,
		'and so is a number of points whose sphere would not fit in memory';
	lives_ok { structure_features($i, points => 1) } 'one sphere point is allowed';
	throws_ok { structure_pi_stacking($i, face_plane_max => 'flat') } qr/face_plane_max must be a number/,
		'and an angle that is not a number';
}

for my $who (qw(structure_features structure_sasa structure_pi_stacking)) {
	no strict 'refs';
	throws_ok { &{"Chem::Structure::Parser::$who"}(undef) } qr/\Q$who\E: expected the hash/,
		"$who refuses undef";
	throws_ok { &{"Chem::Structure::Parser::$who"}({ a => 1 }) } qr/\Q$who\E: expected the hash/,
		"$who refuses a hash that is not a structure";
	throws_ok { &{"Chem::Structure::Parser::$who"}([]) } qr/\Q$who\E: expected the hash/,
		"$who refuses an array reference";
}

# ---- the reader's own options carry through ------------------------------
{
	my $all = structure_sasa(structure_info("$data/mini.pdb"))->{total};
	my $dry = structure_sasa(structure_info("$data/mini.pdb", waters => 0))->{total};
	my $noh = structure_sasa(structure_info("$data/mini.pdb", hydrogens => 0))->{total};
	my $one = structure_sasa(structure_info("$data/mini.pdb", chains => ['B']))->{total};
	isnt($dry, $all, 'dropping the waters changes the surface');
	isnt($noh, $all, 'and so does dropping the hydrogens');
	cmp_ok($one, '<', $all, 'and reading one chain leaves less of it');

	# an NMR ensemble: the properties are of the model the chains were built
	# from, which is the one $info->{model} names.  nmr.pdb's three models are
	# the same tripeptide moved along x, so the surface is the invariant and the
	# centroid is what moved -- which is the pair of assertions worth making.
	my $m1 = structure_features(structure_info("$data/nmr.pdb", model => 1));
	my $m3 = structure_features(structure_info("$data/nmr.pdb", model => 3));
	cmp_ok(abs($m1->{sasa}{total} - $m3->{sasa}{total}), '<', 1e-9,
		'a model that is a translation of another has the same surface');
	cmp_ok(abs($m1->{rg} - $m3->{rg}), '<', 1e-9, 'and the same radius of gyration');
	cmp_ok(abs($m1->{center}[0] - $m3->{center}[0]), '>', 1,
		'and a centroid two angstrom away, which is where it was moved to');
	my $ma = structure_features(structure_info("$data/nmr.pdb", model => 'all'));
	cmp_ok(abs($ma->{sasa}{total} - $m1->{sasa}{total}), '<', 1e-9,
		"model => 'all' reports the model whose chains \$info->{chains} holds");
}

# ---- a gzipped file ------------------------------------------------------
SKIP: {
	# the package is loaded at run time here, so $GzipError is a name perl sees
	# once and warns about; the warning is about this file, not about the module
	no warnings 'once';
	eval { require IO::Compress::Gzip; 1 } or skip 'IO::Compress is not installed', 1;
	my $dir = tempdir(CLEANUP => 1);
	IO::Compress::Gzip::gzip("$data/stack.pdb" => "$dir/stack.pdb.gz")
		or skip "cannot gzip the fixture: $IO::Compress::Gzip::GzipError", 1;
	is_deeply(structure_features(structure_info("$dir/stack.pdb.gz")),
	          structure_features(structure_info("$data/stack.pdb")),
		'a gzipped file gives the same properties as the file it was made from');
}

# ---- the neighbour grid --------------------------------------------------
#
# The surface is computed through a grid of cells one cutoff wide rather than by
# comparing every atom with every other, which is what mdtraj does; t/features.t
# is what proves the two find the same neighbours, because mdtraj's scan is
# exhaustive and the answers match to the last digit.
#
# What is worth checking here is that the answer does not depend on where the
# structure happens to sit in the grid.  Moving it by a fraction of a cell puts
# every atom in a different cell, and a pair whose two atoms end up more than
# one cell apart would be missed.  The offsets below are deliberately not
# multiples of anything: 6.56 A is the cell size for a structure of carbon,
# nitrogen, oxygen and sulphur (twice the largest radius plus the probe).
{
	my $ref;
	for my $shift ([ 0, 0, 0 ], [ 0.37, 0, 0 ], [ 0, 3.13, 0 ], [ 0, 0, -7.91 ],
	               [ 101.3, -55.7, 12.9 ], [ -1000.5, 2000.25, -3000.125 ]) {
		my $i = structure_info("$data/stack.pdb");
		for my $cid (@{ $i->{chain_order} }) {
			my $c = $i->{chains}{$cid};
			for my $rk (@{ $c->{residue_order} }) {
				my $r = $c->{residues}{$rk};
				for my $an (@{ $r->{atom_order} }) {
					my $a = $r->{atoms}{$an};
					$a->{x} += $shift->[0];
					$a->{y} += $shift->[1];
					$a->{z} += $shift->[2];
				}
			}
		}
		my $t = structure_sasa($i)->{total};
		$ref = $t unless defined $ref;
		# the coordinates are three-decimal numbers and the shifts are exact
		# binary fractions or close to them, so the arithmetic is the only
		# source of difference; 1e-6 A^2 of 1736 is 6e-10 relative
		cmp_ok(abs($t - $ref), '<', 1e-6,
			"moving the structure by (@{[ join ', ', @$shift ]}) does not change its surface");
	}
}

# ---- calling twice -------------------------------------------------------
#
# The results are written into $info, so a second call has to overwrite them
# rather than add to them.
{
	my $i = structure_info("$data/stack.pdb");
	my $a = structure_features($i);
	my $b = structure_features($i);
	is_deeply($b, $a, 'asking twice gives the same answer');
	my $c = structure_features($i, probe => 2.0);
	cmp_ok($i->{chains}{A}{sasa}, '>', $a->{sasa}{total} / 2,
		'and asking again with a different probe replaces what was stored');
	cmp_ok(abs($i->{chains}{A}{sasa} - $c->{sasa}{total}), '<', 1e-9,
		'with the new figure, not the old one');
}

done_testing();
