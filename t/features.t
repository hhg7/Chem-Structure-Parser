#!/usr/bin/env perl
# The physical properties, against the implementations they came from.
#
# structure_features(), structure_sasa() and structure_pi_stacking() are
# translations: the Shrake-Rupley surface, the van der Waals radii, the ring
# geometry and the atomic masses are all mdtraj's, and the two sequence numbers
# are Biopython's.  So this compares against those, not against what this module
# currently does.
#
#   mdtraj 1.11.1 -- mdtraj.geometry.shrake_rupley, mdtraj.geometry.pi_stacking,
#     mdtraj.geometry.compute_rg and mdtraj/core/element.py.  What it answered
#     for every structure in t/data is frozen in t/data/features.txt, written by
#     t/data/features.pl from t/data/features.py; the test reads that so the
#     comparison runs on a machine with no python, and re-runs features.py where
#     mdtraj is importable so the frozen answer cannot go stale.
#   Biopython 1.87 -- Bio.SeqUtils.ProtParam.ProteinAnalysis.gravy() and
#     .aromaticity(), and the Kyte-Doolittle scale of
#     Bio.SeqUtils.ProtParamData.kd.  Those two take a sequence rather than a
#     structure, so their answers are written out below as fixtures the way
#     t/foreign.t writes its cases down, with the sequence each came from.
#
# Set STRUCTURE_INFO_PYTHON to a python that can import mdtraj and numpy to run
# the live half as well; /home/con/.pyenv/versions/3.14.2/bin/python3 is one.
require 5.010;
use strict;
use warnings FATAL => 'all';
use Test::More;
use Test::Exception;
use File::Basename 'dirname';
use File::Spec;
use Chem::Structure::Parser qw(
	structure_info structure_features structure_sasa structure_pi_stacking
);

my $dir = File::Spec->catdir(dirname(__FILE__), 'data');

# ---- the frozen answer ---------------------------------------------------

# One section per structure, '== <file>' or '== <file> refused'.  A refused
# file is one mdtraj will not read at all and there is nothing to compare.
sub read_frozen {
	my ($path) = @_;
	open my $fh, '<', $path or die "$path: $!";
	my (%sec, $cur);
	while (my $l = <$fh>) {
		chomp $l;
		if ($l =~ /\A== (\S+)(?:\s+(refused))?\z/) {
			$cur = $1;
			$sec{$cur} = $2 ? undef : [];
			next;
		}
		next unless defined $cur && defined $sec{$cur};
		next if $l =~ /\A#/ || !length $l;
		push @{ $sec{$cur} }, $l;
	}
	close $fh;
	return \%sec;
}

my $frozen = read_frozen(File::Spec->catfile($dir, 'features.txt'));
ok(scalar(keys %$frozen) > 0, 'features.txt holds answers to compare against');

# every structure in t/data has an answer, so that adding one without re-running
# the generator is a failure rather than a silently smaller test
{
	opendir(my $dh, $dir) or die "$dir: $!";
	my @have = sort grep { /\.(?:pdb|ent|cif|mmcif)\z/ } readdir $dh;
	closedir $dh;
	my @missing = grep { !exists $frozen->{$_} } @have;
	is_deeply(\@missing, [], 'features.txt has an answer for every structure in t/data')
		or diag('re-run t/data/features.pl with STRUCTURE_INFO_PYTHON set');
}

# ---- what this module says -----------------------------------------------

# The residues of a structure in the order the file wrote them, which is the
# order mdtraj indexes them in.  Chains are not part of the key: mdtraj takes an
# mmCIF file's chains from label_asym_id rather than auth_asym_id, so mini.cif's
# one protein chain is three chains to it, in the same order.
sub walk {
	my ($info) = @_;
	my @res;
	for my $cid (@{ $info->{chain_order} }) {
		my $c = $info->{chains}{$cid};
		push @res, $c->{residues}{$_} for @{ $c->{residue_order} };
	}
	return \@res;
}

# a stacked pair, named the way the frozen file names it: the two residues'
# positions in that walk, and which of their rings it was
sub pair_key {
	my ($res, $s) = @_;
	my %at = map { $res->[$_]{chain} . '/' . $res->[$_]{key} => $_ } 0 .. $#$res;
	my @e = sort { $a->[0] <=> $b->[0] || $a->[1] cmp $b->[1] }
	        [ $at{"$s->{chain1}/$s->{residue1}"}, $s->{ring1} ],
	        [ $at{"$s->{chain2}/$s->{residue2}"}, $s->{ring2} ];
	return join ' ', map { defined $_->[0] ? "$_->[0]|$_->[1]" : '?' } @e;
}

# ---- the comparison ------------------------------------------------------
#
# Two tolerances, and they are different in kind.
#
# Against the float64 column -- mdtraj's own loop with the dtype changed -- the
# answers are the same calculation and agree to the last digit the generator
# prints.  The observed largest disagreement over every atom and residue of
# every structure in t/data is 4e-9 A^2 absolute, which is the 1e-9 the
# generator rounds to; 1e-7 leaves two orders of magnitude of headroom on that
# and would still catch a single sphere point, the smallest real difference
# there is (0.126 A^2 for a carbon).
#
# Against the float32 column -- mdtraj as it ships -- an atom can differ by a
# whole sphere point, because a point sitting within a float32 ulp of a
# neighbouring atom's surface is accessible at one width and covered at the
# other.  The generator prints what one point is worth for each atom, so the
# bound is one point and does not need a copy of the radius table here.  Over
# t/data that happens to four atoms of 620; the test allows one point everywhere
# rather than recording which four, because which four is a property of float32
# rounding and not of either implementation.
#
# Plus a millionth of the area, because mdtraj's number is itself a float32: an
# atom's area is a count of points times the point value, and both the product
# and the sum that follows it are rounded to 24 bits, so an atom that really is
# one point apart comes out 1.0000025 points apart on mini.pdb.  float32's
# relative epsilon is 1.19e-7, so 1e-6 is most of an order of magnitude of
# headroom and is still a hundredth of the point it is added to.
my $EXACT = 1e-7;

sub compare {
	my ($file, $rows) = @_;
	my $info = structure_info(File::Spec->catfile($dir, $file));
	my $feat = structure_features($info);
	my $res  = walk($info);

	my ($n_atom, $n_res, $worst64, $worst32) = (0, 0, 0, 0);
	my (%want_pi, %seen_res);
	for my $l (@$rows) {
		my @w = split ' ', $l;
		if ($w[0] eq 'A') {
			my ($ri, $name) = split /\|/, $w[1], 2;
			my ($a32, $a64, $pt) = @w[2, 3, 4];
			my $atom = $res->[$ri] ? $res->[$ri]{atoms}{$name} : undef;
			unless ($atom && defined $atom->{sasa}) {
				fail("$file: atom $w[1] is in mdtraj's answer and not in this one");
				next;
			}
			$n_atom++;
			my $d64 = abs($atom->{sasa} - $a64);
			my $d32 = abs($atom->{sasa} - $a32);
			$worst64 = $d64 if $d64 > $worst64;
			my $allow = $pt + 1e-6 * $a32;
			$worst32 = $d32 / $allow if $allow > 0 && $d32 / $allow > $worst32;
		} elsif ($w[0] eq 'R') {
			my ($ri, undef, $num, $name, $a32, $a64) = @w[1 .. 6];
			my $r = $res->[$ri];
			unless ($r) {
				fail("$file: residue $ri is in mdtraj's answer and not in this one");
				next;
			}
			$seen_res{$ri} = 1;
			$n_res++ if $r->{number} == $num && $r->{resname} eq $name;
			my $d64 = abs(($r->{sasa} || 0) - $a64);
			$worst64 = $d64 if $d64 > $worst64;
		} elsif ($w[0] eq 'T') {
			my $name = $w[1];
			if ($name eq 'sasa') {
				cmp_ok(abs($feat->{sasa}{total} - $w[3]), '<', $EXACT,
					"$file: total surface, against the float64 kernel");
				cmp_ok(abs($feat->{sasa}{total} - $w[2]), '<=', $w[4] + 1e-6 * $w[2],
					"$file: total surface, within one point per atom of mdtraj");
			} elsif ($name eq 'mass') {
				cmp_ok(abs($feat->{mass} - $w[2]), '<', 1e-6,
					"$file: mass, against mdtraj's element table");
			} elsif ($name eq 'rg') {
				# mdtraj's compute_rg() default is equal weights, so the centre
				# it measures from is the centroid and the two are the same
				# quantity.  Its xyz is float32 nanometres, which is why this is
				# a relative tolerance and the surfaces above are absolute.
				cmp_ok(abs($feat->{rg} - $w[2]) / $w[2], '<', 1e-6,
					"$file: radius of gyration");
			} elsif ($name eq 'rg_mass') {
				cmp_ok(abs($feat->{rg_mass} - $w[2]) / $w[2], '<', 1e-6,
					"$file: mass-weighted radius of gyration, about the centre of mass");
			} elsif ($name eq 'center' || $name eq 'com') {
				my $got = $feat->{ $name eq 'center' ? 'center' : 'center_of_mass' };
				my $off = 0;
				for my $k (0 .. 2) {
					my $d = abs($got->[$k] - $w[2 + $k]);
					$off = $d if $d > $off;
				}
				cmp_ok($off, '<', 1e-5, "$file: $name");
			}
		} elsif ($w[0] eq 'P') {
			$want_pi{"$w[1] $w[2]"} = 1;
		}
	}

	cmp_ok($worst64, '<', $EXACT,
		sprintf('%s: every atom and residue surface matches the float64 kernel (worst %.2e)',
			$file, $worst64));
	cmp_ok($worst32, '<=', 1.0,
		sprintf('%s: no atom differs from mdtraj by more than one sphere point', $file));
	is($n_res, scalar keys %seen_res,
		"$file: every residue mdtraj found is here, with the same number and name");
	ok($n_atom > 0, "$file: atoms were compared");

	my %got_pi = map { pair_key($res, $_) => 1 } @{ $feat->{pi_stacking} };
	is_deeply([ sort keys %got_pi ], [ sort keys %want_pi ],
		"$file: the same stacked pairs of rings as mdtraj");
}

for my $file (sort keys %$frozen) {
	next unless defined $frozen->{$file};
	compare($file, $frozen->{$file});
}

# a file mdtraj refuses is still a file this module has to survive
for my $file (sort grep { !defined $frozen->{$_} } keys %$frozen) {
	my $info = structure_info(File::Spec->catfile($dir, $file));
	lives_ok { structure_features($info) }
		"$file: mdtraj refuses it; this does not die on it";
}

# ---- the sequence numbers, against Biopython -----------------------------
#
# gravy() sums the Kyte-Doolittle index over the sequence and divides by its
# length; aromaticity() is the relative frequency of F, W and Y.  The sequences
# are the protein chains of t/data's structures as this module reads them, plus
# a few chosen to cover the twenty letters and the two extremes, and the numbers
# beside them are what Biopython 1.87 returns for them:
#
#   python3 -c 'from Bio.SeqUtils.ProtParam import ProteinAnalysis as P;
#               p = P("...."); print(p.gravy(), p.aromaticity())'
#
# Both are exact rational sums, so they are compared to 1e-12 rather than to a
# measured tolerance -- there is no floating point question here beyond the
# order the terms are added in.
my @protparam = (
	# sequence,                    gravy,                aromaticity
	[ 'MAGCMHHSC',                  0.33333333333333326, 0.0    ],
	[ 'MKWVTFISLLFLFSSAYS',         1.2333333333333334,  0.2777777777777778 ],
	[ 'ACDEFGHIKLMNPQRSTVWY',      -0.49000000000000005, 0.15000000000000002 ],
	[ 'FWY',                        0.19999999999999993, 1.0    ],
	[ 'GGGGG',                     -0.4,                 0.0    ],
	[ 'MRNQELARIFEEIGLMSEFLGDNPFRVRAYHQAARTLYDLDTPIEEIAEKGKEALMELPGVGPDLAEKILEFLRTG',
	                               -0.3197368421052632,  0.07894736842105264 ],
);
my %KD = (
	A =>  1.8, R => -4.5, N => -3.5, D => -3.5, C =>  2.5,
	Q => -3.5, E => -3.5, G => -0.4, H => -3.2, I =>  4.5,
	L =>  3.8, K => -3.9, M =>  1.9, F =>  2.8, P => -1.6,
	S => -0.8, T => -0.7, W => -0.9, Y => -1.3, V =>  4.2,
);
for my $c (@protparam) {
	my ($seq, $gravy, $arom) = @$c;
	my ($sum, $n, $a) = (0, 0, 0);
	for my $aa (split //, $seq) {
		$a++ if $aa =~ /[FWY]/;
		next unless exists $KD{$aa};
		$sum += $KD{$aa};
		$n++;
	}
	cmp_ok(abs($sum / $n - $gravy), '<', 1e-12, "Biopython gravy of $seq");
	cmp_ok(abs($a / length($seq) - $arom), '<', 1e-12, "Biopython aromaticity of $seq");
}

# and the same two numbers as structure_features() computes them, over the
# sequence structure_info() read out of the file
{
	my $info = structure_info(File::Spec->catfile($dir, 'mini.pdb'));
	my $f = structure_features($info);
	is($info->{chains}{A}{sequence}, 'MAGCMHHSC', 'mini.pdb chain A reads as expected');
	cmp_ok(abs($f->{hydropathy} - 0.33333333333333326), '<', 1e-12,
		'structure_features hydropathy is gravy over that sequence');
	cmp_ok(abs($f->{aromatic_fraction} - 0), '<', 1e-12,
		'structure_features aromatic_fraction is aromaticity over that sequence');
	is($f->{n_aromatic}, 0, 'no aromatic residues in mini.pdb');
	is($f->{sequence_length}, 9, 'nine residues of protein sequence');
	cmp_ok(abs($info->{chains}{A}{hydropathy} - $f->{hydropathy}), '<', 1e-12,
		'the one protein chain carries the same hydropathy as the structure');

	my $stack = structure_info(File::Spec->catfile($dir, 'stack.pdb'));
	my $sf = structure_features($stack);
	# TRP TRP TRP PHE PHE and two HIS: five of seven are F, W or Y
	is($sf->{n_aromatic}, 5, 'stack.pdb has five F/W/Y residues');
	cmp_ok(abs($sf->{aromatic_fraction} - 5 / 7), '<', 1e-12,
		'aromatic_fraction counts histidine in the denominator and not the numerator');
}

# ---- the live half -------------------------------------------------------
#
# Where mdtraj is importable, run the generator's own dump again and check the
# frozen file still says what mdtraj says.  A fixture regenerated without
# re-running t/data/features.pl would otherwise leave every comparison above
# passing against an answer to an older file.
SKIP: {
	my $py = $ENV{STRUCTURE_INFO_PYTHON};
	skip 'set STRUCTURE_INFO_PYTHON to a python with mdtraj to check the frozen answer', 1
		unless defined $py && length $py;
	my $devnull = File::Spec->devnull;
	skip "$py cannot import mdtraj and numpy", 1
		unless system("$py -c 'import mdtraj, numpy' > $devnull 2> $devnull") == 0;
	my $dump = File::Spec->catfile($dir, 'features.py');
	my @stale;
	for my $file (sort keys %$frozen) {
		next unless defined $frozen->{$file};
		my $path = File::Spec->catfile($dir, $file);
		my @live = grep { length && !/\A#/ }
		           split /\n/, qx($py \Q$dump\E \Q$path\E 2> $devnull);
		push @stale, $file unless "@live" eq "@{ $frozen->{$file} }";
	}
	is_deeply(\@stale, [], 'the frozen answers are still what mdtraj says')
		or diag('re-run t/data/features.pl from t/data/');
}

done_testing();
