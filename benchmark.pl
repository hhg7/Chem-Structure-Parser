#!/usr/bin/env perl
# What the XS is worth, measured rather than asserted.
#
#     perl -Mblib benchmark.pl /path/to/structures
#
# Compares three ways of reading the same files: the module, the module with
# the filters on, and a straightforward pure-Perl reader written the way one
# would write it if there were no XS.  The Perl reader below is not a straw
# man -- it slices the same columns, in the same order, and builds the same
# nested hash -- it is simply doing per line in Perl what the C does per file.
require 5.010;
use strict;
use warnings FATAL => 'all';
use Time::HiRes 'time';
use Chem::Structure::Parser;

my $dir = shift || die "usage: benchmark.pl <directory of .pdb files> [count]\n";
my $want = shift || 200;

opendir(my $dh, $dir) or die "$dir: $!";
my @files = sort grep { /\.(pdb|ent)\z/ } readdir $dh;
closedir $dh;
die "no structures in $dir\n" unless @files;
my $step = @files > $want ? int(@files / $want) : 1;
@files = map { "$dir/$files[$_]" } grep { $_ % $step == 0 } 0 .. $#files;
@files = @files[0 .. $want - 1] if @files > $want;

# --- a pure Perl reader, for the comparison -------------------------------
sub perl_read {
	my ($file) = @_;
	open my $fh, '<', $file or die "$file: $!";
	my %chains;
	my $model = 0;
	while (my $l = <$fh>) {
		if (substr($l, 0, 5) eq 'MODEL') { last if $model++; next }
		next unless substr($l, 0, 6) eq 'ATOM  ' || substr($l, 0, 6) eq 'HETATM';
		my $chain   = substr($l, 21, 1);
		my $resname = substr($l, 17, 3);
		my $resseq  = substr($l, 22, 4);
		my $icode   = substr($l, 26, 1);
		my $name    = substr($l, 12, 4);
		for ($chain, $resname, $resseq, $icode, $name) { s/\A\s+//; s/\s+\z// }
		my $key = "$resseq$icode";
		my $c = $chains{$chain} ||= { residues => {}, order => [] };
		my $r = $c->{residues}{$key};
		unless ($r) {
			$r = $c->{residues}{$key} = {
				resname => $resname, number => $resseq, icode => $icode,
				one => res1($resname), type => res_type($resname),
				atoms => {}, atom_order => [],
			};
			push @{ $c->{order} }, $key;
		}
		next if exists $r->{atoms}{$name};
		$r->{atoms}{$name} = {
			name      => $name,
			serial    => 0 + substr($l, 6, 5),
			x         => 0 + substr($l, 30, 8),
			y         => 0 + substr($l, 38, 8),
			z         => 0 + substr($l, 46, 8),
			occupancy => (length($l) > 59 ? 0 + substr($l, 54, 6) : undef),
			bfactor   => (length($l) > 65 ? 0 + substr($l, 60, 6) : undef),
			element   => (length($l) > 77 ? do { my $e = substr($l, 76, 2); $e =~ s/\s+//g; $e } : ''),
		};
		push @{ $r->{atom_order} }, $name;
	}
	close $fh;
	return \%chains;
}

# Best of three.  A single pass puts whichever variant runs last at a
# disadvantage -- the first draft of this script had reading with the header
# parsing switched off coming out slower than reading with it on, which is not
# a thing that can happen -- and the minimum is the measurement least polluted
# by whatever else the machine was doing.
sub timeit {
	my ($what, $code) = @_;
	my ($best, $atoms);
	for (1 .. 3) {
		my $t = time;
		$atoms = 0;
		$atoms += $code->($_) for @files;
		my $s = time - $t;
		$best = $s if !defined $best || $s < $best;
	}
	printf "  %-38s %6.2f s   %8.0f files/s   %9.0f atoms/s\n",
		$what, $best, scalar(@files) / $best, $atoms / $best;
	return $best;
}

printf "%d files from %s\n\n", scalar @files, $dir;

# features => 0 on every read below: this half of the benchmark is about the
# parse, and comparing a read that computes the solvent-accessible surface with
# a pure-Perl reader that does not would be comparing two different jobs.  What
# the default costs is the structure_info line in the second half.
my $xs = timeit('structure_info, features => 0', sub {
	my $i = structure_info($_[0], features => 0);
	return $i->{stats}{n_atoms};
});
timeit('structure_info, atoms => 0', sub {
	my $i = structure_info($_[0], atoms => 0);
	return $i->{stats}{n_atoms};
});
timeit('structure_info, no water or hydrogen', sub {
	my $i = structure_info($_[0], waters => 0, hydrogens => 0, features => 0);
	return $i->{stats}{n_atoms};
});
timeit('structure_info, meta => 0', sub {
	my $i = structure_info($_[0], meta => 0, features => 0);
	return $i->{stats}{n_atoms};
});
my $raw = timeit('the XS parse alone', sub {
	my $p = Chem::Structure::Parser::_parse_file($_[0], {});
	return $p->{n_atoms};
});
my $pp = timeit('the same reader, in pure Perl', sub {
	my $c = perl_read($_[0]);
	my $n = 0;
	for my $ch (values %$c) {
		$n += scalar @{ $ch->{residues}{$_}{atom_order} } for keys %{ $ch->{residues} };
	}
	return $n;
});
my $pps = timeit('pure Perl, and the same statistics', sub {
	my $c = perl_read($_[0]);
	my (%elements, @bbox, $bmin, $bmax, $bsum, $bn, $n);
	for my $ch (values %$c) {
		for my $rk (keys %{ $ch->{residues} }) {
			my $r = $ch->{residues}{$rk};
			my ($sx, $sy, $sz, $nc) = (0, 0, 0, 0);
			for my $an (@{ $r->{atom_order} }) {
				my $a = $r->{atoms}{$an};
				$n++;
				$elements{ $a->{element} }++ if length $a->{element};
				if (defined $a->{x}) {
					$sx += $a->{x}; $sy += $a->{y}; $sz += $a->{z}; $nc++;
					if (!@bbox) { @bbox = (($a->{x}) x 2, ($a->{y}) x 2, ($a->{z}) x 2) }
					else {
						$bbox[0] = $a->{x} if $a->{x} < $bbox[0];
						$bbox[1] = $a->{x} if $a->{x} > $bbox[1];
						$bbox[2] = $a->{y} if $a->{y} < $bbox[2];
						$bbox[3] = $a->{y} if $a->{y} > $bbox[3];
						$bbox[4] = $a->{z} if $a->{z} < $bbox[4];
						$bbox[5] = $a->{z} if $a->{z} > $bbox[5];
					}
				}
				if (defined $a->{bfactor}) {
					$bsum += $a->{bfactor}; $bn++;
					$bmin = $a->{bfactor} if !defined $bmin || $a->{bfactor} < $bmin;
					$bmax = $a->{bfactor} if !defined $bmax || $a->{bfactor} > $bmax;
				}
			}
			$r->{center} = [ $sx / $nc, $sy / $nc, $sz / $nc ] if $nc;
		}
	}
	return $n;
});

# --- the physical properties ----------------------------------------------
#
# The surface is the expensive one: nine hundred and sixty sphere points per
# atom, each tested against the atom's neighbours.  Timed over the same files so
# that the cost can be read against the cost of reading them in the first place.
my $feat = timeit('structure_info, the default', sub {
	my $i = structure_info($_[0]);
	return $i->{stats}{n_atoms};
});
timeit('structure_sasa alone', sub {
	my $i = structure_info($_[0], meta => 0, features => 0);
	structure_sasa($i);
	return $i->{stats}{n_atoms};
});
timeit('structure_pi_stacking alone', sub {
	my $i = structure_info($_[0], meta => 0, features => 0);
	structure_pi_stacking($i);
	return $i->{stats}{n_atoms};
});
timeit('structure_disulfides alone', sub {
	my $i = structure_info($_[0], meta => 0, features => 0);
	structure_disulfides($i);
	return $i->{stats}{n_atoms};
});
timeit('structure_base_stacks alone', sub {
	my $i = structure_info($_[0], meta => 0, features => 0);
	structure_base_stacks($i, store => 0);
	return $i->{stats}{n_atoms};
});

# --- the same surface, in pure Perl ---------------------------------------
#
# One structure rather than the whole set: this is minutes per file, which is
# the answer the comparison is for.  Written the obvious way, which is also the
# way mdtraj's kernel is written -- every atom against every other for the
# neighbours, then every sphere point against every neighbour.
sub perl_sasa {
	my ($info, $npts) = @_;
	my (@x, @y, @z, @r);
	my %vdw = (C => 1.70, N => 1.55, O => 1.52, S => 1.80, P => 1.80,
	           H => 1.20, SE => 1.90, ZN => 1.39, FE => 2.00, MG => 0.86);
	for my $cid (@{ $info->{chain_order} }) {
		my $c = $info->{chains}{$cid};
		for my $rk (@{ $c->{residue_order} }) {
			my $res = $c->{residues}{$rk};
			for my $an (@{ $res->{atom_order} }) {
				my $a = $res->{atoms}{$an};
				next unless defined $a->{x};
				push @x, $a->{x}; push @y, $a->{y}; push @z, $a->{z};
				my $v = $vdw{ uc($a->{element} || '') };
				push @r, (defined $v ? $v : 2.0) + 1.4;
			}
		}
	}
	my $pi = 4 * atan2(1, 1);
	my $inc = $pi * (3 - sqrt 5);
	my $off = 2 / $npts;
	my (@px, @py, @pz);
	for my $i (0 .. $npts - 1) {
		my $yy = $i * $off - 1 + $off / 2;
		my $rr = sqrt(1 - $yy * $yy);
		push @px, cos($i * $inc) * $rr;
		push @py, $yy;
		push @pz, sin($i * $inc) * $rr;
	}
	my $total = 0;
	for my $i (0 .. $#x) {
		my @nb;
		for my $j (0 .. $#x) {
			next if $i == $j;
			my $s = $r[$i] + $r[$j];
			my $d = ($x[$j] - $x[$i]) ** 2 + ($y[$j] - $y[$i]) ** 2 + ($z[$j] - $z[$i]) ** 2;
			push @nb, $j if $d < $s * $s;
		}
		my $acc = 0;
		POINT: for my $k (0 .. $npts - 1) {
			my $qx = $x[$i] + $r[$i] * $px[$k];
			my $qy = $y[$i] + $r[$i] * $py[$k];
			my $qz = $z[$i] + $r[$i] * $pz[$k];
			for my $j (@nb) {
				my $d = ($qx - $x[$j]) ** 2 + ($qy - $y[$j]) ** 2 + ($qz - $z[$j]) ** 2;
				next POINT if $d < $r[$j] * $r[$j];
			}
			$acc++;
		}
		$total += $acc * 4 * $pi / $npts * $r[$i] ** 2;
	}
	return $total;
}

{
	# the smallest of the set, so the Perl version finishes
	my ($small, $n_small);
	for my $f (@files) {
		my $i = structure_info($f, meta => 0, features => 0);
		next if defined $n_small && $i->{stats}{n_atoms} >= $n_small;
		($small, $n_small) = ($i, $i->{stats}{n_atoms});
	}
	my $t0 = time;
	my $c_area = structure_sasa($small)->{total};
	my $c_time = time - $t0;
	$t0 = time;
	my $p_area = perl_sasa($small, 960);
	my $p_time = time - $t0;
	printf "\n  the smallest structure of the set, %d atoms:\n", $n_small;
	printf "  %-38s %6.2f s   %38.2f A^2\n", 'structure_sasa', $c_time, $c_area;
	printf "  %-38s %6.2f s   %38.2f A^2\n", 'the same surface, in pure Perl', $p_time, $p_area;
	printf "  the C is %.0fx the Perl, and the two areas differ by %.3g A^2\n",
		$p_time / ($c_time || 1e-9), abs($c_area - $p_area);
}

print <<"SUMMARY";

The parse alone is @{[ sprintf '%.1f', $pp / $raw ]}x the same parse written in Perl.

structure_info() comes out at @{[ sprintf '%.1f', $pp / $xs ]}x the plain reader and
@{[ sprintf '%.1f', $pps / $xs ]}x the one that also gathers the statistics -- and it is
doing a good deal more than either: the header records, SEQRES, the gaps
between them and the coordinates, chain types, ligand and ion classification.
What the C buys is the reading; building a hash of hashes out of what was read
costs what it costs, in any language, because the hashes are the answer.

Reading with the properties on -- which is the default -- costs @{[ sprintf '%.1f', $feat / $xs ]}x reading
without them, and is the one part of the module where the C is not competing with a plausible
Perl: nine hundred and sixty sphere points per atom against each of that atom's
neighbours is a loop nobody would write in Perl twice, which is what the
single-structure comparison above is there to show.
SUMMARY
