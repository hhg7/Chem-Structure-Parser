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

# --- the RMSD, and the same RMSD in pure Perl -----------------------------
#
# The job structure_rmsd() does per pair of structures is: walk both into flat
# coordinate arrays, pair the atoms by chain, residue and atom name, and then
# three passes over the pairs for the centroids, the inner-product matrix and
# the deviation.  An NMR ensemble asks for it once per pair of models, which is
# 190 times for the twenty models the archive usually deposits, so this is the
# per-atom work the module's own rule says belongs in C.  perl_rmsd() below
# does the same thing the obvious way, with the same quaternion arithmetic, so
# that the difference is the language and not the method.
sub perl_pairs {
	my ($info) = @_;
	my (@x, @y, @z, @k);
	for my $cid (@{ $info->{chain_order} }) {
		my $c = $info->{chains}{$cid};
		for my $rk (@{ $c->{residue_order} }) {
			my $res = $c->{residues}{$rk};
			for my $an (@{ $res->{atom_order} }) {
				my $a = $res->{atoms}{$an};
				next unless defined $a->{x};
				push @x, $a->{x}; push @y, $a->{y}; push @z, $a->{z};
				push @k, "$cid|$rk|$an";
			}
		}
	}
	return (\@x, \@y, \@z, \@k);
}

sub perl_rmsd {
	my ($a, $b) = @_;
	my ($ax, $ay, $az, $ak) = @$a;
	my ($bx, $by, $bz, $bk) = @$b;
	my %where;
	$where{ $bk->[$_] } = $_ for 0 .. $#$bk;
	my (@i, @j);
	for my $n (0 .. $#$ak) {
		next unless exists $where{ $ak->[$n] };
		push @i, $n;
		push @j, $where{ $ak->[$n] };
	}
	return undef unless @i >= 3;
	my $n = scalar @i;
	my (@ac, @bc) = ((0, 0, 0), (0, 0, 0));
	@ac = (0, 0, 0); @bc = (0, 0, 0);
	for my $p (0 .. $n - 1) {
		$ac[0] += $ax->[$i[$p]]; $ac[1] += $ay->[$i[$p]]; $ac[2] += $az->[$i[$p]];
		$bc[0] += $bx->[$j[$p]]; $bc[1] += $by->[$j[$p]]; $bc[2] += $bz->[$j[$p]];
	}
	$_ /= $n for @ac, @bc;
	my @s = (0) x 9;
	my $e0 = 0;
	for my $p (0 .. $n - 1) {
		my @u = ($ax->[$i[$p]] - $ac[0], $ay->[$i[$p]] - $ac[1], $az->[$i[$p]] - $ac[2]);
		my @v = ($bx->[$j[$p]] - $bc[0], $by->[$j[$p]] - $bc[1], $bz->[$j[$p]] - $bc[2]);
		for my $r (0 .. 2) { $s[3 * $r + $_] += $u[$r] * $v[$_] for 0 .. 2 }
		$e0 += $u[0] ** 2 + $u[1] ** 2 + $u[2] ** 2 + $v[0] ** 2 + $v[1] ** 2 + $v[2] ** 2;
	}
	my @rot = perl_qcp(\@s, $e0 / 2);
	my @t;
	for my $r (0 .. 2) {
		my $v = $bc[$r];
		$v -= $rot[3 * $r + $_] * $ac[$_] for 0 .. 2;
		push @t, $v;
	}
	my $sum = 0;
	for my $p (0 .. $n - 1) {
		my @u = ($ax->[$i[$p]], $ay->[$i[$p]], $az->[$i[$p]]);
		my @w = ($bx->[$j[$p]], $by->[$j[$p]], $bz->[$j[$p]]);
		for my $r (0 .. 2) {
			my $d = $t[$r] - $w[$r];
			$d += $rot[3 * $r + $_] * $u[$_] for 0 .. 2;
			$sum += $d * $d;
		}
	}
	return sqrt($sum / $n);
}

# the quaternion characteristic polynomial, spelled the same way the XS spells
# it, so that the two are the same arithmetic in two languages
sub perl_qcp {
	my ($s, $e0) = @_;
	my ($xx, $xy, $xz, $yx, $yy, $yz, $zx, $zy, $zz) = @$s;
	my $c2 = -2 * ($xx**2 + $yy**2 + $zz**2 + $xy**2 + $yx**2 + $xz**2 + $zx**2 + $yz**2 + $zy**2);
	my $c1 = 8 * ($xx * $yz * $zy + $yy * $zx * $xz + $zz * $xy * $yx
	            - $xx * $yy * $zz - $yz * $zx * $xy - $zy * $yx * $xz);
	my $d1 = $xy**2 + $xz**2 - $yx**2 - $zx**2;
	my $d2 = $yy**2 + $zz**2 - $xx**2 + $yz**2 + $zy**2;
	my $d3 = 2 * ($yz * $zy - $yy * $zz);
	my $c0 = $d1**2 + ($d2 + $d3) * ($d2 - $d3)
	  + (-($xz + $zx) * ($yz - $zy) + ($xy - $yx) * ($xx - $yy - $zz))
	  * (-($xz - $zx) * ($yz + $zy) + ($xy - $yx) * ($xx - $yy + $zz))
	  + (-($xz + $zx) * ($yz + $zy) - ($xy + $yx) * ($xx + $yy - $zz))
	  * (-($xz - $zx) * ($yz - $zy) - ($xy + $yx) * ($xx + $yy + $zz))
	  + (($xy + $yx) * ($yz + $zy) + ($xz + $zx) * ($xx - $yy + $zz))
	  * (-($xy - $yx) * ($yz - $zy) + ($xz + $zx) * ($xx + $yy + $zz))
	  + (($xy + $yx) * ($yz - $zy) + ($xz - $zx) * ($xx - $yy - $zz))
	  * (-($xy - $yx) * ($yz + $zy) + ($xz - $zx) * ($xx + $yy - $zz));
	my $l = $e0;
	for (1 .. 50) {
		my $old = $l;
		my $l2 = $l * $l;
		my $b = ($l2 + $c2) * $l;
		my $a = $b + $c1;
		my $fp = 2 * $l2 * $l + $b + $a;
		last unless abs($fp) > 0;
		$l = abs($l - ($a * $l + $c0) / $fp);
		last if abs($l - $old) <= 4e-16 * abs($l);
	}
	my @a = (
		$xx + $yy + $zz - $l, $yz - $zy,            -($xz - $zx),          $xy - $yx,
		$yz - $zy,            $xx - $yy - $zz - $l, $xy + $yx,             $xz + $zx,
		-($xz - $zx),         $xy + $yx,            $yy - $xx - $zz - $l,  $yz + $zy,
		$xy - $yx,            $xz + $zx,            $yz + $zy,             $zz - $xx - $yy - $l,
	);
	my $m3344 = $a[10] * $a[15] - $a[14] * $a[11];
	my $m3244 = $a[9]  * $a[15] - $a[13] * $a[11];
	my $m3243 = $a[9]  * $a[14] - $a[13] * $a[10];
	my $m3143 = $a[8]  * $a[14] - $a[12] * $a[10];
	my $m3144 = $a[8]  * $a[15] - $a[12] * $a[11];
	my $m3142 = $a[8]  * $a[13] - $a[12] * $a[9];
	my @q = ( $a[5] * $m3344 - $a[6] * $m3244 + $a[7] * $m3243,
	         -$a[4] * $m3344 + $a[6] * $m3144 - $a[7] * $m3143,
	          $a[4] * $m3244 - $a[5] * $m3144 + $a[7] * $m3142,
	         -$a[4] * $m3243 + $a[5] * $m3143 - $a[6] * $m3142 );
	my $qs = $q[0]**2 + $q[1]**2 + $q[2]**2 + $q[3]**2;
	return (1, 0, 0, 0, 1, 0, 0, 0, 1) unless $qs > 0;
	$_ /= sqrt $qs for @q;
	my ($w, $x, $y, $z) = @q;
	return ($w*$w + $x*$x - $y*$y - $z*$z, 2*($x*$y - $w*$z),            2*($z*$x + $w*$y),
	        2*($x*$y + $w*$z),             $w*$w - $x*$x + $y*$y - $z*$z, 2*($y*$z - $w*$x),
	        2*($z*$x - $w*$y),             2*($y*$z + $w*$x),             $w*$w - $x*$x - $y*$y + $z*$z);
}

{
	# the first ensemble of the set, or nothing to say if it holds none
	my ($ens, $nm);
	for my $f (@files) {
		my $i = structure_info($f, meta => 0, features => 0, atoms => 0);
		next unless ($i->{n_models} || 1) > 1;
		$ens = $f;
		$nm = $i->{n_models};
		last;
	}
	if (!defined $ens) {
		print "\n  no NMR ensemble in this set, so no RMSD comparison\n";
	}
	else {
		my $info = structure_info($ens, meta => 0, features => 0, model => 'all');
		my @m = map { $info->{models}{$_} } sort { $a <=> $b } keys %{ $info->{models} };
		my $t0 = time;
		my $c = Chem::Structure::Parser::_rmsd(\@m, {});
		my $c_time = time - $t0;
		$t0 = time;
		my @pre = map { [ perl_pairs($_) ] } @m;
		my @p;
		for my $i (0 .. $#pre) {
			for my $j ($i + 1 .. $#pre) { push @p, perl_rmsd($pre[$i], $pre[$j]) }
		}
		my $p_time = time - $t0;
		my $worst = 0;
		my $k = 0;
		for my $i (0 .. $#pre) {
			for my $j ($i + 1 .. $#pre) {
				my $d = abs($p[$k++] - $c->{rmsd}[$i][$j]);
				$worst = $d if $d > $worst;
			}
		}
		printf "\n  %s, %d models of %d atoms, %d pairs:\n",
			$ens, $nm, $c->{n_atoms}[0], $nm * ($nm - 1) / 2;
		printf "  %-38s %6.2f s\n", 'structure_rmsd', $c_time;
		printf "  %-38s %6.2f s\n", 'the same matrix, in pure Perl', $p_time;
		printf "  the C is %.0fx the Perl, and the worst pair differs by %.3g A\n",
			$p_time / ($c_time || 1e-9), $worst;
	}
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
