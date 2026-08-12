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

my $xs = timeit('structure_info', sub {
	my $i = structure_info($_[0]);
	return $i->{stats}{n_atoms};
});
timeit('structure_info, atoms => 0', sub {
	my $i = structure_info($_[0], atoms => 0);
	return $i->{stats}{n_atoms};
});
timeit('structure_info, no water or hydrogen', sub {
	my $i = structure_info($_[0], waters => 0, hydrogens => 0);
	return $i->{stats}{n_atoms};
});
timeit('structure_info, meta => 0', sub {
	my $i = structure_info($_[0], meta => 0);
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

print <<"SUMMARY";

The parse alone is @{[ sprintf '%.1f', $pp / $raw ]}x the same parse written in Perl.

structure_info() comes out at @{[ sprintf '%.1f', $pp / $xs ]}x the plain reader and
@{[ sprintf '%.1f', $pps / $xs ]}x the one that also gathers the statistics -- and it is
doing a good deal more than either: the header records, SEQRES, the gaps
between them and the coordinates, chain types, ligand and ion classification.
What the C buys is the reading; building a hash of hashes out of what was read
costs what it costs, in any language, because the hashes are the answer.
SUMMARY
