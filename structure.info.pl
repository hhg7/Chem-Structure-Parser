#!/usr/bin/env perl
# Read structure files and print what is in them.
#
#     structure.info.pl 1a22.ent.pdb                  # a summary
#     structure.info.pl --fasta *.pdb                 # the sequences, as FASTA
#     structure.info.pl --tsv   *.cif                 # one row per structure
#     structure.info.pl --ligands *.pdb *.cif         # what is bound to what
#     structure.info.pl --dump 1a22.ent.pdb           # the whole hash of hashes
#     structure.info.pl --features *.pdb              # surface, size, composition
#     structure.info.pl --stacks 1a22.ent.pdb         # stacked aromatic rings
#     structure.info.pl --pairs 1bna.pdb              # Watson-Crick base pairs
#     structure.info.pl --rmsd 2ll7.ent.pdb           # every model against every other
#     structure.info.pl --rmsd a.pdb b.cif c.pdb      # and every file against every other
#
# PDB and mmCIF are read the same way and print the same thing, so a mixed
# directory needs no sorting out first.
require 5.010;
use strict;
use warnings FATAL => 'all';
use Getopt::Long;
use Chem::Structure::Parser;

my %opt = (fasta => 0, tsv => 0, ligands => 0, dump => 0, seqres => 0,
           features => 0, stacks => 0, ss => 0, pairs => 0, rmsd => 0,
           chain => undef, select => 'all');
GetOptions(\%opt, 'fasta', 'tsv', 'ligands', 'dump', 'seqres', 'features', 'stacks',
           'ss', 'pairs', 'rmsd', 'chain=s', 'select=s', 'help')
	or die "see --help\n";

if ($opt{help} || !@ARGV) {
	print <<'USAGE';
usage: structure.info.pl [options] file.pdb|file.cif ...

    --fasta      print every chain's sequence as FASTA
    --seqres     with --fasta, print SEQRES rather than the observed sequence
    --tsv        one tab separated row per structure
    --ligands    one row per bound heterogen
    --dump       the whole hash of hashes, via Data::Dumper
    --features   one row per structure: surface, size, mass, and the numbers
                 the sequence answers -- hydropathy or G+C, as it has
    --stacks     one row per stacked pair of aromatic rings
    --ss         one row per disulfide, with what the file declares beside it
    --pairs      one row per Watson-Crick or wobble base pair
    --rmsd       the RMSD of every structure against every other, as a matrix.
                 A file with more than one model counts as one structure per
                 model, so one NMR entry on its own compares its models
    --select W   with --rmsd, which atoms take part: all, heavy, backbone or ca
    --chain ID   only this chain

With no option, prints a readable summary of each file.
USAGE
	exit(!@ARGV);
}

# --rmsd is the one question that is about the files together rather than about
# each of them in turn, so it is answered here and not in the loop below.
if ($opt{rmsd}) {
	my $r = structure_rmsd(@ARGV, model => 'all', select => $opt{select},
		(defined $opt{chain} ? (chains => [ $opt{chain} ]) : ()));
	# two structures give the number on its own; anything more gives the matrix
	$r = { labels => [ @ARGV[0, 1] ], rmsd => [ [ 0, $r ], [ $r, 0 ] ] }
		unless ref $r;
	my @lab = map { my $l = $_; $l =~ s{.*/}{}; $l } @{ $r->{labels} };
	my $w = 0;
	for my $l (@lab) { $w = length $l if length $l > $w }
	printf "%-*s %s\n", $w, '', join ' ', map { sprintf '%7d', $_ + 1 } 0 .. $#lab;
	for my $i (0 .. $#lab) {
		printf "%-*s %s\n", $w, $lab[$i], join ' ',
			map { defined $_ ? sprintf '%7.3f', $_ : sprintf '%7s', '-' }
			@{ $r->{rmsd}[$i] };
	}
	exit 0;
}

my $header_printed = 0;

for my $file (@ARGV) {
	my $info = eval { structure_info($file) };
	if (!$info) {
		warn "$file: $@";
		next;
	}

	if ($opt{dump}) {
		require Data::Dumper;
		no warnings 'once';    # the package variables are set, never read, here
		local $Data::Dumper::Sortkeys = 1;
		local $Data::Dumper::Indent   = 1;
		print Data::Dumper::Dumper($info);
		next;
	}

	if ($opt{features}) {
		my $f = $info->{features};
		unless ($header_printed++) {
			print join("\t", qw(id file n_atoms mass rg rg_mass sasa apolar polar
			                    hydropathy aromatic_fraction gc_fraction
			                    purine_fraction n_pi_stacking n_disulfides
			                    n_base_pairs n_base_stacks stacked)), "\n";
		}
		print join("\t",
			$info->{id} // '',
			$file,
			$f->{n_atoms},
			sprintf('%.1f', $f->{mass}),
			(defined $f->{rg}      ? sprintf('%.2f', $f->{rg})      : ''),
			(defined $f->{rg_mass} ? sprintf('%.2f', $f->{rg_mass}) : ''),
			sprintf('%.1f', $f->{sasa}{total}),
			sprintf('%.1f', $f->{sasa}{apolar}),
			sprintf('%.1f', $f->{sasa}{polar}),
			(defined $f->{hydropathy}        ? sprintf('%.4f', $f->{hydropathy})        : ''),
			(defined $f->{aromatic_fraction} ? sprintf('%.4f', $f->{aromatic_fraction}) : ''),
			# empty rather than zero for a structure with no nucleic acid in it,
			# the way hydropathy is empty for one with no protein
			(defined $f->{gc_fraction}     ? sprintf('%.4f', $f->{gc_fraction})     : ''),
			(defined $f->{purine_fraction} ? sprintf('%.4f', $f->{purine_fraction}) : ''),
			scalar @{ $f->{pi_stacking} },
			scalar @{ $f->{disulfides} },
			scalar @{ $f->{base_pairs} },
			scalar @{ $f->{base_stacks} },
			# and how many of those cleared the paper's 50%
			scalar(grep { $_->{stacked} } @{ $f->{base_stacks} }),
		), "\n";
		next;
	}

	if ($opt{ss}) {
		# what the coordinates show against what the depositor wrote down; a
		# bond in one column and not the other is a fact about the entry
		my $key = sub {
			my ($c1, $r1, $c2, $r2) = @_;
			return join '|', sort "$c1/$r1", "$c2/$r2";
		};
		my %found = map { $key->(@{$_}{qw(chain1 residue1 chain2 residue2)}) => $_->{distance} }
		            @{ structure_disulfides($info) };
		my %said  = map { $key->(@{$_}{qw(chain1 resseq1 chain2 resseq2)}) => $_->{length} }
		            @{ $info->{ssbond} || [] };
		for my $k (sort keys %{ { %found, %said } }) {
			next if defined $opt{chain} && $k !~ m{(?:\A|\|)\Q$opt{chain}\E/};
			printf "%s\t%s\t%s\t%s\n", $info->{id} // '', $k,
				(defined $found{$k} ? sprintf('%.3f', $found{$k}) : 'not found'),
				(defined $said{$k}  ? $said{$k}                   : 'not declared');
		}
		next;
	}

	if ($opt{pairs}) {
		for my $b (@{ structure_base_pairs($info) }) {
			next if defined $opt{chain}
			        && $b->{chain1} ne $opt{chain} && $b->{chain2} ne $opt{chain};
			printf "%s\t%s\t%d\t%s %s%s\t%s %s%s\t%s\t%.2f\t%.1f\t%.2f\n",
				$info->{id} // '', $b->{type}, $b->{saenger},
				$b->{chain1}, $b->{resname1}, $b->{residue1},
				$b->{chain2}, $b->{resname2}, $b->{residue2},
				join(' ', map { sprintf '%s-%s %.2f', $_->{atom1}, $_->{atom2},
				                        $_->{distance} } @{ $b->{hbonds} }),
				$b->{distance}, $b->{plane_angle}, $b->{stagger};
		}
		next;
	}

	if ($opt{stacks}) {
		for my $s (@{ structure_pi_stacking($info) }) {
			next if defined $opt{chain}
			        && $s->{chain1} ne $opt{chain} && $s->{chain2} ne $opt{chain};
			printf "%s\t%s\t%s %s%s ring %s\t%s %s%s ring %s\t%.2f\t%.1f\n",
				$info->{id} // '', $s->{type},
				$s->{chain1}, $s->{resname1}, $s->{residue1}, $s->{ring1},
				$s->{chain2}, $s->{resname2}, $s->{residue2}, $s->{ring2},
				$s->{distance}, $s->{plane_angle};
		}
		next;
	}

	if ($opt{fasta}) {
		for my $cid (@{ $info->{chain_order} }) {
			next if defined $opt{chain} && $cid ne $opt{chain};
			my $c   = $info->{chains}{$cid};
			my $seq = $opt{seqres} ? $c->{seqres} : $c->{sequence};
			next unless defined $seq && length $seq;
			printf ">%s_%s %s|%s%s\n", $info->{id}, $cid, $c->{type},
				(defined $c->{molecule} ? "$c->{molecule}|" : ''),
				length($seq) . 'aa';
			print "$1\n" while $seq =~ /(.{1,60})/g;
		}
		next;
	}

	if ($opt{tsv}) {
		unless ($header_printed++) {
			print join("\t", qw(id file method resolution r_work r_free n_models
			                    chains n_atoms n_residues sequence_length ligands title)), "\n";
		}
		my $residues = 0;
		my $seqlen   = 0;
		for my $cid (@{ $info->{chain_order} }) {
			$residues += $info->{chains}{$cid}{n_residues};
			$seqlen   += length $info->{chains}{$cid}{sequence};
		}
		my @lig = sort keys %{ structure_ligands($info) };
		print join("\t",
			$info->{id} // '',
			$file,
			join(';', @{ $info->{experiment} }),
			$info->{resolution} // '',
			$info->{r_work}     // '',
			$info->{r_free}     // '',
			$info->{n_models},
			join(';', @{ $info->{chain_order} }),
			$info->{stats}{n_atoms},
			$residues,
			$seqlen,
			join(';', @lig),
			$info->{title} // '',
		), "\n";
		next;
	}

	if ($opt{ligands}) {
		my $lig = structure_ligands($info);
		for my $key (sort keys %$lig) {
			my $r = $lig->{$key};
			printf "%s\t%s\t%s\t%s\t%d\t%s\n", $info->{id} // '', $key, $r->{resname},
				$r->{type}, $r->{n_atoms},
				($info->{het}{ $r->{resname} }{name} // '');
		}
		next;
	}

	print structure_summary($info);
	print "\n" if @ARGV > 1;
}
