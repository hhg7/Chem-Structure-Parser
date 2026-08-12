#!/usr/bin/env perl
# Read structure files and print what is in them.
#
#     structure.info.pl 1a22.ent.pdb                  # a summary
#     structure.info.pl --fasta *.pdb                 # the sequences, as FASTA
#     structure.info.pl --tsv   *.pdb                 # one row per structure
#     structure.info.pl --ligands *.pdb               # what is bound to what
#     structure.info.pl --dump 1a22.ent.pdb           # the whole hash of hashes
require 5.010;
use strict;
use warnings FATAL => 'all';
use Getopt::Long;
use Chem::Structure::Parser;

my %opt = (fasta => 0, tsv => 0, ligands => 0, dump => 0, seqres => 0, chain => undef);
GetOptions(\%opt, 'fasta', 'tsv', 'ligands', 'dump', 'seqres', 'chain=s', 'help')
	or die "see --help\n";

if ($opt{help} || !@ARGV) {
	print <<'USAGE';
usage: structure.info.pl [options] file.pdb ...

    --fasta      print every chain's sequence as FASTA
    --seqres     with --fasta, print SEQRES rather than the observed sequence
    --tsv        one tab separated row per structure
    --ligands    one row per bound heterogen
    --dump       the whole hash of hashes, via Data::Dumper
    --chain ID   only this chain

With no option, prints a readable summary of each file.
USAGE
	exit(!@ARGV);
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
