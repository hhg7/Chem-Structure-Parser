#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use Chem::Structure::Parser;
use Matplotlib::Simple;
use
use Util;

my $r = structure_rmsd("$ENV{HOME}/ui/pepPriML/PPB/PDB/PDBbind.v2020/2ll7.ent.pdb", model => 'all');
p $r;
mkdir 'svg' unless -d 'svg';
imshow(
	data => $r->{rmsd},
	'output.file' => 'svg/2ll7.rmsd.svg',
	cblabel => 'RMSD (Å)',
	title   => '2LL7: RMSD between structures',
	xlabel  => 'Structure A',
	ylabel  => 'Structure B'
);
#qw(json_file_to_ref ref_to_json_file);
#~/.local/share/libgedit-gtksourceview-300/language-specs/perl.lang

#my $f = structure_info('157d.pdb');
#st $f->{chains}{A};

# bare.pdb is the coordinates-and-nothing-else fixture: its atoms are 5.2 A
# apart and all on one straight line, so no residue is peptide-bonded to the
# next and no torsion is defined.  fold.pdb is a real backbone.
=my $s = structure_info('t/data/aform.pdb');
p $s;
=foreach my $chain (sort keys %{ $s->{chains} }) {
	my $c = $s->{chains}{$chain};
	my $t = $c->{torsions};	# one array per angle, by residue_order
	foreach my $i (0 .. $#{ $c->{residue_order} }) {
		my $r = $c->{residue_order}[$i];
		printf "%s %4s %-3s  phi %8s  psi %8s  omega %8s  chi %s\n",
			$chain, $r, $c->{residues}{$r}{resname},
			map({ defined $_ ? sprintf('%.1f', $_) : '-' }
				$t->{phi}[$i], $t->{psi}[$i], $t->{omega}[$i]),
			defined $t->{chi}[$i]
				? join ', ', map { defined $_ ? sprintf '%.1f', $_ : '-' } @{ $t->{chi}[$i] }
				: '-';
	}
}
p $s->{chains};
