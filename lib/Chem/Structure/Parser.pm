#!/usr/bin/env perl
# ABSTRACT: Read a molecular structure file into a hash of hashes, sequences and all, using XS for the coordinate section
require 5.010;
use strict;
package Chem::Structure::Parser;
our $VERSION = 0.03;
require XSLoader;
use autodie ':default';
use warnings FATAL => 'all';
use Exporter 'import';
use Scalar::Util 'reftype';
XSLoader::load('Chem::Structure::Parser', $VERSION);

our @EXPORT_OK = qw(
	structure_info structure_info_string pdb_info cif_info
	structure_atoms structure_residues structure_ligands structure_sequences
	chain_sequence structure_summary is_single_ion
	structure_features structure_sasa structure_pi_stacking structure_disulfides
	structure_contacts structure_hbonds
	aa3to1 aa1to3 res1 res_type formats h
);
our @EXPORT = @EXPORT_OK;

#
# Formats
#
# The module is named for structures, not for PDB, because the file format is
# an accident of history and the thing the caller wants -- chains, residues,
# a single-letter sequence, coordinates -- is the same whichever format it
# arrived in.  So there is one way in, structure_info(), which works out the
# format and hands the file to a reader; the returned hash of hashes has the
# same shape no matter which reader filled it.
#
# Adding a format means writing a reader that returns that shape and adding
# it here.  Nothing else in the module, and nothing in calling code, needs to
# know that a second format exists.
#
my %READER = (
	pdb   => \&_read_pdb,
	mmcif => \&_read_cif,
);

# The XS parser for each format.  Both fill the same hash -- the same column
# arrays, the same residue boundaries, the same counts -- so everything below
# this table is written once and reads either format.  What differs is the
# header: a PDB file hands its header back as raw lines by record name, an
# mmCIF file as tags and loops, and that is the one place the two part company.
my %XS = (
	pdb   => { file => \&_parse_file,     string => \&_parse_string     },
	mmcif => { file => \&_parse_cif_file, string => \&_parse_cif_string },
);

# formats a reader is not written for yet, kept here so that handing one over
# gets a straight answer rather than a puzzling parse of the wrong thing
my %NOT_YET = (
	mol2 => 'Tripos MOL2 (.mol2)',
	sdf  => 'MDL SDF/MOL (.sdf, .mol)',
);

# The names a format goes by, for format => .  A caller writing
# format => 'cif' means the format whose files are called .cif, and dying at
# them over the difference between that and 'mmcif' would be pedantry: there is
# only one thing they could have meant.  The name a format is filed under is
# still the one it reports back as, so $info->{format} has one spelling.
my %ALIAS = (
	cif => 'mmcif', pdbx => 'mmcif', mmcif => 'mmcif',
	ent => 'pdb',   pdb  => 'pdb',
);

# Options
#
# Anything not listed here is a typo, and a typo that is quietly ignored is a
# silent wrong answer later on -- pass 'hydrogen' for 'hydrogens' and you get
# a structure with the hydrogens still in it and no hint of why.
#
my %DEFAULT = (
	model     => 1,       # which MODEL to build chains from; 'all' for every one
	altloc    => 'first', # 'first' or 'highest' occupancy, when an atom has alternates
	hydrogens => 1,       # keep hydrogen/deuterium atoms
	waters    => 1,       # keep waters
	hetatm    => 1,       # keep HETATM records (ligands, ions, modified residues)
	atoms     => 1,       # build the per-atom hashes; 0 stops at the residue level
	features  => 1,       # compute the physical properties; needs atoms => 1
	meta      => 1,       # parse the header records
	anisou    => 0,       # keep ANISOU lines (they double the size of the file)
	chains    => undef,   # arrayref: read only these chains
	format    => undef,   # override format detection
);

# residues that need no explanation.  Everything else that res_type() calls an
# amino acid or a nucleotide is flagged 'modified' in the residue hash.
my %STANDARD = map { $_ => 1 } qw(
	ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL
	DA DC DG DT DU A C G U
);

# names that are a nucleotide in an old file and a free base ligand in a new
# one.  Which they are depends on whether the residue has a sugar; see the note
# where they are re-typed.
my %FREE_BASE = map { $_ => 1 } qw(ADE CYT GUA THY URI);

# single-atom residues that are ions rather than ligands.  The atom-count rule
# where they are typed catches most of them; these are the ones whose residue
# name and element symbol disagree.
my %ION = map { $_ => 1 } qw(
	ZN MG CA MN FE FE2 CU CU1 NA K CL BR IOD CD CO NI HG PT AU AG CS RB SR BA
	LI AL GA IN PB SE4 SO4 PO4 NO3 CO3 NH4 F
);

# The numeric fields of a REMARK are free text, and some of them are not
# numbers: 5m04 writes its pH as "5.4.-5.8", a range with a stray dot in it,
# and [\d.]+ happily matches "5.4." -- which Perl will then refuse to add.
# This matches a number and stops, so a malformed field gives up the part of
# itself that is one, and a field of NULL gives up nothing.
my $NUM = qr/[0-9]*\.?[0-9]+/;

# Public entry points

# structure_info($file, %opt) -- read a structure file into a hash of hashes.
sub structure_info {
	my ($file, %opt) = @_;
	die 'structure_info: no file name given' unless defined $file && length $file;
	die "structure_info: '$file' does not exist"  unless -e $file;
	die "structure_info: '$file' is a directory"  if -d $file;
	my $o   = _options(\%opt, 'structure_info');
	my $fmt = defined $o->{format} ? _alias($o->{format}) : _detect_format($file);
	my $reader = $READER{$fmt}
		or die "structure_info: cannot read '$file': "
		       . (exists $NOT_YET{$fmt}
		          ? "$NOT_YET{$fmt} is not implemented yet; formats read today: " . join(', ', sort keys %READER)
		          : "unrecognized format '$fmt'; formats read today: " . join(', ', sort keys %READER));
	return $reader->($file, $o);
}

# pdb_info($file, %opt) -- structure_info() with the format settled in advance.
sub pdb_info {
	my ($file, %opt) = @_;
	return structure_info($file, %opt, format => 'pdb');
}

# cif_info($file, %opt) -- the same for mmCIF/PDBx.  Both of these exist for a
# caller who already knows what they have; structure_info() works it out and
# returns the same thing either way, so neither is the usual way in.
sub cif_info {
	my ($file, %opt) = @_;
	return structure_info($file, %opt, format => 'mmcif');
}

# structure_info_string($text, %opt) -- the same, from a string already in hand.
sub structure_info_string {
	my ($text, %opt) = @_;
	die 'structure_info_string: text is undefined' unless defined $text;
	my $o = _options(\%opt, 'structure_info_string');
	my $fmt = defined $o->{format} ? _alias($o->{format}) : _sniff_format($text);
	# a string has no name to go on, and the caller has already said this is a
	# structure, so text that looks like nothing in particular is read as PDB.
	# Text that looks like something else still gets a straight answer.
	$fmt = 'pdb' if $fmt eq 'unknown';
	die "structure_info_string: cannot read this text: "
	    . (exists $NOT_YET{$fmt}
	       ? "$NOT_YET{$fmt} is not implemented yet"
	       : "no reader for format '$fmt'")
		unless $READER{$fmt};
	my $parse = $XS{$fmt}{string};
	my $p = $parse->($text, _xs_options($o));
	return _build_structure(_retry_model($p, $o, $parse, $text), $o, undef);
}

# formats() -- the formats that can be read, in list context; in scalar
# context a hashref of every format known, including the unwritten ones.
sub formats {
	return wantarray
		? (sort keys %READER)
		: { (map { $_ => 'supported' } keys %READER), (map { $_ => "not implemented: $NOT_YET{$_}" } keys %NOT_YET) };
}

# Views over a parsed structure
#
# These build what they return.  Nothing in the structure points back up at
# its parent -- a residue does not hold its chain, an atom does not hold its
# residue -- because a hash of hashes with parent links is a cycle, and a
# cycle is a leak that no one notices until the tenth thousand file.

# structure_atoms($info, $chain?) -- every atom as a flat array of hashes,
# each one carrying the chain/residue it came from, in file order.
sub structure_atoms {
	my ($info, $chain) = @_;
	_check_info($info, 'structure_atoms');
	my @out;
	for my $cid (defined $chain ? ($chain) : @{ $info->{chain_order} }) {
		my $c = $info->{chains}{$cid} or die "structure_atoms: no chain '$cid'";
		for my $rk (@{ $c->{residue_order} }) {
			my $r = $c->{residues}{$rk};
			for my $an (@{ $r->{atom_order} }) {
				push @out, {
					%{ $r->{atoms}{$an} },
					chain   => $cid,
					resname => $r->{resname},
					resseq  => $r->{number},
					icode   => $r->{icode},
					reskey  => $rk,
				};
			}
		}
	}
	return \@out;
}

# structure_residues($info, $chain?) -- every residue, in order.
sub structure_residues {
	my ($info, $chain) = @_;
	_check_info($info, 'structure_residues');
	my @out;
	for my $cid (defined $chain ? ($chain) : @{ $info->{chain_order} }) {
		my $c = $info->{chains}{$cid} or die "structure_residues: no chain '$cid'";
		push @out, map { $c->{residues}{$_} } @{ $c->{residue_order} };
	}
	return \@out;
}

# structure_ligands($info) -- the non-water heterogens, keyed NAME_CHAIN_NUM,
# which is what a binding-site table wants as its row label.
sub structure_ligands {
	my ($info) = @_;
	_check_info($info, 'structure_ligands');
	my %lig;
	for my $cid (@{ $info->{chain_order} }) {
		my $c = $info->{chains}{$cid};
		for my $rk (@{ $c->{residue_order} }) {
			my $r = $c->{residues}{$rk};
			next unless $r->{type} eq 'ligand' || $r->{type} eq 'ion';
			$lig{"$r->{resname}_${cid}_$rk"} = $r;
		}
	}
	return \%lig;
}

# structure_sequences($info_or_file, %opt) -- chain id => single-letter
# sequence, observed.  A file name in place of the parsed structure is read
# first, because the sequences are the one view a caller often wants on their
# own, and parsing a file to throw the rest of it away is two lines that read
# as one thought.  %opt is then what structure_info() takes; with a structure
# already in hand there is nothing left for options to affect, so passing them
# there is a mistake and is said to be one rather than quietly ignored.
sub structure_sequences {
	my ($info, %opt) = @_;
	if (ref $info) {
		die 'structure_sequences: options apply to reading a file, not to a structure already parsed: '
		    . join(', ', sort keys %opt) if %opt;
	}
	else {
		die 'structure_sequences: expected a file name or the hash reference from structure_info()'
			unless defined $info && length $info;
		$info = structure_info($info, %opt);
	}
	_check_info($info, 'structure_sequences');
	return { map { $_ => $info->{chains}{$_}{sequence} }
	         grep { length $info->{chains}{$_}{sequence} } @{ $info->{chain_order} } };
}

# chain_sequence($info, $chain, $which?) -- one chain's sequence.  $which is
# 'observed' (default: the residues that have coordinates) or 'seqres' (what
# the crystallographer put in, disordered tails and all).
sub chain_sequence {
	my ($info, $chain, $which) = @_;
	_check_info($info, 'chain_sequence');
	die 'chain_sequence: no chain given' unless defined $chain;
	my $c = $info->{chains}{$chain} or die "chain_sequence: no chain '$chain'";
	$which = 'observed' unless defined $which;
	die "chain_sequence: which must be 'observed' or 'seqres', not '$which'"
		unless $which eq 'observed' || $which eq 'seqres';
	return $which eq 'seqres' ? $c->{seqres} : $c->{sequence};
}

# structure_summary($info) -- a paragraph a human can read, for STDERR or a log.
sub structure_summary {
	my ($info) = @_;
	_check_info($info, 'structure_summary');
	my @l;
	push @l, sprintf('%s  %s', $info->{id} || '????', $info->{title} || '(no title)');
	push @l, sprintf('  file        %s', $info->{file}) if defined $info->{file};
	push @l, sprintf('  method      %s', join(', ', @{ $info->{experiment} })) if @{ $info->{experiment} || [] };
	push @l, sprintf('  resolution  %s A', $info->{resolution}) if defined $info->{resolution};
	push @l, sprintf('  R / R-free  %s / %s',
		defined $info->{r_work} ? $info->{r_work} : '-',
		defined $info->{r_free} ? $info->{r_free} : '-')
		if defined $info->{r_work} || defined $info->{r_free};
	push @l, sprintf('  models      %d%s', $info->{n_models},
		$info->{n_models} > 1 ? " (chains built from model $info->{model})" : '');
	push @l, sprintf('  atoms       %d (%d hetatm, %d water)',
		$info->{stats}{n_atoms}, $info->{stats}{n_hetatm}, $info->{stats}{n_water_atoms});
	for my $cid (@{ $info->{chain_order} }) {
		my $c = $info->{chains}{$cid};
		push @l, sprintf('  chain %-2s    %-11s %4d residues, %5d atoms%s',
			$cid, $c->{type}, $c->{n_residues}, $c->{n_atoms},
			$c->{n_gaps} ? ", $c->{n_gaps} gap" . ($c->{n_gaps} > 1 ? 's' : '') : '');
		push @l, sprintf('              %s', $c->{sequence}) if length $c->{sequence};
		push @l, sprintf('              %s', $c->{molecule}) if defined $c->{molecule};
	}
	my $lig = structure_ligands($info);
	push @l, sprintf('  ligands     %s', join(', ', sort keys %$lig)) if %$lig;
	return join("\n", @l) . "\n";
}

# Physical properties
#
# Everything a structure is asked about after it has been read: how much of it
# the solvent can touch, how big and how heavy it is, which of its aromatic
# rings are stacked, and the two numbers a sequence alone answers -- its mean
# hydropathy and how much of it is aromatic.
#
# The first four are one XS call, because all four have to touch every atom and
# the walk that flattens a hash of chains of residues of atoms into coordinate
# arrays is the expensive part; asking for all of them costs one walk.  The
# last two are here rather than there because they read a sequence string the
# chain already carries, which is a few hundred characters and no per-atom work
# at all.
#
# Every number is checked against the implementation it came from, not against
# this one.  t/features.t compares with mdtraj 1.11 and Biopython 1.87 -- live
# where mdtraj is installed, and against t/data/features.txt where it is not --
# and t/properties.t covers the parts that are this module's own: the options,
# the call forms, and what each function leaves behind in $info.

# Kyte, J; Doolittle, R F (1982) "A simple method for displaying the hydropathic
# character of a protein", J Mol Biol 157(1):105-132, the hydropathy index of
# Table 2.  Spelled the way Bio::SeqUtils' Python cousin spells it --
# Bio.SeqUtils.ProtParamData.kd of Biopython 1.87 -- which is what
# structure_features's mean is compared against.
#
# The twenty are all of it: there is no value for B, J, O, U, X or Z, and a
# residue whose single-letter code is one of those is left out of both the sum
# and the count rather than scored as something it is not.
my %KD = (
	A =>  1.8, R => -4.5, N => -3.5, D => -3.5, C =>  2.5,
	Q => -3.5, E => -3.5, G => -0.4, H => -3.2, I =>  4.5,
	L =>  3.8, K => -3.9, M =>  1.9, F =>  2.8, P => -1.6,
	S => -0.8, T => -0.7, W => -0.9, Y => -1.3, V =>  4.2,
);

# Phenylalanine, tryptophan and tyrosine.  Aromaticity is their relative
# frequency and nothing else -- Lobry, J R; Gautier, C (1994) Nucleic Acids Res
# 22(15):3174-3180, and Biopython's ProteinAnalysis.aromaticity(), which is
# what it is compared against.  Histidine's ring is aromatic and is not in the
# count, because it is not in Lobry's; structure_pi_stacking() does count it,
# because there the question is about the ring rather than about the index.
my %AROMATIC = map { $_ => 1 } qw(F W Y);

# The options the three feature functions take, with what they mean and what
# they default to.  This table is the one place the defaults are written down:
# the XS carries the same numbers as a fallback, and never sees it, because
# every call below passes the whole merged hash.
#
# The angles are degrees and the distances angstrom, which is what the rest of
# the module is in.  mdtraj's are degrees and nanometres; see the note on
# face_distance in structure_pi_stacking's documentation for the one place that
# difference is more than a conversion.
my %FEATURE_DEFAULT = (
	sasa        => 1,     # compute the solvent-accessible surface
	pi_stacking => 1,     # look for stacked aromatic rings
	disulfides  => 1,     # look for SG-SG pairs close enough to be bonded
	interface   => 1,     # also run the surface on each chain alone, for the buried area
	shape       => 1,     # the gyration tensor and the descriptors built from it
	dihedrals   => 1,     # phi, psi, omega and chi1..chi5, onto each residue
	contacts    => 1,     # which residues touch which
	exposure    => 1,     # half-sphere exposure, onto each amino acid residue
	hbonds      => 1,     # backbone hydrogen bonds, by Kabsch and Sander's energy
	secondary   => 1,     # DSSP secondary structure, onto each residue
	# write per-atom, per-residue and per-chain results into $info.  The torsion
	# angles, the half-sphere exposure and the secondary structure have nowhere
	# else to go, so store => 0 does not compute them rather than computing them
	# and dropping them; structure_info() always stores.
	store       => 1,
	probe       => 1.4,   # solvent probe radius, angstrom: mdtraj's 0.14 nm
	points      => 960,   # sphere points per atom; more is more accurate and slower
	face_distance   => 5.5, # face-to-face: centroid separation, angstrom
	face_plane_min  => 0,   # ... angle between the two ring planes, degrees
	face_plane_max  => 35,
	face_normal_min => 0,   # ... a ring's normal against the centroid-to-centroid line
	face_normal_max => 33,
	edge_distance   => 6.5, # edge-to-face: centroid separation, angstrom
	edge_plane_min  => 50,
	edge_plane_max  => 90,
	edge_normal_min => 0,
	edge_normal_max => 30,
	edge_radius     => 1.5, # ... centroid to where the two planes meet, angstrom
	# mdtraj's Topology.create_disulfide_bonds() uses 0.3 nm, and the cutoff is
	# the whole of the rule, so it is an option rather than a constant
	disulfide_distance => 3.0, # the largest SG-SG separation that is a bond, angstrom
	# Bio.PDB.Polypeptide.PPBuilder's `radius': the largest C-to-N separation
	# that still means two residues are joined, and so the largest one across
	# which a phi or a psi means anything
	peptide_bond       => 1.8, # angstrom
	# the closest two heavy atoms of two residues may be and still count as a
	# contact; 4.5 A is the usual convention and the option is there because
	# there is no single right answer
	contact_distance   => 4.5, # angstrom
);

# which of those each function takes.  A geometry threshold passed to
# structure_sasa() is a misunderstanding rather than a harmless extra, and the
# same is true the other way round, so each function only knows its own.
my @SASA_OPT = qw(store probe points interface);
my @PI_OPT   = qw(
	face_distance face_plane_min face_plane_max face_normal_min face_normal_max
	edge_distance edge_plane_min edge_plane_max edge_normal_min edge_normal_max
	edge_radius
);
my @SS_OPT   = qw(store disulfide_distance);
my @TORS_OPT = qw(store peptide_bond);
my @CONT_OPT = qw(store contact_distance);
my @HB_OPT   = qw(peptide_bond);

sub _feature_options {
	my ($opt, $who, $allowed) = @_;
	my %ok = map { $_ => 1 } @$allowed;
	for my $k (sort keys %$opt) {
		die "$who: unknown option '$k'; known options are: " . join(', ', sort @$allowed)
			unless $ok{$k};
	}
	my %o = (%FEATURE_DEFAULT, %$opt);
	die "$who: probe must be a number and must not be negative, not '$o{probe}'"
		unless $o{probe} =~ /\A[0-9]*\.?[0-9]+\z/;
	# the upper bound is memory rather than accuracy: a sphere point is three
	# NVs, so ten million of them is 240 MB on a double perl and four times that
	# on a quadmath one, for an atom whose area is one number
	die "$who: points must be an integer between 1 and 10000000, not '$o{points}'"
		unless $o{points} =~ /\A[0-9]+\z/ && $o{points} > 0 && $o{points} <= 10_000_000;
	for my $k (@PI_OPT) {
		die "$who: $k must be a number, not '$o{$k}'"
			unless $o{$k} =~ /\A-?[0-9]*\.?[0-9]+\z/;
	}
	die "$who: disulfide_distance must be a number and must not be negative, "
	  . "not '$o{disulfide_distance}'"
		unless $o{disulfide_distance} =~ /\A[0-9]*\.?[0-9]+\z/;
	die "$who: peptide_bond must be a positive number, not '$o{peptide_bond}'"
		unless $o{peptide_bond} =~ /\A[0-9]*\.?[0-9]+\z/ && $o{peptide_bond} > 0;
	die "$who: contact_distance must be a positive number, not '$o{contact_distance}'"
		unless $o{contact_distance} =~ /\A[0-9]*\.?[0-9]+\z/ && $o{contact_distance} > 0;
	return \%o;
}

# The properties structure_info() computes on the way past, which is all of
# them.  Kept apart from the wrappers below because it is also what decides what
# $info->{features} holds, and so what the wrappers can hand back without doing
# the work a second time.
sub _all_features {
	my ($info, $who) = @_;
	my $o = _feature_options({}, $who,
		[ qw(sasa pi_stacking disulfides shape dihedrals contacts exposure hbonds
		     secondary),
		  @SASA_OPT, @PI_OPT, @SS_OPT, @TORS_OPT, @CONT_OPT, @HB_OPT ]);
	my $f = _features($info, $o, $who);
	_sequence_features($info, $f, 1);
	return $f;
}

# structure_features($info, %opt) -- the physical properties of a structure.
#
# With no options this is a lookup: structure_info() computed them on the way
# past and left them in $info->{features}, and computing them again would give
# the same answer for the price of the walk.  Name any option and it is
# computed again with that option in force.
sub structure_features {
	my ($info, %opt) = @_;
	_check_info($info, 'structure_features');
	return $info->{features} if !%opt && $info->{features};
	my $o = _feature_options(\%opt, 'structure_features',
		[ qw(sasa pi_stacking disulfides shape dihedrals contacts exposure hbonds
		     secondary),
		  @SASA_OPT, @PI_OPT, @SS_OPT, @TORS_OPT, @CONT_OPT, @HB_OPT ]);
	my $f = _features($info, $o, 'structure_features');
	_sequence_features($info, $f, $o->{store});
	return $f;
}

# structure_sasa($info, %opt) -- the solvent-accessible surface, and nothing else.
sub structure_sasa {
	my ($info, %opt) = @_;
	_check_info($info, 'structure_sasa');
	return $info->{features}{sasa} if !%opt && $info->{features};
	my $o = _feature_options(\%opt, 'structure_sasa', \@SASA_OPT);
	$o->{sasa}        = 1;
	$o->{pi_stacking} = 0;
	$o->{disulfides}  = 0;
	$o->{shape}       = 0;
	$o->{dihedrals}   = 0;
	$o->{contacts}    = 0;
	$o->{exposure}    = 0;
	$o->{hbonds}      = 0;
	$o->{secondary}   = 0;
	return _features($info, $o, 'structure_sasa')->{sasa};
}

# structure_pi_stacking($info, %opt) -- the stacked pairs of aromatic rings.
sub structure_pi_stacking {
	my ($info, %opt) = @_;
	_check_info($info, 'structure_pi_stacking');
	return $info->{features}{pi_stacking} if !%opt && $info->{features};
	my $o = _feature_options(\%opt, 'structure_pi_stacking', \@PI_OPT);
	$o->{sasa}        = 0;
	$o->{pi_stacking} = 1;
	$o->{disulfides}  = 0;
	$o->{shape}       = 0;
	$o->{dihedrals}   = 0;
	$o->{contacts}    = 0;
	$o->{exposure}    = 0;
	$o->{hbonds}      = 0;
	$o->{secondary}   = 0;
	$o->{store}       = 0;
	return _features($info, $o, 'structure_pi_stacking')->{pi_stacking};
}

# structure_contacts($info, %opt) -- which residues touch which.
sub structure_contacts {
	my ($info, %opt) = @_;
	_check_info($info, 'structure_contacts');
	return $info->{features}{contacts} if !%opt && $info->{features};
	my $o = _feature_options(\%opt, 'structure_contacts', \@CONT_OPT);
	$o->{sasa}        = 0;
	$o->{pi_stacking} = 0;
	$o->{disulfides}  = 0;
	$o->{shape}       = 0;
	$o->{dihedrals}   = 0;
	$o->{contacts}    = 1;
	$o->{exposure}    = 0;
	$o->{hbonds}      = 0;
	$o->{secondary}   = 0;
	return _features($info, $o, 'structure_contacts')->{contacts};
}

# structure_hbonds($info, %opt) -- the backbone hydrogen bonds.
sub structure_hbonds {
	my ($info, %opt) = @_;
	_check_info($info, 'structure_hbonds');
	return $info->{features}{hbonds} if !%opt && $info->{features};
	my $o = _feature_options(\%opt, 'structure_hbonds', \@HB_OPT);
	$o->{$_} = 0 for qw(sasa pi_stacking disulfides shape dihedrals contacts
	                    exposure secondary);
	$o->{hbonds} = 1;
	return _features($info, $o, 'structure_hbonds')->{hbonds};
}

# structure_disulfides($info, %opt) -- the SG-SG pairs close enough to be bonded.
sub structure_disulfides {
	my ($info, %opt) = @_;
	_check_info($info, 'structure_disulfides');
	return $info->{features}{disulfides} if !%opt && $info->{features};
	my $o = _feature_options(\%opt, 'structure_disulfides', \@SS_OPT);
	$o->{sasa}        = 0;
	$o->{pi_stacking} = 0;
	$o->{disulfides}  = 1;
	$o->{shape}       = 0;
	$o->{dihedrals}   = 0;
	$o->{contacts}    = 0;
	$o->{exposure}    = 0;
	$o->{hbonds}      = 0;
	$o->{secondary}   = 0;
	return _features($info, $o, 'structure_disulfides')->{disulfides};
}

# The two sequence-level numbers, per protein chain and over the structure.
#
# Over the observed sequence -- the residues that are actually in the file --
# rather than over SEQRES, because these are properties of the model in front of
# the caller and a residue nobody could see has no coordinates to be a property
# of.  chain_sequence($info, $chain, 'seqres') is there for anyone who wants the
# other one.
sub _sequence_features {
	my ($info, $f, $store) = @_;
	my ($sum, $scored, $arom, $len) = (0, 0, 0, 0);
	for my $cid (@{ $info->{chain_order} }) {
		my $c = $info->{chains}{$cid};
		next unless ($c->{type} || '') eq 'protein';
		my $seq = $c->{sequence};
		next unless defined $seq && length $seq;
		my ($csum, $cn, $ca) = (0, 0, 0);
		for my $aa (split //, uc $seq) {
			$ca++ if $AROMATIC{$aa};
			next unless exists $KD{$aa};
			$csum += $KD{$aa};
			$cn++;
		}
		$sum    += $csum;
		$scored += $cn;
		$arom   += $ca;
		$len    += length $seq;
		next unless $store;
		$c->{hydropathy}        = $csum / $cn if $cn;
		$c->{aromatic_fraction} = $ca / length $seq;
	}
	$f->{hydropathy}        = $sum / $scored if $scored;
	$f->{aromatic_fraction} = $arom / $len   if $len;
	$f->{n_aromatic}        = $arom;
	$f->{sequence_length}   = $len;
	return $f;
}

# Options and format detection

sub _options {
	my ($opt, $who) = @_;
	for my $k (sort keys %$opt) {
		die "$who: unknown option '$k'; known options are: " . join(', ', sort keys %DEFAULT)
			unless exists $DEFAULT{$k};
	}
	my %o = (%DEFAULT, %$opt);
	die "$who: altloc must be 'first' or 'highest', not '$o{altloc}'"
		unless $o{altloc} eq 'first' || $o{altloc} eq 'highest';
	if (defined $o{chains}) {
		die "$who: chains must be an array reference"
			unless (reftype($o{chains}) || '') eq 'ARRAY';
		die "$who: chains is empty" unless @{ $o{chains} };
	}
	if (defined $o{model} && $o{model} ne 'all') {
		die "$who: model must be a positive integer or 'all', not '$o{model}'"
			unless $o{model} =~ /\A[0-9]+\z/;
	}
	return \%o;
}

# the options the XS parser understands, which are the ones that let it throw
# a line away before it has built a single SV for it
sub _xs_options {
	my ($o) = @_;
	return {
		model     => ($o->{model} eq 'all' ? -1 : $o->{model}),
		# when atoms are wanted, the parse builds the atom hashes itself; see
		# the note in the XS about not building every atom twice
		atom_hashes => $o->{atoms} ? 1 : 0,
		hydrogens => $o->{hydrogens},
		waters    => $o->{waters},
		hetatm    => $o->{hetatm},
		meta      => $o->{meta},
		anisou    => $o->{anisou},
		(defined $o->{chains} ? (chains => { map { $_ => 1 } @{ $o->{chains} } }) : ()),
	};
}

# the format a caller named, under the name it is filed under here
sub _alias {
	my ($fmt) = @_;
	$fmt = lc $fmt;
	return exists $ALIAS{$fmt} ? $ALIAS{$fmt} : $fmt;
}

sub _detect_format {
	my ($file) = @_;
	my $name = $file;
	$name =~ s/\.(gz|bz2|z)\z//i;
	return 'pdb'   if $name =~ /\.(pdb|ent|pdb\d+)\z/i;
	return 'mmcif' if $name =~ /\.(cif|mmcif|pdbx)\z/i;
	return 'mol2'  if $name =~ /\.mol2\z/i;
	return 'sdf'   if $name =~ /\.(sdf|mol)\z/i;
	return _sniff_format(_head($file));
}

# when the name says nothing, the first few records do
sub _sniff_format {
	my ($text) = @_;
	return 'mmcif' if $text =~ /^(?:data_|loop_|_atom_site\.)/m;
	return 'pdb'   if $text =~ /^(?:HEADER|ATOM  |HETATM|MODEL |REMARK|CRYST1|SEQRES|EXPDTA|TITLE )/m;
	return 'mol2'  if $text =~ /^\@<TRIPOS>/m;
	return 'sdf'   if $text =~ /^\s*M  END\s*$/m;
	return 'unknown';
}

sub _head {
	my ($file) = @_;
	my $text = _slurp_maybe_gzipped($file, 8192);
	return $text;
}

# .gz is worth handling here: a directory of a few thousand structures is
# usually kept compressed, and gunzipping into a temporary file first is both
# slower and something the caller then has to clean up.
sub _slurp_maybe_gzipped {
	my ($file, $limit) = @_;
	if ($file =~ /\.gz\z/i) {
		eval { require IO::Uncompress::Gunzip; 1 }
			or die "Chem::Structure::Parser: '$file' is gzipped but IO::Uncompress::Gunzip is not installed: $@";
		my $z = IO::Uncompress::Gunzip->new($file)
			or die "Chem::Structure::Parser: cannot gunzip '$file': "
			       . do { no warnings 'once'; $IO::Uncompress::Gunzip::GunzipError };
		my ($text, $buf) = ('', '');
		while ($z->read($buf, 65536) > 0) {
			$text .= $buf;
			last if defined $limit && length($text) >= $limit;
		}
		$z->close;
		return $text;
	}
	open my $fh, '<:raw', $file;
	my $text = '';
	if (defined $limit) {
		read $fh, $text, $limit;
	} else {
		local $/;
		$text = <$fh>;
		$text = '' unless defined $text;
	}
	close $fh;
	return $text;
}

#
# The readers
#
# There is one of these per format and they differ only in which XS parser
# they call, because the parsers agree about what they hand back.  Keeping
# them as two named subs rather than one closure is so that %READER reads as a
# list of formats and the die message above can name them.
#

sub _read_pdb { return _read($XS{pdb},   @_) }
sub _read_cif { return _read($XS{mmcif}, @_) }

sub _read {
	my ($xs, $file, $o) = @_;
	my $p;
	if ($file =~ /\.gz\z/i) {
		my $text = _slurp_maybe_gzipped($file, undef);
		$p = $xs->{string}->($text, _xs_options($o));
		$p = _retry_model($p, $o, $xs->{string}, $text);
	} else {
		$p = $xs->{file}->($file, _xs_options($o));
		$p = _retry_model($p, $o, $xs->{file}, $file);
	}
	return _build_structure($p, $o, $file);
}

# An NMR ensemble whose models are numbered from 0, or a file whose only model
# is MODEL 7, would otherwise come back empty for the default model => 1.  The
# parse says which model numbers it saw, so ask again for the first real one
# rather than handing back a structure with no atoms in it.
sub _retry_model {
	my ($p, $o, $parse, $src) = @_;
	return $p if $p->{n_atoms} || $o->{model} eq 'all';
	my $nums = $p->{model_numbers};
	return $p unless @$nums && !grep { $_ == $o->{model} } @$nums;
	my $x = _xs_options($o);
	$x->{model} = $nums->[0];
	my $q = $parse->($src, $x);
	$q->{requested_model} = $nums->[0];
	return $q;
}

sub _build_structure {
	my ($p, $o, $file) = @_;
	my $fmt = $p->{format} || 'pdb';
	my $info = {
		file     => $file,
		format   => $fmt,
		model    => (defined $p->{requested_model} ? $p->{requested_model}
		             : $o->{model} eq 'all' ? 'all' : $o->{model}),
		n_models => $p->{n_models},
		# The counts and extremes come straight from the parse.  They have to
		# touch every atom, and the parse is already reading every atom, so
		# doing them there costs nothing and doing them again here would cost
		# more than the parse itself.
		stats    => {
			n_atoms          => $p->{n_atoms},
			n_hetatm         => 0,
			# Every ATOM/HETATM record the file has, which is not what n_atoms
			# counts: that one is what came back, after the model selection and
			# the hydrogens, waters, hetatm and chains options have had their
			# say, and this one is what was there to be filtered.  The two are
			# equal for a single-model file read with the defaults and are not
			# for anything else -- an NMR ensemble read at its default model => 1
			# returns a twentieth of its atoms and is not a twentieth of a file.
			# total_atoms == n_atoms + n_skipped, always, in both formats.
			total_atoms      => $p->{n_atom_records} + $p->{n_hetatm_records},
			n_hydrogens      => $p->{n_hydrogens},
			n_water_atoms    => $p->{n_water_atoms},
			n_atom_records   => $p->{n_atom_records},
			n_hetatm_records => $p->{n_hetatm_records},
			n_anisou         => $p->{n_anisou},
			n_skipped        => $p->{n_skipped},
			n_lines          => $p->{n_lines},
			elements         => $p->{elements},
			bfactor          => $p->{bfactor_stats},
			bbox             => $p->{bbox},
			center           => $p->{center},
		},
	};

	# the one place the two formats are read differently, and the reason it is
	# the only one: everything below works off $info, which is the same shape
	# whichever of these filled it in
	if ($o->{meta}) {
		$fmt eq 'mmcif' ? _parse_cif_meta($info, $p) : _parse_meta($info, $p->{meta} || {});
	}

	my $by_model = _assemble($p, $o, $info);
	my @models = sort { $a <=> $b } keys %$by_model;
	my $main   = $o->{model} eq 'all' ? (@models ? $models[0] : 1) : $info->{model};
	$main = $models[0] if @models && !exists $by_model->{$main};

	my $sel = $by_model->{$main} || { chains => {}, chain_order => [] };
	$info->{chains}      = $sel->{chains};
	$info->{chain_order} = $sel->{chain_order};
	if ($o->{model} eq 'all') {
		$info->{models} = $by_model;
		$info->{model}  = $main;
	}

	_finish_chains($info);
	_chain_stats($info);
	$info->{id} = _id_from($info, $file);

	# The physical properties, computed on the way past.
	#
	# They are on by default because a structure's surface, size and contacts
	# are as much a part of what it *is* as its sequence, and a caller who has
	# to know to ask for them mostly does not.  What that costs is real and is
	# measured in notes.txt: reading a structure goes from about 214,000 atoms a
	# second to about 19,500, which is eleven times the cost.  Nearly all of it
	# is the solvent-accessible surface, at 960 sphere points per atom -- with
	# sasa => 0 the rest together come to 2.7x the read, and interface => 0
	# alone takes the whole thing from 10.9x to 7.1x.  features => 0 is the way
	# back to the old speed, and is what to reach for when reading a directory
	# for its headers or its sequences.
	#
	# atoms => 0 turns them off on its own rather than dying: it is the
	# documented way to read a file for its residues without its coordinates,
	# structure_sequences() passes it, and a fast path that started dying would
	# be a worse answer than one that quietly has nothing to compute from.
	$info->{features} = _all_features($info, 'structure_info')
		if $o->{features} && $o->{atoms};

	return $info;
}

#  coordinates
#
# The XS parse hands back one array per field plus the index of the first and
# last atom of every residue, so this walks residues, not atoms, and only
# descends into an atom loop when the caller wants atoms at all.
sub _assemble {
	my ($p, $o, $info) = @_;
	my ($rf, $rl) = @{$p}{qw(res_first res_last)};
	my ($chain, $resname, $resseq, $icode, $het, $model)
		= @{$p}{qw(chain resname resseq icode het model)};
	my ($name, $altloc, $serial, $x, $y, $z, $occ, $bf, $elem, $charge)
		= @{$p}{qw(name altloc serial x y z occupancy bfactor element charge)};
	my ($rsx, $rsy, $rsz, $rnxyz, $rsb, $rnb) = @{$p}{qw(sx sy sz n_xyz sb n_b)};
	# the parse tallied the elements of every chain of every model as it read
	# them, keyed model then chain, which is the pair that names a chain here
	my $chain_elem = $p->{chain_elements} || {};
	# the atom hashes, built by the parse when the caller wanted atoms at all
	my $atom_of = ($o->{atoms} && @{ $p->{atoms} || [] }) ? $p->{atoms} : undef;
	my $st = $info->{stats};
	my %by_model;

	for my $r (0 .. $#$rf) {
		my ($i0, $i1) = ($rf->[$r], $rl->[$r]);
		my $m   = $model->[$i0];
		my $cid = $chain->[$i0];
		my $rn  = $resname->[$i0];
		my $num = $resseq->[$i0];
		my $ic  = $icode->[$i0];
		my $key = (defined $num ? $num : '') . $ic;

		my $mm = $by_model{$m} ||= { chains => {}, chain_order => [] };
		my $c  = $mm->{chains}{$cid};
		unless ($c) {
			$c = $mm->{chains}{$cid} = {
				id            => $cid,
				residues      => {},
				residue_order => [],
				n_atoms       => 0,
				n_hetatm      => 0,
				# how many atoms of each element the chain holds, keyed by the
				# IUPAC symbol.  The counts add up to n_atoms, alternate
				# conformers and all, because both count records rather than
				# distinct atoms.
				elements      => (($chain_elem->{$m} || {})->{$cid} || {}),
			};
			push @{ $mm->{chain_order} }, $cid;
		}

		# a residue can be met twice -- altloc groups written apart, or a
		# ligand interleaved with the polymer -- so merge rather than replace
		my $res = $c->{residues}{$key};
		unless ($res) {
			my $type = res_type($rn);
			$res = $c->{residues}{$key} = {
				chain      => $cid,
				resname    => $rn,
				number     => $num,
				icode      => $ic,
				key        => $key,
				one        => res1($rn),
				type       => $type,
				hetero     => $het->[$i0],
				standard   => ($STANDARD{$rn} ? 1 : 0),
				modified   => (($type eq 'amino_acid' || $type eq 'nucleotide') && !$STANDARD{$rn}) ? 1 : 0,
				n_atoms    => 0,
				atoms      => {},
				atom_order => [],
			};
			push @{ $c->{residue_order} }, $key;
		}

		# Counts and sums come out of the parse, which had to read every
		# coordinate anyway.  They are sums rather than means so that a residue
		# met twice can be added up instead of recomputed; _finish_chains
		# divides them and takes the temporaries back out.
		my $n = $i1 - $i0 + 1;
		$res->{n_atoms}  += $n;
		$res->{_sx}      += $rsx->[$r] if defined $rsx->[$r];
		$res->{_sy}      += $rsy->[$r] if defined $rsy->[$r];
		$res->{_sz}      += $rsz->[$r] if defined $rsz->[$r];
		$res->{_nxyz}    += $rnxyz->[$r];
		$res->{_sb}      += $rsb->[$r] if defined $rsb->[$r];
		$res->{_nb}      += $rnb->[$r];
		$c->{n_atoms}    += $n;
		if ($het->[$i0]) {    # the record type is part of a residue's identity
			$c->{n_hetatm} += $n;
			$st->{n_hetatm} += $n;
		}

		# Only the handful of names that are a nucleotide in one file and a
		# free base in another need their atoms looked through; walking every
		# residue's atoms to ask a question about five of them would undo the
		# point of having the parse mark the residues in the first place.
		# Only the handful of names that mean a nucleotide in one file and a
		# free base in another need their atoms looked through.
		my $backbone = 0;
		if ($FREE_BASE{$rn}) {
			for my $i ($i0 .. $i1) {
				my $an = $atom_of ? $atom_of->[$i]{name} : $name->[$i];
				next unless $an eq 'P' || $an =~ /\AC1[*']\z/ || $an =~ /\AO5[*']\z/;
				$backbone = 1;
				last;
			}
		}

		if ($atom_of) {
			my $ra = $res->{atoms};
			my $ro = $res->{atom_order};
			for my $i ($i0 .. $i1) {
				my $atom = $atom_of->[$i];
				my $an   = $atom->{name};
				my $have = $ra->{$an};
				unless ($have) {
					$ra->{$an} = $atom;
					push @$ro, $an;
					# an atom with an altloc keeps the list even when it is the
					# only conformer, so that "was this modelled twice?" is one
					# question rather than two
					$atom->{altlocs} = [ _conformer($atom) ] if length $atom->{altloc};
					next;
				}
				# an alternate conformer: every one is kept on the atom, and
				# the altloc option decides which supplies the coordinates.
				# The one already there is put on the list first if it is not
				# on it yet, which is the case where the first record of the
				# pair had no altloc letter at all -- disordered.pdb writes
				# ARG 27's CZ once with a blank altloc and once as B, and a
				# list holding only the B is a list that has lost a conformer.
				$have->{altlocs} ||= [ _conformer($have) ];
				push @{ $have->{altlocs} }, _conformer($atom);
				if ($o->{altloc} eq 'highest'
				    && defined $atom->{occupancy} && defined $have->{occupancy}
				    && $atom->{occupancy} > $have->{occupancy}) {
					@{$have}{qw(altloc serial x y z occupancy bfactor)}
						= @{$atom}{qw(altloc serial x y z occupancy bfactor)};
				}
			}
		}

		# the type an unknown residue really is, now that its atoms are counted
		if ($res->{type} eq 'other') {
			my $el = $atom_of ? $atom_of->[$i0]{element} : $elem->[$i0];
			$res->{type} = ($ION{$rn} || ($res->{n_atoms} == 1 && uc($el || '') eq uc $rn))
			             ? 'ion' : 'ligand';
		}
		# ADE, CYT, GUA, THY and URI mean two different things depending on
		# how old the file is: in a pre-v3 entry they are the nucleotides of a
		# nucleic acid chain, and in a modern one they are free bases sitting
		# in an active site as ligands.  The sugar tells them apart -- a
		# nucleotide has a C1', a free base has nothing but the base.  Without
		# this, the guanine bound to 1czc is read as a nucleotide and turns up
		# as a G on the end of a 396-residue protein sequence.
		if ($res->{type} eq 'nucleotide' && !$backbone && $FREE_BASE{$rn}) {
			$res->{type}     = 'ligand';
			$res->{one}      = '';
			$res->{modified} = 0;
		}
	}

	return \%by_model;
}

#  per-chain sequence, type and gaps 
sub _finish_chains {
	my ($info) = @_;
	# with model => 'all' the main model's chains are one of the models, so
	# walking the models covers it; walking both would do it twice
	my @all = $info->{models}
	        ? (values %{ $info->{models} })
	        : ({ chains => $info->{chains}, chain_order => $info->{chain_order} });
	for my $set (@all) {
		for my $cid (@{ $set->{chain_order} }) {
			my $c = $set->{chains}{$cid};
			_demote_free_residues($c);
			my (@seq, @poly, %count);
			for my $rk (@{ $c->{residue_order} }) {
				my $r = $c->{residues}{$rk};
				# the sums the parse gathered, turned into the means the
				# residue actually advertises, and then taken back out
				if (my $nc = delete $r->{_nxyz}) {
					$r->{center} = [ $r->{_sx} / $nc, $r->{_sy} / $nc, $r->{_sz} / $nc ];
				}
				if (my $nb = delete $r->{_nb}) {
					$r->{b_mean} = $r->{_sb} / $nb;
				}
				delete @{$r}{qw(_sx _sy _sz _sb)};
				$count{ $r->{type} }++;
				next unless $r->{type} eq 'amino_acid' || $r->{type} eq 'nucleotide';
				push @seq, (length $r->{one} ? $r->{one} : 'X');
				push @poly, $r;
			}
			$c->{sequence}     = join '', @seq;
			$c->{n_residues}   = scalar @{ $c->{residue_order} };
			$c->{n_polymer}    = scalar @poly;
			$c->{n_water}      = $count{water}  || 0;
			$c->{n_ligand}     = ($count{ligand} || 0) + ($count{ion} || 0);
			$c->{residue_types}= \%count;
			$c->{type}         = _chain_type(\%count, \@poly);
			$c->{first}        = @poly ? $poly[0]{key}  : undef;
			$c->{last}         = @poly ? $poly[-1]{key} : undef;

			# gaps: unmodelled stretches, which is where a sequence read off
			# the coordinates quietly differs from the one in SEQRES.
			# missing_residues is the same fact one number at a time: every
			# residue number the polymer skips over, so a caller can ask "is
			# 47 modelled?" without walking the gap list.
			#
			# A jump in the numbering only means missing residues if the chain
			# has that many to be missing, and the numbering says how many it
			# has: a chain running from its first polymer residue to its last
			# covers so many numbers, and the ones it does not use are the
			# missing ones.  That total is the budget.  A jump wider than the
			# whole budget is not a gap but a change of numbering scheme -- an
			# antibody numbered by the Kabat scheme runs 27, 1027, 2027, 28,
			# where the thousands are insertions after 27, and read literally
			# that makes 1a4k a 214-residue light chain missing five thousand
			# residues.  Insertion codes share a number, so it is the distinct
			# numbers that are counted here and not the residues: without
			# that, thrombin's 149A..149E spend budget that its real gap at
			# 217..219 then has none of.  The budget is read off the
			# coordinates alone, so a chain answers the same whether it was
			# read from a PDB file or an mmCIF one and whether or not the
			# headers were parsed.
			my (@gaps, @missing, %numbered);
			my @nums = grep { defined } map { $_->{number} } @poly;
			$numbered{$_} = 1 for @nums;
			my $budget = @nums ? $nums[-1] - $nums[0] + 1 - keys %numbered : 0;
			for my $i (1 .. $#poly) {
				my ($a, $b) = @poly[ $i - 1, $i ];
				next unless defined $a->{number} && defined $b->{number};
				my $n = $b->{number} - $a->{number} - 1;
				next if $n < 1 || $n > $budget;
				push @gaps, { after => $a->{key}, before => $b->{key}, missing => $n };
				push @missing, ($a->{number} + 1) .. ($b->{number} - 1);
			}
			# a residue numbered out of line with the rest can leave the list
			# out of order, and ascending is the whole use of it
			@missing = sort { $a <=> $b } @missing
				if grep { $missing[$_] < $missing[ $_ - 1 ] } 1 .. $#missing;
			$c->{gaps}             = \@gaps;
			$c->{n_gaps}           = scalar @gaps;
			$c->{missing_residues} = \@missing;
		}
	}
	return $info;
}

# one conformer's worth of an atom, for the altlocs list
sub _conformer {
	my ($a) = @_;
	return { map { $_ => $a->{$_} } qw(altloc serial x y z occupancy bfactor) };
}

# A HETATM residue with an amino acid's name is one of two very different
# things.  Numbered among the polymer it is a modified residue -- the MSE that
# replaced a methionine -- and it belongs in the sequence.  Numbered out with
# the ligands it is a free amino acid sitting in a binding site, and it does
# not: 3lms has a glycine at A501, two hundred residues past the end of a
# chain whose SEQRES is 309 long, and counting it makes a 310-residue protein
# out of a 309-residue one.
#
# The numbering is what separates them.  Heterogens are numbered in their own
# range, after the polymer, by long convention; a modified residue takes the
# number of the residue it replaced.  A chain written entirely as HETATM -- a
# synthetic peptide ligand, say -- has no polymer range to compare against, so
# nothing is demoted and the whole thing reads as the peptide it is.
sub _demote_free_residues {
	my ($c) = @_;
	my ($lo, $hi);
	for my $rk (@{ $c->{residue_order} }) {
		my $r = $c->{residues}{$rk};
		next if $r->{hetero};
		next unless $r->{type} eq 'amino_acid' || $r->{type} eq 'nucleotide';
		next unless defined $r->{number};
		$lo = $r->{number} if !defined $lo || $r->{number} < $lo;
		$hi = $r->{number} if !defined $hi || $r->{number} > $hi;
	}
	return $c unless defined $lo;
	for my $rk (@{ $c->{residue_order} }) {
		my $r = $c->{residues}{$rk};
		next unless $r->{hetero};
		next unless $r->{type} eq 'amino_acid' || $r->{type} eq 'nucleotide';
		next unless defined $r->{number};
		# one either side, so that a modified residue capping a terminus is
		# still part of the chain
		next if $r->{number} >= $lo - 1 && $r->{number} <= $hi + 1;
		$r->{type}     = 'ligand';
		$r->{one}      = '';
		$r->{modified} = 0;
		$r->{free}     = 1;    # a free amino acid, not part of the polymer
	}
	return $c;
}

sub _chain_type {
	my ($count, $poly) = @_;
	my $aa  = $count->{amino_acid} || 0;
	my $nuc = $count->{nucleotide} || 0;
	if ($aa || $nuc) {
		return 'protein' if $aa >= $nuc;
		my $deoxy = grep { $_->{resname} =~ /\AD[ACGTUI]\z/ } @$poly;
		return $deoxy * 2 >= $nuc ? 'dna' : 'rna';
	}
	return 'water'  if ($count->{water}  || 0) && !($count->{ligand} || 0) && !($count->{ion} || 0);
	return 'hetero' if ($count->{ligand} || 0) || ($count->{ion} || 0) || ($count->{water} || 0);
	return 'unknown';
}

# SEQRES, COMPND and SOURCE all describe chains; fold them in once the chains
# exist, so that everything about a chain is in one place
sub _chain_stats {
	my ($info) = @_;
	# The free-text COMPND of an old file names no chains, so there was no
	# chain to file it under when the header was read and there is one now: the
	# entry is the one molecule and every chain in it is that molecule.  This is
	# the only place a chain is added to entity_of_chain, because it is the only
	# entity that could not say for itself which chains it means.
	if (($info->{compound}{1} || {})->{free_text} && !%{ $info->{entity_of_chain} }) {
		my $s = $info->{source}{1} || {};
		$info->{entity_of_chain}{$_} = {
			mol_id   => 1,
			molecule => $info->{compound}{1}{molecule},
			organism => $s->{organism_scientific},
		} for @{ $info->{chain_order} };
	}
	for my $cid (@{ $info->{chain_order} }) {
		my $c = $info->{chains}{$cid};
		if (my $s = $info->{seqres}{$cid}) {
			$c->{seqres}        = $s->{sequence};
			$c->{seqres_length} = $s->{length};
			$c->{n_missing}     = $s->{length} - $c->{n_polymer} if defined $s->{length};
		}
		if (my $e = $info->{entity_of_chain}{$cid}) {
			$c->{mol_id}   = $e->{mol_id};
			$c->{molecule} = $e->{molecule} if defined $e->{molecule};
			$c->{organism} = $e->{organism} if defined $e->{organism};
			$c->{fragment} = $e->{fragment} if defined $e->{fragment};
			$c->{ec}       = $e->{ec}       if defined $e->{ec};
		}
		$c->{dbref} = $info->{dbref}{$cid} if $info->{dbref}{$cid};
	}
	return $info;
}

sub _id_from {
	my ($info, $file) = @_;
	return $info->{header}{id_code} if length($info->{header}{id_code} || '');
	return undef unless defined $file;
	my ($base) = $file =~ m{([^/\\]+)\z};
	$base =~ s/\.(gz|bz2|z)\z//i;
	$base =~ s/\.(pdb|ent|cif|mmcif)\z//i;
	$base =~ s/\.ent\z//i;
	$base =~ s/\Apdb//i;
	return uc $base;
}

#
# Header records
#
# Every one of these is a fixed-column record too, but there are only a few
# dozen lines of them in a file, they are irregular, and they are where a new
# quirk turns up every few hundred structures.  That is Perl's job, not C's.
#
# The keys a parsed structure always has, whatever was in the file and
# whichever format it was in.  Set before either reader runs, so that a caller
# can read $info->{resolution} without first asking whether the file was an
# mmCIF, and get undef for "the file does not say" in both.
sub _meta_defaults {
	my ($info) = @_;
	$info->{$_} = undef for qw(title resolution r_work r_free);
	$info->{$_} = []    for qw(keywords experiment authors);
	$info->{$_} = {}    for qw(header compound source seqres het hetnam formul
	                           remarks dbref entity_of_chain cryst1 journal
	                           modres);
	$info->{$_} = []    for qw(helix sheet ssbond link cispep revdat site conect);
	return $info;
}

sub _parse_meta {
	my ($info, $meta) = @_;
	_meta_defaults($info);

	if (my $h = $meta->{HEADER}) {
		my $l = $h->[0];
		$info->{header} = {
			classification => _c($l, 10, 40),
			deposit_date   => _c($l, 50, 9),
			id_code        => _c($l, 62, 4),
		};
	}
	# The entry id, for _untail(): the text records of an old file end in it, and
	# it is the only thing that tells the stationery from the text.
	my $eid = $info->{header}{id_code};

	# a record that is not in the file reads as undef, not as an empty string:
	# "there was no TITLE" and "the TITLE was blank" are different answers
	$info->{title}      = $meta->{TITLE} ? _joined($meta->{TITLE}, 10, $eid) : undef;
	$info->{caveat}     = _joined($meta->{CAVEAT}, 19, $eid) if $meta->{CAVEAT};
	$info->{keywords}   = [ grep { length } map { _t($_) } split /,/, _joined($meta->{KEYWDS}, 10, $eid) ];
	$info->{experiment} = [ grep { length } map { _t($_) } split /;/, _joined($meta->{EXPDTA}, 10, $eid) ];
	$info->{authors}    = [ grep { length } map { _t($_) } split /,/, _joined($meta->{AUTHOR}, 10, $eid) ];
	$info->{model_type} = _joined($meta->{MDLTYP}, 10, $eid) if $meta->{MDLTYP};
	$info->{obsolete}   = _joined($meta->{OBSLTE}, 10, $eid) if $meta->{OBSLTE};
	$info->{split}      = [ split ' ', _joined($meta->{SPLIT}, 10, $eid) ] if $meta->{SPLIT};

	$info->{compound} = _mol_records($meta->{COMPND}, 10, 'molecule', $eid);
	$info->{source}   = _mol_records($meta->{SOURCE}, 10, 'organism_scientific', $eid);
	_entities($info);

	for my $l (@{ $meta->{REVDAT} || [] }) {
		push @{ $info->{revdat} }, {
			num  => _c($l, 7, 3),
			date => _c($l, 13, 9),
			id   => _c($l, 23, 4),
			type => _c($l, 31, 1),
			what => _c($l, 39),
		};
	}

	# JRNL sub-records live in columns 13-16 and continue across lines
	for my $l (@{ $meta->{JRNL} || [] }) {
		my $sub = lc _c($l, 12, 4);
		next unless length $sub;
		my $text = _c(_untail($l, $eid), 19);
		$info->{journal}{$sub} = length($info->{journal}{$sub} || '')
			? _rejoin($info->{journal}{$sub}, $text)
			: $text;
	}
	$info->{journal}{auth} = [ grep { length } map { _t($_) } split /,/, $info->{journal}{auth} ]
		if defined $info->{journal}{auth};

	# REMARKs are kept whole, by number: there are hundreds of kinds and the
	# useful ones are pulled out below.  Anything not pulled out is still there.
	for my $l (@{ $meta->{REMARK} || [] }) {
		my $n = _c($l, 7, 3);
		next unless length $n;
		push @{ $info->{remarks}{$n} }, _c($l, 11);
	}
	for my $l (@{ $info->{remarks}{2} || [] }) {
		$info->{resolution} = $1 + 0 if $l =~ /RESOLUTION\.\s+($NUM)\s+ANGSTROM/;
	}
	# REMARK 3 says it too, as the high resolution limit of the refinement, and
	# that is the one to fall back on: a file written by a refinement program
	# rather than by the archive often has REMARK 3 and no REMARK 2 at all, and
	# 5cvz_final.pdb reads as a structure of no resolution otherwise.  It is
	# also the same number the mmCIF reader already takes from
	# _refine.ls_d_res_high, so the two formats answer alike.  Anchored, because
	# REMARK 3 also carries 'BIN RESOLUTION RANGE HIGH', which is a bin and not
	# the structure.
	if (!defined $info->{resolution}) {
		for my $l (@{ $info->{remarks}{3} || [] }) {
			next unless $l =~ /\ARESOLUTION RANGE HIGH\b[^:]*:\s*($NUM)/;
			$info->{resolution} = $1 + 0;
			last;
		}
	}
	# anchored, because REMARK 3 also carries 'BIN FREE R VALUE' and
	# 'ESTIMATED ERROR OF FREE R VALUE', which are not the R-free.  A value of
	# NULL -- what an unrefined or pre-R-free structure has -- stays undef.
	for my $l (@{ $info->{remarks}{3} || [] }) {
		$info->{r_work} = $1 + 0 if !defined $info->{r_work}
			&& $l =~ /\AR VALUE\s+\(WORKING SET\)\s*:\s*($NUM)/;
		$info->{r_free} = $1 + 0 if !defined $info->{r_free}
			&& $l =~ /\AFREE R VALUE\s*:\s*($NUM)/;
	}
	for my $l (@{ $info->{remarks}{200} || [] }) {
		$info->{temperature} = $1 + 0 if $l =~ /TEMPERATURE\s+\(KELVIN\)\s*:\s*($NUM)/;
		$info->{ph}          = $1 + 0 if $l =~ /\bPH\s*:\s*($NUM)/;
	}
	$info->{biological_assembly} = $info->{remarks}{350} if $info->{remarks}{350};

	# SEQRES -- what was in the crystal, as opposed to what was modelled
	#
	# Thirteen residues to a line, columns 20 to 70, and no further: a file old
	# enough to keep the entry id in columns 73-80 has '1GDR  81' sitting there,
	# and reading to the end of the line makes two more residues out of it.
	# pdb1gdr's 140-residue chain comes back 162 long that way, with an X every
	# thirteenth place, which is a wrong sequence rather than a missing one.
	for my $l (@{ $meta->{SEQRES} || [] }) {
		my $cid = _c($l, 11, 1);
		my $n   = _c($l, 13, 4);
		my @res = split ' ', _c($l, 19, 51);
		my $s = $info->{seqres}{$cid} ||= { chain => $cid, length => ($n =~ /\A\d+\z/ ? $n + 0 : undef), residues => [] };
		push @{ $s->{residues} }, @res;
	}
	for my $cid (keys %{ $info->{seqres} }) {
		my $s = $info->{seqres}{$cid};
		$s->{sequence} = join '', map { my $o = res1($_); length $o ? $o : 'X' } @{ $s->{residues} };
		$s->{length}   = scalar @{ $s->{residues} } unless defined $s->{length};
	}

	for my $l (@{ $meta->{DBREF} || [] }) {
		my $cid = _c($l, 12, 1);
		push @{ $info->{dbref}{$cid} }, {
			chain      => $cid,
			seq_begin  => _c($l, 14, 4),
			seq_end    => _c($l, 20, 4),
			database   => _c($l, 26, 6),
			accession  => _c($l, 33, 8),
			db_id      => _c($l, 42, 12),
			db_begin   => _c($l, 55, 5),
			db_end     => _c($l, 62, 5),
		};
	}
	for my $l (@{ $meta->{SEQADV} || [] }) {
		push @{ $info->{seqadv} }, {
			resname   => _c($l, 12, 3),
			chain     => _c($l, 16, 1),
			resseq    => _c($l, 18, 4),
			database  => _c($l, 24, 4),
			accession => _c($l, 29, 9),
			db_res    => _c($l, 39, 3),
			db_seq    => _c($l, 43, 5),
			comment   => _c($l, 49),
		};
	}
	for my $l (@{ $meta->{MODRES} || [] }) {
		my $r = _c($l, 12, 3);
		$info->{modres}{$r} ||= {
			resname  => $r,
			standard => _c($l, 24, 3),
			comment  => _c($l, 29),
		};
	}

	# heterogens: HET gives the instances, HETNAM/FORMUL name them
	for my $l (@{ $meta->{HET} || [] }) {
		my $id = _c($l, 7, 3);
		push @{ $info->{het}{$id}{instances} }, {
			chain  => _c($l, 12, 1),
			resseq => _c($l, 13, 4),
			icode  => _c($l, 17, 1),
			natoms => _c($l, 20, 5),
		};
		$info->{het}{$id}{het_id} = $id;
	}
	for my $l (@{ $meta->{HETNAM} || [] }) {
		my $id = _c($l, 11, 3);
		my $t  = _c($l, 15);
		$info->{het}{$id}{het_id} = $id;
		$info->{het}{$id}{name} = _rejoin($info->{het}{$id}{name}, $t);
	}
	for my $l (@{ $meta->{HETSYN} || [] }) {
		my $id = _c($l, 11, 3);
		$info->{het}{$id}{synonym} = _rejoin($info->{het}{$id}{synonym}, _c($l, 15));
	}
	for my $l (@{ $meta->{FORMUL} || [] }) {
		my $id = _c($l, 12, 3);
		$info->{het}{$id}{het_id}  = $id;
		$info->{het}{$id}{formula} = _rejoin($info->{het}{$id}{formula}, _c($l, 19));
		$info->{het}{$id}{water}   = 1 if _c($l, 18, 1) eq '*';
	}

	for my $l (@{ $meta->{HELIX} || [] }) {
		# the length is a number or nothing.  A file that keeps its entry id in
		# columns 73-80 puts '1GDR' where the length goes, and a caller adding
		# lengths up has no way to tell that from a length.  Empty rather than
		# undef, because empty is what every other column of a short record
		# gives and a helix should not answer two ways about the same absence.
		my $hlen = _c($l, 71, 5);
		$hlen = '' unless defined _n($hlen);
		push @{ $info->{helix} }, {
			id            => _c($l, 11, 3),
			init_resname  => _c($l, 15, 3),
			init_chain    => _c($l, 19, 1),
			init_resseq   => _c($l, 21, 4),
			end_resname   => _c($l, 27, 3),
			end_chain     => _c($l, 31, 1),
			end_resseq    => _c($l, 33, 4),
			class         => _c($l, 38, 2),
			length        => $hlen,
		};
	}
	for my $l (@{ $meta->{SHEET} || [] }) {
		push @{ $info->{sheet} }, {
			strand        => _c($l, 7, 3),
			id            => _c($l, 11, 3),
			n_strands     => _c($l, 14, 2),
			init_resname  => _c($l, 17, 3),
			init_chain    => _c($l, 21, 1),
			init_resseq   => _c($l, 22, 4),
			end_resname   => _c($l, 28, 3),
			end_chain     => _c($l, 32, 1),
			end_resseq    => _c($l, 33, 4),
			sense         => _c($l, 38, 2),
		};
	}
	for my $l (@{ $meta->{SSBOND} || [] }) {
		push @{ $info->{ssbond} }, {
			chain1  => _c($l, 15, 1),
			resseq1 => _c($l, 17, 4),
			chain2  => _c($l, 29, 1),
			resseq2 => _c($l, 31, 4),
			length  => _c($l, 73, 5),
		};
	}
	for my $l (@{ $meta->{LINK} || [] }) {
		push @{ $info->{link} }, {
			name1    => _c($l, 12, 4), resname1 => _c($l, 17, 3),
			chain1   => _c($l, 21, 1), resseq1  => _c($l, 22, 4),
			name2    => _c($l, 42, 4), resname2 => _c($l, 47, 3),
			chain2   => _c($l, 51, 1), resseq2  => _c($l, 52, 4),
			length   => _c($l, 73, 5),
		};
	}
	for my $l (@{ $meta->{CISPEP} || [] }) {
		push @{ $info->{cispep} }, {
			resname1 => _c($l, 11, 3), chain1 => _c($l, 15, 1),
			resseq1  => _c($l, 17, 4),
			resname2 => _c($l, 25, 3), chain2 => _c($l, 29, 1),
			resseq2  => _c($l, 31, 4),
			angle    => _c($l, 53, 6),
		};
	}
	if (my $c = $meta->{CRYST1}) {
		my $l = $c->[0];
		$info->{cryst1} = {
			a      => _n(_c($l, 6, 9)),  b     => _n(_c($l, 15, 9)),
			c      => _n(_c($l, 24, 9)),  alpha => _n(_c($l, 33, 7)),
			beta   => _n(_c($l, 40, 7)),  gamma => _n(_c($l, 47, 7)),
			sgroup => _c($l, 55, 11),
			z      => _c($l, 66, 4),
		};
	}
	for my $l (@{ $meta->{CONECT} || [] }) {
		my @s = grep { length } map { _t($_) }
		        map { _c($l, $_, 5) } (6, 11, 16, 21, 26);
		push @{ $info->{conect} }, \@s if @s > 1;
	}
	if (my $n = $meta->{NUMMDL}) {
		my $v = _c($n->[0], 10, 4);
		$info->{n_models_declared} = $v + 0 if $v =~ /\A\d+\z/;
	}
	$info->{records} = { map { $_ => scalar @{ $meta->{$_} } } keys %$meta };
	return $info;
}

# COMPND and SOURCE are "TOKEN: value;" lists broken into MOL_ID groups
sub _mol_records {
	my ($lines, $from, $free_key, $entry_id) = @_;
	return {} unless $lines;
	my $text = _joined($lines, $from, $entry_id);
	my %mol;
	my $id = 1;
	for my $piece (split /;/, $text) {
		next unless $piece =~ /\S/;
		my ($k, $v) = $piece =~ /\A\s*([A-Z0-9_ ]+?)\s*:\s*(.*)\z/;
		next unless defined $k;
		$k = lc $k;
		$k =~ s/\s+/_/g;
		$v = _t($v);
		if ($k eq 'mol_id') {
			$id = $v;
			$mol{$id}{mol_id} = $v;
			next;
		}
		$mol{$id}{mol_id} = $id unless exists $mol{$id};
		if ($k eq 'chain') {
			$mol{$id}{chain} = [ grep { length } map { _t($_) } split /,/, $v ];
		} else {
			$mol{$id}{$k} = exists $mol{$id}{$k} ? "$mol{$id}{$k} $v" : $v;
		}
	}
	# A file older than the MOL_ID convention writes the record as free text --
	# 'COMPND    GAMMA DELTA RESOLVASE', 'SOURCE    (ESCHERICHIA COLI)' -- and a
	# reader that knows only about 'MOLECULE:' throws away the one thing the
	# record says.  There is no chain list in that form because there was
	# nothing to distinguish: the whole entry is the one molecule, which is what
	# free_text says and what _chain_stats() does with it.
	if (!%mol && $free_key && $text =~ /\S/) {
		my $v = _t($text);
		$v =~ s/\A\((.*)\)\z/$1/;    # SOURCE used to parenthesise the organism
		%mol = (1 => { mol_id => 1, $free_key => $v, free_text => 1 });
	}
	return \%mol;
}

# one flat record per chain, so a chain hash can say what molecule it is
sub _entities {
	my ($info) = @_;
	my %by_chain;
	for my $id (keys %{ $info->{compound} }) {
		my $c = $info->{compound}{$id};
		my $s = $info->{source}{$id} || {};
		for my $cid (@{ $c->{chain} || [] }) {
			$by_chain{$cid} = {
				mol_id   => $id,
				molecule => $c->{molecule},
				fragment => $c->{fragment},
				ec       => $c->{ec_number} || $c->{ec},
				organism => $s->{organism_scientific},
				taxid    => $s->{organism_taxid},
				expressed_in => $s->{expression_system},
			};
		}
	}
	$info->{entity_of_chain} = \%by_chain;
	return $info;
}

#
# mmCIF header categories
#
# The same facts, filed differently.  A PDB file says the resolution on a
# REMARK 2 line and an mmCIF file says it in _refine.ls_d_res_high, and a
# caller who wants to know the resolution should not have to care which.  So
# this fills in the same $info keys _parse_meta() fills in, from the
# categories that carry the same information.
#
# Where a fact exists in one format and not the other it is left alone rather
# than invented: an mmCIF file has no REMARK records, so $info->{remarks} stays
# empty, and reading it gets the same "nothing there" a PDB file with no
# remarks would give.
#
# Identifiers are the auth_* ones throughout -- pdbx_strand_id, auth_asym_id --
# because those are the chain ids the coordinates were read under and the ones
# the PDB record carried.  Using label_asym_id here would file the annotations
# under chains that the chains hash does not have.
#

sub _parse_cif_meta {
	my ($info, $p) = @_;
	_meta_defaults($info);
	my $cif   = $p->{cif}       || {};
	my $loops = $p->{cif_loops} || {};

	my $id = _cif1($p, '_entry', 'id');
	$info->{header} = {
		classification => _cif1($p, '_struct_keywords', 'pdbx_keywords'),
		deposit_date   => _cif1($p, '_pdbx_database_status', 'recvd_initial_deposition_date'),
		id_code        => defined $id ? uc $id : '',
	};
	# The data_ block name is deliberately not a key of its own.  It is usually
	# the entry id, which $info->{id} already has, and where it is not -- a file
	# written by a simulation program calls its block 'cell' -- it is worse than
	# the file name _id_from() falls back to.  A key only one of the two formats
	# could ever fill in is a key a caller has to test the format for.

	$info->{title}    = _cif1($p, '_struct', 'title');
	$info->{keywords} = [ grep { length } map { _t($_) }
	                      split /,/, (_cif1($p, '_struct_keywords', 'text') || '') ];
	$info->{experiment} = [ grep { defined && length }
	                        map { $_->{method} } @{ _cif_rows($p, '_exptl') } ];
	$info->{authors}    = [ grep { defined && length }
	                        map { $_->{name} } @{ _cif_rows($p, '_audit_author') } ];

	# resolution: refined structures say so in _refine; the others say it
	# wherever their method says it
	for my $where ([ '_refine', 'ls_d_res_high' ],
	               [ '_reflns', 'd_resolution_high' ],
	               [ '_em_3d_reconstruction', 'resolution' ]) {
		my $v = _cifn($p, @$where);
		next unless defined $v;
		$info->{resolution} = $v;
		last;
	}
	$info->{r_work}      = _cifn($p, '_refine', 'ls_r_factor_r_work');
	$info->{r_free}      = _cifn($p, '_refine', 'ls_r_factor_r_free');
	$info->{temperature} = _cifn($p, '_diffrn', 'ambient_temp');
	$info->{ph}          = _cifn($p, '_exptl_crystal_grow', 'ph');
	my $nmr = _cifn($p, '_pdbx_nmr_ensemble', 'conformers_submitted_total_number');
	$info->{n_models_declared} = $nmr if defined $nmr;

	# the paper.  'primary' is the entry's own citation; anything else in the
	# category is a reference it cites.
	my @cites = @{ _cif_rows($p, '_citation') };
	my ($cite) = ((grep { ($_->{id} || '') eq 'primary' } @cites), @cites);
	if ($cite) {
		my %j;
		$j{titl} = $cite->{title}                    if defined $cite->{title};
		$j{pmid} = $cite->{pdbx_database_id_pubmed}  if defined $cite->{pdbx_database_id_pubmed};
		$j{doi}  = $cite->{pdbx_database_id_doi}     if defined $cite->{pdbx_database_id_doi};
		my @ref = grep { defined && length }
		          @{$cite}{qw(journal_abbrev journal_volume page_first year)};
		$j{ref} = join ' ', @ref if @ref;
		my $cid = $cite->{id};
		my @auth = map { $_->{name} }
		           grep { !defined $cid || !defined $_->{citation_id} || $_->{citation_id} eq $cid }
		           @{ _cif_rows($p, '_citation_author') };
		$j{auth} = [ grep { defined && length } @auth ] if @auth;
		$info->{journal} = \%j;
	}

	if (my @cell = grep { defined } map { _cifn($p, '_cell', $_) }
	               qw(length_a length_b length_c angle_alpha angle_beta angle_gamma)) {
		$info->{cryst1} = {
			a     => _cifn($p, '_cell', 'length_a'),
			b     => _cifn($p, '_cell', 'length_b'),
			c     => _cifn($p, '_cell', 'length_c'),
			alpha => _cifn($p, '_cell', 'angle_alpha'),
			beta  => _cifn($p, '_cell', 'angle_beta'),
			gamma => _cifn($p, '_cell', 'angle_gamma'),
			sgroup => _cif1($p, '_symmetry', 'space_group_name_h-m'),
			z      => _cif1($p, '_cell', 'z_pdb'),
		} if @cell;
	}

	_cif_entities($info, $p);
	_cif_seqres($info, $p);
	_cif_het($info, $p);
	_cif_annotations($info, $p);

	# what the file actually contained, by category, which is the mmCIF answer
	# to the question $info->{records} answers for a PDB file
	my %rec = map { $_ => scalar @{ $loops->{$_} } } keys %$loops;
	for my $tag (keys %$cif) {
		my ($cat) = $tag =~ /\A([^.]+)/;
		$rec{$cat} ||= 1;
	}
	$info->{records} = \%rec;
	return $info;
}

#  the entities, and which chains they are 
#
# COMPND and SOURCE in a PDB file are one _entity plus one _entity_src_* here,
# so they are put back into the shape _entities() already knows how to turn
# into a per-chain record.
sub _cif_entities {
	my ($info, $p) = @_;

	# Only the polymer entities claim chains, which is the rule COMPND follows:
	# a chain is the molecule its polymer is, and the ligands, ions and waters
	# sitting in it are not what it is.  They share the chain id -- the zinc in
	# chain A is written as chain A -- so a non-polymer entity that claimed its
	# chains here would take the chain's name over from the protein, and
	# whichever entity happened to be looked at last would win.  Where the
	# ligands are is $info->{het}, which is where a PDB file keeps it too.
	my %chains_of;
	for my $r (@{ _cif_rows($p, '_entity_poly') }) {
		next unless defined $r->{entity_id};
		$chains_of{ $r->{entity_id} } = [ grep { length } map { _t($_) }
		                                  split /,/, ($r->{pdbx_strand_id} || '') ];
	}

	my (%compound, %source);
	for my $r (@{ _cif_rows($p, '_entity') }) {
		my $e = $r->{id};
		next unless defined $e;
		$compound{$e} = {
			mol_id   => $e,
			molecule => $r->{pdbx_description},
			chain    => $chains_of{$e} || [],
			type     => $r->{type},
		};
		$compound{$e}{ec} = $r->{pdbx_ec} if defined $r->{pdbx_ec};
	}
	# an entity that only shows up in _entity_poly still has to have a record,
	# or its chains lose their molecule name
	for my $e (keys %chains_of) {
		$compound{$e} ||= { mol_id => $e, chain => $chains_of{$e} };
	}

	for my $r (@{ _cif_rows($p, '_entity_src_gen') }) {
		my $e = $r->{entity_id};
		next unless defined $e;
		$source{$e} = {
			mol_id              => $e,
			organism_scientific => $r->{pdbx_gene_src_scientific_name},
			organism_taxid      => $r->{pdbx_gene_src_ncbi_taxonomy_id},
			expression_system   => $r->{pdbx_host_org_scientific_name},
		};
	}
	for my $r (@{ _cif_rows($p, '_entity_src_nat') }) {
		my $e = $r->{entity_id};
		next unless defined $e;
		$source{$e} ||= {
			mol_id              => $e,
			organism_scientific => $r->{pdbx_organism_scientific},
			organism_taxid      => $r->{pdbx_ncbi_taxonomy_id},
		};
	}

	$info->{compound} = \%compound;
	$info->{source}   = \%source;
	_entities($info);
	return $info;
}

#  SEQRES 
#
# _entity_poly_seq is the residue list, one row per position, per entity;
# _entity_poly says which chains an entity was crystallised as.  A chain's
# seqres is therefore its entity's list, and two chains of the same entity get
# the same one -- which is what a PDB file writes out twice.
sub _cif_seqres {
	my ($info, $p) = @_;
	my %res_of;
	for my $r (@{ _cif_rows($p, '_entity_poly_seq') }) {
		next unless defined $r->{entity_id} && defined $r->{mon_id};
		push @{ $res_of{ $r->{entity_id} } }, $r->{mon_id};
	}
	for my $r (@{ _cif_rows($p, '_entity_poly') }) {
		my $e = $r->{entity_id};
		next unless defined $e;
		my @chains = grep { length } map { _t($_) } split /,/, ($r->{pdbx_strand_id} || '');
		next unless @chains;
		my $residues = $res_of{$e};
		my $seq;
		if ($residues) {
			$seq = join '', map { my $o = res1($_); length $o ? $o : 'X' } @$residues;
		}
		else {
			# a file with no _entity_poly_seq still carries the sequence as a
			# string; the residue names are what is lost, not the sequence
			my $one = $r->{pdbx_seq_one_letter_code_can} || $r->{pdbx_seq_one_letter_code};
			next unless defined $one;
			$seq = $one;
			$seq =~ s/\s+//g;
			$seq =~ s/\([^)]*\)/X/g;   # a modified residue, spelled out in brackets
		}
		for my $cid (@chains) {
			$info->{seqres}{$cid} = {
				chain    => $cid,
				residues => ($residues ? [ @$residues ] : []),
				sequence => $seq,
				length   => ($residues ? scalar @$residues : length $seq),
			};
		}
	}
	return $info;
}

#  heterogens
#
# _chem_comp describes every residue in the file, standard ones included;
# $info->{het} is what HET/HETNAM/FORMUL describe, which is the rest.
sub _cif_het {
	my ($info, $p) = @_;
	for my $r (@{ _cif_rows($p, '_chem_comp') }) {
		my $cid = $r->{id};
		next unless defined $cid && length $cid;
		next if $STANDARD{$cid};
		my $h = $info->{het}{$cid} ||= { het_id => $cid };
		$h->{name}    = $r->{name}            if defined $r->{name};
		$h->{formula} = $r->{formula}         if defined $r->{formula};
		$h->{synonym} = $r->{pdbx_synonyms}   if defined $r->{pdbx_synonyms};
		$h->{water}   = 1 if $cid eq 'HOH' || $cid eq 'DOD' || $cid eq 'WAT';
	}
	for my $r (@{ _cif_rows($p, '_pdbx_nonpoly_scheme') }) {
		my $cid = $r->{mon_id};
		next unless defined $cid && length $cid;
		my $h = $info->{het}{$cid} ||= { het_id => $cid };
		push @{ $h->{instances} }, {
			chain  => $r->{pdb_strand_id},
			resseq => $r->{pdb_seq_num},
			icode  => (defined $r->{pdb_ins_code} ? $r->{pdb_ins_code} : ''),
		};
	}
	for my $r (@{ _cif_rows($p, '_pdbx_struct_mod_residue') }) {
		my $cid = $r->{auth_comp_id} || $r->{label_comp_id};
		next unless defined $cid && length $cid;
		$info->{modres}{$cid} ||= {
			resname  => $cid,
			standard => $r->{parent_comp_id},
			comment  => $r->{details},
		};
	}
	return $info;
}

#  secondary structure, bonds and database cross-references 
sub _cif_annotations {
	my ($info, $p) = @_;

	for my $r (@{ _cif_rows($p, '_struct_conf') }) {
		next unless ($r->{conf_type_id} || '') =~ /\AHELX/i;
		push @{ $info->{helix} }, {
			id           => $r->{id},
			init_resname => _cif_auth($r, 'beg', 'comp_id'),
			init_chain   => _cif_auth($r, 'beg', 'asym_id'),
			init_resseq  => _cif_auth($r, 'beg', 'seq_id'),
			end_resname  => _cif_auth($r, 'end', 'comp_id'),
			end_chain    => _cif_auth($r, 'end', 'asym_id'),
			end_resseq   => _cif_auth($r, 'end', 'seq_id'),
			class        => $r->{pdbx_pdb_helix_class},
			length       => $r->{pdbx_pdb_helix_length},
		};
	}
	for my $r (@{ _cif_rows($p, '_struct_sheet_range') }) {
		push @{ $info->{sheet} }, {
			strand       => $r->{id},
			id           => $r->{sheet_id},
			init_resname => _cif_auth($r, 'beg', 'comp_id'),
			init_chain   => _cif_auth($r, 'beg', 'asym_id'),
			init_resseq  => _cif_auth($r, 'beg', 'seq_id'),
			end_resname  => _cif_auth($r, 'end', 'comp_id'),
			end_chain    => _cif_auth($r, 'end', 'asym_id'),
			end_resseq   => _cif_auth($r, 'end', 'seq_id'),
		};
	}

	# SSBOND and LINK are one category here, told apart by the bond type
	for my $r (@{ _cif_rows($p, '_struct_conn') }) {
		my $type = lc($r->{conn_type_id} || '');
		if ($type eq 'disulf') {
			push @{ $info->{ssbond} }, {
				chain1  => _cif_ptnr($r, 1, 'asym_id'),
				resseq1 => _cif_ptnr($r, 1, 'seq_id'),
				chain2  => _cif_ptnr($r, 2, 'asym_id'),
				resseq2 => _cif_ptnr($r, 2, 'seq_id'),
				length  => $r->{pdbx_dist_value},
			};
			next;
		}
		next if $type eq 'hydrog';
		push @{ $info->{link} }, {
			name1    => $r->{ptnr1_label_atom_id}, resname1 => _cif_ptnr($r, 1, 'comp_id'),
			chain1   => _cif_ptnr($r, 1, 'asym_id'), resseq1 => _cif_ptnr($r, 1, 'seq_id'),
			name2    => $r->{ptnr2_label_atom_id}, resname2 => _cif_ptnr($r, 2, 'comp_id'),
			chain2   => _cif_ptnr($r, 2, 'asym_id'), resseq2 => _cif_ptnr($r, 2, 'seq_id'),
			length   => $r->{pdbx_dist_value},
		};
	}

	for my $r (@{ _cif_rows($p, '_struct_mon_prot_cis') }) {
		push @{ $info->{cispep} }, {
			resname1 => $r->{auth_comp_id}, chain1 => $r->{auth_asym_id},
			resseq1  => $r->{auth_seq_id},
			resname2 => $r->{pdbx_auth_comp_id_2}, chain2 => $r->{pdbx_auth_asym_id_2},
			resseq2  => $r->{pdbx_auth_seq_id_2},
			angle    => $r->{pdbx_omega_angle},
		};
	}

	my %db;
	for my $r (@{ _cif_rows($p, '_struct_ref') }) {
		next unless defined $r->{id};
		$db{ $r->{id} } = $r;
	}
	for my $r (@{ _cif_rows($p, '_struct_ref_seq') }) {
		my $cid = $r->{pdbx_strand_id};
		next unless defined $cid && length $cid;
		my $ref = (defined $r->{ref_id} && $db{ $r->{ref_id} }) || {};
		push @{ $info->{dbref}{$cid} }, {
			chain     => $cid,
			seq_begin => $r->{pdbx_auth_seq_align_beg},
			seq_end   => $r->{pdbx_auth_seq_align_end},
			database  => $ref->{db_name},
			accession => (defined $r->{pdbx_db_accession} ? $r->{pdbx_db_accession} : $ref->{pdbx_db_accession}),
			db_id     => $ref->{db_code},
			db_begin  => $r->{db_align_beg},
			db_end    => $r->{db_align_end},
		};
	}

	for my $r (@{ _cif_rows($p, '_pdbx_audit_revision_history') }) {
		push @{ $info->{revdat} }, {
			num  => $r->{ordinal},
			date => $r->{revision_date},
			id   => $info->{header}{id_code},
			type => $r->{data_content_type},
		};
	}
	return $info;
}

#  reading the parsed categories 

# _cif_rows($p, $category) -- a category as a list of rows, whether it was
# written as a loop_ or, having only one row, as a run of plain tags.  The
# format allows both for the same category and files use both, so asking for
# the rows has to work either way.
#
# The two are folded together once, on the first call, rather than on each:
# reading the header asks for forty-odd categories and a real entry has
# several hundred plain tags, so doing it per call would walk the lot forty
# times over to answer forty questions.
sub _cif_rows {
	my ($p, $cat) = @_;
	my $by_cat = $p->{_by_category} ||= do {
		my %c;
		while (my ($tag, $val) = each %{ $p->{cif} || {} }) {
			my $dot = index($tag, '.');
			# a core CIF tag has no category half, and is its own category
			my ($k, $item) = $dot < 0 ? ($tag, $tag)
			                          : (substr($tag, 0, $dot), substr($tag, $dot + 1));
			($c{$k} ||= [ {} ])->[0]{$item} = $val;
		}
		# a loop_ is the category, where there is one: a file that wrote both
		# meant the loop, since that is the one that can hold what it holds
		%c = (%c, %{ $p->{cif_loops} || {} });
		\%c;
	};
	return $by_cat->{$cat} || [];
}

# one item from a category that has one row, which is most of them
sub _cif1 {
	my ($p, $cat, $item) = @_;
	my $rows = _cif_rows($p, $cat);
	return undef unless @$rows;
	my $v = $rows->[0]{$item};
	return defined $v && length $v ? $v : undef;
}

# the same, as a number.  '?' and '.' already came back as undef; what is left
# to guard against is a field holding text where a number belongs.
sub _cifn {
	my $v = _cif1(@_);
	return defined $v ? _n($v) : undef;
}

# auth_* first, then label_*: the same rule the coordinates were read under,
# so an annotation and the chain it annotates agree about the chain's name
sub _cif_auth {
	my ($r, $which, $item) = @_;
	for my $k ("${which}_auth_$item", "${which}_label_$item") {
		return $r->{$k} if defined $r->{$k} && length $r->{$k};
	}
	return undef;
}

sub _cif_ptnr {
	my ($r, $n, $item) = @_;
	for my $k ("ptnr${n}_auth_$item", "ptnr${n}_label_$item") {
		return $r->{$k} if defined $r->{$k} && length $r->{$k};
	}
	return undef;
}

#  small helpers 

sub _t {
	my ($s) = @_;
	return '' unless defined $s;
	$s =~ s/\A\s+//;
	$s =~ s/\s+\z//;
	return $s;
}

# _c($line, $from, $length) -- the trimmed contents of a fixed field, clipped
# to what the line actually has.  Records in real files are right-trimmed, so
# any field can begin past the end of its line, and that is not damage: it is
# a field the depositor left empty.  Every fixed-column read below goes
# through this rather than through substr() directly, because a bare substr()
# dies on 15 of the 10,116 entries in PDBbind v2020 -- a SHEET record with no
# sense field on it is enough.
sub _c {
	my ($line, $from, $length) = @_;
	return '' if !defined $line || $from >= length $line;
	my $s = defined $length ? substr($line, $from, $length) : substr($line, $from);
	$s =~ s/\A\s+//;
	$s =~ s/\s+\z//;
	return $s;
}

sub _n {
	my $v = _t($_[0]);
	return $v =~ /\A[-+]?[\d.]+(?:[eE][-+]?\d+)?\z/ ? $v + 0 : undef;
}

# continuation records: text from column $from on, glued back together.  A
# hyphen at the end of a line is a real hyphen in the middle of a word --
# KEYWDS breaks "COMPLEX (HORMONE-" / "RECEPTOR)" across lines -- so it joins
# without a space; anything else takes one.
sub _joined {
	my ($lines, $from, $id) = @_;
	return '' unless $lines && @$lines;
	my $out = '';
	for my $l (@$lines) {
		my $t = _c(_untail($l, $id), $from);
		next unless length $t;
		$out = _rejoin($out, $t);
	}
	return $out;
}

# The text records run to the end of the line, and in a file old enough to keep
# its entry id in columns 73-80 the end of the line is '1GDR   3'.  Given the id
# the HEADER carried, a line whose columns 73-80 hold nothing but that id and a
# line number is cut there: 'GAMMA DELTA RESOLVASE' is the compound and the id
# is the stationery.  Nothing else is cut -- a modern file uses those columns
# for text, and text is not the entry id followed by a number -- so a title that
# runs to column 80 is left alone.
sub _untail {
	my ($line, $id) = @_;
	return $line unless defined $line && defined $id && length $id;
	return $line unless length($line) > 72 && $id =~ /\A\w{1,4}\z/;
	return substr($line, 72) =~ /\A\s*\Q$id\E(?:\s+\d+)?\s*\z/i
		? substr($line, 0, 72)
		: $line;
}

sub _rejoin {
	my ($have, $add) = @_;
	return $add unless defined $have && length $have;
	return $have . $add if $have =~ /-\z/;
	return "$have $add";
}

sub _check_info {
	my ($info, $who) = @_;
	die "$who: expected the hash reference from structure_info()"
		unless defined $info && (reftype($info) || '') eq 'HASH' && exists $info->{chains};
	return 1;
}

#
# Help
#
# h() prints a function's own documentation, in the spirit of R's ?function.
# The text is this file's POD, read at run time, so the help and the shipped
# documentation cannot drift apart.
#
sub h {
	my ($what) = @_;
	my $name = _help_name($what);
	my $sec  = _pod_sections();
	if (defined $name && $sec->{$name}) {
		print STDOUT $sec->{$name};
		return $name;
	}
	print STDOUT "Chem::Structure::Parser $VERSION\n\nDocumented functions:\n";
	print STDOUT "    $_\n" for sort keys %$sec;
	print STDOUT "\nCall h('structure_info') for one of them.\n";
	return undef;
}

sub _help_name {
	my ($what) = @_;
	return undef unless defined $what;
	my $r = ref $what;
	if ($r eq 'CODE') {
		require B;
		my $gv = B::svref_2object($what)->GV;
		return $gv->NAME;
	}
	my $n = "$what";      # a glob stringifies as *Chem::Structure::Parser::res1
	$n =~ s/\A\*//;
	$n =~ s/\A.*:://;
	return $n;
}

sub _pod_sections {
	my %sec;
	# Only the functions.  The POD is generated from README.md, whose '# Changes'
	# section writes each release as '## 0.01 2026-08-21 CDT' -- a =head2 like
	# any other, which listed the version number among the documented functions
	# until this filtered on the export list.
	my %exported = map { $_ => 1 } @EXPORT_OK;
	open my $fh, '<', __FILE__;
	my ($in, $name);
	while (my $l = <$fh>) {
		if ($l =~ /\A=head2\s+(\S+)/) {
			$name = $1;
			$name =~ s/\(.*//;
			if (!$exported{$name}) { $in = 0; $name = undef; next }
			$in = 1;
			$sec{$name} = '';
			next;
		}
		if ($l =~ /\A=head[12]\b/ || $l =~ /\A=cut/) { $in = 0; next }
		# A subsection of a function -- '### Options' in README.md, '=head3
		# Options' here -- prints as its own title and not as the POD command
		# that carries it.
		if ($l =~ /\A=head[3-9]\s+(\S.*?)\s*\z/) {
			$sec{$name} .= "$1\n" if $in && defined $name;
			next;
		}
		$sec{$name} .= $l if $in && defined $name;
	}
	close $fh;
	s/\A\n+//, s/\n+\z/\n/ for values %sec;
	return \%sec;
}

1;

__END__

=encoding utf8

=head1 NAME

Chem::Structure::Parser - Read a molecular structure file into a hash of hashes, sequences and all, using XS for the coordinate section

=head1 Synopsis

Read a molecular structure file and get everything in it back as a hash of
hashes — the header, the annotations, every chain, every residue, every atom,
and the single-letter sequence of each chain — in one call.

 use Chem::Structure::Parser;

 my $info = structure_info('1a22.ent.pdb');

 print $info->{id};                                  # 1A22
 print $info->{resolution};                          # 2.6
 print $info->{chains}{A}{sequence};                 # FPTIPLSRLFDNAMLRAHRLHQL...
 print $info->{chains}{A}{molecule};                 # GROWTH HORMONE
 print $info->{chains}{A}{residues}{54}{resname};    # PHE
 print $info->{chains}{A}{residues}{54}{atoms}{CA}{x};
 print $info->{chains}{A}{elements}{S};              # 7    sulphur atoms in chain A
 print $info->{stats}{elements}{S};                  # 17   and in the structure

 print structure_summary($info);

which prints

 1A22  HUMAN GROWTH HORMONE BOUND TO SINGLE RECEPTOR
   file        1a22.ent.pdb
   method      X-RAY DIFFRACTION
   resolution  2.6 A
   R / R-free  0.187 / -
   models      1
   atoms       3113 (69 hetatm, 69 water)
   chain A     protein      206 residues,  1492 atoms, 2 gaps
               FPTIPLSRLFDNAMLRAHRLHQLAFDTYQEFEEAYIPKEQKYSFLQNPQTSLCFSESIPTP...
               GROWTH HORMONE
   chain B     protein      235 residues,  1621 atoms, 2 gaps
               PKFTKCRSPERETFSCHWTLGPIQLFYTRRNTQEWTQEWKECPDYVSAGENSCYFNSSFTS...
               GROWTH HORMONE RECEPTOR

What the structure I<is> comes back the same way, in one more call:

 my $f = structure_features($info);

 print $f->{sasa}{total};          # 17805.0  solvent-accessible surface, A^2
 print $f->{rg};                   # 22.85    radius of gyration, A
 print $f->{mass};                 # 41307.1  dalton
 print $f->{hydropathy};           # -0.452   mean Kyte-Doolittle hydropathy
 print $f->{aromatic_fraction};    # 0.1263   the F, W and Y share of the sequence
 print scalar @{ $f->{pi_stacking} };  # 2    stacked pairs of aromatic rings

 print $info->{chains}{A}{residues}{54}{rsa};   # 0.103  PHE 54 is mostly buried

which is the Shrake-Rupley surface, the ring geometry, the radii and the masses
as mdtraj implements them and the two sequence indices as Biopython does. All of
it is in C, for the same reason the parse is.

The coordinate section is parsed in C, because across a directory of
structures it is millions of lines: the largest entry in PDBbind v2020 is
33 MB and 411,648 atom records, and it reads in about 1.5 seconds. The header
records are parsed in Perl, because they are irregular and there are only a
few dozen of them in a file.

The module is called C<Chem::Structure::Parser> and not C<PDB::Info> because the
shape of what it hands back has nothing to do with the format it came out of.
It reads PDB and mmCIF/PDBx; C<formats()> says what it reads at any moment, and
a format it knows the name of but cannot read yet says so rather than
misreading it.

=head1 PDB and mmCIF

Reading is the same call either way. C<structure_info()> works out the format
from the file name — C<.pdb>, C<.ent>, C<.cif>, C<.mmcif>, C<.pdbx> — and from the
first records in the file when the name gives nothing away, and the hash that
comes back has the same keys, the same nesting and the same values whichever
it was.

 my $a = structure_info('1a22.pdb');
 my $b = structure_info('1a22.cif');

 $a->{chains}{A}{sequence} eq $b->{chains}{A}{sequence};              # true
 $a->{chains}{A}{residues}{54}{atoms}{CA}{x}
     == $b->{chains}{A}{residues}{54}{atoms}{CA}{x};                  # true

So no calling code branches on the format, and a script written against a
directory of C<.pdb> files works unchanged on a directory of C<.cif> ones.

Equality here means equality rather than approximately: C<t/cif.t> reads
fixture pairs both ways and compares the whole coordinate half of the returned
structure with C<is_deeply>, and C<t/real_cif.t> converts real entries from the
PDB archive into mmCIF and asserts that every chain, residue, atom and count
comes back identical.

Two consequences are worth knowing.

B<The identifiers are the auth_* ones.> An mmCIF file carries two sets: the
C<label_*> identifiers the archive assigns, and the C<auth_*> ones the depositor
used. Only C<auth_*> matches what the PDB record carried, so those are the
chain ids and residue numbers used throughout — in the coordinates and in the
annotations alike. A structure read from a C<.cif> therefore has the same chain
C<A> and the same residue C<54> as the same structure read from a C<.pdb>, not
the C<label_asym_id> lettering that runs on through the waters.

B<Values are converted, not passed through.> Where the two formats spell the
same fact differently, the mmCIF reader produces what the PDB reader would
have: C<_atom_site.pdbx_formal_charge> of C<-1> reads back as C<'1-'>, and C<.>
and C<?> — mmCIF for "not applicable" and "unknown" — read back as the empty
field a PDB record would have had. A charge of C<0> is kept as C<'0'>, because
"the field said zero" and "the field was blank" are different answers and both
formats can say either.

What is I<not> the same is what only one of the formats has. An mmCIF file has
no REMARK records, so C<< $info-E<gt>{remarks} >> is empty for one; a PDB file has no
C<_entity> category, so a chain read from one may not know which entity it
belongs to. Every key is present in both cases, so reading one is a test of
what the file said and never of which format it was.

Where the same fact is filed under different names, it is folded into the same
key:



=begin html

<table>
<thead>
<tr>
  <th><code>$info</code> key</th>
  <th>PDB record</th>
  <th>mmCIF category</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>title</code></td>
  <td><code>TITLE</code></td>
  <td><code>_struct.title</code></td>
</tr>
<tr>
  <td><code>id</code></td>
  <td><code>HEADER</code></td>
  <td><code>_entry.id</code></td>
</tr>
<tr>
  <td><code>experiment</code></td>
  <td><code>EXPDTA</code></td>
  <td><code>_exptl.method</code></td>
</tr>
<tr>
  <td><code>resolution</code></td>
  <td><code>REMARK 2</code>, then <code>REMARK 3</code></td>
  <td><code>_refine.ls_d_res_high</code></td>
</tr>
<tr>
  <td><code>r_work</code>, <code>r_free</code></td>
  <td><code>REMARK 3</code></td>
  <td><code>_refine.ls_R_factor_R_*</code></td>
</tr>
<tr>
  <td><code>keywords</code></td>
  <td><code>KEYWDS</code></td>
  <td><code>_struct_keywords.text</code></td>
</tr>
<tr>
  <td><code>authors</code></td>
  <td><code>AUTHOR</code></td>
  <td><code>_audit_author</code></td>
</tr>
<tr>
  <td><code>journal</code></td>
  <td><code>JRNL</code></td>
  <td><code>_citation</code>, <code>_citation_author</code></td>
</tr>
<tr>
  <td><code>compound</code>, <code>source</code></td>
  <td><code>COMPND</code>, <code>SOURCE</code></td>
  <td><code>_entity</code>, <code>_entity_src_*</code></td>
</tr>
<tr>
  <td><code>seqres</code></td>
  <td><code>SEQRES</code></td>
  <td><code>_entity_poly</code>, <code>_entity_poly_seq</code></td>
</tr>
<tr>
  <td><code>het</code></td>
  <td><code>HET</code>, <code>HETNAM</code>, <code>FORMUL</code></td>
  <td><code>_chem_comp</code>, <code>_pdbx_nonpoly_scheme</code></td>
</tr>
<tr>
  <td><code>helix</code></td>
  <td><code>HELIX</code></td>
  <td><code>_struct_conf</code></td>
</tr>
<tr>
  <td><code>sheet</code></td>
  <td><code>SHEET</code></td>
  <td><code>_struct_sheet_range</code></td>
</tr>
<tr>
  <td><code>ssbond</code>, <code>link</code></td>
  <td><code>SSBOND</code>, <code>LINK</code></td>
  <td><code>_struct_conn</code></td>
</tr>
<tr>
  <td><code>cispep</code></td>
  <td><code>CISPEP</code></td>
  <td><code>_struct_mon_prot_cis</code></td>
</tr>
<tr>
  <td><code>modres</code></td>
  <td><code>MODRES</code></td>
  <td><code>_pdbx_struct_mod_residue</code></td>
</tr>
<tr>
  <td><code>dbref</code></td>
  <td><code>DBREF</code></td>
  <td><code>_struct_ref</code>, <code>_struct_ref_seq</code></td>
</tr>
<tr>
  <td><code>cryst1</code></td>
  <td><code>CRYST1</code></td>
  <td><code>_cell</code>, <code>_symmetry</code></td>
</tr>
<tr>
  <td><code>n_models</code></td>
  <td><code>MODEL</code></td>
  <td><code>_atom_site.pdbx_PDB_model_num</code></td>
</tr>
</tbody>
</table>

=end html



=head1 Installing

 perl Makefile.PL
 make
 make test
 make install

C<make test> reads the fixtures in C<t/data>. If a directory of real structures
is to hand it reads a sample of those too; point it somewhere with

 STRUCTURE_INFO_TEST_DIR=/path/to/pdbs  make test
 STRUCTURE_INFO_TEST_CIF_DIR=/path/to/cifs  make test
 STRUCTURE_INFO_TEST_ALL=1 STRUCTURE_INFO_TEST_DIR=/path make test   # all of them

C<STRUCTURE_INFO_TEST_DIR> is used twice: C<t/real.t> reads those files as PDB
and checks them against a second reader written in plain Perl, and
C<t/real_cif.t> converts each one into mmCIF and asserts that reading it back
gives the same structure to the last digit.
C<STRUCTURE_INFO_TEST_CIF_DIR> takes a directory of real C<.cif> files, either
flat or one subdirectory per structure, and reads those directly.

With no such directory those tests skip, so the distribution builds on a
machine with no structures on it.

=head1 Getting help

C<h> prints any function's section of this document to C<STDOUT> and returns, in
the spirit of R's C<?function> at the prompt. It takes the name three ways:

 h('structure_info');    # by name
 h(*res_type);           # by name, unquoted
 h(\&aa3to1);            # by reference
 h();                    # the list of documented functions

 perl -MChem::Structure::Parser -e 'h(*structure_info)'   # straight from the shell

Note that C<h(res_type)>, with no quotes and no sigil, cannot be made to work:
every function here is exported, so Perl parses the bareword as a call to
C<res_type()> before C<h> is ever reached. Use one of the three forms above.

=head1 Functions/Subroutines

=head2 structure_info

 my $info = structure_info($file, %options);

Reads C<$file> and returns a hash reference. The format is worked out from the
file name — C<.pdb>, C<.ent>, C<.cif>, C<.mmcif>, C<.pdbx> — and from the first
records in the file when the name gives nothing away. C<.gz> files are read as
they are, without unpacking to a temporary file.

=head3 What comes back

Laid out the way C<tree> lays out a directory, this is C<1a22.ent.pdb> — a real
file, real values, the long lists cut short:

 $info
 ├── file            '1a22.ent.pdb'          the path it was read from
 ├── format          'pdb'                   or 'mmcif'
 ├── id              '1A22'                  from HEADER, or from the file name
 ├── title           'HUMAN GROWTH HORMONE BOUND TO SINGLE RECEPTOR'
 ├── header
 │   ├── classification  'COMPLEX (HORMONE/RECEPTOR)'
 │   ├── deposit_date    '15-JAN-98'
 │   └── id_code         '1A22'
 ├── experiment      [ 'X-RAY DIFFRACTION' ]
 ├── resolution      2.6                     REMARK 2
 ├── r_work          0.187                   REMARK 3
 ├── r_free          undef                   this entry does not report one
 ├── temperature     287                     REMARK 200
 ├── ph              6.5
 ├── keywords        [ 'COMPLEX (HORMONE-RECEPTOR)', 'PITUITARY HORMONE', ... ]
 ├── authors         [ 'A.M.DE VOS', 'M.ULTSCH' ]
 ├── journal
 │   ├── auth        [ 'T.CLACKSON', 'M.H.ULTSCH', 'J.A.WELLS', 'A.M.DE VOS' ]
 │   ├── titl        'STRUCTURAL AND FUNCTIONAL ANALYSIS OF THE 1:1 GROWTH...'
 │   ├── ref         'J.MOL.BIOL.                   V. 277  1111 1998'
 │   ├── refn        'ISSN 0022-2836'
 │   ├── pmid        '9571026'
 │   └── doi         '10.1006/JMBI.1998.1669'
 ├── compound                                COMPND, by MOL_ID
 │   ├── 1
 │   │   ├── mol_id      '1'
 │   │   ├── molecule    'GROWTH HORMONE'
 │   │   ├── chain       [ 'A' ]
 │   │   ├── engineered  'YES'
 │   │   └── mutation    'YES'
 │   └── 2           { molecule 'GROWTH HORMONE RECEPTOR', chain [ 'B' ],
 │                     fragment 'EXTRACELLULAR DOMAIN', engineered 'YES' }
 ├── source                                  SOURCE, by MOL_ID
 │   └── 1           { organism_scientific 'HOMO SAPIENS', organism_common
 │                     'HUMAN', organism_taxid '9606', mol_id '1',
 │                     expression_system 'ESCHERICHIA COLI',
 │                     expression_system_taxid '562' }
 ├── entity_of_chain                         COMPND and SOURCE, by chain
 │   ├── A           { mol_id '1', molecule 'GROWTH HORMONE', fragment undef,
 │   │                 ec undef, organism 'HOMO SAPIENS', taxid '9606',
 │   │                 expressed_in 'ESCHERICHIA COLI' }
 │   └── B           { ..., fragment 'EXTRACELLULAR DOMAIN' }
 ├── seqres                                  what SEQRES says was in the crystal
 │   ├── A
 │   │   ├── sequence    'FPTIPLSRLFDNAMLRAHRLHQLAFDTYQEFEEAYIPKEQKYSFLQ...'
 │   │   ├── residues    [ 'PHE', 'PRO', 'THR', 'ILE', ... ]        191 of them
 │   │   └── length      191
 │   └── B               { sequence, residues, length 238 }
 ├── dbref
 │   └── A           [ { database 'UNP', accession 'P01241',
 │                       db_id 'SOMA_HUMAN', seq_begin '1', seq_end '191',
 │                       db_begin '27', db_end '217', chain 'A' } ]
 ├── seqadv          [ { chain 'A', resseq '120', resname 'ARG',
 │                       db_res 'GLY', db_seq '146', comment 'ENGINEERED' } ]
 ├── modres          { }                     no MSE-style residues here
 ├── het
 │   └── HOH         { het_id 'HOH', formula '69(H2 O)', water 1 }
 ├── hetnam          { }
 ├── formul          { }
 ├── helix           [ { id '1', class '1', length '29',
 │                       init_chain 'A', init_resname 'SER', init_resseq '7',
 │                       end_chain 'A', end_resname 'TYR', end_resseq '35' },
 │                     ... ]                                     12 of them
 ├── sheet           [ ... ]                                     12
 ├── ssbond          [ { chain1 'A', resseq1 '53',
 │                       chain2 'A', resseq2 '165', length '2.02' }, ... ]  5
 ├── link            [ ]
 ├── cispep          [ ]
 ├── site            [ ]
 ├── cryst1          { a '67.7', b '67.7', c '228',
 │                     alpha '90', beta '90', gamma '90',
 │                     sgroup 'P 43 21 2', z '8' }
 ├── biological_assembly  [ 32 lines of REMARK 350, verbatim ]
 ├── revdat          [ { num '3', date '18-APR-18', id '1A22',
 │                       type '1', what 'REMARK' }, ... ]
 ├── remarks                                 every REMARK, by number
 │   ├── 2           [ '', 'RESOLUTION.    2.60 ANGSTROMS.' ]
 │   ├── 350         [ ... ]                                     32 lines
 │   └── ...         1, 3, 4, 100, 200, 280, 290, 300, 465, 470, 500
 ├── conect          [ [ 448, 1255 ], ... ]                      10
 ├── records                                 every record type, counted
 │   ├── REMARK      365
 │   ├── SEQRES      34
 │   ├── HELIX       12
 │   └── ...         AUTHOR, COMPND, CONECT, CRYST1, DBREF, SOURCE, SSBOND, ...
 ├── n_models        1                       how many MODEL records the file has
 ├── model           1                       which one the chains below are
 ├── models                                  there only with model => 'all'
 ├── stats
 │   ├── n_atoms         3113    atoms kept: this model, less what was filtered
 │   ├── total_atoms     3113    atoms the file has, every model, unfiltered
 │   ├── n_hetatm        69      of n_atoms, the ones written as HETATM
 │   ├── n_hydrogens     0
 │   ├── n_water_atoms   69
 │   ├── n_lines         3605
 │   ├── n_atom_records  3044    ATOM lines seen, whether kept or not
 │   ├── n_hetatm_records 69     HETATM lines, likewise
 │   ├── n_anisou        0
 │   ├── n_skipped       0       coordinate lines the options threw away
 │   ├── elements        { C 1946, O 643, N 507, S 17 }
 │   │                           every element in the file, keyed by its IUPAC
 │   │                           symbol; the counts add up to n_atoms
 │   ├── bfactor         { min '2.7', max '85.39', mean 30.83, n 3113 }
 │   ├── bbox            { xmin '12.142', xmax '80.34', ymin '2.011', ... }
 │   └── center          [ '46.241', '29.135', '134.559' ]
 ├── chain_order     [ 'A', 'B' ]            the order the file has them in
 └── chains
     ├── A
     │   ├── id              'A'
     │   ├── type            'protein'   protein dna rna water hetero unknown
     │   ├── sequence        'FPTIPLSRLFDNAMLRAHRLHQLAFDTYQEFEEAYIPKEQ...'
     │   │                               single-letter, what has coordinates
     │   ├── seqres          'FPTIPLSRLFDNAMLRAHRLHQLAFDTYQEFEEAYIPKEQ...'
     │   ├── seqres_length   191
     │   ├── n_residues      206
     │   ├── n_polymer       180
     │   ├── n_water         26
     │   ├── n_ligand        0
     │   ├── n_atoms         1492
     │   ├── n_hetatm        26
     │   ├── elements        { C 938, O 301, N 246, S 7 }
     │   │                               the same tally for this chain alone;
     │   │                               adds up to the chain's n_atoms
     │   ├── n_missing       11          SEQRES less what was modelled
     │   ├── gaps            [ { after 129, before 136, missing 6 },
     │   │                     { after 148, before 154, missing 5 } ]
     │   ├── n_gaps          2
     │   ├── missing_residues
     │   │                   [ 130, 131, 132, 133, 134, 135,
     │   │                     149, 150, 151, 152, 153 ]
     │   ├── first           1           the first and last polymer residue keys
     │   ├── last            191
     │   ├── residue_types   { amino_acid 180, water 26 }
     │   ├── molecule        'GROWTH HORMONE'            from COMPND
     │   ├── organism        'HOMO SAPIENS'              from SOURCE
     │   ├── mol_id          '1'
     │   ├── dbref           [ { ... } ]     as in the top-level dbref
     │   │                                   ec and fragment are here too, in a
     │   │                                   chain whose file gives them
     │   ├── residue_order   [ '1', '2', '3', ... '574' ]        file order, 206
     │   └── residues                    keyed number + insertion code
     │       ├── 54
     │       │   ├── resname     'PHE'
     │       │   ├── number      54
     │       │   ├── icode       ''
     │       │   ├── key         '54'
     │       │   ├── chain       'A'
     │       │   ├── one         'F'     '' when there is no letter for it
     │       │   ├── type        'amino_acid'
     │       │   │                       nucleotide water ligand ion
     │       │   ├── standard    1       one of the twenty, or a standard base
     │       │   ├── modified    0       1 for MSE, still an M in the sequence
     │       │   ├── hetero      0       1 when it was written as HETATM
     │       │   ├── free                not here; 1 for a free amino acid
     │       │   │                       bound in a site (see below)
     │       │   ├── n_atoms     11
     │       │   ├── b_mean      22.55
     │       │   ├── center      [ 65.311, 17.127, 140.515 ]
     │       │   ├── atom_order  [ 'N', 'CA', 'C', 'O', 'CB', ... ]
     │       │   └── atoms
     │       │       ├── CA
     │       │       │   ├── name       'CA'
     │       │       │   ├── serial     450
     │       │       │   ├── element    'C'
     │       │       │   ├── charge     ''
     │       │       │   ├── x          '66.446'
     │       │       │   ├── y          '18.25'
     │       │       │   ├── z          '141.982'
     │       │       │   ├── occupancy  '1'
     │       │       │   ├── bfactor    '24.53'
     │       │       │   ├── altloc     ''
     │       │       │   ├── hetero     0
     │       │       │   └── altlocs    [ { altloc, x, y, z, occupancy,
     │       │       │                      bfactor }, ... ]
     │       │       │                  present only when the atom has
     │       │       │                  alternate conformers; every conformer
     │       │       │                  is listed, the chosen one included,
     │       │       │                  and one of them having no letter at
     │       │       │                  all does not take it off the list
     │       │       └── ...    N, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ
     │       └── ...            1 .. 191, then the waters at 512 .. 574
     └── B                      the same again: 235 residues, 1621 atoms

A record that is not in the file reads as C<undef>, and a list that is not in
the file reads as an empty arrayref — C<title> being C<undef> means there was no
TITLE, which is a different thing from a TITLE that was blank.

Everything the module does not take apart is still in C<remarks> and in the
raw record counts, so nothing in the file is lost.

=head3 The two sequences

C<sequence> and C<seqres> are the two different questions people mean by "the
sequence": what was modelled, and what was in the crystal. They differ
wherever a terminus or a loop went unmodelled, which is what C<gaps> counts and
C<n_missing> totals — eleven residues of chain A above, in two stretches.
C<missing_residues> is the same eleven one number at a time, in ascending
order, for asking whether a particular residue was modelled without walking
the gap list.

Both are read off the numbering, so they see the loops a chain skips over and
not the residues that fell off either end — a terminus that went unmodelled
leaves no numbering behind to notice it by, and only C<n_missing> counts those.
Numbering is not always sequential, either: an antibody numbered by the Kabat
scheme runs 27, 1027, 2027, 28, where the thousands are insertions after 27
and not a 999-residue hole. A chain can only be missing as many residues as
the span from its first polymer residue to its last leaves room for, so a jump
wider than that is taken for a change of numbering scheme and not counted.

A scheme that skips a few numbers on purpose is not caught, and cannot be:
a protein numbered by homology to a reference one — chymotrypsin numbering,
and the several conventions like it — leaves unused the numbers its reference
does not need, and the coordinates do not say which of those is a residue that
went unmodelled. 1ahx has all 396 of its SEQRES residues modelled and still
skips nine numbers, which both fields report. About one chain in twenty-five
is like this. Where a file has SEQRES and the two disagree, C<n_missing> is the
one to trust: it counts residues, where C<gaps> and C<missing_residues> count
numbers.

=head3 How residues are keyed

By residue number with the insertion code appended, so C<100>, C<100A> and
C<100B> are three separate keys and nothing is silently overwritten. Waters and
ligands are in C<residues> alongside the polymer, which is why chain A above
has 206 residues to its 191-long SEQRES.

The name is not part of a residue's identity. One position is sometimes
modelled in two chemical states at once, written as complementary altloc
groups — 3zeu has ten methionines that are MSE in altlocs A and B and MET in
C and D, a selenomethionine that only went halfway in — and those are one
residue, not two. It takes the name written first, counts the records of both
states, and keeps the atoms that tell them apart, so an MSE/MET like that has
both an SE and an SD, each carrying the conformers of its own state.

=head3 Counting elements

Two tallies, the same shape: C<< $info-E<gt>{stats}{elements} >> is the whole structure
and C<< $info-E<gt>{chains}{$id}{elements} >> is one chain of it. Both count coordinate
records, which is what C<n_atoms> counts, so both add up:

 use List::Util 'sum0';

 my $info = structure_info('1a22.ent.pdb');

 $info->{stats}{elements};              # { C => 1946, O => 643, N => 507, S => 17 }
 $info->{chains}{A}{elements};          # { C =>  938, O => 301, N => 246, S =>  7 }
 $info->{chains}{B}{elements};          # { C => 1008, O => 342, N => 261, S => 10 }

 sum0(values %{ $info->{chains}{A}{elements} }) == $info->{chains}{A}{n_atoms};  # true
 sum0(values %{ $info->{stats}{elements} })     == $info->{stats}{n_atoms};      # true

Both are tallies of what came back, so the C<model>, C<hydrogens>, C<waters>,
C<hetatm> and C<chains> options are already in them: read an NMR ensemble with
the default C<< model =E<gt> 1 >> and you get one model's worth. With C<< model =E<gt> 'all' >>
each model's chains carry their own tally, under C<< $info-E<gt>{models}{$n}{chains} >>.

The keys are IUPAC symbols — C<Zn>, C<Se>, C<Cl>, C<Fe> — not the shouted spelling
the file uses. A PDB file writes the element in columns 77-78 in capitals, an
mmCIF C<type_symbol> is capitals as often as not, and an element worked out from
the atom name comes out of a table that is capitals throughout, so C<ZN> is what
all three roads arrive with and C<Zn> is what the periodic table calls it. The
correction runs on the symbol once, where it is settled, so the atom's own
C<element>, the chain tally and the structure tally cannot disagree:

 $info->{chains}{A}{residues}{202}{atoms}{ZN}{element};   # 'Zn'
 $info->{chains}{A}{elements}{Zn};                        # 1

The atom is still keyed C<ZN> in C<atoms> there, because that key is the atom's
I<name> out of columns 13-16, not its element.

Only the 118 named elements are corrected. A file whose element column holds
something that spells no element keeps it exactly as written — C<XX> stays C<XX>
rather than becoming a plausible-looking C<Xx> — so a field the module does not
recognise is visibly not an element rather than quietly dressed up as one.

Each count is an unsigned integer. It is counted up from zero and never down,
so there is no sign for it to carry.

=head3 Nothing points back up

A residue does not hold its chain and an atom does not hold its residue — the
names are there, C<< chain =E<gt> 'A' >> on the residue, but not the references. Parent
links would make the whole thing one reference cycle, and a cycle is a leak
that goes unnoticed until the ten-thousandth file.

=head3 Options

 model     => 1          which MODEL to build chains from; 'all' fills in
                         {models} as well.  Default 1, which is also the
                         right answer for a file with no MODEL records
 altloc    => 'first'    which alternate conformer's coordinates win;
                         'highest' takes the highest occupancy instead
 hydrogens => 1          keep hydrogens and deuteriums
 waters    => 1          keep waters
 hetatm    => 1          keep HETATM records
 atoms     => 1          build the atom hashes; 0 stops at the residue
                         level, which is much smaller and faster
 features  => 1          compute the physical properties, into
                         {features} and onto the chains, residues and
                         atoms; see structure_features below
 meta      => 1          parse the header records
 anisou    => 0          keep ANISOU lines
 chains    => ['A','B']  read only these chains
 format    => 'pdb'      skip the format detection: 'pdb' or 'mmcif'
                         ('cif', 'pdbx' and 'ent' name the same two)

Every option is checked. A misspelled one is fatal, because an ignored typo is
a wrong answer that arrives without a word: C<< hydrogen =E<gt> 0 >> that is quietly
dropped gives a structure with the hydrogens still in it and no clue why.

For a very large structure the options are the difference between a hash of
hashes that fits in memory and one that does not. The largest entry in PDBbind
v2020 is 2wy2: 33 MB, 64 models, 411,648 atom records.

 structure_info($f)                # model 1 only    50 MB    0.07 s
 structure_info($f, model => 'all')                 711 MB    1.4 s
 structure_info($f, model => 'all', atoms => 0)     418 MB    0.9 s

C<features> is the expensive one, and it is on by default because a structure's
surface, size and contacts are as much a part of what it is as its sequence, and
a caller who has to know to ask mostly does not. What it costs is measured, over
60 structures of PDBbind:

 structure_info($f, features => 0)         1.50 s   214,000 atoms/s
 structure_info($f)                       16.38 s    19,500 atoms/s   10.9x
 ... with interface => 0                  10.67 s    30,000 atoms/s    7.1x
 ... with sasa => 0                        4.01 s    80,000 atoms/s    2.7x

Nearly all of it is the solvent-accessible surface, at 960 sphere points per
atom; everything else together is 2.7 times the read. C<< interface =E<gt> 0 >> drops the
per-chain surfaces and takes a third off the whole thing.

C<< features =E<gt> 0 >> is what to reach for when reading a directory for its headers or
its sequences. C<< atoms =E<gt> 0 >> turns them off on its own — there are no coordinates
to compute from — rather than dying, so C<< structure_sequences($f, atoms =E<gt> 0) >>
still works.

Filtering happens in the C, before a hydrogen or a water has become a Perl
value, so C<< hydrogens =E<gt> 0 >> is cheaper than reading them and throwing them away.

What was filtered is still counted, so a structure knows how much of its file
it is. C<stats.n_atoms> is what came back and C<stats.total_atoms> is what the
file has -- every ATOM and HETATM record, every model, before any option had a
say -- and C<total_atoms == n_atoms + n_skipped> however the options were set.
2wy2 above, read with the default C<< model =E<gt> 1 >>, gives C<n_atoms> 6,432 and
C<total_atoms> 411,648.

=head2 pdb_info

 my $info = pdb_info($file, %options);

C<structure_info()> with the format settled in advance. Use it when the file is
known to be a PDB whatever it happens to be called.

Telling it wrongly reads no atoms rather than dying, which is what forcing a
format means; C<structure_info()> looks at the file and is the usual way in.

=head2 cif_info

 my $info = cif_info($file, %options);

The same for mmCIF/PDBx.

=head2 structure_info_string

 my $info = structure_info_string($text, %options);

The same, for a structure already in a string. A string has no name to go on,
so text that looks like nothing in particular is read as PDB; text that looks
like another format still gets a straight answer about it.

=head2 structure_atoms

 my $atoms = structure_atoms($info);
 my $atoms = structure_atoms($info, 'A');

Every atom as a flat array of hash references, in file order, each carrying
the C<chain>, C<resname>, C<resseq>, C<icode> and C<reskey> it came from. This is
the shape to hand to a distance calculation or to write out as a table; the
nested form is the shape to look things up in. The hashes are copies, so
writing to them does not scribble on the structure.

=head2 structure_residues

 my $residues = structure_residues($info);
 my $residues = structure_residues($info, 'A');

Every residue in file order. These are the same hash references that are in
the nested structure, not copies, so walking them and looking one up agree.

=head2 structure_ligands

 my $lig = structure_ligands($info);     # { 'NAG_A_301' => { ... } }

The heterogens that are neither water nor part of the polymer, keyed by
residue name, chain and number — which is what a binding-site table wants as
its row label.

=head2 is_single_ion

 is_single_ion($info->{chains}{E});     # 1     a chain that is one zinc
 is_single_ion($info, 'E');             # 1     the same, by chain id
 is_single_ion($info->{chains}{A});     # ''    a chain with a polymer in it

 my @polymers = grep { !is_single_ion($info, $_) } @{ $info->{chain_order} };

True when a chain holds exactly one residue. An ion is often numbered into the
chain it sits in — the zinc of a zinc finger is residue 202 of chain A — and
just as often given a chain of its own, which is a chain with one residue in it
and no sequence to read. This is for the second kind, so that a loop over
C<chain_order> can put them aside before it asks the rest for a sequence.

C<single> counts residues in the chain, and nothing else:

 a chain that is one CL                  # 1
 a chain that is one SO4, five atoms     # 1
 a chain that is one BF4, five atoms     # 1
 a chain of two zincs                    # ''
 a protein chain with a zinc in it       # ''

So the number of atoms in the residue does not come into it, and a sulphate and
a perchlorate answer the same. Neither does the residue's C<type>: that comes off
a table of names, and a table of names cannot be complete — SO4 is on the
module's list and BF4 is not, which is a fact about who wrote the list down and
not about the file. Counting residues asks the table nothing.

The residue is therefore not asked what it is, and a chain that is one sugar,
one buffer molecule, one water or one free amino acid reads true as well. In a
real file those are rare next to the ions and are the same nuisance to a caller
walking chains, but where the difference matters it is a lookup away:

 my $c = $info->{chains}{E};
 my $r = $c->{residues}{ $c->{residue_order}[0] };
 $r->{type} eq 'ion';        # ion, ligand, water, amino_acid, nucleotide
 $r->{resname};              # 'ZN'

C<res_type> is where those types come from, and C<< $c-E<gt>{type} eq 'water' >> is the
narrower question about a chain of nothing but waters.

The argument is either one chain — C<< $info-E<gt>{chains}{$id} >>, or a chain out of
C<< $info-E<gt>{models} >> — or the structure and a chain id, which is the same question
written the way C<chain_sequence()> takes it. Handing it the whole structure
without an id, or a residue, is fatal rather than false: all three are hash
references, and a wrong answer there would be taken at face value.

=head2 structure_sequences

 my $seq  = structure_sequences($info);          # { A => 'FPTIPLSRL...' }
 my $same = structure_sequences('1ubq.pdb');     # read on the spot
 my $fast = structure_sequences('1ubq.pdb', atoms => 0, meta => 0);

The observed single-letter sequence of every chain that has one.

The first argument is either the hash reference from C<structure_info()> or the
name of a file, which is read with the options given. C<< atoms =E<gt> 0 >> is worth
knowing about here, since a sequence needs the residues and not their
coordinates. Options belong with a file name; passing them alongside a
structure that is already parsed is an error, because there is nothing left
for them to change.

=head2 chain_sequence

 my $obs = chain_sequence($info, 'A');
 my $all = chain_sequence($info, 'A', 'seqres');

One chain's sequence: C<observed> is the residues that have coordinates,
C<seqres> is what SEQRES says was in the crystal.

=head2 structure_summary

 print structure_summary($info);

A paragraph a person can read: id, title, method, resolution, models, atom
counts, and a line per chain with its type, size, sequence and molecule. The
example at the top of this document is its output.

=head2 structure_features

 my $f = $info->{features};        # structure_info computed them already
 my $g = structure_features($info);   # the same hash, not a second walk

 $f->{sasa}{total};        # 17805.0   solvent-accessible surface, A^2
 $f->{sasa}{apolar};       # 8211.8    the carbon and sulphur part of it
 $f->{rg};                 # 22.85     radius of gyration, A
 $f->{rg_mass};            # 22.84     the same, weighted by mass
 $f->{mass};               # 41307.1   dalton
 $f->{hydropathy};         # -0.452    mean Kyte-Doolittle over the sequence
 $f->{aromatic_fraction};  # 0.1263    the F, W and Y share of it
 @{ $f->{pi_stacking} };   # the stacked pairs of aromatic rings
 @{ $f->{disulfides} };    # the SG-SG pairs close enough to be bonded

(1a22 again, the structure the summary at the top of this document is of.)

Everything the module can work out about a structure that is a number rather
than a name. One call, because all of it has to touch every atom and the
expensive part is walking the structure into coordinate arrays: asking for all
of it costs one walk.

C<structure_info()> makes that call on the way past and leaves the answer in
C<< $info-E<gt>{features} >>, so most callers never name this function at all. With no
options it is a lookup — it hands back what is already there. Name any option
and it walks the structure again with that option in force. C<< features =E<gt> 0 >> on
the read is how to skip the work; see C<structure_info>'s options above.

The whole-structure figures come back in the hash. What is per-atom, per-residue
or per-chain is written into C<$info> instead, where the atom, residue and chain
already are:

 $info->{chains}{A}{sasa};                          # 8450.6   the chain's surface
 $info->{chains}{A}{hydropathy};                    # -0.3144  and its mean hydropathy
 $info->{chains}{A}{aromatic_fraction};             # 0.1222
 $info->{chains}{A}{residues}{54}{sasa};            # 24.62    one residue's surface
 $info->{chains}{A}{residues}{54}{rsa};             # 0.103    ... as a fraction of its maximum
 $info->{chains}{A}{residues}{54}{atoms}{CZ}{sasa}; # one atom's

 $info->{chains}{A}{buried};                        # 1445.0   what it buries
 $info->{chains}{A}{residues}{54}{phi};             # -127.6   degrees
 $info->{chains}{A}{residues}{54}{psi};             #   20.5
 $info->{chains}{A}{residues}{54}{chi};             # [ -60.6, -82.9 ]
 $info->{chains}{A}{residues}{54}{ss};              # 'E'      or H G I B T S ' '
 $info->{chains}{A}{residues}{54}{ss_simple};       # 'E'      or H or C
 $info->{chains}{A}{residues}{54}{hse_up};          # 12       half-sphere exposure
 $info->{chains}{A}{residues}{54}{n_contacts};      # 13
 $info->{chains}{A}{residues}{23}{disulfide};       # [ { chain, residue, distance } ]

C<< store =E<gt> 0 >> turns that off and leaves C<$info> exactly as it was; the totals
still come back. Asking twice replaces what was stored rather than adding to it,
so a second call with a different probe radius leaves the second answer behind.

Three of the properties are I<only> per-residue — the torsion angles, the
half-sphere exposure and the secondary structure — so C<< store =E<gt> 0 >> does not
compute them at all rather than computing them and dropping them on the floor.
C<structure_info()> always stores, so this is only reachable by calling
C<structure_features()> yourself.

=head3 What comes back



=begin html

<table>
<thead>
<tr>
  <th>key</th>
  <th>what it is</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>n_atoms</code></td>
  <td>atoms the walk found: the ones in <code>$info</code>, after whatever <code>structure_info()</code> was told to leave out</td>
</tr>
<tr>
  <td><code>n_residues</code>, <code>n_chains</code></td>
  <td>and how they were grouped</td>
</tr>
<tr>
  <td><code>n_no_element</code></td>
  <td>atoms whose element field spells no element this module knows; they get a 2.0 A radius and no mass</td>
</tr>
<tr>
  <td><code>sasa</code></td>
  <td><code>total</code>, <code>apolar</code>, <code>polar</code>, and the <code>probe</code> and <code>points</code> used</td>
</tr>
<tr>
  <td><code>mass</code></td>
  <td>the sum of the atoms' standard atomic weights, in dalton</td>
</tr>
<tr>
  <td><code>rg</code></td>
  <td>radius of gyration about the centroid, in angstrom</td>
</tr>
<tr>
  <td><code>rg_mass</code></td>
  <td>the same, weighted by mass and taken about the centre of mass</td>
</tr>
<tr>
  <td><code>center</code>, <code>center_of_mass</code></td>
  <td><code>[x, y, z]</code>, in angstrom</td>
</tr>
<tr>
  <td><code>hydropathy</code></td>
  <td>the mean Kyte-Doolittle index over the protein chains' observed sequences</td>
</tr>
<tr>
  <td><code>aromatic_fraction</code>, <code>n_aromatic</code></td>
  <td>how much of that sequence is phenylalanine, tryptophan or tyrosine</td>
</tr>
<tr>
  <td><code>sequence_length</code></td>
  <td>how long the sequence those two are over is</td>
</tr>
<tr>
  <td><code>pi_stacking</code></td>
  <td>the arrayref <code>structure_pi_stacking()</code> returns</td>
</tr>
<tr>
  <td><code>disulfides</code></td>
  <td>the arrayref <code>structure_disulfides()</code> returns</td>
</tr>
<tr>
  <td><code>contacts</code></td>
  <td>the arrayref <code>structure_contacts()</code> returns</td>
</tr>
<tr>
  <td><code>hbonds</code></td>
  <td>the arrayref <code>structure_hbonds()</code> returns</td>
</tr>
<tr>
  <td><code>shape</code></td>
  <td><code>gyration_tensor</code>, <code>principal_moments</code>, <code>asphericity</code>, <code>acylindricity</code>, <code>anisotropy</code></td>
</tr>
</tbody>
</table>

=end html



C<rg>, C<center> and C<center_of_mass> are absent from a structure with no atoms in
it, and C<rg_mass> and C<center_of_mass> from one whose atoms have no mass between
them, because there is no such number rather than because it is zero. The same
goes for C<hydropathy> and C<aromatic_fraction> when there is no protein.

The surface is of the structure as C<$info> holds it. Reading with C<< waters =E<gt> 0 >>
and asking for the surface afterwards gives the surface of a protein with no
water in the way, which is a different — and usually more useful — number than
the surface of the file. C<< hydrogens =E<gt> 0 >> likewise: crystallographic structures
mostly have no hydrogens to begin with, and one that does will give a smaller
surface than its neighbours in the archive unless they are taken out.

=head3 Options



=begin html

<table>
<thead>
<tr>
  <th>option</th>
  <th>default</th>
  <th>what it does</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>sasa</code></td>
  <td>1</td>
  <td>compute the solvent-accessible surface</td>
</tr>
<tr>
  <td><code>pi_stacking</code></td>
  <td>1</td>
  <td>look for stacked aromatic rings</td>
</tr>
<tr>
  <td><code>disulfides</code></td>
  <td>1</td>
  <td>look for SG-SG pairs close enough to be bonded</td>
</tr>
<tr>
  <td><code>interface</code></td>
  <td>1</td>
  <td>also run the surface on each chain alone, for the buried area</td>
</tr>
<tr>
  <td><code>shape</code></td>
  <td>1</td>
  <td>the gyration tensor and the descriptors built from it</td>
</tr>
<tr>
  <td><code>dihedrals</code></td>
  <td>1</td>
  <td>phi, psi, omega and chi1-chi5, onto each residue</td>
</tr>
<tr>
  <td><code>contacts</code></td>
  <td>1</td>
  <td>which residues touch which</td>
</tr>
<tr>
  <td><code>exposure</code></td>
  <td>1</td>
  <td>half-sphere exposure, onto each amino acid residue</td>
</tr>
<tr>
  <td><code>hbonds</code></td>
  <td>1</td>
  <td>backbone hydrogen bonds, by Kabsch and Sander's energy</td>
</tr>
<tr>
  <td><code>secondary</code></td>
  <td>1</td>
  <td>secondary structure, onto each residue</td>
</tr>
<tr>
  <td><code>store</code></td>
  <td>1</td>
  <td>write the per-atom, per-residue and per-chain figures into <code>$info</code></td>
</tr>
<tr>
  <td><code>probe</code></td>
  <td>1.4</td>
  <td>solvent probe radius, angstrom</td>
</tr>
<tr>
  <td><code>points</code></td>
  <td>960</td>
  <td>sphere points per atom, from 1 to 10,000,000</td>
</tr>
<tr>
  <td><code>face_distance</code></td>
  <td>5.5</td>
  <td>face-to-face: the largest centroid separation, angstrom</td>
</tr>
<tr>
  <td><code>face_plane_min</code>, <code>face_plane_max</code></td>
  <td>0, 35</td>
  <td>... and the angle between the ring planes, degrees</td>
</tr>
<tr>
  <td><code>face_normal_min</code>, <code>face_normal_max</code></td>
  <td>0, 33</td>
  <td>... and between a ring's normal and the line joining the centroids</td>
</tr>
<tr>
  <td><code>edge_distance</code></td>
  <td>6.5</td>
  <td>edge-to-face: the largest centroid separation, angstrom</td>
</tr>
<tr>
  <td><code>edge_plane_min</code>, <code>edge_plane_max</code></td>
  <td>50, 90</td>
  <td>... and the angle between the ring planes</td>
</tr>
<tr>
  <td><code>edge_normal_min</code>, <code>edge_normal_max</code></td>
  <td>0, 30</td>
  <td>... and between a normal and the centroid line</td>
</tr>
<tr>
  <td><code>edge_radius</code></td>
  <td>1.5</td>
  <td>... and how close to a centroid the two planes' shared line must pass</td>
</tr>
<tr>
  <td><code>disulfide_distance</code></td>
  <td>3.0</td>
  <td>the largest SG-SG separation that counts as a bond, angstrom</td>
</tr>
<tr>
  <td><code>peptide_bond</code></td>
  <td>1.8</td>
  <td>the largest C-to-N separation that still joins two residues, angstrom</td>
</tr>
<tr>
  <td><code>contact_distance</code></td>
  <td>4.5</td>
  <td>the largest heavy-atom separation that counts as a contact, angstrom</td>
</tr>
</tbody>
</table>

=end html



A structure read with C<< atoms =E<gt> 0 >> has no coordinates to work from, and saying
so is more use than reporting no surface:

 my $info = structure_info('1ubq.pdb', atoms => 0);
 structure_features($info);
 # dies: this structure has no atom hashes to work from;
 #       read it again without atoms => 0

=head3 Shape, and what the chains bury

C<< $f-E<gt>{shape} >> is the gyration tensor and the three numbers built from its
eigenvalues, which are mdtraj's C<geometry/shape.py>:

 $f->{shape}{principal_moments};   # [ 74.1, 92.6, 132.1 ]  A^2, ascending
 $f->{shape}{asphericity};         # 41.84   how far from a sphere
 $f->{shape}{acylindricity};       # 21.10   how far from a cylinder
 $f->{shape}{anisotropy};          # 0.0233  0 for a sphere, 1 for a line

The three moments sum to C<rg> squared, which is the same sum read two ways.

C<< $f-E<gt>{sasa}{buried} >> is what the chains bury against each other — the surface
they have apart, less the surface they have together — and each chain carries
C<sasa_alone> and C<buried> of its own. A two-body interface is usually quoted as
half of the total, because the area is counted once on each side of it.

 $f->{sasa}{buried} / 2;              # 1471.2 A^2 of interface
 $info->{chains}{A}{buried};          # 1445.0 A^2, chain A's side of it

It costs one more surface calculation per chain, and a chain's calculation
touches only that chain's atoms, so all of them together cost about what the
first one cost rather than the number of chains times it. C<< interface =E<gt> 0 >> turns
it off.

=head3 Torsion angles

C<phi>, C<psi> and C<omega> on each residue, in degrees, and C<chi> as a list of
chi1 upwards — mdtraj's C<compute_phi>, C<compute_psi>, C<compute_omega> and
C<compute_chi1> through C<compute_chi5>.

An angle that would be measured across a chain break is not reported: the two
residues must be peptide-bonded first, at Biopython's C<PPBuilder> radius of
1.8 Å. mdtraj takes the residue before this one to be whichever came before it
in the file and computes a phi across whatever gap is there, which is a number
rather than an answer. Four collinear atoms get no angle either, for the same
reason.

C<omega> near zero is a cis peptide bond, which C<< $info-E<gt>{cispep} >> is the
depositor's own record of — two answers to one question, as with the disulfides.

=head3 Half-sphere exposure

C<hse_up> and C<hse_down> on each residue: the CA atoms of other residues within
12 Å, split by which side of the plane through this residue's CA they fall —
C<hse_up> towards the side chain. It says something the accessible surface does
not, because a residue can be buried and still have its side chain pointing into
a cavity. This is Biopython's C<Bio.PDB.HSExposure.HSExposureCB>, and the two
agree exactly.

Only the twenty standard amino acids, because Biopython's is built on
C<CaPPBuilder> with C<aa_only>, so a selenomethionine is invisible to it — it
neither gets a figure nor counts towards anybody else's. Glycine gets the
virtual CB Biopython builds for it.

=head2 structure_sasa

 my $s = structure_sasa($info);
 $s->{total};                                # 17805.0 A^2
 $info->{chains}{A}{residues}{54}{rsa};      # 0.103 -- mostly buried

 my $vdw = structure_sasa($info, probe => 0);   # the van der Waals surface
 my $fast = structure_sasa($info, points => 100);

The solvent-accessible surface and nothing else: the same calculation
C<structure_features()> runs, without the ring geometry. It returns the C<sasa>
hash and writes the per-atom, per-residue and per-chain surfaces into C<$info>
the same way. C<store>, C<probe> and C<points> are the options it takes.

A water molecule is a sphere about 1.4 A across, which is where the default
probe comes from; rolling a larger one gives a larger surface, because it cannot
reach into the dips. C<< probe =E<gt> 0 >> gives the van der Waals surface, which is the
smallest of them.

C<points> is accuracy against time. The area of an atom is 4*pi*r^2 times the
share of its sphere points no neighbour covers, so one point is worth about
0.13 A^2 for a carbon at 960 points and ten times that at 96; the whole surface
of a small protein moves by well under a percent between the two, and any single
atom can move by rather more.

=head3 Relative accessibility

C<rsa> is a residue's surface as a fraction of the most it could have — the
number the buried-or-exposed question is actually asked of, since 130 A^2 is
most of an alanine and a sliver of a tryptophan. It is on every amino acid
residue and on nothing else: the single-letter codes of the nucleotides are
amino acid codes too, and a guanine divided by glycine's maximum would be a
number rather than an answer.

A residue can come out above 1. The maxima are of a Gly-X-Gly tripeptide
stretched out, and a residue at the end of a chain with nothing next to it can
beat that.

C<apolar> is the part of the surface belonging to carbon and sulphur atoms and
C<polar> is everything else, which is the split Chothia made when he first added
a protein's buried surface up. Only the element symbol decides it, so a sulphur
in a sulphate counts as apolar; the per-atom figures are there for anyone who
wants a chemistry-aware split.

=head2 structure_pi_stacking

 for my $s (@{ structure_pi_stacking($info) }) {
     printf "%s %s%s %s and %s %s%s %s are %s stacked, %.2f A apart\n",
         $s->{chain1}, $s->{resname1}, $s->{residue1}, $s->{ring1},
         $s->{chain2}, $s->{resname2}, $s->{residue2}, $s->{ring2},
         $s->{type}, $s->{distance};
 }
 # A TRP5 6 and A HIS64 5 are face stacked, 3.70 A apart
 # A PHE66 6 and A PHE226 6 are edge stacked, 5.36 A apart

Every pair of aromatic rings in the structure that is stacked on the other, as
an arrayref of hashes. Two arrangements count, and they are the two ProLIF and
mdtraj look for: I<face>, two rings lying flat on each other, and I<edge>, one
ring pointing its edge at the other's face.



=begin html

<table>
<thead>
<tr>
  <th>key</th>
  <th>what it is</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>type</code></td>
  <td><code>face</code> or <code>edge</code></td>
</tr>
<tr>
  <td><code>chain1</code>, <code>residue1</code>, <code>resname1</code>, <code>ring1</code></td>
  <td>the first ring: its chain, the residue key it is keyed by in <code>$info</code>, the residue name, and <code>6</code> or <code>5</code> for the ring's size</td>
</tr>
<tr>
  <td><code>chain2</code>, <code>residue2</code>, <code>resname2</code>, <code>ring2</code></td>
  <td>the second</td>
</tr>
<tr>
  <td><code>distance</code></td>
  <td>between the two ring centroids, angstrom</td>
</tr>
<tr>
  <td><code>plane_angle</code></td>
  <td>between the two ring planes, degrees, folded into 0 to 90</td>
</tr>
<tr>
  <td><code>normal_angle1</code>, <code>normal_angle2</code></td>
  <td>between each ring's normal and the line joining the centroids, likewise</td>
</tr>
<tr>
  <td><code>intersect_distance</code></td>
  <td>edge stacks only: how far the line where the two planes meet passes from the nearer centroid</td>
</tr>
</tbody>
</table>

=end html



It writes nothing into C<$info>, and takes the eleven geometry options in the
table above and none of the others.

=head3 Which rings

The rings are the ones the format's own atom naming fixes: phenylalanine and
tyrosine's six-membered ring, tryptophan's six and five, histidine's five
(including the HID/HIE/HIP and HSD/HSE/HSP spellings a force field leaves
behind), the six- and five-membered rings of adenine and guanine, and the
six-membered ring of cytosine, thymine and uracil, in DNA and in RNA.

A ligand contributes none. Working out that a ligand has an aromatic ring means
perceiving its bonds, and this module never reads a CONECT record or guesses a
bond — so a stack between a drug and a tyrosine is not something it can find,
and it says so here rather than quietly finding nothing.

A ring missing any of its atoms is skipped, and the two rings of one tryptophan
or one purine are never paired with each other: they are fused into one aromatic
system, not two systems stacked.

=head3 On the face-to-face distance

C<face_distance> defaults to 5.5 B<angstrom>. mdtraj's C<pi_stacking()>, which
this is a translation of, has C<max_face_to_face_centroid_distance=5.5> in a
function whose other three distances are nanometres — 55 A, far enough that any
two aromatic rings in a small protein would qualify on distance alone. ProLIF,
whose geometry mdtraj's is taken from, has 5.5 A, and mdtraj's other three
distances are ProLIF's converted. So 5.5 A is what was meant, and it is what
this uses. C<< face_distance =E<gt> 55 >> gets the number mdtraj ships.

=head2 structure_contacts

 for my $c (@{ structure_contacts($info) }) {
     printf "%s%s - %s%s  %.2f A\n",
         $c->{chain1}, $c->{residue1}, $c->{chain2}, $c->{residue2}, $c->{distance};
 }

Which residues touch which: pairs whose closest heavy atoms are within
C<contact_distance>, with that distance. C<< $residue-E<gt>{n_contacts} >> counts them per
residue. Hydrogens are left out, which is what makes the number comparable
between a structure that has them and one that does not.

This is mdtraj's C<compute_contacts()> with its default C<closest-heavy> scheme.
mdtraj's C<all> pairs up residues in the same chain that are three or more apart
in it; this reports those and the neighbouring and cross-chain pairs too,
because a caller looking at a complex wants the interface and dropping it
silently would be strange. C<t/features.t> compares the subset mdtraj has an
opinion about, and finds the same distances.

=head2 structure_hbonds

 for my $b (@{ structure_hbonds($info) }) {
     printf "%s%s N-H ... O=C %s%s  %.2f kcal/mol\n",
         $b->{donor_chain}, $b->{donor_residue},
         $b->{acceptor_chain}, $b->{acceptor_residue}, $b->{energy};
 }

The backbone hydrogen bonds, by Kabsch and Sander's electrostatic definition:
a charge of -0.42 e on the carbonyl oxygen and +0.20 e on the amide hydrogen,
their opposites on the carbon and the nitrogen, and the four-way Coulomb sum as
the energy. Anything below -0.5 kcal/mol is a bond.

The amide hydrogen is not read from the file, it is placed — one angstrom from
N along the previous residue's C=O — which is what lets the definition be used
on a crystal structure that has no hydrogens in it. Proline has no amide
hydrogen and never donates. Each nitrogen keeps its two best acceptors.

This is mdtraj's C<kabsch_sander()>, and the two agree exactly: over 1A42, 1A22,
1AHW and 3AU6 they find the same bonds and the same energies to 4e-5 kcal/mol,
which is float32 rounding on mdtraj's side.

=head3 Secondary structure

The same bonds, read for the patterns they fall into, give every residue a
letter — the Kabsch–Sander dictionary:



=begin html

<table>
<thead>
<tr>
  <th>letter</th>
  <th>what it is</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>H</code></td>
  <td>alpha helix</td>
</tr>
<tr>
  <td><code>G</code></td>
  <td>3-10 helix</td>
</tr>
<tr>
  <td><code>I</code></td>
  <td>pi helix</td>
</tr>
<tr>
  <td><code>E</code></td>
  <td>extended strand</td>
</tr>
<tr>
  <td><code>B</code></td>
  <td>isolated beta bridge</td>
</tr>
<tr>
  <td><code>T</code></td>
  <td>hydrogen-bonded turn</td>
</tr>
<tr>
  <td><code>S</code></td>
  <td>bend</td>
</tr>
<tr>
  <td><code> </code></td>
  <td>none of them</td>
</tr>
</tbody>
</table>

=end html



C<< $residue-E<gt>{ss} >> is that letter and C<< $residue-E<gt>{ss_simple} >> is the three-state
reduction of it — C<H> for the three helices, C<E> for the two sheet letters, C<C>
for everything else. A residue with no backbone gets neither, because it is not
coil, it is not protein.

B<This is the one number in the module that does not reproduce its reference
exactly.> Against mdtraj's C<compute_dssp()> it agrees on 96.9% of residues
eight-state and 98.2% three-state, measured over 1A42, 1A22, 1AHW, 3AU6 and
1A4K — 3,314 residues. The residual is almost all beta sheet extension: mdtraj's
implementation joins ladders across bulges, a rule past the 1983 paper's core
definitions, and this does not. C<t/features.t> bounds the disagreement rather
than asserting equality, and says so where it does. If you need DSSP's exact
answer, run DSSP.

=head2 structure_disulfides

 for my $b (@{ structure_disulfides($info) }) {
     printf "%s%s - %s%s  %.2f A\n",
         $b->{chain1}, $b->{residue1}, $b->{chain2}, $b->{residue2}, $b->{distance};
 }
 # A23 - A88    2.01 A
 # A134 - A194  2.04 A

 $info->{chains}{A}{residues}{23}{disulfide};
 # [ { chain => 'A', residue => '88', distance => 2.011 } ]

The disulfide bonds the coordinates show: pairs of cysteines whose SG atoms are
close enough to be bonded. Each bond is also written onto both of its residues,
as a list naming the partner — a list because a cysteine that appears to hold
two bonds means something is wrong with the entry, and reporting both is more
use than dropping one.

The rule is a cysteine with an C<SG> and no C<HG>, paired with another under 3.0 Å.
The C<HG> test is what separates a cysteine whose thiol hydrogen was modelled — so
it is reduced, and holds no bond — from one that was not; a crystal structure
with no hydrogens has no C<HG> anywhere and every cysteine is a candidate, which
is right. C<disulfide_distance> moves the cutoff.

Only residues named C<CYS>. AMBER and CHARMM rename a bonded cysteine to C<CYX>,
and a structure that has been through a force field needs its residues named the
way the archive names them.

=head3 Against what the file says

C<< $info-E<gt>{ssbond} >> is the other answer: what the depositor wrote in an SSBOND
record, or in an mmCIF C<_struct_conn> row of type C<disulf>. Neither is the
authority, and comparing them is worth doing:

 my $key = sub {
     my ($c1, $r1, $c2, $r2) = @_;
     return join '|', sort "$c1/$r1", "$c2/$r2";
 };
 my %found = map { $key->(@{$_}{qw(chain1 residue1 chain2 residue2)}) => $_->{distance} }
             @{ structure_disulfides($info) };
 my %said  = map { $key->(@{$_}{qw(chain1 resseq1 chain2 resseq2)}) => $_->{length} }
             @{ $info->{ssbond} };
 my @undeclared = grep { !exists $said{$_}  } keys %found;
 my @unseen     = grep { !exists $found{$_} } keys %said;

Over a 60-entry spread of PDBbind the two agree on every entry that has
disulfides, and where both name the same bond their lengths agree to the two
decimals SSBOND is written in. Over a wider sweep they agree on 21 of 22; the
one that does not is 1A4K, a Fab that is in the file twice, whose SSBOND records
cover one copy and whose coordinates show both. A disagreement is a fact about
the entry, not about either answer.

=head2 aa3to1

 aa3to1('ALA');    # 'A'
 aa3to1('MSE');    # 'M'   selenomethionine is still a methionine
 aa3to1('HOH');    # ''    water is not an amino acid
 aa3to1('NAG');    # ''

The single-letter code of an amino acid, and the empty string for anything
that is not one. Leading and trailing blanks and case do not matter, because
the name usually arrives straight out of columns 18 to 20.

Modified residues map to the residue they were made from — C<MSE> to C<M>, C<SEP>
to C<S>, C<HYP> to C<P>, the D-amino acids to their L partners — because a
structure that soaked in selenomethionine has the same sequence as one that
did not, and a sequence with an C<X> every seventh position is no use to
anyone.

=head2 aa1to3

 aa1to3('A');      # 'ALA'
 aa1to3('X');      # 'UNK'
 aa1to3('B');      # 'ASX'   ASP or ASN, as the format spells it
 aa1to3('*');      # ''      not a single-letter code

C<aa3to1> backwards: the three-letter name a single-letter code stands for, and
the empty string for anything that is not one of the twenty-six. Blanks and
case do not matter, since the letter usually comes out of a sequence string
rather than out of a file.

Every letter of the alphabet has a name, because the ambiguity codes have one
of their own — C<B> is ASX, C<Z> is GLX, C<J> is XLE, C<X> is UNK. Going this way
there is only ever one answer: C<aa3to1> maps sixty-odd names onto C<C>, and
only CYS comes back.

Amino acids only, as the name says. C<aa1to3('A')> is ALA and not adenine, and
C<aa1to3('T')> is THR and not thymine — a caller who wants C<' DA'> already
knows the chain is DNA, and a function that guessed from a bare letter would be
wrong half the time.

It is in the XS rather than in Perl because it is both faster and smaller
there: the table is 104 bytes of read-only memory in the shared object, shared
between every process that loads the module, against 3,350 bytes of hash per
interpreter, and the lookup is one bounds check and one array index instead of
a hash lookup — about 4.5× the throughput measured a letter at a time.

=head2 res1

 res1('ALA');      # 'A'
 res1(' DA');      # 'A'   deoxyadenosine
 res1('PSU');      # 'U'   pseudouridine
 res1('HOH');      # ''

C<aa3to1()> widened to nucleotides, which is what building a sequence wants
when the chain might be DNA or RNA.

=head2 res_type

 res_type('ALA');  # 'amino_acid'
 res_type('DA');   # 'nucleotide'
 res_type('HOH');  # 'water'
 res_type('NAG');  # 'other'

What kind of residue a name is. C<other> covers ligands, ions and sugars;
C<structure_info()> narrows those to C<ligand> or C<ion> once it can see how many
atoms the residue has and what they are.

=head2 formats

 my @can = formats(); # ('mmcif', 'pdb')
 my $all = formats(); # every format known, supported or not

=head2 h

Prints a function's documentation to STDOUT. See I<Getting help> above.

=head1 Two residues that are not what they look like

Both of these were found by running the module over PDBbind v2020 and asking
where the sequence it read disagreed with SEQRES. Both are in the test suite.

B<A free base is not a nucleotide.> C<ADE>, C<CYT>, C<GUA>, C<THY> and C<URI> mean
one thing in a file written before 2007 — the nucleotides of a nucleic acid
chain — and another in a file written since: a free base sitting in an active
site as a ligand. The sugar tells them apart, since a nucleotide has a C<C1'>
and a free base has nothing but the base. Without that check the guanine bound
to 1czc reads as a nucleotide and turns up as a C<G> on the end of a
396-residue protein sequence.

B<A free amino acid is not part of the chain.> A HETATM residue with an amino
acid's name is a modified residue when it is numbered among the polymer — the
MSE that replaced a methionine belongs in the sequence — and a free amino acid
bound in a site when it is numbered out with the ligands, in which case it
does not. 3lms has a glycine at A501, two hundred residues past the end of a
chain whose SEQRES is 309 long. Those are flagged C<< free =E<gt> 1 >> and typed as
ligands.

Neither is a rule the format states; both are what the format means.

=head1 Files that keep their entry id in columns 73-80

An entry deposited before about 1996 carried its id and a line number in the
last eight columns of every record, and the archive still distributes those
files as they were deposited. Every field a reader takes to the end of the line
is wrong on one of them, and wrong in a way nothing downstream can see: the
SEQRES of a 140-residue chain comes back 162 long with an C<X> every thirteenth
place, the compound is the compound with C<1GDR   3> after it, C<HELIX> reports a
length of C<1GDR>, and columns 77-78 make 105 atoms of element C<1> — which also
stops C<< hydrogens =E<gt> 0 >> from finding any hydrogens, since it is the element that
says which atoms those are.

So the columns are read as columns. SEQRES takes 20-70 and no more; an element
field that is not letters is not an element and the atom name is used instead;
a charge field that is not a digit and a sign reads as the empty string a blank
one would have given; a C<HELIX> length that is not a number reads as empty; and
a text record whose columns 73-80 hold nothing but the entry id and a line
number is cut there — text that is not the entry id is left alone, so a title
that really does run to column 80 is not truncated.

C<COMPND> and C<SOURCE> predate the C<MOL_ID> convention in a file like this and
are free text: C<COMPND    GAMMA DELTA RESOLVASE>. There is no chain list in that
form because there was nothing to distinguish, so the entry is the one molecule
and every chain in it gets it, and C<< $info-E<gt>{compound}{1}{free_text} >> is 1 to say
the record was read that way rather than parsed into tokens.

C<t/data/pdb1gdr.ent> is one such file, a 1993 entry, and C<t/foreign.t> reads it.

=head1 Where the physical properties come from

None of the arithmetic in C<structure_features()> is this module's own. Each
piece is a translation of a published method as somebody else implemented it,
and the tests compare against those implementations rather than against what
this module currently does — C<t/features.t> reads what mdtraj answered for every
structure in C<t/data>, frozen into C<t/data/features.txt>, and re-runs mdtraj
where mdtraj is installed so the frozen answer cannot go stale.



=begin html

<table>
<thead>
<tr>
  <th>what</th>
  <th>from</th>
  <th>as implemented in</th>
</tr>
</thead>
<tbody>
<tr>
  <td>solvent-accessible surface</td>
  <td>Shrake, A; Rupley, J A (1973) <i>J Mol Biol</i> 79(2):351-71</td>
  <td>mdtraj 1.11's <code>mdtraj.geometry.shrake_rupley</code></td>
</tr>
<tr>
  <td>gyration tensor and shape</td>
  <td></td>
  <td>mdtraj's <code>geometry/shape.py</code></td>
</tr>
<tr>
  <td>torsion angles</td>
  <td></td>
  <td><code>mdtraj.compute_phi</code>, <code>compute_psi</code>, <code>compute_omega</code>, <code>compute_chi1</code>-<code>chi5</code></td>
</tr>
<tr>
  <td>residue contacts</td>
  <td></td>
  <td><code>mdtraj.compute_contacts</code>, <code>closest-heavy</code></td>
</tr>
<tr>
  <td>backbone hydrogen bonds</td>
  <td>Kabsch, W; Sander, C (1983) <i>Biopolymers</i> 22(12):2577-637</td>
  <td><code>mdtraj.geometry.kabsch_sander</code></td>
</tr>
<tr>
  <td>secondary structure</td>
  <td>the same paper</td>
  <td><code>mdtraj.compute_dssp</code> — the one place the agreement is not exact; see <code>structure_hbonds</code></td>
</tr>
<tr>
  <td>half-sphere exposure</td>
  <td>Hamelryck, T (2005) <i>Proteins</i> 59(1):38-48</td>
  <td>Biopython's <code>Bio.PDB.HSExposure.HSExposureCB</code></td>
</tr>
<tr>
  <td>the peptide-bond cutoff</td>
  <td></td>
  <td>Biopython's <code>Bio.PDB.Polypeptide.PPBuilder</code> <code>radius</code></td>
</tr>
<tr>
  <td>disulfides</td>
  <td></td>
  <td>mdtraj's <code>Topology.create_disulfide_bonds</code> rule</td>
</tr>
<tr>
  <td>van der Waals radii</td>
  <td>Bondi, A (1964) <i>J Phys Chem</i> 68:441, extended by Mantina, M <i>et al.</i> (2009) <i>J Phys Chem A</i> 113:5806, with Shannon, R D (1976) <i>Acta Cryst</i> A32:751 ionic radii for the ions that are always ionised</td>
  <td>mdtraj's <code>_ATOMIC_RADII</code></td>
</tr>
<tr>
  <td>atomic masses</td>
  <td></td>
  <td>mdtraj's <code>mdtraj/core/element.py</code></td>
</tr>
<tr>
  <td>pi-stacking geometry</td>
  <td>ProLIF's FaceToFace and EdgeToFace</td>
  <td><code>mdtraj.geometry.pi_stacking</code></td>
</tr>
<tr>
  <td>radius of gyration</td>
  <td></td>
  <td><code>mdtraj.geometry.compute_rg</code></td>
</tr>
<tr>
  <td>maximum accessible surface, for <code>rsa</code></td>
  <td>Tien, M Z <i>et al.</i> (2013) <i>PLoS ONE</i> 8(11):e80635, Table 1, the theoretical column</td>
  <td></td>
</tr>
<tr>
  <td>hydropathy</td>
  <td>Kyte, J; Doolittle, R F (1982) <i>J Mol Biol</i> 157(1):105-132</td>
  <td>Biopython's <code>Bio.SeqUtils.ProtParamData.kd</code> and <code>ProteinAnalysis.gravy()</code></td>
</tr>
<tr>
  <td>aromaticity</td>
  <td>Lobry, J R; Gautier, C (1994) <i>Nucleic Acids Res</i> 22(15):3174-3180</td>
  <td>Biopython's <code>ProteinAnalysis.aromaticity()</code></td>
</tr>
</tbody>
</table>

=end html



mdtraj works in nanometres and float32; this module works in angstrom and NV,
which is what the two file formats are written in and what the rest of the
module already returns. The formulae are the same ones, so the answers agree to
the width of a float32: run in float64, mdtraj's own Shrake-Rupley loop and this
one give the same surface for every atom of every structure in C<t/data> to nine
digits, and mdtraj as it ships differs on four atoms of 620, by one sphere point
each — a point sitting within a float32 ulp of a neighbouring atom's surface is
accessible at one width and covered at the other.

Three of these are deliberately not what the reference does, and each is argued
where it is written down. mdtraj computes a torsion angle, and places a Kabsch–
Sander amide hydrogen, between whichever residues are next to each other in the
file — across a chain break, where the chemistry never joined them; here the two
must be peptide-bonded first. mdtraj's C<Topology.create_disulfide_bonds()>
compares angstrom coordinates against a nanometre cutoff and so finds no
disulfide in any file; the rule it documents is the one implemented here. And
the face-to-face pi-stacking distance, above.

Two more are deliberately not mdtraj's: C<rg_mass>, and
C<compute_rg(traj, masses=m)> weights the distances by mass but still measures
them from the geometric centroid; C<rg_mass> measures from the centre of mass,
which is what the quantity means. C<rg> takes mdtraj's default of equal weights,
where the two centres are the same point and the two answers agree exactly.

The last is the secondary structure, which is the only number here that does not
reproduce its reference: 96.9% agreement eight-state, 98.2% three-state, over
3,314 residues. It is written up under C<structure_hbonds> above, and
C<t/features.t> bounds the disagreement instead of pretending there is none.

=head1 What is parsed in C, and why

The C side does one pass over the bytes. It splits ATOM/HETATM records into
their fields, marks where each residue begins and ends, sums what has to be
summed over every atom — the element tally, the bounding box, the B-factors,
each residue's centre — and groups every other record by record name for Perl
to take apart. Residue name lookup, three letters to one letter and the amino
acid/nucleotide/water question, is a switch on three packed bytes, and one
table serves C<aa3to1()>, C<res1()> and C<res_type()> so the three can never
disagree.

When atoms are wanted the parse builds the atom hashes itself, rather than
handing back columns for Perl to rebuild them from; building every atom twice
cost more than everything else in the read put together. When they are not
wanted — C<< atoms =E<gt> 0 >> — it builds none, and the Perl that follows walks
residues rather than atoms.

Everything else is Perl. The header records are irregular, they are a few
dozen lines per file rather than hundreds of thousands, and they are where the
next surprise will turn up; none of that is worth writing in C.

The mmCIF reader is a second pass written to the same division. A PDB file is
fixed columns and an mmCIF file is tag/value pairs and C<loop_> tables, so none
of the column arithmetic carries over and the tokenizer — quoting, semicolon
text fields, comments, the two spellings of null — is its own code. What it is
not is a second answer: it fills in the same output, the same column arrays and
residue boundaries and counts, so everything downstream of it, in C and in
Perl, is written once. C<_atom_site> goes through that path; every other
category is handed to Perl as tags and loops, which is the same place the line
between the two languages falls for PDB.

On 200 structures from PDBbind v2020, the parse runs at about 2.8 times the
speed of the same parse written in Perl. C<structure_info()> as a whole comes
out close to a pure-Perl reader that gathers the same statistics, while also
reading the headers, SEQRES, the gaps, the chain types and the ligands; see
C<benchmark.pl>, which measures all of it rather than asserting any of it.

=head1 Author

David E. Condon L<mailto:dec986@gmail.com>

=head1 COPYRIGHT AND LICENSE

This library is free software; you can redistribute it and/or modify it under
the same terms as Perl itself.
