#!/usr/bin/env perl
# ABSTRACT: Read a molecular structure file into a hash of hashes, sequences and all, using XS for the coordinate section
require 5.010;
use strict;
package Structure::Info;
our $VERSION = '0.01';
require XSLoader;
use autodie ':default';
use warnings FATAL => 'all';
use Carp 'croak';
use Exporter 'import';
use Scalar::Util 'reftype';
XSLoader::load('Structure::Info', $VERSION);

our @EXPORT_OK = qw(
	structure_info structure_info_string pdb_info
	structure_atoms structure_residues structure_ligands structure_sequences
	chain_sequence structure_summary
	aa3to1 res1 res_type formats h
);
our @EXPORT = @EXPORT_OK;

# ===========================================================================
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
# ===========================================================================
my %READER = (
	pdb => \&_read_pdb,
);

# formats a reader is not written for yet, kept here so that handing one over
# gets a straight answer rather than a puzzling parse of the wrong thing
my %NOT_YET = (
	mmcif => 'mmCIF/PDBx (.cif)',
	mol2  => 'Tripos MOL2 (.mol2)',
	sdf   => 'MDL SDF/MOL (.sdf, .mol)',
);

# ===========================================================================
# Options
#
# Anything not listed here is a typo, and a typo that is quietly ignored is a
# silent wrong answer later on -- pass 'hydrogen' for 'hydrogens' and you get
# a structure with the hydrogens still in it and no hint of why.
# ===========================================================================
my %DEFAULT = (
	model     => 1,       # which MODEL to build chains from; 'all' for every one
	altloc    => 'first', # 'first' or 'highest' occupancy, when an atom has alternates
	hydrogens => 1,       # keep hydrogen/deuterium atoms
	waters    => 1,       # keep waters
	hetatm    => 1,       # keep HETATM records (ligands, ions, modified residues)
	atoms     => 1,       # build the per-atom hashes; 0 stops at the residue level
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

# ===========================================================================
# Public entry points
# ===========================================================================

# structure_info($file, %opt) -- read a structure file into a hash of hashes.
sub structure_info {
	my ($file, %opt) = @_;
	croak 'structure_info: no file name given' unless defined $file && length $file;
	croak "structure_info: '$file' does not exist"  unless -e $file;
	croak "structure_info: '$file' is a directory"  if -d $file;
	my $o   = _options(\%opt, 'structure_info');
	my $fmt = defined $o->{format} ? lc $o->{format} : _detect_format($file);
	my $reader = $READER{$fmt}
		or croak "structure_info: cannot read '$file': "
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

# structure_info_string($text, %opt) -- the same, from a string already in hand.
sub structure_info_string {
	my ($text, %opt) = @_;
	croak 'structure_info_string: text is undefined' unless defined $text;
	my $o = _options(\%opt, 'structure_info_string');
	my $fmt = defined $o->{format} ? lc $o->{format} : _sniff_format($text);
	# a string has no name to go on, and the caller has already said this is a
	# structure, so text that looks like nothing in particular is read as PDB.
	# Text that looks like something else still gets a straight answer.
	$fmt = 'pdb' if $fmt eq 'unknown';
	croak "structure_info_string: cannot read this text: "
	    . (exists $NOT_YET{$fmt}
	       ? "$NOT_YET{$fmt} is not implemented yet"
	       : "no reader for format '$fmt'")
		unless $READER{$fmt};
	my $p = _parse_string($text, _xs_options($o));
	return _build_pdb(_retry_model($p, $o, \&_parse_string, $text), $o, undef);
}

# formats() -- the formats that can be read, in list context; in scalar
# context a hashref of every format known, including the unwritten ones.
sub formats {
	return wantarray
		? (sort keys %READER)
		: { (map { $_ => 'supported' } keys %READER), (map { $_ => "not implemented: $NOT_YET{$_}" } keys %NOT_YET) };
}

# ===========================================================================
# Views over a parsed structure
#
# These build what they return.  Nothing in the structure points back up at
# its parent -- a residue does not hold its chain, an atom does not hold its
# residue -- because a hash of hashes with parent links is a cycle, and a
# cycle is a leak that no one notices until the tenth thousand file.
# ===========================================================================

# structure_atoms($info, $chain?) -- every atom as a flat array of hashes,
# each one carrying the chain/residue it came from, in file order.
sub structure_atoms {
	my ($info, $chain) = @_;
	_check_info($info, 'structure_atoms');
	my @out;
	for my $cid (defined $chain ? ($chain) : @{ $info->{chain_order} }) {
		my $c = $info->{chains}{$cid} or croak "structure_atoms: no chain '$cid'";
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
		my $c = $info->{chains}{$cid} or croak "structure_residues: no chain '$cid'";
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

# structure_sequences($info) -- chain id => single-letter sequence, observed.
sub structure_sequences {
	my ($info) = @_;
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
	croak 'chain_sequence: no chain given' unless defined $chain;
	my $c = $info->{chains}{$chain} or croak "chain_sequence: no chain '$chain'";
	$which = 'observed' unless defined $which;
	croak "chain_sequence: which must be 'observed' or 'seqres', not '$which'"
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

# ===========================================================================
# Options and format detection
# ===========================================================================

sub _options {
	my ($opt, $who) = @_;
	for my $k (sort keys %$opt) {
		croak "$who: unknown option '$k'; known options are: " . join(', ', sort keys %DEFAULT)
			unless exists $DEFAULT{$k};
	}
	my %o = (%DEFAULT, %$opt);
	croak "$who: altloc must be 'first' or 'highest', not '$o{altloc}'"
		unless $o{altloc} eq 'first' || $o{altloc} eq 'highest';
	if (defined $o{chains}) {
		croak "$who: chains must be an array reference"
			unless (reftype($o{chains}) || '') eq 'ARRAY';
		croak "$who: chains is empty" unless @{ $o{chains} };
	}
	if (defined $o{model} && $o{model} ne 'all') {
		croak "$who: model must be a positive integer or 'all', not '$o{model}'"
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
			or croak "Structure::Info: '$file' is gzipped but IO::Uncompress::Gunzip is not installed: $@";
		my $z = IO::Uncompress::Gunzip->new($file)
			or croak "Structure::Info: cannot gunzip '$file': "
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

# ===========================================================================
# The PDB reader
# ===========================================================================

sub _read_pdb {
	my ($file, $o) = @_;
	my $p;
	if ($file =~ /\.gz\z/i) {
		my $text = _slurp_maybe_gzipped($file, undef);
		$p = _parse_string($text, _xs_options($o));
		$p = _retry_model($p, $o, \&_parse_string, $text);
	} else {
		$p = _parse_file($file, _xs_options($o));
		$p = _retry_model($p, $o, \&_parse_file, $file);
	}
	return _build_pdb($p, $o, $file);
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

sub _build_pdb {
	my ($p, $o, $file) = @_;
	my $meta = $p->{meta} || {};
	my $info = {
		file     => $file,
		format   => 'pdb',
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

	_parse_meta($info, $meta) if $o->{meta};

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
	return $info;
}

# --- coordinates -----------------------------------------------------------
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
				# the altloc option decides which supplies the coordinates
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

# --- per-chain sequence, type and gaps -------------------------------------
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
			# the coordinates quietly differs from the one in SEQRES
			my @gaps;
			for my $i (1 .. $#poly) {
				my ($a, $b) = @poly[ $i - 1, $i ];
				next unless defined $a->{number} && defined $b->{number};
				next if $b->{number} == $a->{number} + 1;
				next if $b->{number} == $a->{number}; # insertion code, not a gap
				push @gaps, { after => $a->{key}, before => $b->{key},
				              missing => $b->{number} - $a->{number} - 1 };
			}
			$c->{gaps}   = \@gaps;
			$c->{n_gaps} = scalar @gaps;
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

# ===========================================================================
# Header records
#
# Every one of these is a fixed-column record too, but there are only a few
# dozen lines of them in a file, they are irregular, and they are where a new
# quirk turns up every few hundred structures.  That is Perl's job, not C's.
# ===========================================================================
sub _parse_meta {
	my ($info, $meta) = @_;

	# defaults, so that callers can read a key without testing for it first
	$info->{$_} = undef for qw(title resolution r_work r_free);
	$info->{$_} = []    for qw(keywords experiment authors);
	$info->{$_} = {}    for qw(header compound source seqres het hetnam formul
	                           remarks dbref entity_of_chain cryst1 journal
	                           modres);
	$info->{$_} = []    for qw(helix sheet ssbond link cispep revdat site conect);

	if (my $h = $meta->{HEADER}) {
		my $l = $h->[0];
		$info->{header} = {
			classification => _c($l, 10, 40),
			deposit_date   => _c($l, 50, 9),
			id_code        => _c($l, 62, 4),
		};
	}
	# a record that is not in the file reads as undef, not as an empty string:
	# "there was no TITLE" and "the TITLE was blank" are different answers
	$info->{title}      = $meta->{TITLE} ? _joined($meta->{TITLE}, 10) : undef;
	$info->{caveat}     = _joined($meta->{CAVEAT}, 19) if $meta->{CAVEAT};
	$info->{keywords}   = [ grep { length } map { _t($_) } split /,/, _joined($meta->{KEYWDS}, 10) ];
	$info->{experiment} = [ grep { length } map { _t($_) } split /;/, _joined($meta->{EXPDTA}, 10) ];
	$info->{authors}    = [ grep { length } map { _t($_) } split /,/, _joined($meta->{AUTHOR}, 10) ];
	$info->{model_type} = _joined($meta->{MDLTYP}, 10) if $meta->{MDLTYP};
	$info->{obsolete}   = _joined($meta->{OBSLTE}, 10) if $meta->{OBSLTE};
	$info->{split}      = [ split ' ', _joined($meta->{SPLIT}, 10) ] if $meta->{SPLIT};

	$info->{compound} = _mol_records($meta->{COMPND}, 10);
	$info->{source}   = _mol_records($meta->{SOURCE}, 10);
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
		my $text = _c($l, 19);
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
	for my $l (@{ $meta->{SEQRES} || [] }) {
		my $cid = _c($l, 11, 1);
		my $n   = _c($l, 13, 4);
		my @res = split ' ', _c($l, 19);
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
		push @{ $info->{helix} }, {
			id            => _c($l, 11, 3),
			init_resname  => _c($l, 15, 3),
			init_chain    => _c($l, 19, 1),
			init_resseq   => _c($l, 21, 4),
			end_resname   => _c($l, 27, 3),
			end_chain     => _c($l, 31, 1),
			end_resseq    => _c($l, 33, 4),
			class         => _c($l, 38, 2),
			length        => _c($l, 71, 5),
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
	my ($lines, $from) = @_;
	return {} unless $lines;
	my $text = _joined($lines, $from);
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

# --- small helpers ---------------------------------------------------------

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
	my ($lines, $from) = @_;
	return '' unless $lines && @$lines;
	my $out = '';
	for my $l (@$lines) {
		my $t = _c($l, $from);
		next unless length $t;
		$out = _rejoin($out, $t);
	}
	return $out;
}

sub _rejoin {
	my ($have, $add) = @_;
	return $add unless defined $have && length $have;
	return $have . $add if $have =~ /-\z/;
	return "$have $add";
}

sub _check_info {
	my ($info, $who) = @_;
	croak "$who: expected the hash reference from structure_info()"
		unless defined $info && (reftype($info) || '') eq 'HASH' && exists $info->{chains};
	return 1;
}

# ===========================================================================
# Help
#
# h() prints a function's own documentation, in the spirit of R's ?function.
# The text is this file's POD, read at run time, so the help and the shipped
# documentation cannot drift apart.
# ===========================================================================
sub h {
	my ($what) = @_;
	my $name = _help_name($what);
	my $sec  = _pod_sections();
	if (defined $name && $sec->{$name}) {
		print STDOUT $sec->{$name};
		return $name;
	}
	print STDOUT "Structure::Info $VERSION\n\nDocumented functions:\n";
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
	my $n = "$what";      # a glob stringifies as *Structure::Info::res1
	$n =~ s/\A\*//;
	$n =~ s/\A.*:://;
	return $n;
}

sub _pod_sections {
	my %sec;
	open my $fh, '<', __FILE__;
	my ($in, $name);
	while (my $l = <$fh>) {
		if ($l =~ /\A=head2\s+(\S+)/) {
			$name = $1;
			$name =~ s/\(.*//;
			$in = 1;
			$sec{$name} = '';
			next;
		}
		if ($l =~ /\A=head[12]\b/ || $l =~ /\A=cut/) { $in = 0; next }
		$sec{$name} .= $l if $in && defined $name;
	}
	close $fh;
	s/\A\n+//, s/\n+\z/\n/ for values %sec;
	return \%sec;
}

1;

__END__

=head1 NAME

Structure::Info - read a molecular structure file into a hash of hashes

=head1 SYNOPSIS

    use Structure::Info;

    my $info = structure_info('1a22.ent.pdb');

    print $info->{id};                          # 1A22
    print $info->{resolution};                  # 2.6
    print $info->{chains}{A}{sequence};         # FPTIPLSRLFDNAMLRAHRLHQLAFD...
    print $info->{chains}{A}{residues}{54}{resname};              # LEU
    print $info->{chains}{A}{residues}{54}{atoms}{CA}{x};         # 23.918

    print structure_summary($info);

=head1 DESCRIPTION

One call reads a structure file and returns everything in it as a hash of
hashes: the header, the annotations, every chain, every residue, every atom,
and the single-letter sequence of each chain.  The coordinate section is
parsed in C, because it is the part that is millions of lines across a
directory of structures; the header records are parsed in Perl, because they
are irregular and there are only a few dozen of them per file.

The module is named for structures rather than for PDB because the shape of
what it returns has nothing to do with the file format it came from.  Today
it reads PDB; C<formats()> says what it reads at any moment.

=head2 structure_info

    my $info = structure_info($file, %options);

Reads C<$file> and returns a hash reference.  The format is worked out from
the file name, and from the first records in the file when the name gives
nothing away.  C<.gz> files are read as they are.

The returned hash:

    file        the path it was read from
    format      'pdb'
    id          the four-character PDB id, from HEADER or the file name
    title       TITLE, continuation lines joined
    header      { classification, deposit_date, id_code }
    experiment  [ 'X-RAY DIFFRACTION' ]
    resolution  2.6                         (from REMARK 2)
    r_work      0.196                       (from REMARK 3)
    r_free      0.278
    keywords    [ ... ]                     (KEYWDS)
    authors     [ ... ]                     (AUTHOR)
    journal     { auth => [...], titl, ref, pmid, doi }
    compound    { 1 => { molecule, chain => [...], engineered, ... } }
    source      { 1 => { organism_scientific, organism_taxid, ... } }
    seqres      { A => { sequence, residues => [...], length } }
    het         { NAG => { name, formula, instances => [...] } }
    helix sheet ssbond link cispep          [ { ... }, ... ]
    cryst1      { a, b, c, alpha, beta, gamma, sgroup, z }
    remarks     { 2 => [ lines ], 350 => [ lines ], ... }
    conect      [ [ serial, serial, ... ], ... ]
    n_models    how many MODEL records the file has
    model       which one the chains below were built from
    models      every model, when called with model => 'all'
    chains      { A => { ... } }
    chain_order [ 'A', 'B' ]                the order they appear in the file
    stats       { n_atoms, n_hetatm, n_hydrogens, elements => { C => 1234 },
                  bfactor => { min, max, mean }, bbox, center }

A chain:

    id            'A'
    type          protein | dna | rna | water | hetero | unknown
    sequence      single-letter, the residues that have coordinates
    seqres        single-letter, what SEQRES says was in the crystal
    n_residues n_polymer n_water n_ligand n_atoms n_hetatm
    n_missing     seqres length less the residues actually modelled
    gaps          [ { after, before, missing } ]   unmodelled stretches
    n_gaps
    first last    the first and last polymer residue keys
    molecule      from COMPND, e.g. 'GROWTH HORMONE'
    organism      from SOURCE
    dbref         [ { database, accession, db_id, ... } ]
    residue_order [ '1', '2', '3', '3A', ... ]     file order
    residues      { '54' => { ... } }

A residue, keyed by residue number with its insertion code appended, so that
C<100>, C<100A> and C<100B> are three separate keys:

    resname   'LEU'          number 54       icode ''
    one       'L'            the single-letter code, '' when there is none
    type      amino_acid | nucleotide | water | ligand | ion
    standard  1 for the twenty and the standard nucleotides
    modified  1 for a residue like MSE, which is still an M in the sequence
    hetero    1 when it was written as HETATM
    n_atoms   center [ x, y, z ]   b_mean
    atom_order [ 'N', 'CA', 'C', 'O', ... ]
    atoms      { CA => { ... } }

An atom:

    name 'CA'   serial 384   element 'C'   charge ''
    x y z       occupancy    bfactor       altloc   hetero
    altlocs   [ { altloc, x, y, z, occupancy, bfactor }, ... ]
              present only when the atom has alternate conformers; every
              conformer is listed, including the one chosen above

Options:

    model     => 1          which MODEL to build chains from; 'all' fills in
                            {models} as well.  Default 1, which is also the
                            right answer for a file with no MODEL records
    altloc    => 'first'    which alternate conformer's coordinates win;
                            'highest' picks the highest occupancy instead
    hydrogens => 1          keep hydrogens and deuteriums
    waters    => 1          keep waters
    hetatm    => 1          keep HETATM records
    atoms     => 1          build the atom hashes; 0 stops at the residue
                            level, which is much smaller and rather faster
    meta      => 1          parse the header records
    anisou    => 0          keep ANISOU lines
    chains    => ['A','B']  read only these chains
    format    => 'pdb'      skip format detection

Every option is checked; an unknown one is fatal, because an ignored typo is
a wrong answer that arrives silently.

For a very large structure the options are the difference between a hash of
hashes that fits comfortably in memory and one that does not.  The largest
entry in PDBbind v2020 is 2wy2: 33 MB, 64 models, 411,648 atom records.

    structure_info($f)                # model 1 only    50 MB    0.07 s
    structure_info($f, model => 'all')                 711 MB    1.4 s
    structure_info($f, model => 'all', atoms => 0)     418 MB    0.9 s

The filtering happens in the C, before a hydrogen or a water has become a Perl
value, so C<< hydrogens => 0 >> is cheaper than reading them and throwing them
away afterwards.

=head2 pdb_info

    my $info = pdb_info($file, %options);

C<structure_info()> with the format settled in advance.  Use it when the file
is known to be PDB whatever it is called.

=head2 structure_info_string

    my $info = structure_info_string($text, %options);

The same, for a structure already in a string.

=head2 structure_atoms

    my $atoms = structure_atoms($info);
    my $atoms = structure_atoms($info, 'A');

Every atom as a flat array of hash references, in file order, each one
carrying the C<chain>, C<resname>, C<resseq>, C<icode> and C<reskey> it came
from.  This is the shape to hand to a distance calculation or to write out as
a table; the nested form is the shape to look things up in.

=head2 structure_residues

    my $residues = structure_residues($info);
    my $residues = structure_residues($info, 'A');

Every residue in file order, as an array reference of the same hash
references that are in the nested structure.

=head2 structure_ligands

    my $lig = structure_ligands($info);     # { 'NAG_A_301' => { ... } }

The heterogens that are neither water nor part of the polymer, keyed by
residue name, chain and number.

=head2 structure_sequences

    my $seq = structure_sequences($info);   # { A => 'FPTIPLSRL...' }

The observed single-letter sequence of every chain that has one.

=head2 chain_sequence

    my $obs = chain_sequence($info, 'A');
    my $all = chain_sequence($info, 'A', 'seqres');

One chain's sequence: C<observed> is the residues with coordinates, C<seqres>
is what SEQRES says was in the crystal.  The two differ wherever the ends or
a loop went unmodelled, which is exactly the difference C<gaps> counts.

=head2 structure_summary

    print structure_summary($info);

A paragraph a person can read: id, title, method, resolution, models, atom
counts, and a line per chain with its type, size and sequence.

=head2 aa3to1

    aa3to1('ALA');    # 'A'
    aa3to1('MSE');    # 'M'   selenomethionine is still a methionine
    aa3to1('HOH');    # ''    water is not an amino acid
    aa3to1('NAG');    # ''

The single-letter code of an amino acid, and the empty string for anything
that is not one.  Modified residues map to the residue they were made from,
because a sequence with an X every seventh position is no use to anyone.
Leading and trailing blanks, and case, do not matter.

=head2 res1

    res1('ALA');      # 'A'
    res1(' DA');      # 'A'   deoxyadenosine
    res1('HOH');      # ''

C<aa3to1()> widened to nucleotides, which is what building a sequence wants
when the chain might be DNA or RNA.

=head2 res_type

    res_type('ALA');  # 'amino_acid'
    res_type('DA');   # 'nucleotide'
    res_type('HOH');  # 'water'
    res_type('NAG');  # 'other'

What kind of residue a name is.  C<other> covers ligands, ions and sugars;
C<structure_info()> narrows those to C<ligand> or C<ion> once it can see how
many atoms the residue has.

=head2 formats

    my @can = formats();          # ('pdb')
    my $all = formats();          # every format known, supported or not

=head2 h

    h('structure_info');    # by name
    h(*res_type);           # by name, unquoted
    h(\&aa3to1);            # by reference
    h();                    # the list of documented functions

Prints a function's documentation to STDOUT, in the spirit of R's C<?fn>.

Note that C<h(res_type)>, with no quotes and no sigil, cannot work: every
function here is exported, so Perl parses the bareword as a call to
C<res_type()> before C<h> ever sees it.

=head1 TWO RESIDUES THAT ARE NOT WHAT THEY LOOK LIKE

Both of these were found by running the module over PDBbind v2020 and asking
where the sequence it read disagreed with the SEQRES the file declared.  Both
are in the test suite.

B<A free base is not a nucleotide.>  C<ADE>, C<CYT>, C<GUA>, C<THY> and C<URI>
mean one thing in a file written before 2007 -- the nucleotides of a nucleic
acid chain -- and another in a file written since: a free base sitting in an
active site as a ligand.  The sugar tells them apart, since a nucleotide has a
C1' and a free base has nothing but the base.  Without that check the guanine
bound to 1czc reads as a nucleotide and turns up as a C<G> on the end of a
396-residue protein sequence.

B<A free amino acid is not part of the chain.>  A HETATM residue with an amino
acid's name is a modified residue when it is numbered among the polymer -- the
MSE that replaced a methionine belongs in the sequence -- and a free amino acid
bound in a site when it is numbered out with the ligands, in which case it does
not.  3lms has a glycine at A501, two hundred residues past the end of a chain
whose SEQRES is 309 long.  Those are flagged C<< free => 1 >> and typed as
ligands.

Neither is a rule the format states; both are what the format means.

=head1 WHAT IS PARSED IN C, AND WHY

The C side does one pass over the bytes.  It splits ATOM/HETATM records into
their fields, marks where each residue begins and ends, sums what has to be
summed over every atom -- the element tally, the bounding box, the B-factors,
each residue's centre -- and groups every other record by record name for Perl
to take apart.  Residue name lookup, three letters to one letter and the amino
acid/nucleotide/water question, is a switch on three packed bytes, and one
table serves C<aa3to1()>, C<res1()> and C<res_type()> so the three can never
disagree.

When atoms are wanted the parse builds the atom hashes itself, rather than
handing back columns for Perl to rebuild them from; building every atom twice
cost more than everything else in the read put together.  When they are not
wanted -- C<< atoms => 0 >> -- it builds none, and the Perl that follows walks
residues rather than atoms.

Everything else is Perl.  The header records are irregular, they are a few
dozen lines per file rather than hundreds of thousands, and they are where the
next surprise will turn up; none of that is worth writing in C.

On 200 structures from PDBbind v2020, the parse runs at about 2.8 times the
speed of the same parse written in Perl.  C<structure_info()> as a whole comes
out close to a pure-Perl reader that gathers the same statistics, while also
reading the headers, SEQRES, the gaps, the chain types and the ligands; see
C<benchmark.pl>, which measures all of it rather than asserting any of it.

=head1 AUTHOR

David E. Condon <dec986@gmail.com>

=head1 LICENSE

This library is free software; you can redistribute it and/or modify it under
the same terms as Perl itself.

=cut
