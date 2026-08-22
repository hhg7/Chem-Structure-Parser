#!/usr/bin/env perl
# Against another reader.
#
# t/real.t checks the C against a second reader written in Perl, which catches a
# column read wrongly but cannot catch a column read wrongly by both, because
# one person wrote both.  This checks it against a reader nobody here wrote:
# gemmi, which reads PDB and mmCIF in C++, is what the crystallography
# programs use, and has its own thirty years of quirks to answer for.
#
# Every atom of model 1 is compared -- chain, residue number, insertion code,
# atom name, alternate location, and the three coordinates -- as a multiset, so
# a file the two disagree about says which atoms and not merely how many.
#
# Two differences are expected and are not failures, because both readings are
# defensible and this module documents which it takes:
#
#   * a residue modelled in two chemical states at once (3jqh writes residue 1
#     as PRO in altloc A and SER in altloc B) is one residue here and two in
#     gemmi, so the residue name is not part of what is compared.  What is
#     compared is that the name this module chose is one of the names gemmi
#     found at that number.
#
#   * a file with no CRYST1 has no cell here and a 1x1x1 cell in gemmi, which
#     is gemmi's default rather than the file's answer.
#
# It runs when python3 can import gemmi -- 'pip install gemmi' -- and skips
# when it cannot, so the distribution builds on a machine without it.  The
# structures are t/data, plus a spread of STRUCTURE_INFO_TEST_DIR and
# STRUCTURE_INFO_TEST_CIF_DIR when those are set, the same directories the rest
# of the suite reads.
require 5.010;
use strict;
use warnings FATAL => 'all';
use Cwd 'abs_path';
use File::Basename 'dirname';
use File::Spec;
use File::Temp 'tempdir';
use Chem::Structure::Parser;
use Test::More;

my $py = $ENV{STRUCTURE_INFO_PYTHON} || 'python3';
my $has = do {
	no warnings 'exec';
	system("$py -c 'import gemmi' >/dev/null 2>&1") == 0;
};
plan skip_all => "$py cannot import gemmi" unless $has;

my $data = dirname(abs_path(__FILE__)) . '/data';
my $tmp  = tempdir(CLEANUP => 1);

# --- the other reader ------------------------------------------------------
#
# One line per atom of model 1, sorted, so that the comparison does not depend
# on either reader's idea of the order.  The residue name is on the line too,
# but as its own field, for the microheterogeneity check below.
my $script = "$tmp/dump.py";
open my $out, '>', $script or die "$script: $!";
print $out <<'PY';
import sys
import gemmi

st = gemmi.read_structure(sys.argv[1])
rows = []
if len(st):
    for chain in st[0]:
        for res in chain:
            for at in res:
                alt = '' if at.altloc in ('', '\x00') else at.altloc
                rows.append('%s|%d|%s|%s|%s|%.3f|%.3f|%.3f|%s' % (
                    chain.name, res.seqid.num, res.seqid.icode.strip(),
                    at.name, alt, at.pos.x, at.pos.y, at.pos.z, res.name))
print('\n'.join(sorted(rows)))
print('#models %d' % len(st))
print('#resolution %s' % (st.resolution or ''))
print('#spacegroup %s' % (st.spacegroup_hm or ''))
PY
close $out;

sub gemmi_read {
	my ($file) = @_;
	# gemmi refuses a file it cannot make sense of, loudly, and a file it
	# refuses is a fact about the file rather than a failure of this test: its
	# complaint goes to /dev/null and the caller is told there is nothing to
	# compare against.  t/data/pdb1gdr.ent is one -- gemmi will not have the
	# '90' that the entry id leaves in the charge columns.
	open my $fh, '-|' or do {
		open STDERR, '>', File::Spec->devnull or exit 1;
		exec $py, $script, $file or exit 1;
	};
	my @lines = <$fh>;
	close $fh;
	return undef if $?;
	chomp @lines;
	my %info = (rows => []);
	for my $l (@lines) {
		if ($l =~ /\A\#(\w+) ?(.*)\z/) { $info{$1} = $2; next }
		push @{ $info{rows} }, $l if length $l;
	}
	return \%info;
}

# --- this one --------------------------------------------------------------
sub ours {
	my ($info) = @_;
	my @rows;
	for my $cid (@{ $info->{chain_order} }) {
		my $c = $info->{chains}{$cid};
		for my $rk (@{ $c->{residue_order} }) {
			my $r = $c->{residues}{$rk};
			for my $an (@{ $r->{atom_order} }) {
				my $a = $r->{atoms}{$an};
				# every record the atom stands for: the conformer that supplied
				# the coordinates, and any other one on its list
				my @conf = ({ map { $_ => $a->{$_} } qw(altloc serial x y z) });
				my %seen = ($a->{serial} => 1);
				for my $cf (@{ $a->{altlocs} || [] }) {
					push @conf, $cf unless $seen{ $cf->{serial} }++;
				}
				for my $cf (@conf) {
					next unless defined $cf->{x};
					push @rows, sprintf '%s|%s|%s|%s|%s|%.3f|%.3f|%.3f|%s',
						$cid, $r->{number}, $r->{icode}, $an,
						$cf->{altloc}, $cf->{x}, $cf->{y}, $cf->{z}, $r->{resname};
				}
			}
		}
	}
	return [ sort @rows ];
}

# the fields both readers must agree about, which is everything on the line
# except the residue name
sub key { my $r = shift; $r =~ s/\|[^|]*\z//; return $r }

# enough of a difference to see what it is, and not a screenful of it
sub _few { return @_ > 5 ? (@_[0 .. 4], sprintf('... and %d more', @_ - 5)) : @_ }

sub compare {
	my ($file) = @_;
	my $them = gemmi_read($file);
	if (!$them) {
		note("gemmi would not read $file; skipped");
		return;
	}
	my $info = structure_info($file);
	my $us   = ours($info);
	my $name = $file;
	$name =~ s{.*/}{};
	# gemmi reads the atom name out of _atom_site.label_atom_id, so an mmCIF
	# written with the auth_* names alone -- which t/data/quirks.cif is, on
	# purpose -- gives it no atoms at all.  Nothing to compare against is not a
	# disagreement about anything.
	if (!@{ $them->{rows} } && @$us) {
		note("gemmi found no atoms in $name; skipped");
		return;
	}

	my (%ours, %theirs, %named);
	$ours{ key($_) }++   for @$us;
	for my $r (@{ $them->{rows} }) {
		$theirs{ key($r) }++;
		my ($cid, $num, $ic, undef, undef, undef, undef, undef, $rn) = split /\|/, $r, 9;
		$named{"$cid|$num|$ic"}{$rn} = 1;
	}

	is(scalar @$us, scalar @{ $them->{rows} }, "$name: the same number of atoms");
	my @only_ours   = grep { ($theirs{$_} || 0) != $ours{$_} } sort keys %ours;
	my @only_theirs = grep { ($ours{$_}   || 0) != $theirs{$_} } sort keys %theirs;
	ok(!@only_ours && !@only_theirs, "$name: every atom agrees") or do {
		diag("  here and not in gemmi: $_") for _few(@only_ours);
		diag("  in gemmi and not here: $_") for _few(@only_theirs);
	};

	# the residue name: one residue in two chemical states is one residue here
	# and two there, so the test is that the name chosen is one of the names
	my @wrong_name;
	for my $cid (@{ $info->{chain_order} }) {
		my $c = $info->{chains}{$cid};
		for my $rk (@{ $c->{residue_order} }) {
			my $r     = $c->{residues}{$rk};
			my $names = $named{"$cid|$r->{number}|$r->{icode}"} or next;
			push @wrong_name, "$cid $rk $r->{resname} vs " . join('/', sort keys %$names)
				unless $names->{ $r->{resname} };
		}
	}
	ok(!@wrong_name, "$name: every residue name is one gemmi found there")
		or do { diag("  $_") for _few(@wrong_name) };

	if (length $them->{resolution} && defined $info->{resolution}) {
		# REMARK 2 writes the resolution to two decimals, so half of the last
		# place is as far apart as two readings of the same field can be; more
		# than that is a different number rather than a rounding.
		ok(abs($info->{resolution} - $them->{resolution}) < 0.005,
			"$name: the same resolution ($info->{resolution})");
	}
	if (length($them->{spacegroup} || '') && length($info->{cryst1}{sgroup} || '')) {
		(my $a = $info->{cryst1}{sgroup}) =~ s/\s+//g;
		(my $b = $them->{spacegroup})     =~ s/\s+//g;
		is($a, $b, "$name: the same space group");
	}
}

# --- what to read ----------------------------------------------------------
my @files = grep { -f } map { "$data/$_" } qw(mini.pdb mini.cif nmr.pdb nmr.cif
                                              quirks.cif pdb1gdr.ent);

sub spread {
	my ($dir, $want, $re) = @_;
	return () unless defined $dir && -d $dir;
	opendir(my $dh, $dir) or return ();
	my @all = sort grep { /$re/ } readdir $dh;
	closedir $dh;
	return () unless @all;
	my $step = @all > $want ? int(@all / $want) : 1;
	my @f = map { "$dir/$all[$_]" } grep { $_ % $step == 0 } 0 .. $#all;
	return @f > $want ? @f[0 .. $want - 1] : @f;
}

my $want = $ENV{STRUCTURE_INFO_TEST_ALL} ? 1e6 : 15;
push @files, spread($ENV{STRUCTURE_INFO_TEST_DIR}, $want, qr/\.(pdb|ent)(\.gz)?\z/);
push @files, spread($ENV{STRUCTURE_INFO_TEST_CIF_DIR}, $want, qr/\.(cif|mmcif)(\.gz)?\z/);

diag(sprintf 'comparing %d structures against gemmi', scalar @files);
compare($_) for @files;

done_testing();
