#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use XS::Check;

# Run XS::Check over Parser.xs, then the house rules below; print one line per
# finding and exit 1 if there were any.  Adapted from ~/Scripts/stats/xs.check.pl,
# which is where the history below happened.
#
# XS::Check reports through warn() unless given a reporter, and Devel::Confess
# (which this script used to load) turned each of those into a five-line stack
# trace: 196 findings came out as 1177 lines, and a real one -- merge()'s
# segfault on a suffixes array with a hole -- sat in the middle of them
# unnoticed.  Two of its checks also need correcting before their findings mean
# anything here; see decl_type() and %svpv_audited below.

open my $xs, '<', 'Parser.xs';
my @src = (undef, <$xs>);	# $src[$n] is line $n
close $xs;

my @found;	# [line, message]
my $check = XS::Check->new(reporter => sub {
	my %r = @_;
	push @found, [ $r{line}, $r{message} ];
});
$check->check_file('Parser.xs');

# A source line as a key: indentation and any trailing // comment removed, so
# that re-indenting or re-wording the comment keeps an audit, and changing the
# code itself does not.
sub line_key {
	my ($n) = @_;
	(my $l = $src[$n]) =~ s{\s*//.*$}{};
	$l =~ s/^\s+|\s+$//g;
	return $l;
}

# The type of the declaration of $var nearest above the call on line $n, or
# undef.  XS::Check keeps one type per variable *name* for the whole file --
# whichever declaration it read last -- so in a file with a `len` or an `l` in
# most functions its "not a constant type" and "is not a STRLEN variable" were
# reporting some other function's variable: 43 findings, of which two were
# real.  The nearest declaration above the use is the one in scope, unless a
# closer one sits in a sibling block that has already closed, which this does
# not try to tell apart.  Statements are split at ; ( and , so that one line of
# `STRLEN l; const char *s = SvPV(sv, l);`, and a parameter list, are each
# seen as their parts.
my $kw = qr/^(?:return|else|case|goto|sizeof|do|if|while|for|switch)$/;
sub decl_type {
	my ($var, $n) = @_;
	for (my $i = $n; $i >= 1; $i--) {
		my $text = $src[$i];
		$text =~ s{//.*$}{};
		$text = substr($text, 0, index($text, 'SvPV')) if $i == $n && index($text, 'SvPV') >= 0;
		for my $part (reverse split /[;(,]/, $text) {
			next unless $part =~ /^\s*((?:\w+\s+|\w+\s*\*+\s*)+?)\Q$var\E\s*(?:[=\[)]|$)/;
			my $type = $1;
			my ($first) = $type =~ /^(\w+)/;
			next if $first =~ $kw;
			$type =~ s/\s+$//;
			return $type;
		}
	}
	return undef;
}

# SvPV() calls whose use of the buffer has been read and found right for the
# encoding question XS::Check raises ("Specify either SvPVbyte or SvPVutf8").
# Neither of those is a safe default: SvPVbyte croaks on a wide character and
# SvPVutf8 upgrades the caller's own SV in place.  What matters is whether the
# code goes on to honour SvUTF8 -- or has no need to -- and that is what each
# reason says.  Keyed by line_key(), so editing an audited line brings its
# finding back until it is looked at again.
#
# None is audited here yet: every SvPV in this Parser.xs is SvPV_const,
# SvPV_nolen or SvPVX, which XS::Check does not ask about.  The table that came
# over from Stats::LikeR named that module's lines, and was emptied rather than
# left to match an unrelated line here by accident.
my %svpv_audited = ();

my @report;
for my $f (@found) {
	my ($n, $msg) = @$f;
	if ($msg =~ /^(\S+) not a constant type$/) {
		my $type = decl_type($1, $n);
		push @report, [$n, defined $type ? "$1 is declared '$type', not const" : "$1: no declaration found above the call"]
			unless defined $type && $type =~ /\bconst\b/;
	} elsif ($msg =~ /^(\S+) is not a STRLEN variable/) {
		my $type = decl_type($1, $n);
		push @report, [$n, defined $type ? "$1 is declared '$type', not STRLEN" : "$1: no declaration found above the call"]
			unless defined $type && $type =~ /\bSTRLEN\b/;
	} elsif ($msg =~ /^Specify either SvPVbyte or SvPVutf8/) {
		push @report, [$n, $msg] unless exists $svpv_audited{ line_key($n) };
	} elsif ($msg eq q{Remove the 'Perl_' prefix from Perl_sv_cmp}
	         && $src[$n] =~ /\bsortsv\w*\s*\(.*\bPerl_sv_cmp\s*\)/
	         && $src[$n] !~ /\bPerl_sv_cmp\s*\(/) {
		# sortsv() takes a function pointer, and sv_cmp is a function-like
		# macro, so it has no address: the Perl_ name is the only spelling
		next;
	} else {
		push @report, [$n, $msg];
	}
}

# Perl_isnan/Perl_isinf/Perl_isfinite must never be called from Parser.xs: on
# every perl before 5.22 they can expand to perl.h's Perl_fp_class() block,
# which is written with an empty parameter list and compares against FP_CLASS_*
# names <ieeefp.h> does not define.  Configure only reaches that block where it
# fails to find isinf(), so the breakage is invisible here and fatal on
# illumos/Solaris -- it is what stopped Stats::LikeR 0.298 building there.
# Parser.xs needs none of them today (grid_build() tests finiteness as x - x ==
# 0); anything that does gets a width-correct wrapper, as CLAUDE.md says under
# "Numbers must be exact".  Comments may still name the macros; only a call is
# an error, so comments are blanked before the match -- string literals are
# not, and none here names one.
my $in_block = 0;	# 1 = inside a /* */ comment carried over from an earlier line
for my $n (1 .. $#src) {
	my $code = '';
	my $l = $src[$n];
	while (length $l) {
		if ($in_block) {
			my $e = index($l, '*/');
			last if $e < 0;
			$l = substr($l, $e + 2);
			$in_block = 0;
		} elsif ($l =~ m{^(.*?)(/\*|//)}s) {
			$code .= $1;
			last if $2 eq '//';
			$l = substr($l, length($1) + 2);
			$in_block = 1;
		} else {
			$code .= $l;
			last;
		}
	}
	push @report, [$n, 'Perl_is*() must not be called -- see "Numbers must be exact" in CLAUDE.md']
		if $code =~ m/\bPerl_is(?:nan|inf|finite)\s*\(/;
}

for my $r (sort { $a->[0] <=> $b->[0] } @report) {
	say "Parser.xs:$r->[0]: $r->[1]";
}
say scalar(@report), ' finding', (@report == 1 ? '' : 's'), ' in Parser.xs';
exit(@report ? 1 : 0);
