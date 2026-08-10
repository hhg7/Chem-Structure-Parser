#!/usr/bin/env perl
require 5.010;
use strict;
use warnings FATAL => 'all';
use Test::More;

BEGIN { use_ok('Structure::Info') or BAIL_OUT('Structure::Info will not load') }

ok(defined $Structure::Info::VERSION, "VERSION is set ($Structure::Info::VERSION)");

# the XS half has to be there: without it every other test is testing nothing
ok(defined &Structure::Info::_parse_file,   'XS _parse_file is bootstrapped');
ok(defined &Structure::Info::_parse_string, 'XS _parse_string is bootstrapped');

for my $f (qw(
	structure_info structure_info_string pdb_info
	structure_atoms structure_residues structure_ligands structure_sequences
	chain_sequence structure_summary aa3to1 res1 res_type formats h
)) {
	can_ok('Structure::Info', $f);
	ok(defined &{"main::$f"}, "$f is exported into the caller");
}

is_deeply([ formats() ], ['pdb'], 'formats() lists what can be read');
my $all = formats();
is(ref $all, 'HASH', 'formats() in scalar context is a hashref');
like($all->{mmcif}, qr/not implemented/, 'formats() names the formats not written yet');

done_testing();
