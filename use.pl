#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Chem::Structure::Parser;
#use Util;
#qw(json_file_to_ref ref_to_json_file);
#~/.local/share/libgedit-gtksourceview-300/language-specs/perl.lang

my $f = structure_info('157d.pdb');
p $f;
