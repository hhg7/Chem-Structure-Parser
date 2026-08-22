use strict;
use warnings;
use ExtUtils::MakeMaker;
use Config;

# -std=gnu99 for the declarations-after-statements and // comments the XS uses;
# only when the compiler is one that understands the flag.
my $std = ($Config{gccversion} || $Config{cc} =~ /\b(?:g?cc|clang)\b/) ? ' -std=gnu99' : '';

# MakeMaker's directory scan classifies .pm, .pl and .pod alike, so every
# author-only helper in the distribution root -- test.all.perls.pl,
# structure.info.pl, benchmark.pl, and read.me.pod, which md2pod.pl writes --
# landed in TO_INST_PM next to lib/Chem/Structure/Parser.pm and was copied into
# $(INST_LIB)/Chem/Structure/ by `make install` (read.me.pod also got a
# Chem::Structure::read.me man page of its own).  They are not part of the
# module (and PruneFiles keeps them out of the release), so drop them here:
# libscan is the documented hook for rejecting a file the scan picked up.  Only
# root-level names are filtered; anything with a directory in its path falls
# through to the default, which keeps its README.pod and version-control
# filtering.
{
	package MY;
	sub libscan {
		my ($self, $path) = @_;
		return '' if $path =~ /\.(?:pl|pod)\z/ && $path !~ m{/};
		return $self->SUPER::libscan($path);
	}
}

WriteMakefile(
	NAME             => 'Chem::Structure::Parser',
	VERSION_FROM     => 'lib/Chem/Structure/Parser.pm',
	ABSTRACT_FROM    => 'lib/Chem/Structure/Parser.pm',
	AUTHOR           => 'David E. Condon <dec986@gmail.com>',
	LICENSE          => 'perl',
	MIN_PERL_VERSION => '5.010',
	CONFIGURE_REQUIRES => {
		'ExtUtils::MakeMaker' => 6.52,
	},
	PREREQ_PM => {
		'Carp'         => 0,
		'Exporter'     => 0,
		'Scalar::Util' => 0,
		'XSLoader'     => 0,
		'autodie'      => 0,
	},
	TEST_REQUIRES => {
		'Test::More'      => 0,
		'Test::Exception' => 0,
	},
	CCFLAGS => $Config{ccflags} . $std,
	OBJECT  => '$(O_FILES)',
	test    => { TESTS => 't/*.t' },
	clean   => { FILES => 'Chem-Structure-Parser-*' },
	META_MERGE => {
		'meta-spec' => { version => 2 },
		resources   => {
			repository => {
				type => 'git',
				url  => 'https://github.com/hhg7/Chem-Structure-Parser.git',
				web  => 'https://github.com/hhg7/Chem-Structure-Parser',
			},
			bugtracker => { web => 'https://github.com/hhg7/Chem-Structure-Parser/issues' },
		},
	},
);
