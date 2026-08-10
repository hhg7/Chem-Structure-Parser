use strict;
use warnings;
use ExtUtils::MakeMaker;
use Config;

# -std=gnu99 for the declarations-after-statements and // comments the XS uses;
# only when the compiler is one that understands the flag.
my $std = ($Config{gccversion} || $Config{cc} =~ /\b(?:g?cc|clang)\b/) ? ' -std=gnu99' : '';

WriteMakefile(
	NAME             => 'Structure::Info',
	VERSION_FROM     => 'lib/Structure/Info.pm',
	ABSTRACT_FROM    => 'lib/Structure/Info.pm',
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
	clean   => { FILES => 'Structure-Info-*' },
	META_MERGE => {
		'meta-spec' => { version => 2 },
		resources   => {
			repository => {
				type => 'git',
				url  => 'https://github.com/hhg7/pdb.info.git',
				web  => 'https://github.com/hhg7/pdb.info',
			},
			bugtracker => { web => 'https://github.com/hhg7/pdb.info/issues' },
		},
	},
);
