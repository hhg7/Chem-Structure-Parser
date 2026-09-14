use strict;
use warnings;
use ExtUtils::MakeMaker;
use Config;
use File::Temp ();
use File::Spec;

# -std=gnu99 for the declarations-after-statements and // comments the XS uses;
# only when the compiler is one that understands the flag.
my $std = ($Config{gccversion} || $Config{cc} =~ /\b(?:g?cc|clang)\b/) ? ' -std=gnu99' : '';

# Parser.xs wants the `l' suffixed libm functions on a long-double perl, so that
# an NV is not silently computed at double precision and converted back -- see
# the nv_* block at the head of the file.  They are C99, but some libms (the
# thinner BSD ones especially) do not ship the whole set, so the set is
# link-tested as a unit: all of it, or the module falls back to the double
# functions, which costs accuracy on that perl and nothing else.  __float128
# needs no probe -- quadmath.h and -lquadmath come with the quadmath perl
# itself, and a double perl uses the unsuffixed functions it always had.
#
# Keep this list in step with the nv_* macros in Parser.xs: a macro with no
# entry here is a function the probe never tested for, and the build breaks on
# a machine nobody here owns.
my $define = '';
if ($Config{nvtype} eq 'long double') {
	my $src = "#include <math.h>\n"
	        . "int main(void) { long double x = 1.0L;\n"
	        . "  volatile long double r = fabsl(x) + sqrtl(x) + acosl(x) + cosl(x) + sinl(x)\n"
	        . "    + atan2l(x, x);\n"
	        . "  return (int)r * 0; }\n";
	my $dir = File::Temp::tempdir(CLEANUP => 1);
	my $c   = File::Spec->catfile($dir, 'ldmath.c');
	my $bin = File::Spec->catfile($dir, 'ldmath');
	if (open my $fh, '>', $c) {
		print {$fh} $src;
		close $fh;
		my $quiet = File::Spec->devnull;
		$define = '-DCSP_HAVE_LONG_DOUBLE_MATH'
			if system("$Config{cc}$std -o $bin $c -lm > $quiet 2> $quiet") == 0;
	}
	warn "Chem::Structure::Parser: this libm is missing some long double math\n"
	   . "  functions; the feature calculations fall back to double precision.\n"
		unless $define;
}

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
	# "autodie" is deliberately absent: the module checks the calls the pragma
	# would have checked, in its own words, rather than asking every installer
	# for it.  Kept in step with the [Prereqs] block in dist.ini.
	PREREQ_PM => {
		'Carp'         => 0,
		'Exporter'     => 0,
		'Scalar::Util' => 0,
		'XSLoader'     => 0,
	},
	TEST_REQUIRES => {
		'Test::More'      => 0,
		'Test::Exception' => 0,
	},
	CCFLAGS => $Config{ccflags} . $std,
	DEFINE  => $define,
	LIBS    => ['-lm'],
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
