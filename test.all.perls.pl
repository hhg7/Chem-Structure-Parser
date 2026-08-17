#!/usr/bin/env perl
# test.all.perls.pl - build & test Chem::Structure::Parser against every perlbrew perl.
#
# Automates the manual loop of:
#     perlbrew use perl-5.10.1 && ./compile.sh
#     perlbrew use perl-5.12.5 && ./compile.sh
#     perlbrew use perl-5.44.0 && ./compile.sh
#
# Each version is built in its own environment (no `perlbrew use` needed: the
# perl binary is invoked directly and PATH is rewritten for the child), the
# tree is cleaned between versions, and a pass/fail summary is printed at the
# end.  Exit status is non-zero if any version failed.
#
# This is an XS distribution, so a stale Parser.o or Parser.c compiled against a
# different perl will happily link and then crash: the cleaning here removes
# the XS products explicitly rather than trusting a Makefile to still exist.
#
# By default four perls are built at once with a four-job make each (-P 4 -j 4).
# They cannot share the distribution root for exactly the reason above -- one
# Makefile, one Parser.c, one blib -- so each parallel child builds in a private
# copy of the tree under the log directory and reports its result back to the
# parent, which leaves the distribution root untouched.  Nothing forces
# Parallel::ForkManager here: the fork/throttle/reap loop is a dozen lines and
# the script already forks in run_cmd, so this keeps the dev tooling free of
# non-core prerequisites.
#
# -P 1 is the old behaviour: one perl at a time, in the distribution root, which
# is then left built and installed against the perl that runs last -- the
# newest plain one, see the run order below.

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use Getopt::Long 'GetOptions';
# -P (parallel perls) and -p (pick a perl) are different options, so the
# default case-folding of single-letter aliases has to go.
Getopt::Long::Configure('no_ignore_case');
use File::Spec;
use File::Copy 'copy';
use File::Path 'remove_tree';
use IO::Handle;
use Cwd 'getcwd';
use Data::Dumper;
use POSIX 'strftime';
use Time::HiRes 'time';

my $PERLBREW_ROOT = $ENV{PERLBREW_ROOT} || File::Spec->catdir($ENV{HOME}, 'perl5', 'perlbrew');

# XS build products.  `make clean` handles these when a Makefile survived from
# the previous perl; when it did not, they have to go by hand.
my @XS_PRODUCTS = qw(Parser.c Parser.o Parser.bs Parser.so Parser.def);

my ($help, $list, $install, $deps, $clean, $stop, $quiet, $jobs, $log_dir, $optimize, $skip_old, @only);
my ($par, $keep_work);
$install  = 1; # make install, as compile.sh does
$clean    = 1;
$skip_old = 1; # perls older than MIN_PERL_VERSION cannot pass; do not pretend
$optimize = '-O2 -Wall -Wextra'; # same as compile.sh
$log_dir  = File::Spec->catdir('.build', 'multiperl');
# Defaults are "./test.all.perls.pl -P 4 -j 4": four perls at once, four make
# jobs and four test jobs each.  That is up to 16 compilers, which is the point
# -- the sweep drops from ~2 minutes to ~25 seconds.  -P 1 restores the old
# serial build in the distribution root; -j 0 the old single-job make.
$par      = 4;
$jobs     = 4;

GetOptions(
    'perl|p=s@'     => \@only,
    'install!'      => \$install,
    'deps!'         => \$deps,
    'clean!'        => \$clean,
    'stop-on-fail!' => \$stop,
    'skip-old!'     => \$skip_old,
    'jobs|j=i'      => \$jobs,
    'parallel|P:i'  => \$par,
    'keep-work!'    => \$keep_work,
    'optimize=s'    => \$optimize,
    'log-dir=s'     => \$log_dir,
    'quiet|q'       => \$quiet,
    'list|l'        => \$list,
    'help|h'        => \$help,
) or usage(1);
usage(0) if $help;

sub usage {
    my $rc = shift;
    print STDERR <<"END";
usage: $0 [options]

Builds and tests the distribution in the current directory against each
perl installed under $PERLBREW_ROOT.

Run order: the quadmath builds first, since their test suite is the slowest and
must not sit in the queue behind faster perls; then the rest oldest first; then
the newest plain perl (double NV, no ithreads) last, so a serial run leaves the
tree built and installed against the reference build rather than an outlier.

options:
  -p, --perl VERSION   only this perl (repeatable); accepts "5.10.1",
                       "perl-5.10.1" or an exact directory name such as
                       "5.44.0-quadmath".  default: every installed perl
  -l, --list           list the perls that would be tested, in run order, with
                       each one's NV width and threading, then exit
      --no-install     skip "make install" (build + test only)
      --no-clean       skip cleaning before each version
      --deps           cpanm any missing prerequisite for that perl first
                       (PREREQ_PM, TEST_REQUIRES, and the modules the test
                       suite loads optionally, e.g. Test::LeakTrace)
      --no-skip-old    attempt perls older than MIN_PERL_VERSION anyway
      --stop-on-fail   abort at the first version that fails (in parallel mode:
                       launch no more perls; the running ones finish)
  -j, --jobs N         parallel make, and HARNESS_OPTIONS=j<N> for the tests.
                       default: $jobs; -j 0 for a single-job make
  -P, --parallel [N]   build/test N perls at once, each in a private copy of the
                       tree under --log-dir/work (the distribution root is left
                       untouched).  bare -P, or -P 0, uses one child per perl,
                       capped at the CPU count; -P 1 builds serially in the
                       distribution root and leaves it built against the newest
                       plain perl.  note that -P and -j multiply: the defaults
                       can run 16 compilers.  default: $par
      --keep-work      keep the private build trees of perls that passed
                       (failed ones are always kept for inspection)
      --optimize STR   OPTIMIZE= passed to Makefile.PL (default: $optimize)
      --log-dir DIR    where per-version logs go (default: $log_dir)
  -q, --quiet          only write logs; do not echo build output.  implied for
                       -P > 1, where interleaved output would be unreadable
  -h, --help           this message

exit status: 0 if every perl built, tested and installed cleanly, else 1.
END
    exit $rc;
}

# ---------------------------------------------------------------- discovery --

my $perls_dir = File::Spec->catdir($PERLBREW_ROOT, 'perls');
die "$0: no perlbrew perls directory at $perls_dir\n" unless -d $perls_dir;

opendir my $dh, $perls_dir or die "$0: cannot read $perls_dir: $!\n";
my @installed = grep { -x File::Spec->catfile($perls_dir, $_, 'bin', 'perl') }
                grep { !/^\.\.?$/ } readdir $dh;
closedir $dh;

# numeric sort, oldest first, which is the order the run order below starts from.
sub vkey {
    my $v = shift;
    $v =~ s/^perl-//;
    my @p = ($v =~ /(\d+)/g);
    push @p, 0 while @p < 3;
    return sprintf '%05d%05d%05d', @p[0 .. 2];
}
@installed = sort { vkey($a) cmp vkey($b) } @installed;

my @targets = @installed;
if (@only) {
    my %have = map { $_ => 1 } @installed;
    my (@want, @missing);
    for my $arg (map { split /,/ } @only) {
        # an exact directory name wins, so builds like "5.44.0-quadmath" that
        # carry no perl- prefix are still selectable.
        my ($match) = grep { $have{$_} } $arg, "perl-$arg";
        if (defined $match) { push @want, $match }
        else                { push @missing, $arg }
    }
    die "$0: not installed under perlbrew: @missing\n(installed: @installed)\n" if @missing;
    my %seen;
    @targets = sort { vkey($a) cmp vkey($b) } grep { !$seen{$_}++ } @want;
}
die "$0: no perls found in $perls_dir\n" unless @targets;

# ------------------------------------------------------------------ run order --

# What each perl actually is.  The NV width and the threading model decide the
# order below, and the answer has to come from the interpreter: "5.44.0-quadmath"
# is a local naming habit rather than a promise, and a threaded or long-double
# build usually says nothing at all in its directory name.  A perl that will not
# answer is described as unknown, never as plain, so a broken interpreter cannot
# become the reference build the tree is left standing on.
my %facts;
sub facts {
    my $version = shift;
    return $facts{$version} if $facts{$version};
    my $perl = File::Spec->catfile($perls_dir, $version, 'bin', 'perl');
    my %f = (nv => '?', threads => 0, known => 0);
    if (open my $fh, '-|', $perl, '-MConfig', '-e',
            'print "$Config{nvtype}\t", ($Config{useithreads} ? 1 : 0)') {
        my $line = <$fh>;
        close $fh;
        if (defined $line && $line =~ /^(\S[^\t]*)\t([01])/) {
            %f = (nv => $1, threads => $2, known => 1);
        }
    }
    $f{quadmath} = $f{nv} eq '__float128';
    # plain: the double-NV, unthreaded build every other configuration is a
    # variation on.  Unknown does not qualify.
    $f{plain} = $f{known} && $f{nv} eq 'double' && !$f{threads};
    return $facts{$version} = \%f;
}

sub nv_label {
    my $f = shift;
    my %short = ('double' => 'double', 'long double' => 'long-double',
                 '__float128' => 'quadmath');
    return ($short{ $f->{nv} } || $f->{nv}) . ($f->{threads} ? '-thr' : '');
}

# Run order, given @targets oldest first:
#
#   1. the quadmath builds, because their test suite is far and away the
#      slowest (minutes against seconds here), and the whole matrix is only as
#      fast as the slowest perl that had to wait for a free slot;
#   2. everything else, still oldest first;
#   3. the newest plain perl last, because a serial run leaves the distribution
#      root built and installed against whichever perl went last, and that
#      should be the ordinary reference build rather than a long-double,
#      threaded or quadmath outlier.
sub order_targets {
    my @t = @_;
    my (@quad, @rest);
    push @{ facts($_)->{quadmath} ? \@quad : \@rest }, $_ for @t;
    my ($plain) = grep { facts($_)->{plain} } reverse @rest;   # newest first
    @rest = grep { $_ ne $plain } @rest if defined $plain;
    return (@quad, @rest, defined $plain ? $plain : ());
}
@targets = order_targets(@targets);

if ($list) {
    printf "%-20s %-16s %s\n", $_, nv_label(facts($_)),
        File::Spec->catfile($perls_dir, $_, 'bin', 'perl') for @targets;
    exit 0;
}

die "$0: no Makefile.PL in " . getcwd() . " - run this from the distribution root\n"
    unless -f 'Makefile.PL';
die "$0: no Parser.xs in " . getcwd() . " - is this the Chem-Structure-Parser tree?\n"
    unless -f 'Parser.xs';

# ------------------------------------------------------------------ logging --

sub mkdirp {
    my @parts = File::Spec->splitdir(shift);
    my $path;   # undef, not '': splitdir gives an absolute path a leading ''
                # and catdir('', 'home') is '/home', where 'home' would be a
                # directory of that name in the cwd.
    for my $p (@parts) {
        $path = defined $path ? File::Spec->catdir($path, $p) : $p;
        next if !length($path) || -d $path;
        mkdir $path or die "$0: mkdir $path: $!\n";
    }
}

mkdirp($log_dir);
my $stamp = strftime '%Y%m%d-%H%M%S', localtime;

# ------------------------------------------------------------- parallelism --

sub cpu_count {
    open my $fh, '<', '/proc/cpuinfo' or return 4;
    my $n = grep { /^processor\s*:/ } <$fh>;
    close $fh;
    return $n || 4;
}

# bare -P (or -P 0) means "all of them", within reason: more compilers than
# cores only makes every perl slower.
$par = 0 if $par < 0;
$par ||= @targets < cpu_count() ? scalar @targets : cpu_count();
$par = @targets if $par > @targets;

my $work_root    = File::Spec->catdir($log_dir, 'work');
my $in_parallel  = 0; # set once the dispatch below knows there is more than one
if ($par > 1) {
    # children chdir into their own tree, so the log paths must not be relative
    $log_dir   = File::Spec->rel2abs($log_dir);
    $work_root = File::Spec->rel2abs($work_root);
}

# ------------------------------------------------------------- child runner --

# Run @cmd with STDERR folded into STDOUT, echoing to the terminal and to
# $logfh.  Returns ($exit_code, \@lines).
sub run_cmd {
	my ($cmd, $logfh) = @_;
	print $logfh "\n\$ @$cmd\n";
	print "\$ @$cmd\n" unless $quiet;

	my $pid = open my $fh, '-|';
	die "$0: fork failed: $!\n" unless defined $pid;
	if (!$pid) {                              # child
	  open STDERR, '>&', \*STDOUT or die "$0: dup STDERR: $!\n";
	  $| = 1;
	  { exec { $cmd->[0] } @$cmd; }
	  print "exec @$cmd failed: $!\n";
	  exit 127;
	}

	my @lines;
	while (defined(my $line = <$fh>)) {
	  push @lines, $line;
	  print $logfh $line;
	  print $line unless $quiet;
	}
	close $fh;
	my $status = $?;
	my $code = $status == -1              ? -1
		    : ($status & 127)            ? 128 + ($status & 127)
		    :                              ($status >> 8);
	return ($code, \@lines);
}

# ------------------------------------------------------- distribution facts --

sub slurp_makefile_pl {
    open my $fh, '<', 'Makefile.PL' or return '';
    local $/;
    my $src = <$fh>;
    close $fh;
    return $src;
}
my $mpl = slurp_makefile_pl();

# Prerequisite names, scraped from Makefile.PL so this stays in sync with
# dist.ini however the Makefile.PL was written: hand-rolled (PREREQ_PM => {})
# or generated by Dist::Zilla ("PREREQ_PM" => {}).
sub prereqs {
    my %seen;
    my @mods;
    for my $key (qw(PREREQ_PM TEST_REQUIRES BUILD_REQUIRES)) {
        next unless $mpl =~ /"?\Q$key\E"?\s*=>\s*\{(.*?)\}/s;
        my $block = $1;
        push @mods, grep { !$seen{$_}++ && $_ ne 'perl' }
                    ($block =~ /"?([A-Za-z][\w:]*)"?\s*=>/g);
    }
    return @mods;
}

# Modules the suite loads only if they are there (t/leaks.t needs
# Test::LeakTrace, t/gzip.t needs IO::Compress::Gzip).  Without them those
# files skip, which looks like a pass while testing nothing, so --deps
# installs them too.
sub optional_test_deps {
    my (%seen, @mods);
    for my $t (sort glob 't/*.t') {
        open my $fh, '<', $t or next;
        while (my $line = <$fh>) {
            next unless $line =~ /eval\s*\{/;
            push @mods, grep { !$seen{$_}++ } ($line =~ /require\s+([A-Z][\w:]*)/g);
        }
        close $fh;
    }
    return @mods;
}

# MIN_PERL_VERSION, so an ancient perl is reported as skipped rather than as a
# failure the distribution never claimed to support.  Both the numeric form
# ('5.010', '5.010001') and the dotted form ('5.10.1') turn up in the wild.
sub min_perl_version {
    return undef unless $mpl =~ /"?MIN_PERL_VERSION"?\s*=>\s*['"]?v?([\d._]+)/;
    my $v = $1;
    $v =~ s/_//g;
    return $v =~ /^(\d+)\.(\d{3})(\d{0,3})$/
        ? sprintf('%d.%d.%d', $1, $2, $3 || 0)   # 5.010001 -> 5.10.1
        : $v;
}

my @prereqs  = (prereqs(), optional_test_deps());
my $min_perl = min_perl_version();
my $min_key  = defined $min_perl ? vkey($min_perl) : undef;

# ------------------------------------------------------------ build one perl --

# Build, test and install the distribution in the current directory with one
# perl.  Returns the result hashref; prints its own progress unless $silent
# (parallel children are silent and the parent reports for them).
sub build_one {
    my ($version, $silent) = @_;
    my $root = File::Spec->catdir($perls_dir, $version);
    my $bin  = File::Spec->catdir($root, 'bin');
    my $perl = File::Spec->catfile($bin, 'perl');

    my $log = File::Spec->catfile($log_dir, "$version.$stamp.log");
    open my $logfh, '>', $log or die "$0: cannot write $log: $!\n";
    $logfh->autoflush(1);
    STDOUT->autoflush(1);

    unless ($silent) {
        print "\n", '=' x 72, "\n";
        printf "== %s   (log: %s)\n", $version, $log;
        print '=' x 72, "\n";
    }
    print $logfh "== $version at " . strftime('%F %T', localtime)
               . ' in ' . getcwd() . "\n";

    # Emulate `perlbrew use $version` for the children: this perl's bin first,
    # every other perlbrew perl stripped out, and local::lib / PERL5LIB
    # leftovers from the calling shell removed so nothing bleeds across
    # versions.
    local %ENV = %ENV;
    my @path = grep { index($_, File::Spec->catdir($perls_dir, '')) != 0 }
               split /:/, ($ENV{PATH} || '/usr/bin:/bin');
    $ENV{PATH}          = join ':', $bin, @path;
    $ENV{PERLBREW_ROOT} = $PERLBREW_ROOT;
    $ENV{PERLBREW_PERL} = $version;
    $ENV{PERLBREW_PATH} = $bin;
    delete @ENV{qw(PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MM_OPT PERL_MB_OPT
                   PERLBREW_LIB PERL_MM_USE_DEFAULT)};
    $ENV{HARNESS_OPTIONS} = "j$jobs" if $jobs;

    my $t0 = time;
    my %r = (version => $version, log => $log, steps => [], warnings => 0);

    # sanity: the interpreter really is the version we think it is, and say so
    # when it is a threaded build -- an XS bug can be invisible on an
    # unthreaded perl and a hard compile error under MULTIPLICITY, so the
    # summary has to make that coverage visible rather than implied by a name.
    my ($vc, $vout) = run_cmd([$perl, '-MConfig', '-e',
        'printf "%vd%s\n", $^V, $Config{useithreads} ? "-thr" : ""'], $logfh);
    $r{reported} = $vc == 0 && @$vout ? do { my $s = $vout->[0]; chomp $s; $s } : '?';

    my @steps;
    if ($deps && @prereqs) {
        my @missing;
        for my $mod (@prereqs) {
            my ($c) = run_cmd([$perl, "-M$mod", '-e', '1'], $logfh);
            push @missing, $mod if $c != 0;
        }
        if (@missing) {
            my $cpanm = -x File::Spec->catfile($bin, 'cpanm')
                      ? File::Spec->catfile($bin, 'cpanm')
                      : File::Spec->catfile($PERLBREW_ROOT, 'bin', 'cpanm');
            push @steps, ['deps', [$perl, $cpanm, '--notest', @missing]];
        }
    }

    if ($clean) {
        push @steps, ['clean', ['make', 'clean'], 1] if -f 'Makefile'; # 1 = failure tolerated
        clean_xs($logfh);
    }

    push @steps, ['Makefile.PL', [$perl, 'Makefile.PL', "OPTIMIZE=$optimize"]];
    push @steps, ['make',        ['make', $jobs ? ("-j$jobs") : ()]];
    push @steps, ['make test',   ['make', 'test']];
    push @steps, ['make install',['make', 'install']] if $install;

    my $failed;
    for my $step (@steps) {
        my ($label, $cmd, $soft) = @$step;
        my $t_step = time;
        my ($code, $lines) = run_cmd($cmd, $logfh);
        my $secs = time - $t_step;
        printf $logfh "-- step '%s' exited %d after %.1fs\n", $label, $code, $secs;

        if ($label eq 'make test') {
            for my $l (@$lines) {
                $r{files}  = $1 if $l =~ /^(Files=\d+.*)/;
                $r{result} = $1 if $l =~ /^Result:\s*(\S+)/;
                $r{passed} = 1  if $l =~ /^All tests successful/;
                # a skipped file tests nothing; XS coverage that quietly
                # vanished on one perl is worth seeing in the summary.
                $r{skips}++     if $l =~ m{^t/\S+\.t\s+.*\bskipped\b}i;
            }
        }
        # residual build warnings are worth surfacing even on success
        $r{warnings} += grep { /: warning:/ } @$lines if $label eq 'make';

        push @{ $r{steps} }, { label => $label, code => $code, seconds => $secs };
        next if $code == 0 || $soft;
        $failed = $label;
        last;
    }

    if ($clean && !$failed && !$in_parallel && $version ne $targets[-1]) {
        # leave a clean tree behind before the next perl takes over.  in
        # parallel mode the whole private tree goes instead, so this would
        # only cost a `make clean` nobody reads.
        run_cmd(['make', 'clean'], $logfh);
        clean_xs($logfh);
    }

    $r{seconds} = time - $t0;
    $r{failed}  = $failed;
    close $logfh;

    report_one(\%r) unless $silent;
    return \%r;
}

# the two progress lines a finished perl prints, from the parent in either mode
sub report_one {
    my $r = shift;
    printf "-- %s: %s in %.1fs%s\n", $r->{version},
        ($r->{failed} ? "FAILED at '$r->{failed}'" : 'ok'),
        $r->{seconds},
        (defined $r->{files} ? " ($r->{files})" : '');
    print '-- ', join('  ', map { sprintf '%s %.1fs', $_->{label}, $_->{seconds} }
                            @{ $r->{steps} }), "\n";
}

# Remove the XS products by hand.  `make clean` cannot be relied on: the
# Makefile it needs is itself one of the things being replaced, and a Parser.o
# built by another perl links without complaint and then dies at load.
sub clean_xs {
    my $logfh = shift;
    my @gone = grep { -e $_ && unlink $_ } @XS_PRODUCTS;
    print $logfh "-- removed stale XS products: @gone\n" if @gone;
    print "-- removed stale XS products: @gone\n" if @gone && !$quiet;
}

# ------------------------------------------------------------ private trees --

# Directories and files that must not be copied into a private build tree:
# the build products of whichever perl last used the source tree (copying them
# would recreate exactly the stale-object hazard the cleaning avoids), the
# repository, and the log/work directory itself.
my %SKIP_DIR  = map { $_ => 1 } qw(.git blib .build);
my @SKIP_FILE = (qr/\.(?:o|a|so|bs|dylib|dll)$/, qr/^Parser\.(?:c|def)$/,
                 qr/^Makefile(?:\.old)?$/, qr/^pm_to_blib$/, qr/^MYMETA\./,
                 qr/^Chem-Structure-Parser-.*\.tar\.gz$/);

sub copy_tree {
    my ($src, $dst) = @_;
    mkdirp($dst);
    opendir my $dh, $src or die "$0: cannot read $src: $!\n";
    my @entries = grep { !/^\.\.?$/ } readdir $dh;
    closedir $dh;
    for my $e (@entries) {
        my $s = File::Spec->catfile($src, $e);
        my $d = File::Spec->catfile($dst, $e);
        if (-d $s) {
            next if $SKIP_DIR{$e};
            # a --log-dir inside the tree would otherwise copy itself forever
            my $abs = File::Spec->rel2abs($s);
            next if $abs eq $work_root || $abs eq $log_dir;
            copy_tree($s, $d);
            next;
        }
        next if grep { $e =~ $_ } @SKIP_FILE;
        copy($s, $d) or die "$0: copy $s -> $d: $!\n";
        chmod((stat $s)[2] & 07777, $d);
    }
}

# The child's result hashref has to cross a fork, and %r is plain data, so
# Data::Dumper out / eval in beats any IPC here: it also leaves the file behind
# next to the log when something needs explaining.
sub write_result {
    my ($file, $r) = @_;
    open my $fh, '>', $file or die "$0: cannot write $file: $!\n";
    local $Data::Dumper::Indent   = 0;
    local $Data::Dumper::Sortkeys = 1;
    print $fh Data::Dumper->Dump([$r], ['R']);
    close $fh or die "$0: close $file: $!\n";
}

sub read_result {
    my $file = shift;
    open my $fh, '<', $file or return undef;
    local $/;
    my $src = <$fh>;
    close $fh;
    my $R;
    eval $src;                      ## no critic -- our own Dumper output
    return ref $R eq 'HASH' ? $R : undef;
}

# ------------------------------------------------------------------ drivers --

sub run_serial {
    my @queue = @_;
    my @out;
    for my $version (@queue) {
        my $r = build_one($version);
        push @out, $r;
        next unless $r->{failed} && $stop;
        my %done = map { $_->{version} => 1 } @out;
        my @rest = grep { !$done{$_} } @queue;
        print "-- --stop-on-fail: skipping @rest\n" if @rest;
        last;
    }
    return @out;
}

# Fork up to $par children, each in its own copy of the tree, reaping as they
# finish and starting the next perl in the freed slot.
sub run_parallel {
    my @order   = @_;
    my @queue   = @order;
    my $src     = getcwd();
    my %kid;                        # pid => { version, work, result }
    my (@out, $halt);

    my $reaper = sub {
        my $pid    = shift;
        my $status = shift;
        my $kid    = delete $kid{$pid} or return;
        my $r      = read_result($kid->{result});
        if (!$r) {
            # the child died without reporting: exec failure, signal, OOM
            $r = { version => $kid->{version}, log => $kid->{log},
                   reported => '?', steps => [], warnings => 0, seconds => 0,
                   failed => sprintf('child exited %d%s', $status >> 8,
                                     ($status & 127) ? ' on signal ' . ($status & 127) : '') };
        }
        push @out, $r;
        report_one($r);
        $halt = 1 if $r->{failed} && $stop;
        if ($r->{failed} || $keep_work) { print "-- $kid->{version}: build tree kept at $kid->{work}\n" }
        else                            { remove_tree($kid->{work}) }
    };

    local $SIG{INT} = local $SIG{TERM} = sub {
        my $sig = shift;
        print "\n-- $sig: stopping " . keys(%kid) . " running build(s)\n";
        # each child leads its own process group (see the fork below), so one
        # signal per group takes its make, its compilers and its prove with it;
        # signalling the child perl alone would orphan those.
        kill 'TERM', map { -$_ } keys %kid;
        sleep 1;
        kill 'KILL', map { -$_ } keys %kid;
        exit 130;
    };

    while (@queue || %kid) {
        while (@queue && keys(%kid) < $par && !$halt) {
            my $version = shift @queue;
            my $work    = File::Spec->catdir($work_root, "$version.$stamp");
            my $result  = File::Spec->catfile($log_dir, "$version.$stamp.result");
            my $log     = File::Spec->catfile($log_dir, "$version.$stamp.log");

            remove_tree($work) if -d $work;
            printf "-- %-20s starting (tree: %s)\n", $version, $work;
            copy_tree($src, $work);

            my $pid = fork;
            die "$0: fork failed: $!\n" unless defined $pid;
            if (!$pid) {                                  # child
                $SIG{$_} = 'DEFAULT' for qw(INT TERM);
                setpgrp 0, 0;   # so ^C reaches this whole build, once, via the
                                # parent's handler rather than the terminal
                chdir $work or die "$0: chdir $work: $!\n";
                my $r = build_one($version, 1);
                write_result($result, $r);
                exit 0;
            }
            $kid{$pid} = { version => $version, work => $work,
                           result => $result, log => $log };
        }
        last unless %kid;
        my $pid = waitpid -1, 0;
        last if $pid <= 0;
        $reaper->($pid, $?);
    }

    if ($halt && @queue) {
        print "-- --stop-on-fail: skipping @queue\n";
    }
    # report in version order, not completion order
    my %by = map { $_->{version} => $_ } @out;
    return map { $by{$_} } grep { $by{$_} } @order;
}

# --------------------------------------------------------------- main loop --

my (@results, @skipped);
my $t_all = time;

my @queue;
for my $version (@targets) {
    if ($skip_old && defined $min_key && vkey($version) lt $min_key) {
        printf "-- %s: SKIP (older than MIN_PERL_VERSION %s)\n", $version, $min_perl;
        push @skipped, { version => $version, why => "< $min_perl" };
        next;
    }
    push @queue, $version;
}

if ($par > 1 && @queue > 1) {
    $in_parallel = 1;
    $quiet       = 1; # N interleaved build logs on one terminal is noise
    printf "-- %d perl(s), %d at a time%s; per-version output goes to the logs\n",
        scalar @queue, $par, ($jobs ? " with make -j$jobs each" : '');
    print "-- note: --deps shares one ~/.cpanm across the children; if a "
        . "prerequisite install misbehaves, run once with -P 1 --deps first\n"
        if $deps;
    @results = run_parallel(@queue);
}
else {
    @results = run_serial(@queue);
}

# ----------------------------------------------------------------- summary --

my $bad = grep { $_->{failed} } @results;
print "\n", '=' x 72, "\n";
printf "%-20s %-11s %-8s %-7s %-6s %s\n",
    qw(PERL REPORTED STATUS TIME WARN TESTS);
print '-' x 79, "\n";
for my $r (@results) {
    printf "%-20s %-11s %-8s %6.1fs %-6s %s\n",
        $r->{version},
        $r->{reported},
        ($r->{failed} ? 'FAIL' : 'PASS'),
        $r->{seconds},
        ($r->{warnings} || 0),
        ($r->{failed} ? "failed at '$r->{failed}' - see $r->{log}"
                      : ($r->{result}
                            ? "Result: $r->{result}"
                              . ($r->{skips} ? " ($r->{skips} file(s) skipped)" : '')
                            : 'no test summary')),
        ;
}
for my $s (@skipped) {
    printf "%-20s %-11s %-8s %7s %-6s %s\n", $s->{version}, '-', 'SKIP', '-', '-',
        "not attempted: perl $s->{why}";
}
print '-' x 79, "\n";
my $not_run = @targets - @results - @skipped;
printf "%d/%d perl(s) passed%s in %.1fs.  Logs in %s\n",
    scalar(@results) - $bad, scalar(@targets) - scalar(@skipped),
    ($not_run ? " ($not_run not run)" : ''),
    time - $t_all, $log_dir;
# every perl installed into its own site_perl, but nothing built here: say so,
# because after a serial run the root held a tree built against the newest plain
# perl.
print "-- built in private trees; this directory is unbuilt (-P 1 to build here)\n"
    if $in_parallel;

exit($bad || $not_run ? 1 : 0);
