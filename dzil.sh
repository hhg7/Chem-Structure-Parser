#!/bin/sh
# dzil.sh -- build the release tarball, and commit what building it took.
#
# Adapted from ~/Scripts/SimpleFlow/dzil.sh, which is where the ordering below
# is argued; the two places this distribution differs from it are marked.
# md2pod.pl runs first because the POD half of lib/Chem/Structure/Parser.pm is
# generated from README.md: building without it would ship the documentation as
# it was before the last edit to the README.  It also runs changes_file_ok over
# Changes, which is written by hand, so a release whose notes do not satisfy
# CPAN::Changes::Spec stops here.
#
# The tarball is listed rather than released.  What to look for is what
# PruneFiles in dist.ini and MANIFEST.SKIP are there to keep out -- the XS build
# products (Parser.c, Parser.o, Parser.bs, Parser.so), blib/, the author-only
# root helpers (test.all.perls.pl, benchmark.pl, structure.info.pl, md2pod.pl,
# read.me.pod, compile.sh, dzil.sh, CLAUDE.md, notes.txt), and any
# Chem-Structure-Parser-*/ directory left over from an earlier build.  Parser.xs
# and ppport.h *do* belong in it: they are what the module is compiled from.
set -e

# The version dzil is about to build, read from the one place it is set:
# dist.ini's [VersionFromModule] takes it from this same line.  Needed below to
# name the tarball to list and the build to record in the commit message.  This
# module writes it unquoted, where SimpleFlow quotes it, so the pattern is not
# the same one.
VERSION=$(sed -n 's/^our $VERSION = \([0-9][0-9.]*\);.*$/\1/p' lib/Chem/Structure/Parser.pm)
if [ -z "$VERSION" ]; then
	echo "dzil.sh: no \$VERSION found in lib/Chem/Structure/Parser.pm" >&2
	exit 1
fi

perl md2pod.pl

# Staged before the build rather than committed after it, which is the half of
# the job that used to want a second run: [Git::GatherDir] gathers what git
# knows about, so a file that exists only in the working tree -- a new
# SECURITY.md, a new test, a new fixture -- is left out of the tarball however
# the PruneFiles allowlist reads, and only the run after the commit shipped it.
# The index is what `git ls-files` reports, so staging is enough for the build
# two lines below to see it.
#
# The pathspec is the distribution itself, which is the same set PruneFiles
# allows through.  A bare `git add -A` would also sweep up whatever else is
# sitting untracked in the root -- notes.txt, a stray cover_db, the editor's
# droppings.
git add -A -- lib t Parser.xs ppport.h Makefile.PL \
              README.md Changes CONTRIBUTING.md SECURITY.md

dzil clean
dzil build

echo "==== tarball contents (verify: no .c/.o/.so/.bs/blib, no Chem-Structure-Parser/) ===="
tar tzf "Chem-Structure-Parser-$VERSION.tar.gz"

# The commit comes last so that one run finishes the job.  SimpleFlow's stages
# the tarball and the unpacked directory into it as well, because there they
# are tracked; here .gitignore keeps both out of the history, so the build
# itself leaves nothing to commit.  -a takes the tracked files md2pod.pl just
# rewrote, on top of what was staged above.
git commit -am "Update generated docs; build $VERSION" || true # an empty commit is not an error

echo "If that looks clean, run: dzil release"
