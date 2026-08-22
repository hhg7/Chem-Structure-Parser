#!/bin/sh
# dzil.sh -- build the release tarball.
#
# Adapted from ~/Scripts/stats/dzil.sh.  md2pod.pl runs first because the POD
# half of lib/Chem/Structure/Parser.pm and the whole of Changes are generated
# from README.md: building without it would ship the documentation as it was
# before the last edit to the README.
#
# The tarball is listed rather than released.  What to look for is what
# PruneFiles in dist.ini and MANIFEST.SKIP are there to keep out -- the XS build
# products (Parser.c, Parser.o, Parser.bs, Parser.so), blib/, the author-only
# root helpers (test.all.perls.pl, benchmark.pl, structure.info.pl, md2pod.pl,
# read.me.pod, compile.sh, dzil.sh, CLAUDE.md, notes.txt), and any
# Chem-Structure-Parser-*/ directory left over from an earlier build.  Parser.xs
# and ppport.h *do* belong in it: they are what the module is compiled from.
set -e
perl md2pod.pl
git commit -am "Update generated docs" || true   # -a = tracked files only; NOT -A
dzil clean
dzil build
echo "==== tarball contents (verify: no .c/.o/.so/.bs/blib, no Chem-Structure-Parser/) ===="
tar tzf Chem-Structure-Parser-*.tar.gz
echo "If that looks clean, run: dzil release"
