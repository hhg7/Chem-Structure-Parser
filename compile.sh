#!/bin/sh
make clean 2>/dev/null; perl Makefile.PL OPTIMIZE='-O2 -Wall -Wextra' && make && make test && make install
