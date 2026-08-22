# Chem::Structure::Parser — instructions for Claude

These are `~/Scripts/stats/CLAUDE.md` read against this distribution: the same
doctrine, spelled with this module's types, macros, tests and reference
implementations. Where the two files differ, this one is right *here*. The
comment doctrine in `~/.claude/CLAUDE.md` applies unchanged, with the local
specifics at the foot of this file.

## Never write release notes

The release notes are the maintainer's, and are not Claude's to write — not
even to add an entry for work Claude just did, and not even when asked to
"update the changelog" as part of a larger task. When a change would normally
warrant one, say so in the reply and leave the wording to the maintainer.

There are two files to keep away from, because there are two copies of the
notes:

- **`Changes` is generated.** `md2pod.pl` writes it from README.md's
  `# Changes` section, in the form `CPAN::Changes::Spec` asks for, and the next
  run overwrites whatever is in it. So a hand edit there is lost work, not just
  a rule broken. Do not create, edit, append to, or revert it: no `Edit`, no
  `Write`, no `sed -i`/`perl -pi`, no `git checkout`/`git revert` that touches
  it, no patch that includes it.
- **README.md's `# Changes` section is the source those notes are written in**,
  and the same rule covers it. Editing the README elsewhere — a function's
  documentation, an example, a new section — is ordinary work; editing that one
  section is not.

## The documentation is generated from README.md

README.md is the single source. `md2pod.pl` (run it directly, or through
`./dzil.sh`, which runs it first) produces three things from it:

- the POD half of `lib/Chem/Structure/Parser.pm` — everything after the `1;`
  line is replaced, so **never hand-edit the POD**; edit README.md and
  regenerate. A section that exists only in the POD is a section the next run
  deletes.
- `read.me.pod`, the same POD on its own, for inspection. It is generated, so
  `.gitignore` has it, and `MY::libscan` in `Makefile.PL` drops root-level
  `*.pod` and `*.pl` so `make install` does not scatter it (and a
  `Chem::Structure::read.me` man page) into `site_perl`.
- `Changes`, from the `# Changes` section, as above.

Both `pod_file_ok` on the module and `changes_file_ok` on `Changes` run before
`md2pod.pl` exits; a release that fails either is not built. The script is
adapted from `~/Scripts/stats/md2pod.pl`, which is where its workarounds are
explained — take a fix to a conversion bug back there too.

README.md follows the conventions of `~/Scripts/stats/README.md`, and
`md2pod.pl` is written against them:

- Code and output are **4-space indented blocks**, never ``` fences. The
  Changes generator reads a fence, but the POD converter reformats prose it
  mistakes for a list, so indented is the form that survives both.
- `#` for a top-level section, `##` for one function, `###` for a subsection
  inside it. `h()` shows a function's section by looking its `=head2` up in the
  module's own POD, and `_pod_sections()` keeps only names in `@EXPORT_OK` — so
  a `##` heading under `# Functions/Subroutines` is a function, and a `##`
  heading anywhere else must not be spelled like one.
- Options and record mappings are GFM tables; `md2pod.pl` turns each into an
  `=begin html` block, because POD has no table.
- A release in `# Changes` is `## <version> <date> TZ` (`## 0.01 2026-08-21
  CDT`), then paragraphs at column 0 — one paragraph becomes one ` - ` bullet.
  `###` inside a release becomes a `[group]` and `####` a `- group:`. A version
  with no date makes `changes_file_ok` fail, and `md2pod.pl` says which heading
  to fix.
- `# COPYRIGHT AND LICENSE` is the last section: it is where the Changes
  generator stops reading.

## The two formats must return the same structure

The one claim the module makes above all others is that which format a
structure arrived in does not change the structure: `structure_info()` on
`1cka.pdb` and on `1cka.cif` returns hashes that compare equal by `is_deeply`
over the whole coordinate half, not merely equal in the fields someone thought
to check. `t/cif.t` is where that is asserted.

So a field added to one reader is not finished until it is in the other:

- `parse_buf()` (PDB) and `parse_cif_buf()` (mmCIF) both fill the same
  parallel column arrays, the same residue boundaries, the same counts. Add a
  column to `col_name[]`/the `C_*` enum and both readers must fill it.
- Where the formats disagree about how to write a fact, mmCIF converts to what
  the PDB reader would have produced — a formal charge of `-1` comes back as
  `"1-"`, `auth_*` is preferred over `label_*`, and `.`/`?` come back as the
  empty field a PDB record would have had. Keep new fields on that side of the
  line and Perl never learns there are two formats.
- The fixture pair goes in `t/data` (written by `t/data/generate.pl`, not
  typed) and the `is_deeply` pair goes in `t/cif.t`.

The header is the one place the two part company: a PDB file hands its header
back as raw lines grouped by record name, an mmCIF file as tags and loops.

## What is written in C, and what is not

`notes.txt` records the measurements this division came out of; do not undo it
on taste.

- C makes one pass over the bytes and does three things: splits ATOM/HETATM
  into parallel per-field arrays, marks residue boundaries, and groups every
  other record by record name. Residue name lookup is a switch on three packed
  bytes, and one table serves `aa3to1()`, `res1()` and `res_type()` so the
  three can never disagree.
- Header records (HEADER, COMPND, REMARK, SEQRES, …) are parsed in Perl. They
  are a few dozen, irregular, and the part most likely to need fixing later.
  Do not move them into C.
- Per-atom work belongs in the parse. The first round of benchmarking had the
  whole read coming out *slower* than a naive pure-Perl reader because every
  per-atom statistic was gathered in Perl over data C had just read; the second
  showed each atom being built twice. Anything that must touch every atom goes
  in the C pass.
- Atoms are parallel column arrays, not a hash per atom: that is roughly an
  order of magnitude off peak memory on a 33 MB entry, and it is why
  `atom_hashes` is an option rather than the default.

## C types must match the value's real domain

Every variable and parameter in `Parser.xs` gets the type that describes what
the value can actually be, not the type that happens to be convenient. Plain
`int` is not a default — reach for it only when the domain genuinely is
"signed, machine-word-ish, and not one of the cases below" (`errno`, and the
`int` that `toupper()` returns, are the honest cases in this file).

- **Two states only** → `bool`, assigned `TRUE`/`FALSE`, returned as
  `TRUE`/`FALSE` (perl's own, already in scope via `perl.h`; do not add
  `<stdbool.h>`). Never `bool flag = 0` and never `return 1` from a
  `bool` function.
- **A length, an offset, or a column range in the file buffer** → `STRLEN`.
  Every accessor in the column-field block (`fld()`, `trim_ptr()`, `fld_iv()`,
  `fld_nv()`) is written that way, and so is `pos`, because these lengths cross
  the perl string API.
- **A count or index that does not cross the perl API** → `size_t`
  (`ntags`, `cap`, `col_i`, and the `i` that indexes them in
  `parse_cif_buf()`). Do not "save" a `size_t` with `unsigned int` when the
  value tracks the size of something the file controls.
- **An integer handed to or taken from perl** → `UV`, and `IV` when it must be
  signed. The tallies (`n_atom`, `n_models`, `n_skipped`, `lineno`, …) are `UV`
  because they only ever count up; the model and residue numbers next to them
  stay `IV` because they are read out of the file's own columns, where a
  negative number is legal — residue numbering routinely starts before 1 — and
  is better carried as the negative it was written as than wrapped to a huge
  unsigned one.
- **A packed fixed-width key** → `U32`. `K3()` shifts a byte 16 places and
  `K2()` eight; plain `unsigned` is only guaranteed 16 bits, and the switches
  in `res_lookup()` and `elem_case()` are the whole residue and element tables.
- **A small, enumerated set of possibilities** → `short int` when `-1` is a
  sentinel (`cif_atom_field()` and the `fld_of[]` it fills: "a column this
  reader ignores is -1"), `unsigned short int` when none is needed
  (`cif_tok.kind`, one of the seven `CT_*`). Comment what each value means at
  its declaration, as `RT_*` and `CT_*` do.
- **A loop counter bounded by a small literal** → `unsigned short int`, as in
  `for (unsigned short int i = 0; i < NCOL; i++)`. A counter bounded by a
  runtime count is `size_t` instead.
- **A byte** → `unsigned char` (`res_info.type`, and every `ctype.h` argument:
  `isdigit((unsigned char)s[i])`, never a bare `char`).
- **Floating point** → `NV` throughout, never bare `double`, so the module
  keeps working on long-double and quadmath perls.

Widening or narrowing must never silently change a perl-visible signature or
conversion (`SvUV` vs `SvIV` vs `SvNV`, `%zu` vs `%" UVuf "`); fix the format
strings and the `Sv*` calls in the same edit. Mind the mixed-sign comparison
warnings this can introduce: the build is expected to stay clean under
`OPTIMIZE='-O2 -Wall -Wextra'`, which is what `compile.sh` and
`test.all.perls.pl` use.

## Use `CSP_RESTRICT`, never bare `restrict`

Pointer parameters and pointer locals get `CSP_RESTRICT` unless there is a real
reason they cannot. Apply it even when the compiler emits byte-identical code:
it documents the non-aliasing contract for the next reader, and perl's
`ccflags` carry `-fno-strict-aliasing`, so `restrict` is one of the few
aliasing facts a compiler here still has.

It is spelled through the macro, never as the C99 keyword. The block below the
includes in `Parser.xs` maps it to `restrict`, `__restrict__` or `__restrict`
per compiler and defines it away where there is none; bare `restrict` compiles
on this machine only because `Makefile.PL` appends `-std=gnu99` for gcc and
clang, and says nothing for MSVC (`__restrict`), an older gcc (`__restrict__`)
or a C89 vendor compiler. `-DCSP_RESTRICT=` is the escape hatch, and is how the
before/after timings in `notes.txt` were taken.

Pair it with `const` when the pointee is not written, which is the existing
style: `const char *CSP_RESTRICT line`, `const cif_tok *CSP_RESTRICT t`.

Leave it off — and say why in a short comment — when:

- two parameters may legitimately point at the same object or into the same
  buffer (in-place transforms, overlapping slices of one allocation);
- pointers are swapped or re-pointed at each other's targets during the
  function's lifetime (the lexer's `save`/`lx.p` pair is close to this: it is
  written as one buffer read forwards, and that is what makes it safe);
- the pointee can be reached by another route inside the same call — through a
  struct field also passed separately, a recursive call, or a callback;
- it points at perl-managed memory the perl API may alias behind your back: an
  `SV`/`AV`/`HV` internal, a `PV` that may be COW or a shared hash key, or a
  buffer an intervening perl call can realloc.

The whole-file buffer is the documented exception to that last point: the
comment at the head of `parse_buf()` argues it, and the argument is that the
buffer is only ever read, never written through a second pointer, while
everything the parse allocates is a fresh SV. Do not extend the exception past
that without writing down why.

Do not churn existing declarations solely to add or remove the qualifier; apply
this when writing new code or when already editing the function.

## Numbers must be exact, and must stay exact at every NV width

The coordinate reader is the module's one piece of numerical work, and its
correctness argument is written at the head of `str2nv_fixed()`. Read it before
changing anything it touches.

- `str2nv_fixed()` does not approximate `strtod()`: it returns the same NV, bit
  for bit, for `[+-]?digits[.digits]`, and the reasoning is load-bearing. The
  mantissa is an exact integer in an NV, the fifteen-digit cap keeps it inside
  `2**53` (exact even on the narrowest NV perl configures), `10**frac` is exact
  to `10**22`, and an IEEE 754 division of two exact operands is correctly
  rounded — which is what `strtod()` is required to return. Dividing by
  `10**frac` rather than multiplying by `10**-frac` is not a style choice:
  `1e-3` is not representable and multiplying would round twice.
- Changing `STR2NV_MAX_DIGITS`, the `str2nv_pow10[]` table, the leading-zero
  rule, or the divide is a change to that argument. `t/numbers.t` is what holds
  you to it, through `_str2nv_paths()`, which calls both C paths on one string.
  It cannot be checked against perl's own conversion: the two paths are never
  both reached for a real field, and `Perl_my_atof` on an older
  `-Duselongdouble` perl is a hand-rolled decimal accumulator rather than
  `strtod()`.
- Anything `str2nv_fixed()` declines falls through to `str2nv_slow()`, and that
  goes through `STR2NV_STRTOD()`, which picks `Strtod`, `strtoflt128`,
  `strtold` or `strtod` for the build. Never call `strtod()` directly: on a
  wide-NV perl it silently parses into 53 bits and `0.60` comes back as
  `0.599999999999999978`.
- **Never call libm bare.** This file calls none today, and that is worth
  keeping: `sqrt(x)`, `log(x)` and friends take a `double`, so on a
  long-double or `__float128` build they round the NV down, compute there, and
  convert back — no warning, just a less accurate answer than the perl running
  it. Anything new that needs one gets a width-correct wrapper block
  (`nv_sqrt`, `nv_log`, … as `Stats::LikeR` has) plus a link probe in
  `Makefile.PL`, not a bare call. The same goes for classification: not
  `isnan()`/`isfinite()` bare, and not `Perl_isnan`/`Perl_isinf`, which on
  perls before 5.22 can route through a `Perl_fp_class()` block in `perl.h`
  that has never compiled.
- Format NVs and IVs with perl's `my_snprintf` and `NVgf`/`NVff`/`IVdf`/`UVuf`,
  as `chain_tally()` does — never the C library `snprintf` with a hand-written
  `%g`, and never `%zu`/`%lld`, which MSVC's older CRT does not take.
- Do not assume a literal is exactly representable and do not hardcode a
  `double`-sized epsilon; scale a tolerance off `NV_EPSILON`.

## Memory, and the perl API before libc

- Allocate with `Newx`/`Newxz`/`Renew`/`Safefree`, never `malloc`/`free`, so
  the module uses whatever allocator its perl was built with. Free on the
  croak path too: `t/leaks.t` checks that a failed open does not leak the
  buffer it had already allocated.
- `t/leaks.t` is not optional coverage. XS that builds Perl data structures is
  exactly where reference counts go wrong, and a leak of a few hundred SVs per
  file is invisible on one structure and fatal on a directory of twenty
  thousand. It also checks for cycles: a residue that pointed back at its
  chain would never be freed and no single-file test would notice.
- Reach for the perl API before libc whenever perl has an equivalent
  (`my_snprintf`, `strEQ`/`foldEQ`, `Zero`, `PerlIO`, `PerlProc_*`): it is the
  layer already ported to every target below. `slurp()` is the one place that
  is not honoured — it reads with `fopen`/`fread`/`ferror` — so if you touch
  that function, move it to `PerlIO_open`/`PerlIO_read`/`PerlIO_error` rather
  than adding more stdio around it.

## Tests come from other people's readers, not from this one

If a claim about the format has an equivalent in a reader someone else wrote,
the test for it comes from **that reader's suite or from its answer on a real
file**, not from what this parser currently does. Invented cases confirm the
code does what it does; another reader's cases pin it to what the format means.

The two readers used here, and where they are on this machine:

- **gemmi** (reads PDB and mmCIF in C++, is what the crystallography programs
  use): gemmi 0.7.5 is importable from
  `/home/con/.pyenv/versions/3.14.2/bin/python3`. The default `python3` cannot
  import it, so point `STRUCTURE_INFO_PYTHON` at that interpreter to compare
  against gemmi live; without it `t/oracle.t` compares against gemmi's frozen
  answer instead. Its `tests/` directory is thirty years of the format's bad
  behaviour collected by people who had to read it too; `t/data/pdb1gdr.ent`
  came from there.
- **Biopython** 1.85, `/home/con/.local/lib/python3.12/site-packages/Bio` —
  `parse_pdb_header` and `pdb-seqres` are the header-side comparison, and
  `Tests/PDB/` is the file corpus (not checked out here: cases taken from it
  are written into `t/foreign.t` as text).
- **Real structures**: PDBbind v2020, 10,116 entries, at
  `/home/con/ui/pepPriML/PPB/PDB/PDBbind.v2020`, which `t/real.t` and
  `t/real_cif.t` find by default; `STRUCTURE_INFO_TEST_DIR` and
  `STRUCTURE_INFO_TEST_CIF_DIR` point them elsewhere, and
  `STRUCTURE_INFO_TEST_ALL` reads the whole set rather than a spread.

### Neither half may skip

- `t/oracle.t` **runs** the comparison — every atom of model 1 against gemmi,
  as a multiset, so a disagreement says which atoms. gemmi is a tool the
  distribution cannot require, so the comparison does not depend on having it:
  what gemmi answered for every structure in `t/data` is frozen in
  `t/data/oracle.txt`, and that is what the test reads. It must not `skip_all`.
  It is also the only file allowed to shell out, and it does so only where
  gemmi is importable, to two further ends: checking that the frozen answer is
  still what gemmi says — a fixture regenerated without re-running the
  generator would otherwise leave the comparison passing against an answer to
  an older file — and extending the comparison to a spread of
  `STRUCTURE_INFO_TEST_DIR` and `STRUCTURE_INFO_TEST_CIF_DIR`.
- `t/data/oracle.txt` is generated: `t/data/oracle.pl` writes it by running
  `t/data/oracle.py`, the one copy of the dump, which `t/oracle.t` runs too so
  that the frozen answer and the live one cannot become answers to different
  questions. Re-run it from `t/data/` after adding or changing a structure
  there, with `STRUCTURE_INFO_PYTHON` set, and commit what it wrote; the test
  fails if the directory holds a structure the file has no answer for.
- `t/foreign.t` is where what that comparison *found* is written down: a
  hardcoded fixture, the assertion, and the file the case came out of, named so
  it can be looked up in the suite that first thought of it. It must never
  skip, never call `python3`, and never need a structure that is not in
  `t/data`. A cross-validation test that skips is a test that silently never
  runs on a CPAN smoker.
- Two known differences are documented in `t/oracle.t` and are not failures: a
  residue modelled in two chemical states at once is one residue here and two
  in gemmi, and a file with no CRYST1 has no cell here and gemmi's default
  1×1×1 cell there. Add to that list only with the reasoning, never to quiet a
  disagreement.
- `t/real.t`/`t/real_cif.t` may skip for want of structures — that is data, not
  a reference implementation — but what they check must be self-evident from
  the file (the chains account for every atom, the residues for every chain
  atom, the sequence is as long as the polymer).
- Fixtures are generated, not typed: `t/data/generate.pl` writes them, because
  a PDB record is fixed-column and hand editing shifts a column silently. Run
  it from `t/data/` and commit both the script and what it wrote.
- Record provenance in the file's header comment: which reader, which version,
  which directory or function, and which real entry the case came out of.
  `t/foreign.t`, `t/classify.t` and `t/oracle.t` are the models.

### Maximum coverage

Aim to exercise every path, not one representative call:

- Cross the argument space: every option (`model`, `altloc`, `hydrogens`,
  `waters`, `hetatm`, `anisou`, `meta`, `lineno`, `atom_hashes`, `chains`) and
  combinations of them, not the default; both formats for each; the file, the
  string and the gzipped forms of every entry point.
- Cover the format's edges the way the archive does: an empty file, a file with
  no CRYST1, the pre-1996 entries that keep the id in columns 73–80, a
  right-trimmed line that stops before the element or B-factor column, an
  insertion code, a negative or zero residue number, a `*****` serial that
  overflowed its five columns, a blank chain id, microheterogeneity in altlocs,
  an NMR ensemble, a HETATM-only entry, a single-ion chain, a residue whose
  three letters mean two different things.
- Cover the Perl surface too: each accepted call form, each documented return
  field, every `croak` message and argument-validation path, and the leak
  checks in `t/leaks.t`.
- Justify every tolerance in a comment against what was actually observed —
  `t/structure.t`'s coordinate comparison says the observed difference is 0 on
  all three NV widths, and what the `1e-9` is really there for. Never widen one
  to make a failure go away.
- Gate genuinely slow cases behind `STRUCTURE_INFO_TEST_ALL`,
  `EXTENDED_TESTING` or `AUTHOR_TESTING` rather than dropping them.

## Every change must hold across the whole support matrix

A modification is not finished when it works on the default perl.

### All installed perls

`./test.all.perls.pl` builds and tests against every perlbrew perl; run it for
anything touching `Parser.xs` (`--no-install` to leave the installed copies
alone, `-p` for a subset). The local matrix under
`/home/con/perl5/perlbrew/perls/`, by `$Config{nvtype}`:

- `perl-5.44.0` (default), `perl-5.42.3` (threaded), `perl-5.10.1` — `double`
- `perl-5.12.5` — `long double` (archname `x86_64-linux-ld`)
- `5.44.0-quadmath` — `__float128`

The fast check, which does not clobber the current `Makefile`: generate the
`.c` with each perl's `ExtUtils::ParseXS` and compile it against that perl's
`CORE` directory. It is not a substitute for running the suite.

`-std=c99` (not `gnu99`) is the closest local proxy for a vendor compiler:
compile the generated `.c` with it when in doubt. It must stay clean.

### Windows, Solaris and every BSD

There is no local Windows, Solaris, illumos or BSD perl, so this is discipline
applied while writing, not something a run here will catch.

- Keep the `_GNU_SOURCE` / `__EXTENSIONS__` block at the top of `Parser.xs`
  intact and add to it rather than around it: `__EXTENSIONS__` is what exposes
  the POSIX and XPG symbols on Solaris once `-std=c99`-style strictness is on.
- No GNU extensions: no VLAs, no nested functions, no statement expressions,
  no `typeof`, no `__builtin_*` without a fallback, no zero-length arrays, no
  case ranges, no `void *` arithmetic.
- No POSIX-only headers or functions: no `<unistd.h>`, `<sys/time.h>`,
  `<strings.h>`; no `strcasecmp`/`strncasecmp` (`cif_iskw()` is the
  case-insensitive compare this file uses), no `fork`, no hardcoded `/` in a
  path.
- No glibc-only or BSD-only libc: no `qsort_r`, `memmem`, `strchrnul`,
  `strndup`, `asprintf`, `getline`, `fmemopen`, `reallocarray`, `arc4random`.
  `memchr` and `memcmp` are C89 and are what the line and keyword scans use.
- Printf lengths must be perl's: `%" UVuf "`, `%" IVdf "`, `%" NVgf "`.
- Assume nothing about `char`'s signedness, struct field order, or alignment.
  SPARC and some ARM/BSD combinations fault on unaligned access, so never cast
  a `char *` into the file buffer to a wider pointer type and dereference it —
  the packed keys are built by shifting bytes for exactly this reason.

### Back-compatible to perl 5.10

`dist.ini` declares `perl = 5.010`, `Makefile.PL` sets
`MIN_PERL_VERSION => '5.010'`, and `lib/Chem/Structure/Parser.pm` says
`require 5.010`. Both sides must honour it.

- Perl code (module and tests) is limited to 5.10 syntax: no signatures, no
  postfix dereference, no `say` without the feature, no `isa`, no chained
  comparisons. Tests follow the existing header — `require 5.010; use strict;
  use warnings FATAL => 'all';` — not `use 5.044`. The author-only helpers in
  the distribution root (`test.all.perls.pl`, `structure.info.pl`,
  `benchmark.pl`) are exempt: they are not shipped (`PruneFiles` in `dist.ini`
  and the `MY::libscan` in `Makefile.PL` keep them out) and use modern perl.
- Any perl C API newer than 5.10 comes through `ppport.h`, and `ppport.h` is
  regenerated when a new one is used. If it does not cover the call, write the
  fallback rather than raising the minimum.
- Do not rely on behaviour that changed after 5.10 (hash iteration order,
  `SvPV` return guarantees, COW string semantics) without checking on
  `perl-5.10.1`.
- Raising `MIN_PERL_VERSION` is a maintainer decision, not a way to make a
  build error go away.

### `Makefile.PL` and `dist.ini` are two copies of one thing

The shipped `Makefile.PL` is checked in, and `dist.ini`'s
`[MakeMaker::Awesome]` header carries the same logic for the released
`Makefile.PL`. Change both in the same edit — the `-std=gnu99` probe and the
`MY::libscan` that keeps the root `*.pl` helpers out of `$(INST_LIB)` exist in
both places.

## Comments

`~/.claude/CLAUDE.md` is the doctrine and it applies as written. What it means
here:

- `Parser.xs` uses the house block form — `/*` with text on the same line,
  continuation lines at column 0 with no `*` gutter, `*/` closing the last line
  of text — and `//` for a single fact about one line, which is the right form
  for annotating an enumerated value at its declaration (`RT_*`, `CT_*`, the
  `C_*` column names). C99 line comments are fine: `Makefile.PL` compiles with
  `-std=gnu99`.
- Write `--` for an em dash in code comments and in POD. `README.md` is
  markdown prose and uses real em dashes; leave it that way.
- Indent with tabs, in both the XS and the Perl. Spaces appear only to align a
  continuation line of a parameter list, which is alignment and not
  indentation. No trailing whitespace.
- Provenance belongs in the header comment of the file that depends on it: a
  case from another reader's suite, a number from a benchmark, an
  externally-imposed column range. "gemmi's `tests/pdb1gdr.ent`, a 1993 entry"
  is provenance; "test file" is not.
- Keep them true. A comment that has drifted is worse than none because it is
  trusted; update it in the same edit as the code or delete it. Do not leave
  commented-out code behind.
