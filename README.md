# Synopsis

Read a molecular structure file and get everything in it back as a hash of
hashes — the header, the annotations, every chain, every residue, every atom,
and the single-letter sequence of each chain — in one call.

```perl
use Structure::Info;

my $info = structure_info('1a22.ent.pdb');

print $info->{id};                                  # 1A22
print $info->{resolution};                          # 2.6
print $info->{chains}{A}{sequence};                 # FPTIPLSRLFDNAMLRAHRLHQL...
print $info->{chains}{A}{molecule};                 # GROWTH HORMONE
print $info->{chains}{A}{residues}{54}{resname};    # LEU
print $info->{chains}{A}{residues}{54}{atoms}{CA}{x};

print structure_summary($info);
```

```
1A22  HUMAN GROWTH HORMONE BOUND TO SINGLE RECEPTOR
  file        1a22.ent.pdb
  method      X-RAY DIFFRACTION
  resolution  2.6 A
  R / R-free  0.187 / -
  models      1
  atoms       3113 (69 hetatm, 69 water)
  chain A     protein      206 residues,  1492 atoms, 2 gaps
              FPTIPLSRLFDNAMLRAHRLHQLAFDTYQEFEEAYIPKEQKYSFLQNPQTSLCFSESIPTP...
              GROWTH HORMONE
  chain B     protein      235 residues,  1621 atoms, 2 gaps
              PKFTKCRSPERETFSCHWTLGPIQLFYTRRNTQEWTQEWKECPDYVSAGENSCYFNSSFTS...
              GROWTH HORMONE RECEPTOR
```

The coordinate section is parsed in C, because across a directory of
structures it is millions of lines: the largest entry in PDBbind v2020 is
33 MB and 411,648 atom records, and it reads in about 1.5 seconds. The header
records are parsed in Perl, because they are irregular and there are only a
few dozen of them in a file.

The module is called `Structure::Info` and not `PDB::Info` because the shape
of what it hands back has nothing to do with the format it came out of. Today
it reads PDB; `formats()` says what it reads at any moment, and a format it
knows the name of but cannot read yet says so rather than misreading it.

# Installing

    perl Makefile.PL
    make
    make test
    make install

`make test` reads the fixtures in `t/data`. If a directory of real structures
is to hand it reads a sample of those too; point it somewhere with

    STRUCTURE_INFO_TEST_DIR=/path/to/structures make test
    STRUCTURE_INFO_TEST_ALL=1 STRUCTURE_INFO_TEST_DIR=/path make test   # all of them

With no such directory those tests skip, so the distribution builds on a
machine with no structures on it.

# Getting help

`h` prints any function's section of this document to `STDOUT` and returns, in
the spirit of R's `?function` at the prompt. It takes the name three ways:

    h('structure_info');    # by name
    h(*res_type);           # by name, unquoted
    h(\&aa3to1);            # by reference
    h();                    # the list of documented functions

    perl -MStructure::Info -e 'h(*structure_info)'   # straight from the shell

Note that `h(res_type)`, with no quotes and no sigil, cannot be made to work:
every function here is exported, so Perl parses the bareword as a call to
`res_type()` before `h` is ever reached. Use one of the three forms above.

# Functions/Subroutines

## structure_info

    my $info = structure_info($file, %options);

Reads `$file` and returns a hash reference. The format is worked out from the
file name, and from the first records in the file when the name gives nothing
away. `.gz` files are read as they are, without unpacking to a temporary file.

### What comes back

    file        the path it was read from
    format      'pdb'
    id          the four-character id, from HEADER or from the file name
    title       TITLE, continuation lines joined
    header      { classification, deposit_date, id_code }
    experiment  [ 'X-RAY DIFFRACTION' ]
    resolution  2.6                        from REMARK 2
    r_work      0.187                      from REMARK 3
    r_free      0.278
    temperature 100                        from REMARK 200
    ph          7.5
    keywords    [ ... ]                    KEYWDS
    authors     [ ... ]                    AUTHOR
    journal     { auth => [...], titl, ref, pmid, doi }
    compound    { 1 => { molecule, chain => [...], engineered, ... } }
    source      { 1 => { organism_scientific, organism_taxid, ... } }
    seqres      { A => { sequence, residues => [...], length } }
    dbref       { A => [ { database, accession, db_id, ... } ] }
    seqadv      [ ... ]                    where the entry differs from the database
    modres      { MSE => { standard => 'MET', comment } }
    het         { NAG => { name, formula, instances => [...] } }
    helix sheet ssbond link cispep         [ { ... }, ... ]
    cryst1      { a, b, c, alpha, beta, gamma, sgroup, z }
    remarks     { 2 => [ lines ], 350 => [ lines ], ... }
    conect      [ [ serial, serial, ... ], ... ]
    records     { REMARK => 340, SEQRES => 24, ... }   every record, counted
    n_models    how many MODEL records the file has
    model       which one the chains below were built from
    models      every model, when called with model => 'all'
    chains      { A => { ... } }
    chain_order [ 'A', 'B' ]               the order the file has them in
    stats       { n_atoms, n_hetatm, n_hydrogens, n_water_atoms,
                  elements => { C => 1946, N => 507, ... },
                  bfactor => { min, max, mean, n }, bbox, center }

A record that is not in the file reads as `undef`, and a list that is not in
the file reads as an empty arrayref — `title` being `undef` means there was no
TITLE, which is a different thing from a TITLE that was blank.

Everything the module does not take apart is still in `remarks` and in the
raw record counts, so nothing in the file is lost.

### A chain

    id            'A'
    type          protein | dna | rna | water | hetero | unknown
    sequence      single-letter, the residues that have coordinates
    seqres        single-letter, what SEQRES says was in the crystal
    n_residues n_polymer n_water n_ligand n_atoms n_hetatm
    n_missing     the SEQRES length less the residues actually modelled
    gaps          [ { after, before, missing } ]   unmodelled stretches
    n_gaps
    first last    the first and last polymer residue keys
    residue_types { amino_acid => 191, water => 69, ligand => 1 }
    molecule      from COMPND, e.g. 'GROWTH HORMONE'
    organism      from SOURCE
    ec fragment mol_id dbref
    residue_order [ '1', '2', '3', '3A', ... ]      file order
    residues      { '54' => { ... } }

`sequence` and `seqres` are the two different questions people mean by "the
sequence": what was modelled, and what was in the crystal. They differ
wherever a terminus or a loop went unmodelled, which is what `gaps` counts and
`n_missing` totals.

### A residue

Keyed by residue number with its insertion code appended, so `100`, `100A` and
`100B` are three separate keys and nothing is silently overwritten.

    resname   'LEU'          number 54       icode ''      key '54'
    one       'L'            the single-letter code, '' when there is none
    type      amino_acid | nucleotide | water | ligand | ion
    standard  1 for the twenty and the standard nucleotides
    modified  1 for a residue like MSE, which is still an M in the sequence
    hetero    1 when it was written as HETATM
    free      1 for a free amino acid bound in a site (see below)
    n_atoms   center [ x, y, z ]   b_mean
    atom_order [ 'N', 'CA', 'C', 'O', ... ]
    atoms      { CA => { ... } }

### An atom

    name 'CA'   serial 384   element 'C'   charge ''
    x y z       occupancy    bfactor       altloc   hetero
    altlocs   [ { altloc, x, y, z, occupancy, bfactor }, ... ]
              present only when the atom has alternate conformers; every
              conformer is listed, including the one chosen above

Nothing in the structure points back up at its parent — a residue does not
hold its chain, an atom does not hold its residue. Parent links would make the
whole thing one reference cycle, and a cycle is a leak that goes unnoticed
until the ten-thousandth file.

### Options

    model     => 1          which MODEL to build chains from; 'all' fills in
                            {models} as well.  Default 1, which is also the
                            right answer for a file with no MODEL records
    altloc    => 'first'    which alternate conformer's coordinates win;
                            'highest' takes the highest occupancy instead
    hydrogens => 1          keep hydrogens and deuteriums
    waters    => 1          keep waters
    hetatm    => 1          keep HETATM records
    atoms     => 1          build the atom hashes; 0 stops at the residue
                            level, which is much smaller and faster
    meta      => 1          parse the header records
    anisou    => 0          keep ANISOU lines
    chains    => ['A','B']  read only these chains
    format    => 'pdb'      skip the format detection

Every option is checked. A misspelled one is fatal, because an ignored typo is
a wrong answer that arrives without a word: `hydrogen => 0` that is quietly
dropped gives a structure with the hydrogens still in it and no clue why.

For a very large structure the options are the difference between a hash of
hashes that fits in memory and one that does not. The largest entry in PDBbind
v2020 is 2wy2: 33 MB, 64 models, 411,648 atom records.

    structure_info($f)                # model 1 only    50 MB    0.07 s
    structure_info($f, model => 'all')                 711 MB    1.4 s
    structure_info($f, model => 'all', atoms => 0)     418 MB    0.9 s

Filtering happens in the C, before a hydrogen or a water has become a Perl
value, so `hydrogens => 0` is cheaper than reading them and throwing them away.

## pdb_info

    my $info = pdb_info($file, %options);

`structure_info()` with the format settled in advance. Use it when the file is
known to be a PDB whatever it happens to be called.

## structure_info_string

    my $info = structure_info_string($text, %options);

The same, for a structure already in a string. A string has no name to go on,
so text that looks like nothing in particular is read as PDB; text that looks
like another format still gets a straight answer about it.

## structure_atoms

    my $atoms = structure_atoms($info);
    my $atoms = structure_atoms($info, 'A');

Every atom as a flat array of hash references, in file order, each carrying
the `chain`, `resname`, `resseq`, `icode` and `reskey` it came from. This is
the shape to hand to a distance calculation or to write out as a table; the
nested form is the shape to look things up in. The hashes are copies, so
writing to them does not scribble on the structure.

## structure_residues

    my $residues = structure_residues($info);
    my $residues = structure_residues($info, 'A');

Every residue in file order. These are the same hash references that are in
the nested structure, not copies, so walking them and looking one up agree.

## structure_ligands

    my $lig = structure_ligands($info);     # { 'NAG_A_301' => { ... } }

The heterogens that are neither water nor part of the polymer, keyed by
residue name, chain and number — which is what a binding-site table wants as
its row label.

## structure_sequences

    my $seq = structure_sequences($info);   # { A => 'FPTIPLSRL...' }

The observed single-letter sequence of every chain that has one.

## chain_sequence

    my $obs = chain_sequence($info, 'A');
    my $all = chain_sequence($info, 'A', 'seqres');

One chain's sequence: `observed` is the residues that have coordinates,
`seqres` is what SEQRES says was in the crystal.

## structure_summary

    print structure_summary($info);

A paragraph a person can read: id, title, method, resolution, models, atom
counts, and a line per chain with its type, size, sequence and molecule. The
example at the top of this document is its output.

## aa3to1

    aa3to1('ALA');    # 'A'
    aa3to1('MSE');    # 'M'   selenomethionine is still a methionine
    aa3to1('HOH');    # ''    water is not an amino acid
    aa3to1('NAG');    # ''

The single-letter code of an amino acid, and the empty string for anything
that is not one. Leading and trailing blanks and case do not matter, because
the name usually arrives straight out of columns 18 to 20.

Modified residues map to the residue they were made from — `MSE` to `M`, `SEP`
to `S`, `HYP` to `P`, the D-amino acids to their L partners — because a
structure that soaked in selenomethionine has the same sequence as one that
did not, and a sequence with an `X` every seventh position is no use to
anyone.

## res1

    res1('ALA');      # 'A'
    res1(' DA');      # 'A'   deoxyadenosine
    res1('PSU');      # 'U'   pseudouridine
    res1('HOH');      # ''

`aa3to1()` widened to nucleotides, which is what building a sequence wants
when the chain might be DNA or RNA.

## res_type

    res_type('ALA');  # 'amino_acid'
    res_type('DA');   # 'nucleotide'
    res_type('HOH');  # 'water'
    res_type('NAG');  # 'other'

What kind of residue a name is. `other` covers ligands, ions and sugars;
`structure_info()` narrows those to `ligand` or `ion` once it can see how many
atoms the residue has and what they are.

## formats

    my @can = formats();          # ('pdb')
    my $all = formats();          # every format known, supported or not

## h

Prints a function's documentation to STDOUT. See *Getting help* above.

# Two residues that are not what they look like

Both of these were found by running the module over PDBbind v2020 and asking
where the sequence it read disagreed with SEQRES. Both are in the test suite.

**A free base is not a nucleotide.** `ADE`, `CYT`, `GUA`, `THY` and `URI` mean
one thing in a file written before 2007 — the nucleotides of a nucleic acid
chain — and another in a file written since: a free base sitting in an active
site as a ligand. The sugar tells them apart, since a nucleotide has a `C1'`
and a free base has nothing but the base. Without that check the guanine bound
to 1czc reads as a nucleotide and turns up as a `G` on the end of a
396-residue protein sequence.

**A free amino acid is not part of the chain.** A HETATM residue with an amino
acid's name is a modified residue when it is numbered among the polymer — the
MSE that replaced a methionine belongs in the sequence — and a free amino acid
bound in a site when it is numbered out with the ligands, in which case it
does not. 3lms has a glycine at A501, two hundred residues past the end of a
chain whose SEQRES is 309 long. Those are flagged `free => 1` and typed as
ligands.

Neither is a rule the format states; both are what the format means.

# What the C is worth

`benchmark.pl` measures it rather than asserting it. Against a pure-Perl
reader that slices the same columns in the same order and builds the same
nested hash — not a straw man, just the same job done per line in Perl —
reading 200 structures from PDBbind v2020, best of three:

    structure_info                           6.30 s     168,000 atoms/s
    structure_info, atoms => 0               4.05 s     262,000 atoms/s
    the XS parse alone                       1.97 s     538,000 atoms/s
    the same reader, in pure Perl            5.48 s     193,000 atoms/s
    pure Perl, and the same statistics       7.14 s     148,000 atoms/s

**The parse is 2.8x the same parse written in Perl.** That is what the C buys.

`structure_info()` as a whole lands at 0.9x the plain Perl reader and 1.1x the
one that also gathers the statistics — while doing a good deal more than
either: the header records, SEQRES, the gaps between the declared sequence and
the modelled one, chain typing, ligand and ion classification. Building a hash
of hashes out of what was read costs what it costs in any language, because
the hashes are the answer. Reading is the part worth writing in C, and reading
is the part that got faster.

Two rounds of that measurement changed the code. The first showed the whole
read coming out *slower* than the naive Perl reader, because every per-atom
statistic was being gathered in Perl over data C had just finished reading;
those moved into the parse. The second showed each atom being built twice —
once as columns in C, once as a hash in Perl — and the parse now builds the
hash itself, which took another 25% off both the time and the memory.

# What is written in C, and what is not

The C side makes one pass over the bytes and does three things: it splits
ATOM/HETATM records into parallel per-field arrays, it marks where each
residue begins and ends, and it groups every other record by record name for
Perl to take apart. Residue name lookup — three letters to one letter, and the
amino acid/nucleotide/water question — is a switch on three packed bytes, and
one table serves `aa3to1()`, `res1()` and `res_type()` so the three can never
disagree.

The per-atom result is columnar — one array per field rather than one hash per
atom — which is both faster to build and about an order of magnitude smaller
while it is being built. Perl turns it into the nested hash afterwards, and
only for the atoms that were kept.

Everything else is Perl. The header records are irregular, they are a few
dozen lines per file rather than hundreds of thousands, and they are where the
next surprise will turn up; none of that is worth writing in C.

Two details in there are worth knowing about, because they are where a PDB
reader is usually wrong:

**Coordinates.** `x`, `y` and `z` are adjacent eight-column fields with no
separator between them, so a reader that scans for a number instead of slicing
the columns runs `-12.500 111.000` together. They are sliced, and converted
with `strtod` on the slice.

**Elements.** Columns 77-78 hold the element symbol, and when they are missing
— older files, and files written by programs that should know better — it has
to come from the atom name. A one-letter element is right-justified from
column 14, so ` CA ` is a carbon alpha, and a two-letter element starts in
column 13, so `CA  ` is a calcium. A hydrogen with a long name like `HG11`
also starts in column 13, and taking its first two letters makes it mercury;
what separates them is that it fills all four columns, which no two-letter
element's name does.

# Author

David E. Condon <dec986@gmail.com>

# License

This library is free software; you can redistribute it and/or modify it under
the same terms as Perl itself.
