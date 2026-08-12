# Synopsis

Read a molecular structure file and get everything in it back as a hash of
hashes — the header, the annotations, every chain, every residue, every atom,
and the single-letter sequence of each chain — in one call.

```perl
use Chem::Structure::Parser;

my $info = structure_info('1a22.ent.pdb');

print $info->{id};                                  # 1A22
print $info->{resolution};                          # 2.6
print $info->{chains}{A}{sequence};                 # FPTIPLSRLFDNAMLRAHRLHQL...
print $info->{chains}{A}{molecule};                 # GROWTH HORMONE
print $info->{chains}{A}{residues}{54}{resname};    # PHE
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

The module is called `Chem::Structure::Parser` and not `PDB::Info` because the
shape of what it hands back has nothing to do with the format it came out of.
It reads PDB and mmCIF/PDBx; `formats()` says what it reads at any moment, and
a format it knows the name of but cannot read yet says so rather than
misreading it.

# PDB and mmCIF

Reading is the same call either way. `structure_info()` works out the format
from the file name — `.pdb`, `.ent`, `.cif`, `.mmcif`, `.pdbx` — and from the
first records in the file when the name gives nothing away, and the hash that
comes back has the same keys, the same nesting and the same values whichever
it was.

```perl
my $a = structure_info('1a22.pdb');
my $b = structure_info('1a22.cif');

$a->{chains}{A}{sequence} eq $b->{chains}{A}{sequence};              # true
$a->{chains}{A}{residues}{54}{atoms}{CA}{x}
	== $b->{chains}{A}{residues}{54}{atoms}{CA}{x};                  # true
```

So no calling code branches on the format, and a script written against a
directory of `.pdb` files works unchanged on a directory of `.cif` ones.

Equality here means equality rather than approximately: `t/cif.t` reads
fixture pairs both ways and compares the whole coordinate half of the returned
structure with `is_deeply`, and `t/real_cif.t` converts real entries from the
PDB archive into mmCIF and asserts that every chain, residue, atom and count
comes back identical.

Two consequences are worth knowing.

**The identifiers are the auth_\* ones.** An mmCIF file carries two sets: the
`label_*` identifiers the archive assigns, and the `auth_*` ones the depositor
used. Only `auth_*` matches what the PDB record carried, so those are the
chain ids and residue numbers used throughout — in the coordinates and in the
annotations alike. A structure read from a `.cif` therefore has the same chain
`A` and the same residue `54` as the same structure read from a `.pdb`, not
the `label_asym_id` lettering that runs on through the waters.

**Values are converted, not passed through.** Where the two formats spell the
same fact differently, the mmCIF reader produces what the PDB reader would
have: `_atom_site.pdbx_formal_charge` of `-1` reads back as `'1-'`, and `.`
and `?` — mmCIF for "not applicable" and "unknown" — read back as the empty
field a PDB record would have had. A charge of `0` is kept as `'0'`, because
"the field said zero" and "the field was blank" are different answers and both
formats can say either.

What is *not* the same is what only one of the formats has. An mmCIF file has
no REMARK records, so `$info->{remarks}` is empty for one; a PDB file has no
`_entity` category, so a chain read from one may not know which entity it
belongs to. Every key is present in both cases, so reading one is a test of
what the file said and never of which format it was.

Where the same fact is filed under different names, it is folded into the same
key:

| `$info` key | PDB record | mmCIF category |
| --- | --- | --- |
| `title` | `TITLE` | `_struct.title` |
| `id` | `HEADER` | `_entry.id` |
| `experiment` | `EXPDTA` | `_exptl.method` |
| `resolution` | `REMARK 2` | `_refine.ls_d_res_high` |
| `r_work`, `r_free` | `REMARK 3` | `_refine.ls_R_factor_R_*` |
| `keywords` | `KEYWDS` | `_struct_keywords.text` |
| `authors` | `AUTHOR` | `_audit_author` |
| `journal` | `JRNL` | `_citation`, `_citation_author` |
| `compound`, `source` | `COMPND`, `SOURCE` | `_entity`, `_entity_src_*` |
| `seqres` | `SEQRES` | `_entity_poly`, `_entity_poly_seq` |
| `het` | `HET`, `HETNAM`, `FORMUL` | `_chem_comp`, `_pdbx_nonpoly_scheme` |
| `helix` | `HELIX` | `_struct_conf` |
| `sheet` | `SHEET` | `_struct_sheet_range` |
| `ssbond`, `link` | `SSBOND`, `LINK` | `_struct_conn` |
| `cispep` | `CISPEP` | `_struct_mon_prot_cis` |
| `modres` | `MODRES` | `_pdbx_struct_mod_residue` |
| `dbref` | `DBREF` | `_struct_ref`, `_struct_ref_seq` |
| `cryst1` | `CRYST1` | `_cell`, `_symmetry` |
| `n_models` | `MODEL` | `_atom_site.pdbx_PDB_model_num` |

# Installing

    perl Makefile.PL
    make
    make test
    make install

`make test` reads the fixtures in `t/data`. If a directory of real structures
is to hand it reads a sample of those too; point it somewhere with

    STRUCTURE_INFO_TEST_DIR=/path/to/pdbs  make test
    STRUCTURE_INFO_TEST_CIF_DIR=/path/to/cifs  make test
    STRUCTURE_INFO_TEST_ALL=1 STRUCTURE_INFO_TEST_DIR=/path make test   # all of them

`STRUCTURE_INFO_TEST_DIR` is used twice: `t/real.t` reads those files as PDB
and checks them against a second reader written in plain Perl, and
`t/real_cif.t` converts each one into mmCIF and asserts that reading it back
gives the same structure to the last digit.
`STRUCTURE_INFO_TEST_CIF_DIR` takes a directory of real `.cif` files, either
flat or one subdirectory per structure, and reads those directly.

With no such directory those tests skip, so the distribution builds on a
machine with no structures on it.

# Getting help

`h` prints any function's section of this document to `STDOUT` and returns, in
the spirit of R's `?function` at the prompt. It takes the name three ways:

    h('structure_info');    # by name
    h(*res_type);           # by name, unquoted
    h(\&aa3to1);            # by reference
    h();                    # the list of documented functions

    perl -MChem::Structure::Parser -e 'h(*structure_info)'   # straight from the shell

Note that `h(res_type)`, with no quotes and no sigil, cannot be made to work:
every function here is exported, so Perl parses the bareword as a call to
`res_type()` before `h` is ever reached. Use one of the three forms above.

# Functions/Subroutines

## structure_info

    my $info = structure_info($file, %options);

Reads `$file` and returns a hash reference. The format is worked out from the
file name — `.pdb`, `.ent`, `.cif`, `.mmcif`, `.pdbx` — and from the first
records in the file when the name gives nothing away. `.gz` files are read as
they are, without unpacking to a temporary file.

### What comes back

Laid out the way `tree` lays out a directory, this is `1a22.ent.pdb` — a real
file, real values, the long lists cut short:

```
$info
├── file            '1a22.ent.pdb'          the path it was read from
├── format          'pdb'                   or 'mmcif'
├── id              '1A22'                  from HEADER, or from the file name
├── title           'HUMAN GROWTH HORMONE BOUND TO SINGLE RECEPTOR'
├── header
│   ├── classification  'COMPLEX (HORMONE/RECEPTOR)'
│   ├── deposit_date    '15-JAN-98'
│   └── id_code         '1A22'
├── experiment      [ 'X-RAY DIFFRACTION' ]
├── resolution      2.6                     REMARK 2
├── r_work          0.187                   REMARK 3
├── r_free          undef                   this entry does not report one
├── temperature     287                     REMARK 200
├── ph              6.5
├── keywords        [ 'COMPLEX (HORMONE-RECEPTOR)', 'PITUITARY HORMONE', ... ]
├── authors         [ 'A.M.DE VOS', 'M.ULTSCH' ]
├── journal
│   ├── auth        [ 'T.CLACKSON', 'M.H.ULTSCH', 'J.A.WELLS', 'A.M.DE VOS' ]
│   ├── titl        'STRUCTURAL AND FUNCTIONAL ANALYSIS OF THE 1:1 GROWTH...'
│   ├── ref         'J.MOL.BIOL.                   V. 277  1111 1998'
│   ├── refn        'ISSN 0022-2836'
│   ├── pmid        '9571026'
│   └── doi         '10.1006/JMBI.1998.1669'
├── compound                                COMPND, by MOL_ID
│   ├── 1
│   │   ├── mol_id      '1'
│   │   ├── molecule    'GROWTH HORMONE'
│   │   ├── chain       [ 'A' ]
│   │   ├── engineered  'YES'
│   │   └── mutation    'YES'
│   └── 2           { molecule 'GROWTH HORMONE RECEPTOR', chain [ 'B' ],
│                     fragment 'EXTRACELLULAR DOMAIN', engineered 'YES' }
├── source                                  SOURCE, by MOL_ID
│   └── 1           { organism_scientific 'HOMO SAPIENS', organism_common
│                     'HUMAN', organism_taxid '9606', mol_id '1',
│                     expression_system 'ESCHERICHIA COLI',
│                     expression_system_taxid '562' }
├── entity_of_chain                         COMPND and SOURCE, by chain
│   ├── A           { mol_id '1', molecule 'GROWTH HORMONE', fragment undef,
│   │                 ec undef, organism 'HOMO SAPIENS', taxid '9606',
│   │                 expressed_in 'ESCHERICHIA COLI' }
│   └── B           { ..., fragment 'EXTRACELLULAR DOMAIN' }
├── seqres                                  what SEQRES says was in the crystal
│   ├── A
│   │   ├── sequence    'FPTIPLSRLFDNAMLRAHRLHQLAFDTYQEFEEAYIPKEQKYSFLQ...'
│   │   ├── residues    [ 'PHE', 'PRO', 'THR', 'ILE', ... ]        191 of them
│   │   └── length      191
│   └── B               { sequence, residues, length 238 }
├── dbref
│   └── A           [ { database 'UNP', accession 'P01241',
│                       db_id 'SOMA_HUMAN', seq_begin '1', seq_end '191',
│                       db_begin '27', db_end '217', chain 'A' } ]
├── seqadv          [ { chain 'A', resseq '120', resname 'ARG',
│                       db_res 'GLY', db_seq '146', comment 'ENGINEERED' } ]
├── modres          { }                     no MSE-style residues here
├── het
│   └── HOH         { het_id 'HOH', formula '69(H2 O)', water 1 }
├── hetnam          { }
├── formul          { }
├── helix           [ { id '1', class '1', length '29',
│                       init_chain 'A', init_resname 'SER', init_resseq '7',
│                       end_chain 'A', end_resname 'TYR', end_resseq '35' },
│                     ... ]                                     12 of them
├── sheet           [ ... ]                                     12
├── ssbond          [ { chain1 'A', resseq1 '53',
│                       chain2 'A', resseq2 '165', length '2.02' }, ... ]  5
├── link            [ ]
├── cispep          [ ]
├── site            [ ]
├── cryst1          { a '67.7', b '67.7', c '228',
│                     alpha '90', beta '90', gamma '90',
│                     sgroup 'P 43 21 2', z '8' }
├── biological_assembly  [ 32 lines of REMARK 350, verbatim ]
├── revdat          [ { num '3', date '18-APR-18', id '1A22',
│                       type '1', what 'REMARK' }, ... ]
├── remarks                                 every REMARK, by number
│   ├── 2           [ '', 'RESOLUTION.    2.60 ANGSTROMS.' ]
│   ├── 350         [ ... ]                                     32 lines
│   └── ...         1, 3, 4, 100, 200, 280, 290, 300, 465, 470, 500
├── conect          [ [ 448, 1255 ], ... ]                      10
├── records                                 every record type, counted
│   ├── REMARK      365
│   ├── SEQRES      34
│   ├── HELIX       12
│   └── ...         AUTHOR, COMPND, CONECT, CRYST1, DBREF, SOURCE, SSBOND, ...
├── n_models        1                       how many MODEL records the file has
├── model           1                       which one the chains below are
├── models                                  there only with model => 'all'
├── stats
│   ├── n_atoms         3113
│   ├── n_hetatm        69
│   ├── n_hydrogens     0
│   ├── n_water_atoms   69
│   ├── n_lines         3605
│   ├── n_atom_records  3044    ATOM lines seen, whether kept or not
│   ├── n_hetatm_records 69     HETATM lines, likewise
│   ├── n_anisou        0
│   ├── n_skipped       0       coordinate lines the options threw away
│   ├── elements        { C 1946, O 643, N 507, S 17 }
│   ├── bfactor         { min '2.7', max '85.39', mean 30.83, n 3113 }
│   ├── bbox            { xmin '12.142', xmax '80.34', ymin '2.011', ... }
│   └── center          [ '46.241', '29.135', '134.559' ]
├── chain_order     [ 'A', 'B' ]            the order the file has them in
└── chains
    ├── A
    │   ├── id              'A'
    │   ├── type            'protein'   protein dna rna water hetero unknown
    │   ├── sequence        'FPTIPLSRLFDNAMLRAHRLHQLAFDTYQEFEEAYIPKEQ...'
    │   │                               single-letter, what has coordinates
    │   ├── seqres          'FPTIPLSRLFDNAMLRAHRLHQLAFDTYQEFEEAYIPKEQ...'
    │   ├── seqres_length   191
    │   ├── n_residues      206
    │   ├── n_polymer       180
    │   ├── n_water         26
    │   ├── n_ligand        0
    │   ├── n_atoms         1492
    │   ├── n_hetatm        26
    │   ├── n_missing       11          SEQRES less what was modelled
    │   ├── gaps            [ { after 129, before 136, missing 6 },
    │   │                     { after 148, before 154, missing 5 } ]
    │   ├── n_gaps          2
    │   ├── first           1           the first and last polymer residue keys
    │   ├── last            191
    │   ├── residue_types   { amino_acid 180, water 26 }
    │   ├── molecule        'GROWTH HORMONE'            from COMPND
    │   ├── organism        'HOMO SAPIENS'              from SOURCE
    │   ├── mol_id          '1'
    │   ├── dbref           [ { ... } ]     as in the top-level dbref
    │   │                                   ec and fragment are here too, in a
    │   │                                   chain whose file gives them
    │   ├── residue_order   [ '1', '2', '3', ... '574' ]        file order, 206
    │   └── residues                    keyed number + insertion code
    │       ├── 54
    │       │   ├── resname     'PHE'
    │       │   ├── number      54
    │       │   ├── icode       ''
    │       │   ├── key         '54'
    │       │   ├── chain       'A'
    │       │   ├── one         'F'     '' when there is no letter for it
    │       │   ├── type        'amino_acid'
    │       │   │                       nucleotide water ligand ion
    │       │   ├── standard    1       one of the twenty, or a standard base
    │       │   ├── modified    0       1 for MSE, still an M in the sequence
    │       │   ├── hetero      0       1 when it was written as HETATM
    │       │   ├── free                not here; 1 for a free amino acid
    │       │   │                       bound in a site (see below)
    │       │   ├── n_atoms     11
    │       │   ├── b_mean      22.55
    │       │   ├── center      [ 65.311, 17.127, 140.515 ]
    │       │   ├── atom_order  [ 'N', 'CA', 'C', 'O', 'CB', ... ]
    │       │   └── atoms
    │       │       ├── CA
    │       │       │   ├── name       'CA'
    │       │       │   ├── serial     450
    │       │       │   ├── element    'C'
    │       │       │   ├── charge     ''
    │       │       │   ├── x          '66.446'
    │       │       │   ├── y          '18.25'
    │       │       │   ├── z          '141.982'
    │       │       │   ├── occupancy  '1'
    │       │       │   ├── bfactor    '24.53'
    │       │       │   ├── altloc     ''
    │       │       │   ├── hetero     0
    │       │       │   └── altlocs    [ { altloc, x, y, z, occupancy,
    │       │       │                      bfactor }, ... ]
    │       │       │                  present only when the atom has
    │       │       │                  alternate conformers; every conformer
    │       │       │                  is listed, the chosen one included
    │       │       └── ...    N, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ
    │       └── ...            1 .. 191, then the waters at 512 .. 574
    └── B                      the same again: 235 residues, 1621 atoms
```

A record that is not in the file reads as `undef`, and a list that is not in
the file reads as an empty arrayref — `title` being `undef` means there was no
TITLE, which is a different thing from a TITLE that was blank.

Everything the module does not take apart is still in `remarks` and in the
raw record counts, so nothing in the file is lost.

### The two sequences

`sequence` and `seqres` are the two different questions people mean by "the
sequence": what was modelled, and what was in the crystal. They differ
wherever a terminus or a loop went unmodelled, which is what `gaps` counts and
`n_missing` totals — eleven residues of chain A above, in two stretches.

### How residues are keyed

By residue number with the insertion code appended, so `100`, `100A` and
`100B` are three separate keys and nothing is silently overwritten. Waters and
ligands are in `residues` alongside the polymer, which is why chain A above
has 206 residues to its 191-long SEQRES.

### Nothing points back up

A residue does not hold its chain and an atom does not hold its residue — the
names are there, `chain => 'A'` on the residue, but not the references. Parent
links would make the whole thing one reference cycle, and a cycle is a leak
that goes unnoticed until the ten-thousandth file.

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
    format    => 'pdb'      skip the format detection: 'pdb' or 'mmcif'
                            ('cif', 'pdbx' and 'ent' name the same two)

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

Telling it wrongly reads no atoms rather than dying, which is what forcing a
format means; `structure_info()` looks at the file and is the usual way in.

## cif_info

    my $info = cif_info($file, %options);

The same for mmCIF/PDBx.

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

    my $seq  = structure_sequences($info);          # { A => 'FPTIPLSRL...' }
    my $same = structure_sequences('1ubq.pdb');     # read on the spot
    my $fast = structure_sequences('1ubq.pdb', atoms => 0, meta => 0);

The observed single-letter sequence of every chain that has one.

The first argument is either the hash reference from `structure_info()` or the
name of a file, which is read with the options given. `atoms => 0` is worth
knowing about here, since a sequence needs the residues and not their
coordinates. Options belong with a file name; passing them alongside a
structure that is already parsed is an error, because there is nothing left
for them to change.

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

## aa1to3

    aa1to3('A');      # 'ALA'
    aa1to3('X');      # 'UNK'
    aa1to3('B');      # 'ASX'   ASP or ASN, as the format spells it
    aa1to3('*');      # ''      not a single-letter code

`aa3to1` backwards: the three-letter name a single-letter code stands for, and
the empty string for anything that is not one of the twenty-six. Blanks and
case do not matter, since the letter usually comes out of a sequence string
rather than out of a file.

Every letter of the alphabet has a name, because the ambiguity codes have one
of their own — `B` is ASX, `Z` is GLX, `J` is XLE, `X` is UNK. Going this way
there is only ever one answer: `aa3to1` maps sixty-odd names onto `C`, and
only CYS comes back.

Amino acids only, as the name says. `aa1to3('A')` is ALA and not adenine, and
`aa1to3('T')` is THR and not thymine — a caller who wants `' DA'` already
knows the chain is DNA, and a function that guessed from a bare letter would be
wrong half the time.

It is in the XS rather than in Perl because it is both faster and smaller
there: the table is 104 bytes of read-only memory in the shared object, shared
between every process that loads the module, against 3,350 bytes of hash per
interpreter, and the lookup is one bounds check and one array index instead of
a hash lookup — about 4.5× the throughput measured a letter at a time.

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

    my @can = formats(); # ('mmcif', 'pdb')
    my $all = formats(); # every format known, supported or not

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

# Author

David E. Condon <dec986@gmail.com>

# License

This library is free software; you can redistribute it and/or modify it under
the same terms as Perl itself.
