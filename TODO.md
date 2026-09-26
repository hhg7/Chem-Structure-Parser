Found but not fixed (they need your judgement or larger changes)

- Free-residue demotion (_demote_free_residues) leaves some ligands in the polymer. 1hsl, 2rok, 3aya and 3lms are examples; 10 chains end up with a negative n_missing. The comment's claim about 3lms is also false.
- Header readings are missed or differ between the formats:
  - DBREF1/DBREF2 are ignored (232 entries).
  - SHELXL-style R values are missed (81 entries).
  - Each REVDAT continuation line becomes its own revision (1,789 entries with a blank date).
  - mmCIF helix ids come back as HELX_P1 rather than AA.
  - mmCIF sheets have no strand count or sense.
  - mmCIF never fills a chain's fragment, and synthetic entities get no organism.
- Minor:
  - mmCIF chain ids longer than 7 characters are cut short.
  - An undef option value dies with a bare Perl warning.
  - dssp => 1 with features => 0 silently returns nothing.
  - A malformed number such as 0-200.00 reads as 0.
  - A few leaks or missing checks are reachable only through the internal XS entry points.

Where speed and memory could still be saved (measured on 4fqr, 90,792 atoms)

- Surface calculation: about 91% of a default read. It already uses a neighbour grid, and neither reviewer found an O(n²) loop anywhere.
- With features => 0: the per-atom loop in _assemble takes about 20% of the time and freeing the structure 22%. Moving that loop into the C pass is the largest remaining gain.
- atoms => 0 still builds every per-atom column. That is an estimated 57 MB of the 96 MB peak.
- Contacts: 26.9 MB of the 158 MB structure.

My scratch harness (mdcmp.t plus the dumps) is in the session scratchpad if you want to rerun the comparison.
