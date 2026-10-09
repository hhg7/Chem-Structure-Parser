Found but not fixed (they need your judgement or larger changes)

- Free-residue demotion (_demote_free_residues) leaves some ligands in the polymer. 1hsl, 2rok, 3aya and 3lms are examples; 10 chains end up with a negative n_missing. The comment's claim about 3lms is also false.
- mmCIF header values are folded into the same keys as the PDB ones but keep mmCIF's spelling: deposit_date is 2020-01-01 where the PDB reader gives 01-JAN-20, authors are 'Condon, D.E.' against 'D.E.CONDON', a revision's type is 'Structure model' where REVDAT has 0 or 1, and het{HOH}{instances} lists every water from _pdbx_nonpoly_scheme where the PDB reader lists none. Converting them is a decision about which spelling a caller should get.
- A numeric field in a REMARK that is not a number reads as its first number: a pH of '0-200.00' is 0.
- mmCIF: a closing quote followed directly by '#' ('abc'#c) is not taken as closing. gemmi takes it; CIF 1.1's grammar wants whitespace before the comment.
- mmCIF without type_symbol: MSE's SE reads as sulphur, the name having no columns to say it is a two-letter element. The ion rule (an atom named for its residue) does not reach it.
- structure_rmsd() pairs atoms by chain, residue key and atom name and nothing else. Two files whose numbering is offset pair the wrong residues without a word wherever the numbers overlap: 185l renumbered by +100 against 184l, with the chain renamed, gave 10.18 A over 62 CA atoms. A check that paired residues have the same name, or a sequence-alignment pairing, would catch it.
- The order of structure_contacts() and the other pair lists within a residue is the order the neighbour grid's cells are walked, which depends on the box; nothing documents an order, and nothing sorts them.
- A few leaks or missing checks are reachable only through the internal XS entry points.

Where speed and memory could still be saved

Measured 2026-10-09 on this machine, perl 5.44.0, gcc -O2, while other jobs were running, so the figures are a guide rather than a benchmark.

- Decompression: every entry of PDBbind v2020 is now a .bz2. Over 200 of them, bunzip2 alone took 2.84 s of a 9.10 s structure_info(features => 0) read (31%); the same files gzipped (-9) took 0.81 s of 7.87 s, at 35% more disk. On 4fqr it is 0.36 s of 0.55 s. The module's read loop is as fast as IO::Uncompress's one-shot functions, so the cost is bzip2's, not the code's.
- Holding the unpacked text adds to the peak: 4fqr with features => 0 peaks at 175 MB from its .bz2 and 143 MB from the plain file.
- The physical properties are about three times the read on 4fqr (1.66 s against 0.55 s, peak 247 MB). The surface is most of it; structure_features() with only the wanted options is the way round it.
- The contacts list is 25.7 MB of 4fqr's 152 MB structure: 55,117 hashes of five keys, about 490 bytes each. A more compact shape would change the API.
- An atom hash has 32 buckets for its 11 keys on perl 5.44, because hv_ksplit(a, 16) asks for room for 24. Asking for 11 would leave about half the atoms at 16 buckets: roughly 5 MB of 4fqr's 107 MB.
- With features => 0, freeing the structure is about a fifth of the time again on top of the read.
