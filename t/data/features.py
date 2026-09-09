# What mdtraj says about the structures in this directory: solvent-accessible
# surface, radius of gyration, mass, and which aromatic rings are stacked.
#
# Written for t/features.t and t/data/features.pl, which are the only two things
# that run it: the test compares this module against mdtraj live where mdtraj is
# installed, and the generator freezes this output into features.txt so that the
# comparison still runs where it is not.  One copy of the script, so the frozen
# answer and the live one cannot be answers to different questions.
#
# Usage: python3 features.py FILE
#
# mdtraj 1.11.1, importable from /home/con/.pyenv/versions/3.14.2/bin/python3.
#
# Two things here are not a straight call into mdtraj, and both are marked at
# the point they happen:
#
# Every surface is printed twice, as mdtraj computes it and as the same
# algorithm computes it in float64.  mdtraj's kernel is float32 throughout --
# coordinates, radii, the centred sphere points and the distances between them
# -- so a sphere point sitting within a float32 ulp of a neighbouring atom's
# surface can be called accessible by one width and covered by the other.  On
# mini.pdb two of 62,400 points flip, worth 0.26 A^2 of 1551.  The float64
# column is mdtraj's own loop with the dtype changed and nothing else, and it is
# what pins the answer exactly; the float32 column is what mdtraj ships and is
# compared against with a tolerance of one sphere point per atom.
#
# The disulfides are mdtraj's rule and not mdtraj's function.  Its
# Topology.create_disulfide_bonds() is handed the positions the PDB reader
# parsed, which are angstrom, and compares them with 0.3 -- "this is supposed to
# be nm. I think we're good", says the comment beside it -- so it tests for an
# SG-SG separation under 0.3 A and finds nothing on any file.  The SG-SG bonds
# that do appear in topology.bonds come from the file's own CONECT records,
# which is the declared answer rather than a computed one.  So the rule is
# transcribed below with its units made consistent, and both answers are printed
# so the difference stays visible.
#
# pi_stacking() is handed unitcell_vectors of None.  A PDB CRYST1 record
# describes the crystal, not a simulation box, and mdtraj's centroid distances
# apply the minimum image convention when a box is set -- which would measure
# from a ring to a symmetry copy of another ring that is not in the file at all.
import sys

import numpy as np
import mdtraj as md
from mdtraj.geometry.sasa import _ATOMIC_RADII

# Biopython, for the half-sphere exposure only: HSExposureCB is Biopython's and
# mdtraj has no equivalent.  It reads the file itself, so the H lines are
# written for PDB files alone -- t/cif.t is what says the two formats agree.
from Bio.PDB import PDBParser
from Bio.PDB.HSExposure import HSExposureCB

PROBE = 0.14      # nm, mdtraj's default
NPTS = 960        # sphere points per atom, mdtraj's default

# The aromatic rings this module knows how to find, in the order it lists their
# atoms.  mdtraj's compute_ring_normal() takes the plane from the first two
# atoms of the group, so the order is part of the answer and has to be the same
# on both sides; this is the table in Parser.xs, spelled in Python.
RINGS = {
    'PHE': [('6', ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'])],
    'TYR': [('6', ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'])],
    'TRP': [('6', ['CD2', 'CE2', 'CZ2', 'CH2', 'CZ3', 'CE3']),
            ('5', ['CG', 'CD1', 'NE1', 'CE2', 'CD2'])],
    'A':   [('6', ['N1', 'C2', 'N3', 'C4', 'C5', 'C6']),
            ('5', ['C4', 'C5', 'N7', 'C8', 'N9'])],
    'C':   [('6', ['N1', 'C2', 'N3', 'C4', 'C5', 'C6'])],
    'T':   [('6', ['N1', 'C2', 'N3', 'C4', 'C5', 'C6'])],
    'U':   [('6', ['N1', 'C2', 'N3', 'C4', 'C5', 'C6'])],
}
for _h in ('HIS', 'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP'):
    RINGS[_h] = [('5', ['CG', 'ND1', 'CE1', 'NE2', 'CD2'])]
for _p in ('DA', 'DG'):
    RINGS[_p] = RINGS['A']
for _p in ('DC', 'DT', 'DU'):
    RINGS[_p] = RINGS['C']
RINGS['G'] = RINGS['A']

# structure_pi_stacking()'s defaults, in the units mdtraj takes them in.  The
# face-to-face distance is 5.5 A, not the 5.5 nm mdtraj's signature carries: see
# the note in the Perl.  The other three distances are mdtraj's own.
PI_KW = dict(
    max_face_to_face_centroid_distance=0.55,
    face_plane_angle_range=(0.0, 35.0),
    face_normal_to_centroid_angle_range=(0.0, 33.0),
    max_edge_to_face_centroid_distance=0.65,
    edge_plane_angle_range=(50.0, 90.0),
    edge_normal_to_centroid_angle_range=(0.0, 30.0),
    edge_intersection_radius=0.15,
)


def sasa_float64(xyz_nm, radii_nm, npts):
    """mdtraj's Shrake-Rupley kernel, in float64.  Same sphere points, same
    neighbour rule, same constant; only the width differs."""
    xyz = np.asarray(xyz_nm, dtype=np.float64)
    rad = np.asarray(radii_nm, dtype=np.float64)
    n = len(rad)
    i = np.arange(npts)
    inc = np.pi * (3.0 - np.sqrt(5.0))
    off = 2.0 / npts
    y = i * off - 1.0 + off / 2.0
    pts = np.stack([np.cos(i * inc) * np.sqrt(1.0 - y * y), y,
                    np.sin(i * inc) * np.sqrt(1.0 - y * y)], axis=1)
    out = np.zeros(n)
    idx = np.arange(n)
    for a in range(n):
        d = np.linalg.norm(xyz - xyz[a], axis=1)
        nb = np.where((d < rad + rad[a]) & (idx != a))[0]
        p = xyz[a] + rad[a] * pts
        if len(nb):
            dd = np.linalg.norm(p[:, None, :] - xyz[nb][None, :, :], axis=2)
            acc = int(np.all(dd >= rad[nb][None, :], axis=1).sum())
        else:
            acc = npts
        out[a] = acc * 4.0 * np.pi / npts * rad[a] ** 2
    return out


# Residues are named by their position in the file, not by chain and number.
#
# Two things make a name the wrong key here.  mdtraj keeps no insertion code, so
# HIS 8 and HIS 8A of mini.pdb are both "8" to it; and it takes an mmCIF file's
# chains from label_asym_id, which is the archive's lettering rather than the
# depositor's, so mini.cif's chain A comes back as B, C, D, E, F and G.  Both
# readers walk the file in the order it was written, so the nth residue is the
# nth residue in either -- and a reader that split them differently would show
# up as a mismatch, which is the point.
def res_index(res):
    return res.index


def disulfides(t):
    """mdtraj's Topology.create_disulfide_bonds() rule, in nanometres: a CYS
    with an SG and no HG, paired with another under 0.3 nm."""
    def is_cyx(res):
        names = [a.name for a in res.atoms]
        return "SG" in names and "HG" not in names

    cyx = [r for r in t.topology.residues if r.name == "CYS" and is_cyx(r)]
    sg = {r.index: [a.index for a in r.atoms if a.name == "SG"][0] for r in cyx}
    pos = t.xyz[0]
    out = []
    for i in range(len(cyx)):
        for j in range(i):
            a, b = sg[cyx[i].index], sg[cyx[j].index]
            d = float(np.linalg.norm(pos[a] - pos[b]))
            if d < 0.3:
                lo, hi = sorted((cyx[i].index, cyx[j].index))
                out.append((lo, hi, d * 10.0))
    return sorted(out)


# Biopython's Bio.PDB.Polypeptide.PPBuilder `radius': the largest C-to-N
# separation that still means two residues are joined.
PEPTIDE_BOND = 0.18   # nm


def linked(t, prev, nxt):
    """Are these two residues actually peptide-bonded?  mdtraj does not ask --
    it takes the residue next to this one in the file and computes a phi across
    whatever gap is there -- so the answer is printed beside every angle and the
    test uses it."""
    if prev is None or nxt is None or prev.chain is not nxt.chain:
        return False
    c = [a.index for a in prev.atoms if a.name == "C"]
    n = [a.index for a in nxt.atoms if a.name == "N"]
    if not c or not n:
        return False
    return bool(np.linalg.norm(t.xyz[0][c[0]] - t.xyz[0][n[0]]) < PEPTIDE_BOND)


# The shorter projected length against the bond it came from, which is the sine
# of the angle between them: below this the four atoms are collinear and there
# is no torsion to measure.  The same rule and the same constant as Parser.xs,
# which is where the two measurements behind the constant are written down --
# 0.736 at worst for real geometry, 1.1e-6 at best for mini.pdb's collinear one.
COLLINEAR = 1e-3


def torsion_defined(t, idx):
    """Is there an angle here at all?  mini.pdb lays a residue's atoms out along
    a straight line on purpose, and four collinear atoms have no torsion --
    atan2(0, 0) gives zero, which is a number where there is no answer."""
    p = t.xyz[0][list(idx)].astype(np.float64)
    b0, b1, b2 = p[0] - p[1], p[2] - p[1], p[3] - p[2]
    n1 = np.linalg.norm(b1)
    l0, l2 = np.linalg.norm(b0), np.linalg.norm(b2)
    if n1 == 0 or l0 == 0 or l2 == 0:
        return False
    u = b1 / n1
    v = b0 - np.dot(b0, u) * u
    w = b2 - np.dot(b2, u) * u
    return bool(np.linalg.norm(v) > l0 * COLLINEAR and np.linalg.norm(w) > l2 * COLLINEAR)


def torsions(t):
    """mdtraj's phi, psi, omega and chi1..chi5, by residue index, each with
    whether the residues it was measured across are bonded and whether the four
    atoms define an angle at all."""
    top = t.topology
    by_index = {r.index: r for r in top.residues}
    rows = []

    def neighbour(res, step):
        other = by_index.get(res.index + step)
        return other if other is not None and other.chain is res.chain else None

    for name, fn, owner in (("phi", md.compute_phi, 1),
                            ("psi", md.compute_psi, 0),
                            ("omega", md.compute_omega, 0)):
        idx, ang = fn(t)
        for k in range(len(idx)):
            res = top.atom(int(idx[k][owner])).residue
            if name == "phi":
                ok = linked(t, neighbour(res, -1), res)
            else:
                ok = linked(t, res, neighbour(res, 1))
            rows.append((res.index, name, float(np.degrees(ang[0][k])), int(ok),
                         int(torsion_defined(t, idx[k]))))

    for k, fn in enumerate((md.compute_chi1, md.compute_chi2, md.compute_chi3,
                            md.compute_chi4, md.compute_chi5), start=1):
        idx, ang = fn(t)
        for j in range(len(idx)):
            res = top.atom(int(idx[j][0])).residue
            rows.append((res.index, "chi%d" % k, float(np.degrees(ang[0][j])), 1,
                         int(torsion_defined(t, idx[j]))))
    return sorted(rows)


def recorded_sg_bonds(t):
    """What mdtraj's topology actually holds, which on a PDB file is whatever
    its CONECT records said.  Printed for the contrast, not compared against."""
    n = 0
    for a, b in t.topology.bonds:
        if a.name == "SG" and b.name == "SG":
            n += 1
    return n


def load(path):
    # mdtraj dispatches on the file name, and a 1993 entry is called .ent
    if path.endswith(('.ent', '.pdb')):
        return md.load_pdb(path)
    return md.load(path)


def rings_of(top):
    """Every aromatic ring in the topology, as (label, residue, atom indices)."""
    out = []
    for res in top.residues:
        for label, names in RINGS.get(res.name.strip(), []):
            by = {a.name: a.index for a in res.atoms}
            if not all(n in by for n in names):
                continue
            out.append((label, res, tuple(by[n] for n in names)))
    return out


def main(path):
    t = load(path)
    top = t.topology
    n = top.n_atoms
    print('#atoms %d' % n)
    if n == 0:
        return
    radii = np.array([_ATOMIC_RADII[a.element.symbol] for a in top.atoms]) + PROBE
    a32 = md.shrake_rupley(t, probe_radius=PROBE, n_sphere_points=NPTS)[0] * 100.0
    a64 = sasa_float64(t.xyz[0], radii, NPTS) * 100.0
    print('#residues %d' % top.n_residues)
    # The area one sphere point is worth, per atom.  It is the whole of the
    # difference between the two widths: an atom's area is a count of points
    # times this, so a point that one width calls accessible and the other calls
    # covered moves the answer by exactly this much and nothing can move it by
    # less.  Printed so that t/features.t can hold the float32 comparison to one
    # point per atom without carrying a copy of the radius table.
    point = 4.0 * np.pi * (radii * 10.0) ** 2 / NPTS
    per_res = {}
    rows = []
    for i, at in enumerate(top.atoms):
        k = res_index(at.residue)
        rows.append((k, at.name, a32[i], a64[i], point[i]))
        s = per_res.setdefault(k, [0.0, 0.0, 0.0])
        s[0] += a32[i]
        s[1] += a64[i]
        s[2] += point[i]
    for k, name, x32, x64, p in sorted(rows):
        print('A %d|%s %.9f %.9f %.9f' % (k, name, x32, x64, p))
    for res in top.residues:
        k = res_index(res)
        got = per_res.get(k, [0.0, 0.0, 0.0])
        # '-' for a chain with no id: a pre-1996 entry can leave column 22
        # blank, and a blank field would shift every column after it
        cid = (res.chain.chain_id or '').strip() or '-'
        print('R %d %s %d %s %.9f %.9f %.9f' % (
            k, cid, res.resSeq, res.name, got[0], got[1], got[2]))
    print('T sasa %.9f %.9f %.9f' % (a32.sum(), a64.sum(), point.sum()))

    # Each chain's surface with the other chains taken away, which with the
    # surface it has in the structure is the area the chains bury between them.
    #
    # Grouped by chain_id and not by mdtraj's own chains: mdtraj starts a new
    # chain at every TER record, so one deposited chain becomes three or four of
    # them -- the polymer, its heterogens, its waters -- while this module keeps
    # a chain whole.  Grouping by the author's chain id puts mdtraj's back
    # together and the two agree atom for atom.
    #
    # Only for PDB files.  mdtraj takes an mmCIF file's chains from
    # label_asym_id, the archive's lettering, so there is no author id to group
    # by and no comparison to make; t/cif.t is what says the two formats give
    # the same answer here.
    if path.endswith(('.pdb', '.ent')):
        by_chain = {}
        for a in top.atoms:
            by_chain.setdefault((a.residue.chain.chain_id or '').strip() or '-', []).append(a.index)
        for cid in sorted(by_chain):
            idx = by_chain[cid]
            sub = t.atom_slice(idx)
            r = np.array([_ATOMIC_RADII[a.element.symbol] for a in sub.topology.atoms]) + PROBE
            i32 = md.shrake_rupley(sub, probe_radius=PROBE, n_sphere_points=NPTS)[0].sum() * 100.0
            i64 = sasa_float64(sub.xyz[0], r, NPTS).sum() * 100.0
            pt = (4.0 * np.pi * (r * 10.0) ** 2 / NPTS).sum()
            print('I %s %d %.9f %.9f %.9f' % (cid, len(idx), i32, i64, pt))

    xyz = np.asarray(t.xyz[0], dtype=np.float64) * 10.0
    mass = np.array([a.element.mass for a in top.atoms], dtype=np.float64)
    mu = xyz.mean(0)
    print('T rg %.9f' % np.sqrt(((xyz - mu) ** 2).sum(1).mean()))
    # mass-weighted, about the centre of mass.  compute_rg(traj, masses=m)
    # weights by mass but measures from the geometric centroid, which is not
    # what the quantity means; this is that call with the centre corrected.
    com = (mass[:, None] * xyz).sum(0) / mass.sum()
    print('T rg_mass %.9f' % np.sqrt((mass * ((xyz - com) ** 2).sum(1)).sum() / mass.sum()))
    print('T mass %.9f' % mass.sum())
    print('T center %.9f %.9f %.9f' % tuple(mu))
    print('T com %.9f %.9f %.9f' % tuple(com))

    for lo, hi, d in disulfides(t):
        print('S %d %d %.9f' % (lo, hi, d))
    print('#sg_bonds_recorded %d' % recorded_sg_bonds(t))

    for ri, name, value, ok, defined in torsions(t):
        print('D %d %s %.9f %d %d' % (ri, name, value, ok, defined))

    # Residue contacts, mdtraj's closest-heavy scheme.  Only the pairs mdtraj's
    # `all' considers -- same chain, three or more apart in it -- because those
    # are the ones it has an opinion about; this module reports the neighbouring
    # and the cross-chain pairs too, and the test does not hold those to mdtraj.
    # Capped at 8 A so the frozen file stays a fixture rather than a matrix.
    try:
        cd, cpairs = md.compute_contacts(t, contacts='all', scheme='closest-heavy')
        for (a, b), dist in zip(cpairs, cd[0]):
            if dist * 10.0 <= 8.0:
                print('C %d %d %.9f' % (a, b, dist * 10.0))
    except ValueError:
        pass    # no residue pairs three apart in one chain: nothing to compare

    if path.endswith(('.pdb', '.ent')):
        st = PDBParser(QUIET=True).get_structure('x', path)
        model = next(iter(st))
        HSExposureCB(model)
        for chain in model:
            for res in chain:
                u = res.xtra.get('EXP_HSE_B_U')
                if u is None:
                    continue
                cid = (chain.id or '').strip() or '-'
                _, num, icode = res.id
                print('H %s %d %s %d %d' % (cid, num, (icode or '').strip() or '-',
                                            u, res.xtra.get('EXP_HSE_B_D')))

    # Backbone hydrogen bonds, mdtraj's kabsch_sander.  Its H is placed from
    # the previous residue in the file whether or not the two are bonded, so the
    # linked flag rides along the way it does for the torsions.
    ks = md.geometry.kabsch_sander(t)[0].tocoo()
    by_index = {r.index: r for r in top.residues}
    for a, d, e in sorted(zip(ks.row, ks.col, ks.data)):
        prev = by_index.get(int(d) - 1)
        cur = by_index.get(int(d))
        ok = linked(t, prev, cur) if (prev is not None and cur is not None) else False
        print('B %d %d %.9f %d' % (int(a), int(d), float(e), int(ok)))

    # Secondary structure, mdtraj's compute_dssp.  See the head of the block in
    # Parser.xs: this module's assignment is the Kabsch-Sander dictionary and
    # agrees with mdtraj on most residues but not all, so the test bounds the
    # disagreement rather than demanding equality.
    dssp = md.compute_dssp(t, simplified=False)[0]
    for i, v in enumerate(dssp):
        v = v.strip() or '_'
        print('X %d %s' % (i, v))

    rings = rings_of(top)
    print('#rings %d' % len(rings))
    if not rings:
        return
    t.unitcell_vectors = None    # a crystal cell is not a periodic box; see above
    groups = [g for (_, _, g) in rings]
    named = {}
    for label, res, g in rings:
        named[g] = '%d|%s' % (res_index(res), label)
    res_of = {g: res for (_, res, g) in rings}
    hits = set()
    for pair in md.geometry.pi_stacking(t, groups, groups, **PI_KW)[0]:
        a, b = pair
        if res_of[a] is res_of[b]:
            continue        # the two rings of one tryptophan are fused, not stacked
        hits.add(tuple(sorted((named[a], named[b]),
                              key=lambda t: (int(t.split('|')[0]), t))))
    for h in sorted(hits, key=lambda p: [(int(t.split('|')[0]), t) for t in p]):
        print('P %s %s' % h)


if __name__ == '__main__':
    main(sys.argv[1])
