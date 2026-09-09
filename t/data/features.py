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
# pi_stacking() is handed unitcell_vectors of None.  A PDB CRYST1 record
# describes the crystal, not a simulation box, and mdtraj's centroid distances
# apply the minimum image convention when a box is set -- which would measure
# from a ring to a symmetry copy of another ring that is not in the file at all.
import sys

import numpy as np
import mdtraj as md
from mdtraj.geometry.sasa import _ATOMIC_RADII

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
