#!/usr/bin/env python3
"""Quantify how much each fragment's internal geometry changes between two
xyz files with identical atom order (e.g. each fragment's own isolated-
optimized geometry, simply placed together, vs. the geometry after a joint
"docking" optimization) -- isolating real conformational change (bond/
angle/dihedral flexing) from simple rigid-body translation/rotation of the
fragment as a whole, via the Kabsch algorithm.

For each fragment (atom index range, 1-based inclusive):
  1. Center both the "before" and "after" coordinate sets on their own
     centroid.
  2. Find the optimal rotation (SVD of the cross-covariance matrix, with a
     determinant check to reject reflections) that best superimposes
     "before" onto "after".
  3. RMSD after that optimal alignment = internal deformation only.
  4. Per-atom displacement after alignment identifies which specific atoms
     moved the most.
  5. Optionally, specific bond lengths (atom index pairs) are compared
     before/after.

Usage:
    python3 fragment_geometry_change.py BEFORE.xyz AFTER.xyz \
        --frag 1-61 --frag-name rubpy \
        --frag 62-77 --frag-name mos \
        --bond 61,1 --bond 61,21 --bond 62,65 [...]
"""
import argparse
import numpy as np


def read_xyz(path):
    with open(path) as f:
        lines = f.readlines()
    n = int(lines[0])
    elems, coords = [], []
    for line in lines[2:2 + n]:
        parts = line.split()
        elems.append(parts[0])
        coords.append([float(x) for x in parts[1:4]])
    return elems, np.array(coords)


def parse_range(spec):
    atoms = []
    for part in spec.split(","):
        if "-" in part:
            lo, hi = part.split("-")
            atoms.extend(range(int(lo), int(hi) + 1))
        else:
            atoms.append(int(part))
    return atoms


def kabsch_rmsd(P, Q):
    """P, Q: (N,3) arrays, same atom order. Returns (rmsd, per_atom_dist)."""
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    P_aligned = (R @ Pc.T).T
    diff = P_aligned - Qc
    per_atom = np.linalg.norm(diff, axis=1)
    rmsd = np.sqrt(np.mean(per_atom**2))
    return rmsd, per_atom


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("before_xyz")
    ap.add_argument("after_xyz")
    ap.add_argument("--frag", action="append", required=True,
                     help="atom range, 1-based, e.g. 1-61 (repeatable)")
    ap.add_argument("--frag-name", action="append", required=True,
                     help="name for the matching --frag (repeatable, same order)")
    ap.add_argument("--bond", action="append", default=[],
                     help="atom1,atom2 (1-based) to report bond length before/after (repeatable)")
    ap.add_argument("--top-movers", type=int, default=5,
                     help="how many highest-displacement atoms to list per fragment")
    args = ap.parse_args()

    elems_b, coords_b = read_xyz(args.before_xyz)
    elems_a, coords_a = read_xyz(args.after_xyz)
    assert elems_b == elems_a, "atom order must match between the two xyz files"

    for frag_spec, frag_name in zip(args.frag, args.frag_name):
        idx = [i - 1 for i in parse_range(frag_spec)]
        P = coords_b[idx]
        Q = coords_a[idx]
        rmsd, per_atom = kabsch_rmsd(P, Q)
        print(f"\n=== {frag_name} (atoms {frag_spec}, {len(idx)} atoms) ===")
        print(f"Internal RMSD (Kabsch-aligned): {rmsd:.4f} Angstrom")
        order = np.argsort(-per_atom)[:args.top_movers]
        print(f"Top {args.top_movers} movers (post-alignment displacement):")
        for k in order:
            atom_num = idx[k] + 1
            print(f"  atom {atom_num:>4} ({elems_b[idx[k]]}): {per_atom[k]:.4f} Angstrom")

    if args.bond:
        print("\n=== Bond length changes ===")
        for b in args.bond:
            i, j = [int(x) - 1 for x in b.split(",")]
            d_before = np.linalg.norm(coords_b[i] - coords_b[j])
            d_after = np.linalg.norm(coords_a[i] - coords_a[j])
            print(f"  {i+1}({elems_b[i]})-{j+1}({elems_b[j]}): "
                  f"{d_before:.4f} -> {d_after:.4f} Angstrom "
                  f"(delta {d_after-d_before:+.4f})")


if __name__ == "__main__":
    main()
