#!/usr/bin/env python3
"""
Combine two independently-optimized molecular fragments into a single xyz
file, placed side-by-side along x with a safe non-clashing gap between them
(each fragment keeps its own internal geometry; only a rigid translation is
applied to the second one).

Usage:
    python3 make_combined_xyz.py FRAGMENT1.xyz FRAGMENT2.xyz [--gap 6.0]
                                  [--out combined.xyz]
"""
import argparse
import numpy as np


def read_xyz(path):
    with open(path) as f:
        lines = f.readlines()
    n = int(lines[0])
    comment = lines[1].rstrip("\n")
    atoms = []
    coords = []
    for line in lines[2:2 + n]:
        parts = line.split()
        atoms.append(parts[0])
        coords.append([float(x) for x in parts[1:4]])
    return atoms, np.array(coords), comment


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fragment1", help="first fragment's xyz file (kept at the origin)")
    p.add_argument("fragment2", help="second fragment's xyz file (translated along +x)")
    p.add_argument("--gap", type=float, default=6.0,
                    help="closest-approach gap in Angstrom between the two fragments' "
                         "bounding spheres (default: 6.0)")
    p.add_argument("--out", default="combined.xyz", help="output xyz path (default: combined.xyz)")
    args = p.parse_args()

    atoms1, coords1, _ = read_xyz(args.fragment1)
    atoms2, coords2, _ = read_xyz(args.fragment2)

    centroid1 = coords1.mean(axis=0)
    centroid2 = coords2.mean(axis=0)
    radius1 = np.max(np.linalg.norm(coords1 - centroid1, axis=1))
    radius2 = np.max(np.linalg.norm(coords2 - centroid2, axis=1))

    # re-center each fragment on its own centroid, then place fragment 2
    # along +x at a distance = radius1 + gap + radius2 from fragment 1's
    # centroid (now at the origin)
    coords1c = coords1 - centroid1
    offset_x = radius1 + args.gap + radius2
    coords2c = (coords2 - centroid2) + np.array([offset_x, 0.0, 0.0])

    all_atoms = atoms1 + atoms2
    all_coords = np.vstack([coords1c, coords2c])

    with open(args.out, "w") as f:
        f.write(f"{len(all_atoms)}\n")
        f.write(f"{args.fragment1} + {args.fragment2}, placed {offset_x:.2f} A apart "
                f"center-to-center along x (gap={args.gap:.1f} A); "
                f"no interaction optimized\n")
        for a, (x, y, z) in zip(all_atoms, all_coords):
            f.write(f"{a:<3s}{x:15.8f}{y:15.8f}{z:15.8f}\n")

    print(f"fragment1 radius = {radius1:.2f} A")
    print(f"fragment2 radius = {radius2:.2f} A")
    print(f"center-to-center offset = {offset_x:.2f} A (gap = {args.gap:.1f} A)")
    print(f"wrote {args.out} ({len(all_atoms)} atoms)")


if __name__ == "__main__":
    main()
