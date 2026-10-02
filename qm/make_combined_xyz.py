#!/usr/bin/env python3
"""
Combine the optimized Ru(bpy)3^2+ and [Mo3S13]2- geometries into a single
xyz file, placed side-by-side along x with a safe non-clashing gap between
them (both fragments keep their own optimized internal geometry; only a
rigid translation is applied to the Mo3S13 cluster).
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUBPY_XYZ = "/gpfs/data/stippner/Luca/04_uvvis/rubpy/qm/rubpy3_opt_uvvis.xyz"
MOS_XYZ = "/gpfs/data/stippner/Luca/04_uvvis/mos/qm/mo3s13_opt_uvvis.xyz"
OUT_XYZ = os.path.join(HERE, "rubpy_mos.xyz")

GAP = 6.0  # Angstrom, closest-approach gap between the two fragments' bounding spheres


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
    atoms1, coords1, c1 = read_xyz(RUBPY_XYZ)
    atoms2, coords2, c2 = read_xyz(MOS_XYZ)

    centroid1 = coords1.mean(axis=0)
    centroid2 = coords2.mean(axis=0)
    radius1 = np.max(np.linalg.norm(coords1 - centroid1, axis=1))
    radius2 = np.max(np.linalg.norm(coords2 - centroid2, axis=1))

    # re-center each fragment on its own centroid, then place fragment 2
    # along +x at a distance = radius1 + GAP + radius2 from fragment 1's
    # centroid (now at the origin)
    coords1c = coords1 - centroid1
    offset_x = radius1 + GAP + radius2
    coords2c = (coords2 - centroid2) + np.array([offset_x, 0.0, 0.0])

    all_atoms = atoms1 + atoms2
    all_coords = np.vstack([coords1c, coords2c])

    with open(OUT_XYZ, "w") as f:
        f.write(f"{len(all_atoms)}\n")
        f.write(f"Ru(bpy)3^2+ (opt. B3LYP-D3/def2-SVP/CPCM(MeCN)) "
                f"+ [Mo3S13]2- (same level), placed {offset_x:.2f} A apart "
                f"center-to-center along x; charge sum = 0, no interaction optimized\n")
        for a, (x, y, z) in zip(all_atoms, all_coords):
            f.write(f"{a:<3s}{x:15.8f}{y:15.8f}{z:15.8f}\n")

    print(f"rubpy radius  = {radius1:.2f} A")
    print(f"mo3s13 radius = {radius2:.2f} A")
    print(f"center-to-center offset = {offset_x:.2f} A (gap = {GAP:.1f} A)")
    print(f"wrote {OUT_XYZ} ({len(all_atoms)} atoms)")


if __name__ == "__main__":
    main()
