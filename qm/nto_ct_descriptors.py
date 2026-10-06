#!/usr/bin/env python3
"""Compute the Lambda (Peach/Tozer) overlap index and D_CT (Le Bahers
charge-transfer distance) for an NTO hole/particle pair, from the two
Gaussian cube files make_nto_cubes.sh generates.

Lambda = normalized spatial overlap integral of |hole| and |particle|
orbital densities:
    Lambda = ( Int |psi_h(r)| |psi_p(r)| dV ) / sqrt( Int psi_h^2 dV * Int psi_p^2 dV )
Lambda -> 1: hole and particle occupy the same region (local excitation).
Lambda -> 0: no spatial overlap (strong charge transfer). Peach et al.
(J. Chem. Phys. 2008) flag Lambda < ~0.3-0.4 as diagnostic of a TD-DFT
state whose energy a standard (non-range-separated) functional is likely
to get wrong.

D_CT = distance between the centroids of the hole and particle densities
(the NTO-pair analogue of Le Bahers' density-difference-based D_CT;
Angstrom J. Chem. Theory Comput. 2011), i.e. how far the "hole center" and
"particle center" sit apart in space.

Usage:
    python3 nto_ct_descriptors.py HOLE.cube PARTICLE.cube [--csv row.csv --label "state44"]
"""
import argparse
import csv as csvmod
import numpy as np

BOHR_TO_ANG = 0.529177210903


def read_cube(path):
    with open(path) as f:
        lines = f.readlines()
    natoms_raw = int(lines[2].split()[0])
    # A negative atom count (ORCA does this) means there's one extra header
    # line right after the atom block: "NVal  orbital_index..." -- skip it.
    multi_dataset = natoms_raw < 0
    natoms = abs(natoms_raw)
    origin = np.array([float(x) for x in lines[2].split()[1:4]])
    nx, dx = int(lines[3].split()[0]), np.array([float(x) for x in lines[3].split()[1:4]])
    ny, dy = int(lines[4].split()[0]), np.array([float(x) for x in lines[4].split()[1:4]])
    nz, dz = int(lines[5].split()[0]), np.array([float(x) for x in lines[5].split()[1:4]])
    data_start = 6 + natoms + (1 if multi_dataset else 0)
    vals = []
    for line in lines[data_start:]:
        vals.extend(float(x) for x in line.split())
    data = np.array(vals).reshape(nx, ny, nz)
    voxel_vol = abs(np.dot(dx, np.cross(dy, dz)))  # bohr^3
    # grid point coordinates (bohr)
    ix, iy, iz = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    coords = (origin[None, None, None, :]
              + ix[..., None] * dx[None, None, None, :]
              + iy[..., None] * dy[None, None, None, :]
              + iz[..., None] * dz[None, None, None, :])
    return data, coords, voxel_vol


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("hole_cube")
    ap.add_argument("particle_cube")
    ap.add_argument("--label", default="")
    ap.add_argument("--csv", default=None, help="append a row to this CSV")
    args = ap.parse_args()

    psi_h, coords_h, dV_h = read_cube(args.hole_cube)
    psi_p, coords_p, dV_p = read_cube(args.particle_cube)
    assert psi_h.shape == psi_p.shape, "hole/particle cubes must share the same grid"
    dV = dV_h

    norm_h = np.sum(psi_h**2) * dV
    norm_p = np.sum(psi_p**2) * dV
    overlap = np.sum(np.abs(psi_h) * np.abs(psi_p)) * dV
    lam = overlap / np.sqrt(norm_h * norm_p)

    centroid_h = np.sum(coords_h * (psi_h**2)[..., None], axis=(0, 1, 2)) * dV / norm_h
    centroid_p = np.sum(coords_p * (psi_p**2)[..., None], axis=(0, 1, 2)) * dV / norm_p
    d_ct_bohr = np.linalg.norm(centroid_p - centroid_h)
    d_ct_ang = d_ct_bohr * BOHR_TO_ANG

    print(f"{args.label or args.hole_cube}")
    print(f"  Lambda (hole/particle overlap) = {lam:.4f}")
    print(f"  D_CT (centroid separation)     = {d_ct_ang:.3f} Angstrom")
    if lam < 0.3:
        print("  -> Lambda < 0.3: strong CT character; B3LYP-class functional energy is suspect")
    elif lam < 0.5:
        print("  -> Lambda 0.3-0.5: moderate/partial CT character")
    else:
        print("  -> Lambda > 0.5: predominantly local excitation, not CT")

    if args.csv:
        row = {"label": args.label or args.hole_cube, "lambda": round(lam, 4),
               "d_ct_angstrom": round(d_ct_ang, 3)}
        import os
        write_header = not os.path.exists(args.csv)
        with open(args.csv, "a", newline="") as f:
            w = csvmod.DictWriter(f, fieldnames=list(row.keys()))
            if write_header:
                w.writeheader()
            w.writerow(row)


if __name__ == "__main__":
    main()
