#!/usr/bin/env python3
"""3D scatter view of solute body positions at a few snapshots across a
trajectory -- the 3D companion to com_density_timeseries.py's x-y
projection. matplotlib has no good volumetric/density rendering, so rather
than a pooled 3D density this plots a handful of representative single
frames (e.g. early/mid/late) side by side in 3D, colored by species, which
is far more legible than a dense 3D point cloud.

Requires cpptraj (AmberTools) on PATH.

Usage:
    python3 com_3d_snapshots.py PRMTOP TRAJ \\
        --body-mask CAT ":1-10" --body-mask CAT ":11-20" ... \\
        --body-mask PS ":101" --body-mask PS ":102" ... \\
        --species-color CAT "#386641" --species-color PS "#7798ab" \\
        --anchor-mask ":1-10" --n-snapshots 3 --dt-ps 10 \\
        --out com_3d_snapshots.png
"""
import argparse
import subprocess
import tempfile
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def extract_com(prmtop, traj, all_masks, anchor_mask, workdir):
    out = os.path.join(workdir, "com.dat")
    lines = [f"parm {prmtop}", f"trajin {traj}", f"autoimage anchor {anchor_mask} origin"]
    for i, mask in enumerate(all_masks):
        lines.append(f"vector B{i} {mask} out {out} center")
    lines.append("run")
    script = os.path.join(workdir, "extract.in")
    with open(script, "w") as f:
        f.write("\n".join(lines) + "\n")
    result = subprocess.run(["cpptraj", "-i", script], capture_output=True, text=True)
    if result.returncode != 0 or not os.path.exists(out):
        print(result.stdout)
        print(result.stderr)
        raise RuntimeError(f"cpptraj COM extraction failed for {traj}")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prmtop")
    ap.add_argument("traj")
    ap.add_argument("--body-mask", nargs=2, action="append", required=True,
                     metavar=("SPECIES", "MASK"))
    ap.add_argument("--species-color", nargs=2, action="append", required=True,
                     metavar=("SPECIES", "HEXCOLOR"))
    ap.add_argument("--anchor-mask", required=True)
    ap.add_argument("--n-snapshots", type=int, default=3)
    ap.add_argument("--dt-ps", type=float, default=10.0)
    ap.add_argument("--box-half", type=float, default=30.0)
    ap.add_argument("--marker-size", type=float, default=80)
    ap.add_argument("--title", default="Solute body positions in 3D\n"
                     "(single frames, relative to the solute assembly's own collective center)")
    ap.add_argument("--out", default="com_3d_snapshots.png")
    args = ap.parse_args()

    species_colors = dict(args.species_color)
    species_of_body = [sp for sp, _ in args.body_mask]
    all_masks = [mask for _, mask in args.body_mask]
    n_bodies = len(all_masks)

    with tempfile.TemporaryDirectory() as workdir:
        com_file = extract_com(args.prmtop, args.traj, all_masks, args.anchor_mask, workdir)
        data = np.loadtxt(com_file, comments="#")

    frames = data[:, 0]
    bodies = data[:, 1:].reshape(len(frames), n_bodies, 6)[:, :, :3]
    bodies = bodies - bodies.mean(axis=1, keepdims=True)

    n = len(frames)
    snap_idx = np.linspace(0, n - 1, args.n_snapshots).astype(int)

    fig = plt.figure(figsize=(5 * args.n_snapshots, 5.5))
    for col, idx in enumerate(snap_idx):
        ax = fig.add_subplot(1, args.n_snapshots, col + 1, projection="3d")
        pts = bodies[idx]
        for sp in sorted(set(species_of_body)):
            mask = [s == sp for s in species_of_body]
            p = pts[mask]
            ax.scatter(p[:, 0], p[:, 1], p[:, 2], c=species_colors[sp], s=args.marker_size,
                       edgecolors="white", linewidths=0.6, label=sp, depthshade=True)
        b = args.box_half
        ax.set_xlim(-b, b); ax.set_ylim(-b, b); ax.set_zlim(-b, b)
        ax.set_xlabel("Δx (Å)", fontsize=8)
        ax.set_ylabel("Δy (Å)", fontsize=8)
        ax.set_zlabel("Δz (Å)", fontsize=8)
        ax.set_title(f"{frames[idx] * args.dt_ps / 1000.0:.0f} ns", fontsize=11)
        ax.tick_params(labelsize=6)
        if col == 0:
            ax.legend(loc="upper left", fontsize=9)

    fig.suptitle(args.title, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    fig.savefig(args.out, dpi=150)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
