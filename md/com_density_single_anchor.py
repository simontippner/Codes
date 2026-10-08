#!/usr/bin/env python3
"""Spatial density heatmap of solute body positions, centered on ONE SPECIFIC
body (e.g. one chosen PS molecule) rather than the solute assembly's
collective centroid -- the single-colormap style from the first version of
this analysis (see com_density_timeseries.py for the species-colored,
group-centroid-recentered version this complements).

Unlike com_density_timeseries.py, which deliberately avoids pinning any one
body at the origin (that was diagnosed as a reference-choice artifact when
it happened implicitly via cpptraj's autoimage anchor), this script does
exactly that ON PURPOSE: cpptraj's `autoimage anchor <mask> origin` is used
directly with NO centroid correction afterward, so the chosen anchor body
sits fixed at the origin (marked with a star) in every panel by
construction, and what's plotted is genuinely "what does the local
environment around this one specific molecule look like." That is a
different, equally legitimate question from "where does the solute
assembly condense relative to its own center" -- just make sure it's the
question being asked before using this script instead of
com_density_timeseries.py.

Defaults to short (e.g. 0.1 ns) windows for near-instantaneous snapshots
rather than the long (e.g. 40 ns) pooled windows used for the condensation
story -- note short windows are necessarily sparse (few frames x few
bodies), which shows up as a speckled rather than smooth density; that is
expected, not a bug, and the tradeoff for temporal sharpness.

Requires cpptraj (AmberTools) on PATH.

Usage:
    python3 com_density_single_anchor.py PRMTOP TRAJ \\
        --body-mask ":1-10" --body-mask ":11-20" ... \\
        --anchor-mask ":101" \\
        --n-panels 6 --window-ns 0.1 --dt-ps 10 --box-half 32 \\
        --out com_density_single_anchor.png
"""
import argparse
import subprocess
import tempfile
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def extract_com(prmtop, traj, body_masks, anchor_mask, workdir):
    out = os.path.join(workdir, "com.dat")
    lines = [f"parm {prmtop}", f"trajin {traj}", f"autoimage anchor {anchor_mask} origin"]
    for i, mask in enumerate(body_masks):
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
    ap.add_argument("--body-mask", action="append", required=True,
                     help="AMBER mask for one solute body's COM; repeat per body")
    ap.add_argument("--anchor-mask", required=True,
                     help="the ONE body to pin at the origin, e.g. a single PS molecule's mask")
    ap.add_argument("--n-panels", type=int, default=6)
    ap.add_argument("--window-ns", type=float, default=0.1,
                     help="width of each panel's pooled window, in ns (short -> sharp/sparse, long -> smooth/averaged)")
    ap.add_argument("--dt-ps", type=float, default=10.0, help="ps between saved frames")
    ap.add_argument("--box-half", type=float, default=32.0)
    ap.add_argument("--bins", type=int, default=40)
    ap.add_argument("--cmap", default="viridis")
    ap.add_argument("--title", default=None)
    ap.add_argument("--out", default="com_density_single_anchor.png")
    args = ap.parse_args()

    n_bodies = len(args.body_mask)
    with tempfile.TemporaryDirectory() as workdir:
        com_file = extract_com(args.prmtop, args.traj, args.body_mask, args.anchor_mask, workdir)
        data = np.loadtxt(com_file, comments="#")

    frames = data[:, 0]
    xyz = data[:, 1:].reshape(len(frames), n_bodies, 6)[:, :, :3]
    # NO centroid subtraction: the anchor body stays pinned at the origin by
    # construction, intentionally (see module docstring).

    n = len(frames)
    bin_frames = max(1, int(round(args.window_ns * 1000.0 / args.dt_ps)))
    centers = np.linspace(bin_frames, max(n - bin_frames, bin_frames), args.n_panels).astype(int)
    centers = np.clip(centers, bin_frames // 2, n - bin_frames // 2 - 1)

    hists, labels = [], []
    for c in centers:
        lo, hi = max(0, c - bin_frames // 2), min(n, c + bin_frames // 2)
        sl = xyz[lo:hi]
        xs, ys = sl[:, :, 0].ravel(), sl[:, :, 1].ravel()
        h, _, _ = np.histogram2d(xs, ys, bins=args.bins,
                                  range=[[-args.box_half, args.box_half], [-args.box_half, args.box_half]])
        hists.append(h)
        labels.append(f"{frames[c] * args.dt_ps / 1000.0:.2f} ns")

    vmax = max(h.max() for h in hists)
    norm = LogNorm(vmin=1, vmax=vmax)

    fig, axes = plt.subplots(1, args.n_panels, figsize=(3.0 * args.n_panels, 3.6))
    if args.n_panels == 1:
        axes = [axes]
    im = None
    for ax, h, lab in zip(axes, hists, labels):
        im = ax.imshow(h.T, origin="lower",
                        extent=[-args.box_half, args.box_half, -args.box_half, args.box_half],
                        cmap=args.cmap, aspect="equal", norm=norm)
        ax.plot(0, 0, marker="*", color="red", markersize=10, markeredgecolor="white", markeredgewidth=0.5)
        ax.set_title(lab, fontsize=10)
        ax.set_xticks([-args.box_half + 2, 0, args.box_half - 2])
        ax.set_yticks([-args.box_half + 2, 0, args.box_half - 2])
        ax.tick_params(labelsize=7)
    axes[0].set_ylabel("Δy (Å)")
    axes[0].set_xlabel("Δx (Å)")
    fig.colorbar(im, ax=axes, shrink=0.8, pad=0.01,
                 label=f"solute COM density ({args.window_ns:g} ns window)")
    title = args.title or (f"Solute density centered on one fixed body (red star), "
                            f"{args.window_ns:g} ns windows")
    fig.suptitle(title, fontsize=12)
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print("wrote", args.out)


if __name__ == "__main__":
    main()
