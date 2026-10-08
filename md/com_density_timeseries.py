#!/usr/bin/env python3
"""Spatial density heatmap of multiple solute bodies' centers of mass, binned
into time windows, to visualize aggregation/condensation over the course of
an MD trajectory.

For each frame, extracts the COM of every given body (via cpptraj `vector
... center`), re-centers all bodies on their own collective centroid each
frame (so what's plotted is genuine spread/condensation relative to the
solute assembly's own center -- NOT relative to one fixed/pinned reference
body, which would be a reference-choice artifact), then pools the recentered
positions within each time window into a 2D (x-y) histogram and plots it on
a shared log color scale.

One or more systems (trajectories) can be given as separate rows, e.g. to
contrast a weakly-aggregating vs. a strongly-aggregating condition side by
side on the same time-window and color axes.

Requires cpptraj (AmberTools) on PATH.

Usage:
    python3 com_density_timeseries.py \\
        --system "0% water" PRMTOP1 TRAJ1 \\
        --system "100% water" PRMTOP2 TRAJ2 \\
        --body-mask ":1-10" --body-mask ":11-20" ... \\
        --anchor-mask ":1-10" \\
        --windows 6 --dt-ps 10 --box-half 32 \\
        --out com_density_timeseries.png

--body-mask: one AMBER mask per solute body whose COM should be tracked
(repeat for each body -- e.g. 10 CAT clusters + 10 RUB complexes = 20 flags).
--anchor-mask: a single-molecule mask used only to resolve periodic-boundary
imaging consistently across all bodies each frame (passed to cpptraj
`autoimage anchor ... origin`) -- does NOT bias the plotted density, since
every body's position is re-centered on the collective centroid afterward.
--dt-ps: trajectory save interval, for converting frame number to ns in
panel titles (default 10 ps/frame, this project's convention).
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
    lines = [f"parm {prmtop}", f"trajin {traj}",
             f"autoimage anchor {anchor_mask} origin"]
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


def load(path, n_bodies):
    data = np.loadtxt(path, comments="#")
    frames = data[:, 0]
    rest = data[:, 1:]
    bodies = rest.reshape(len(frames), n_bodies, 6)[:, :, :3]
    centroid = bodies.mean(axis=1, keepdims=True)
    return frames, bodies - centroid


def histograms(frames, bodies, n_windows, box_half, bins=40):
    n_frames = len(frames)
    edges = np.linspace(0, n_frames, n_windows + 1).astype(int)
    hs, labels = [], []
    for i in range(n_windows):
        lo, hi = edges[i], edges[i + 1]
        sl = bodies[lo:hi]
        h, _, _ = np.histogram2d(sl[:, :, 0].ravel(), sl[:, :, 1].ravel(), bins=bins,
                                  range=[[-box_half, box_half], [-box_half, box_half]])
        hs.append(h)
        labels.append((lo, hi))
    return hs, labels, edges


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--system", nargs=3, action="append", required=True,
                     metavar=("LABEL", "PRMTOP", "TRAJ"),
                     help="one row of the figure; repeat for multiple systems")
    ap.add_argument("--body-mask", action="append", required=True,
                     help="AMBER mask for one solute body's COM; repeat per body")
    ap.add_argument("--anchor-mask", required=True,
                     help="single-molecule mask for periodic-image resolution (does not bias the plot)")
    ap.add_argument("--windows", type=int, default=6)
    ap.add_argument("--dt-ps", type=float, default=10.0, help="ps between saved frames")
    ap.add_argument("--box-half", type=float, default=32.0, help="plot range, +/- Angstrom")
    ap.add_argument("--bins", type=int, default=40)
    ap.add_argument("--title", default="Solute-body COM density over time\n"
                     "(relative to the solutes' own collective center)")
    ap.add_argument("--out", default="com_density_timeseries.png")
    args = ap.parse_args()

    n_bodies = len(args.body_mask)
    all_hs, all_labels, all_frames = [], [], []
    with tempfile.TemporaryDirectory() as workdir:
        for label, prmtop, traj in args.system:
            com_file = extract_com(prmtop, traj, args.body_mask, args.anchor_mask, workdir)
            frames, bodies = load(com_file, n_bodies)
            hs, win_frames, _ = histograms(frames, bodies, args.windows, args.box_half, args.bins)
            all_hs.append(hs)
            all_labels.append(label)
            all_frames.append((frames, win_frames))

    vmax = max(h.max() for hs in all_hs for h in hs)
    norm = LogNorm(vmin=1, vmax=vmax)

    n_sys = len(all_hs)
    fig, axes = plt.subplots(n_sys, args.windows, figsize=(2.6 * args.windows, 3.1 * n_sys),
                              squeeze=False)
    fig.subplots_adjust(top=1 - 0.16 / n_sys, left=0.07, right=0.88,
                         wspace=0.08, hspace=0.25)
    im = None
    for row, (hs, label, (frames, win_frames)) in enumerate(zip(all_hs, all_labels, all_frames)):
        for i, (h, (lo, hi)) in enumerate(zip(hs, win_frames)):
            ax = axes[row][i]
            im = ax.imshow(h.T, origin="lower",
                            extent=[-args.box_half, args.box_half, -args.box_half, args.box_half],
                            cmap="viridis", aspect="equal", norm=norm)
            t_lo = frames[lo] * args.dt_ps / 1000.0
            t_hi = frames[min(hi, len(frames) - 1)] * args.dt_ps / 1000.0
            ax.set_title(f"{t_lo:.0f}–{t_hi:.0f} ns", fontsize=9)
            if i == 0:
                ax.set_xticks([-args.box_half + 2, 0, args.box_half - 2])
                ax.set_yticks([-args.box_half + 2, 0, args.box_half - 2])
                ax.tick_params(labelsize=7)
            else:
                ax.set_xticks([]); ax.set_yticks([])
        axes[row][0].set_ylabel(label + "\nΔy (Å)", fontsize=10, fontweight="bold")
    axes[-1][0].set_xlabel("Δx (Å)", fontsize=9)

    cax = fig.add_axes([0.90, 0.20, 0.015, 0.55])
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(f"COM density, log scale\n(counts per bin, pooled over window)", fontsize=9)

    fig.suptitle(args.title, fontsize=13)
    fig.savefig(args.out, dpi=160)
    print("wrote", args.out, "(shared color scale vmax =", vmax, ")")


if __name__ == "__main__":
    main()
