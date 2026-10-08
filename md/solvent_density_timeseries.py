#!/usr/bin/env python3
"""Spatial density heatmap of SOLVENT positions (not solute), binned into
time windows, relative to the solute assembly's own collective center --
the solvent-side companion to com_density_timeseries.py. Shows whether
solvent is depleted from the region where solutes condense (classic
hydrophobic-collapse signature) or stays roughly uniform.

Pipeline: one cpptraj pass per system computes (a) the solute bodies' COM
vectors (same `vector ... center` calls as com_density_timeseries.py, used
only to get each frame's solute centroid for recentering) and (b) strips
the trajectory down to just the solvent's representative atom (e.g. water
O, methanol O) and writes it out as a small NetCDF trajectory -- far
cheaper than extracting a per-molecule COM for every one of the thousands
of solvent molecules the way com_density_timeseries.py does for the (few)
solute bodies. The stripped trajectory is read directly with scipy, each
frame's solvent coordinates are recentered on that frame's solute
centroid, and pooled into per-time-window 2D histograms exactly like the
solute-side script.

Requires cpptraj (AmberTools) and scipy on PATH.

Usage:
    python3 solvent_density_timeseries.py \\
        --system "0% water" PRMTOP1 TRAJ1 \\
        --system "100% water" PRMTOP2 TRAJ2 \\
        --solute-body-mask ":1-10" --solute-body-mask ":11-20" ... \\
        --solvent-mask ":WAT@O" --solvent-mask ":MOH@O1" \\
        --anchor-mask ":1-10" \\
        --windows 5 --dt-ps 10 --box-half 32 \\
        --out solvent_density_timeseries.png

--solute-body-mask: one AMBER mask per solute body (repeat per body) --
only used to compute each frame's solute centroid, not plotted.
--solvent-mask: one representative-atom AMBER mask per solvent species to
include (repeat for multiple solvents -- all pooled into one density).
--anchor-mask: a single-molecule mask used only to resolve periodic-boundary
imaging consistently each frame (passed to cpptraj `autoimage anchor ...
origin`) -- does NOT bias the plotted density, since density is recentered
on the solute centroid afterward.
"""
import argparse
import subprocess
import tempfile
import os
import numpy as np
import scipy.io
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize


def extract(prmtop, traj, body_masks, solvent_masks, anchor_mask, workdir):
    com_out = os.path.join(workdir, "com.dat")
    nc_out = os.path.join(workdir, "solvent.nc")
    # NOT (:WAT@O|:MOH@O1) -- NOT (:WAT@O)|(:MOH@O1): cpptraj's strip mask
    # parser throws "Mask::ToRPN: unbalanced parentheses" when a `!(...)`
    # wraps an inner expression that itself has parens around each term,
    # even though the parens are balanced by any normal count. Confirmed via
    # a minimal reproduction (fails for the real masks below, succeeds with
    # the same masks and no inner parens) -- same class of cpptraj mask-
    # parser bug as the comma/distance-operator ones found earlier in this
    # project (see make_qmmm_frame.py).
    solvent_expr = "|".join(solvent_masks)
    lines = [f"parm {prmtop}", f"trajin {traj}", f"autoimage anchor {anchor_mask} origin"]
    for i, mask in enumerate(body_masks):
        lines.append(f"vector B{i} {mask} out {com_out} center")
    lines.append(f"strip !({solvent_expr})")
    lines.append(f"trajout {nc_out} netcdf")
    lines.append("run")
    script = os.path.join(workdir, "extract.in")
    with open(script, "w") as f:
        f.write("\n".join(lines) + "\n")
    result = subprocess.run(["cpptraj", "-i", script], capture_output=True, text=True)
    if result.returncode != 0 or not os.path.exists(nc_out):
        print(result.stdout)
        print(result.stderr)
        raise RuntimeError(f"cpptraj solvent extraction failed for {traj}")
    return com_out, nc_out


def load_centroid(com_path, n_bodies):
    data = np.loadtxt(com_path, comments="#")
    frames = data[:, 0]
    rest = data[:, 1:]
    bodies = rest.reshape(len(frames), n_bodies, 6)[:, :, :3]
    return frames, bodies.mean(axis=1)  # (n_frames, 3)


def load_solvent(nc_path):
    nc = scipy.io.netcdf_file(nc_path, "r", mmap=False)
    return nc.variables["coordinates"].data.copy()  # (n_frames, n_atoms, 3)


def window_hist(frames, solvent_xyz, centroid, n_windows, box_half, bins):
    n_frames = len(frames)
    edges = np.linspace(0, n_frames, n_windows + 1).astype(int)
    out, win_ranges = [], []
    for i in range(n_windows):
        lo, hi = edges[i], edges[i + 1]
        rel = solvent_xyz[lo:hi] - centroid[lo:hi, None, :]
        xs, ys = rel[:, :, 0].ravel(), rel[:, :, 1].ravel()
        h, _, _ = np.histogram2d(xs, ys, bins=bins, range=[[-box_half, box_half], [-box_half, box_half]])
        out.append(h)
        win_ranges.append((lo, hi))
    return out, win_ranges


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--system", nargs=3, action="append", required=True,
                     metavar=("LABEL", "PRMTOP", "TRAJ"))
    ap.add_argument("--solute-body-mask", action="append", required=True,
                     help="AMBER mask for one solute body (for centroid only); repeat per body")
    ap.add_argument("--solvent-mask", action="append", required=True,
                     help="AMBER mask for a solvent's representative atom; repeat per solvent type")
    ap.add_argument("--anchor-mask", required=True)
    ap.add_argument("--cmap", default="viridis")
    ap.add_argument("--scale", choices=["linear", "log"], default="linear",
                     help="solvent density varies only modestly (bulk vs. a depleted pocket near a "
                     "solute cluster, not orders of magnitude like the sparse solute case) -- linear "
                     "shows that contrast; log washes it out under the huge uniform bulk-density value")
    ap.add_argument("--windows", type=int, default=5)
    ap.add_argument("--dt-ps", type=float, default=10.0)
    ap.add_argument("--box-half", type=float, default=32.0)
    ap.add_argument("--bins", type=int, default=40)
    ap.add_argument("--title", default="Solvent COM density over time\n(relative to the solutes' own collective center)")
    ap.add_argument("--out", default="solvent_density_timeseries.png")
    args = ap.parse_args()

    n_bodies = len(args.solute_body_mask)
    all_windows, all_labels, all_frames = [], [], []
    with tempfile.TemporaryDirectory() as workdir:
        for label, prmtop, traj in args.system:
            com_out, nc_out = extract(prmtop, traj, args.solute_body_mask, args.solvent_mask,
                                       args.anchor_mask, workdir)
            frames, centroid = load_centroid(com_out, n_bodies)
            solvent_xyz = load_solvent(nc_out)
            hs, win_ranges = window_hist(frames, solvent_xyz, centroid, args.windows, args.box_half, args.bins)
            all_windows.append(hs)
            all_labels.append(label)
            all_frames.append((frames, win_ranges))
            os.remove(nc_out)

    global_max = max(h.max() for hs in all_windows for h in hs)
    if args.scale == "log":
        norm = LogNorm(vmin=1, vmax=global_max)
    else:
        norm = Normalize(vmin=0, vmax=global_max)

    n_sys = len(all_windows)
    row_h = min(3.1, 2.6 * 1.15 * args.windows / max(n_sys, 1))
    fig_w, fig_h = 2.6 * args.windows, row_h * n_sys
    fig, axes = plt.subplots(n_sys, args.windows, figsize=(fig_w, fig_h), squeeze=False)
    top_margin_in = 0.85
    fig.subplots_adjust(top=1 - top_margin_in / fig_h, left=0.07, right=0.86,
                         wspace=0.08, hspace=0.25)

    im = None
    for row, (hs, label, (frames, win_ranges)) in enumerate(zip(all_windows, all_labels, all_frames)):
        for i, (h, (lo, hi)) in enumerate(zip(hs, win_ranges)):
            ax = axes[row][i]
            im = ax.imshow(h.T, origin="lower",
                            extent=[-args.box_half, args.box_half, -args.box_half, args.box_half],
                            cmap=args.cmap, aspect="equal", norm=norm)
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

    cax = fig.add_axes([0.875, 0.20, 0.013, 0.55])
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(f"solvent density, {args.scale} scale\n(counts per bin, pooled over window)", fontsize=8.5)

    fig.suptitle(args.title, fontsize=13, x=0.44)
    fig.savefig(args.out, dpi=160)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
