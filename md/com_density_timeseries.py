#!/usr/bin/env python3
"""Spatial density heatmap of solute bodies' centers of mass, binned into
time windows, to visualize aggregation/condensation over the course of an MD
trajectory, distinguishing species WITHOUT blending colors (a hue blend of
two densities is hard to read at a glance -- flat regions of a single solid
blended tone don't obviously parse as "a mix of two things").

Instead: one species is the FIELD -- a standard single-colormap density
heatmap with a normal colorbar -- and the other species is MARKED --
its individual body positions plotted directly as scatter points on top.
This is also the natural choice whenever one species has very few bodies
(e.g. a single catalyst cluster against ten photosensitizer copies): a
"density" of one body is a weak concept, whereas marking its actual
position(s) is direct and unambiguous.

For each frame, extracts the COM of every given body (via cpptraj `vector
... center`), re-centers all bodies on their own collective centroid each
frame (so what's plotted is genuine spread/condensation relative to the
solute assembly's own center -- NOT relative to one fixed/pinned reference
body, which would be a reference-choice artifact), then pools the recentered
positions within each time window.

One or more systems (trajectories) can be given as separate rows, e.g. to
compare several compositions side by side on the same time-window and
color axes.

Requires cpptraj (AmberTools) on PATH.

Usage (PS density as the field, CAT positions marked -- this project's
convention: CAT=forest-green marker, PS=viridis field):
    python3 com_density_timeseries.py \\
        --system "0% water" PRMTOP1 TRAJ1 \\
        --system "100% water" PRMTOP2 TRAJ2 \\
        --field-species PS --mark-species CAT --mark-color "#386641" \\
        --body-mask CAT ":1-10" --body-mask CAT ":11-20" ... \\
        --body-mask PS ":101" --body-mask PS ":102" ... \\
        --anchor-mask ":1-10" \\
        --windows 5 --dt-ps 10 --box-half 32 \\
        --out com_density_timeseries.png

--body-mask: SPECIES followed by one AMBER mask for one solute body's COM
(repeat per body per species).
--field-species / --mark-species: which species gets the colormap density
vs. which gets scatter-marked. Exactly one field species; one or more mark
species (repeat --mark-species/--mark-color together per marked species).
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
from matplotlib.lines import Line2D


def extract_com(prmtop, traj, all_masks, anchor_mask, workdir):
    out = os.path.join(workdir, "com.dat")
    lines = [f"parm {prmtop}", f"trajin {traj}",
             f"autoimage anchor {anchor_mask} origin"]
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


def load(path, n_bodies):
    data = np.loadtxt(path, comments="#")
    frames = data[:, 0]
    rest = data[:, 1:]
    bodies = rest.reshape(len(frames), n_bodies, 6)[:, :, :3]
    centroid = bodies.mean(axis=1, keepdims=True)
    return frames, bodies - centroid


def window_data(frames, bodies, species_of_body, field_species, mark_species_list,
                 n_windows, box_half, bins):
    """Per time window: field species' 2D histogram, plus each marked
    species' raw (x, y) scatter points (strided for readability)."""
    n_frames = len(frames)
    edges = np.linspace(0, n_frames, n_windows + 1).astype(int)
    field_idx = [j for j, s in enumerate(species_of_body) if s == field_species]
    mark_idx = {sp: [j for j, s in enumerate(species_of_body) if s == sp] for sp in mark_species_list}
    out, win_ranges = [], []
    for i in range(n_windows):
        lo, hi = edges[i], edges[i + 1]
        sl = bodies[lo:hi]  # (nf, n_bodies, 3)
        xs = sl[:, field_idx, 0].ravel()
        ys = sl[:, field_idx, 1].ravel()
        h, _, _ = np.histogram2d(xs, ys, bins=bins, range=[[-box_half, box_half], [-box_half, box_half]])
        marks = {}
        n_window_frames = hi - lo
        for sp in mark_species_list:
            n_bodies_sp = max(len(mark_idx[sp]), 1)
            stride = max(1, (n_window_frames * len(mark_idx[sp])) // 400)
            pts = sl[::stride][:, mark_idx[sp], :2].reshape(-1, 2)
            marks[sp] = pts
        out.append((h, marks))
        win_ranges.append((lo, hi))
    return out, win_ranges


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--system", nargs=3, action="append", required=True,
                     metavar=("LABEL", "PRMTOP", "TRAJ"),
                     help="one row of the figure; repeat for multiple systems")
    ap.add_argument("--body-mask", nargs=2, action="append", required=True,
                     metavar=("SPECIES", "MASK"),
                     help="AMBER mask for one solute body's COM, tagged with its species; repeat per body")
    ap.add_argument("--field-species", required=True,
                     help="species shown as the colormap density (should have enough bodies for a real density)")
    ap.add_argument("--mark-species", action="append", default=[],
                     help="species shown as scatter-marked points; repeat per marked species (optional)")
    ap.add_argument("--mark-color", action="append", default=[],
                     help="hex color for the corresponding --mark-species, in the same order")
    ap.add_argument("--cmap", default="viridis")
    ap.add_argument("--anchor-mask", required=True,
                     help="single-molecule mask for periodic-image resolution (does not bias the plot)")
    ap.add_argument("--windows", type=int, default=5)
    ap.add_argument("--dt-ps", type=float, default=10.0, help="ps between saved frames")
    ap.add_argument("--box-half", type=float, default=32.0, help="plot range, +/- Angstrom")
    ap.add_argument("--bins", type=int, default=40)
    ap.add_argument("--title", default="Solute-body COM density over time\n"
                     "(relative to the solutes' own collective center)")
    ap.add_argument("--out", default="com_density_timeseries.png")
    args = ap.parse_args()

    if len(args.mark_species) != len(args.mark_color):
        raise ValueError("--mark-species and --mark-color must be given the same number of times")
    mark_colors = dict(zip(args.mark_species, args.mark_color))

    species_of_body = [sp for sp, _ in args.body_mask]
    all_masks = [mask for _, mask in args.body_mask]
    n_bodies = len(all_masks)

    all_windows, all_labels, all_frames = [], [], []
    with tempfile.TemporaryDirectory() as workdir:
        for label, prmtop, traj in args.system:
            com_file = extract_com(prmtop, traj, all_masks, args.anchor_mask, workdir)
            frames, bodies = load(com_file, n_bodies)
            wins, win_ranges = window_data(frames, bodies, species_of_body, args.field_species,
                                            args.mark_species, args.windows, args.box_half, args.bins)
            all_windows.append(wins)
            all_labels.append(label)
            all_frames.append((frames, win_ranges))

    global_max = max(h.max() for wins in all_windows for h, _ in wins)
    norm = LogNorm(vmin=1, vmax=global_max)

    n_sys = len(all_windows)
    # Row height: weasyprint (and most report pipelines) stretch the image to
    # the page's full content WIDTH, so only the image's ASPECT RATIO (not
    # its absolute size) controls the rendered height -- a tall aspect ratio
    # (many rows) can render taller than one page even though it looked fine
    # at this script's own native resolution, pushing the figure caption
    # onto a separate page. Cap the row height so height/width stays <=~1.15
    # once there are more than a few rows, trading a little panel padding
    # for a figure (+ its caption) that reliably fits one page together.
    row_h = min(3.1, 2.6 * 1.15 * args.windows / max(n_sys, 1))
    fig_w, fig_h = 2.6 * args.windows, row_h * n_sys
    fig, axes = plt.subplots(n_sys, args.windows, figsize=(fig_w, fig_h), squeeze=False)
    top_margin_in = 0.85
    fig.subplots_adjust(top=1 - top_margin_in / fig_h, left=0.07, right=0.86,
                         wspace=0.08, hspace=0.25)

    im = None
    marker_size = 9
    for row, (wins, label, (frames, win_ranges)) in enumerate(zip(all_windows, all_labels, all_frames)):
        for i, ((h, marks), (lo, hi)) in enumerate(zip(wins, win_ranges)):
            ax = axes[row][i]
            im = ax.imshow(h.T, origin="lower",
                            extent=[-args.box_half, args.box_half, -args.box_half, args.box_half],
                            cmap=args.cmap, aspect="equal", norm=norm)
            for sp in args.mark_species:
                pts = marks[sp]
                if len(pts):
                    ax.scatter(pts[:, 0], pts[:, 1], s=marker_size, c=mark_colors[sp],
                               edgecolors="white", linewidths=0.3, alpha=0.75, zorder=3)
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
    cbar.set_label(f"{args.field_species} COM density, log scale\n(counts per bin, pooled over window)", fontsize=8.5)

    legend_handles = [Line2D([0], [0], marker="o", linestyle="none", markersize=7,
                              markerfacecolor=mark_colors[sp], markeredgecolor="white",
                              markeredgewidth=0.4, label=f"{sp} position")
                       for sp in args.mark_species]
    fig.legend(handles=legend_handles, loc="upper right", bbox_to_anchor=(0.995, 0.995),
               ncol=len(args.mark_species), fontsize=10, frameon=False)

    fig.suptitle(args.title, fontsize=13, x=0.44)
    fig.savefig(args.out, dpi=160)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
