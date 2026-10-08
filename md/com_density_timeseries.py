#!/usr/bin/env python3
"""Spatial density heatmap of multiple solute bodies' centers of mass, binned
into time windows, to visualize aggregation/condensation over the course of
an MD trajectory -- with each species shown in its own color, so mixed vs.
segregated aggregation is visually obvious, not just inferable from a single
pooled density.

For each frame, extracts the COM of every given body (via cpptraj `vector
... center`), re-centers all bodies on their own collective centroid each
frame (so what's plotted is genuine spread/condensation relative to the
solute assembly's own center -- NOT relative to one fixed/pinned reference
body, which would be a reference-choice artifact), then pools the recentered
positions within each time window into per-species 2D (x-y) histograms.

Each pixel's color is a bivariate encoding: HUE is the local species mix
(linear interpolation between each species' color, weighted by its share of
the local density -- pure color = that species dominates locally, a blend =
genuinely mixed), and BRIGHTNESS/alpha is the total density on a shared log
scale (faint/white = little data, vivid = a lot) -- so composition and
magnitude are both visible in one image without needing two separate plots.

One or more systems (trajectories) can be given as separate rows, e.g. to
compare several compositions side by side on the same time-window and
color axes.

Requires cpptraj (AmberTools) on PATH.

Usage (two species, matching this project's CAT=forest-green / PS=slate-blue
convention):
    python3 com_density_timeseries.py \\
        --system "0% water" PRMTOP1 TRAJ1 \\
        --system "100% water" PRMTOP2 TRAJ2 \\
        --species CAT "#386641" --body-mask CAT ":1-10" --body-mask CAT ":11-20" ... \\
        --species PS  "#7798ab" --body-mask PS  ":101" --body-mask PS ":102" ... \\
        --anchor-mask ":1-10" \\
        --windows 6 --dt-ps 10 --box-half 32 \\
        --out com_density_timeseries.png

--body-mask: SPECIES followed by one AMBER mask for one solute body's COM
(repeat per body per species).
--species: SPECIES label followed by its hex color (repeat per species;
2 is the common case but any number works).
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
from matplotlib.colors import LogNorm, to_rgb, rgb_to_hsv, hsv_to_rgb
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter


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


def species_histograms(frames, bodies, species_of_body, species_list, n_windows, box_half, bins):
    """Returns, per time window: dict species -> 2D histogram, plus the (lo,hi) frame range."""
    n_frames = len(frames)
    edges = np.linspace(0, n_frames, n_windows + 1).astype(int)
    out, win_ranges = [], []
    for i in range(n_windows):
        lo, hi = edges[i], edges[i + 1]
        sl = bodies[lo:hi]  # (nf, n_bodies, 3)
        hs = {}
        for sp in species_list:
            idx = [j for j, s in enumerate(species_of_body) if s == sp]
            xs = sl[:, idx, 0].ravel()
            ys = sl[:, idx, 1].ravel()
            h, _, _ = np.histogram2d(xs, ys, bins=bins,
                                      range=[[-box_half, box_half], [-box_half, box_half]])
            hs[sp] = h
        out.append(hs)
        win_ranges.append((lo, hi))
    return out, win_ranges


def composite_rgba(hs_by_species, species_colors, species_list, norm, smooth_sigma=0.8, saturation_boost=1.9):
    """Bivariate color: hue = local species mix (density-weighted blend of
    each species' flat color), brightness/alpha = total density (log scale,
    shared across the whole figure via `norm`). Lightly Gaussian-smoothed
    (the raw per-bin species fraction is "salt-and-pepper" speckled --
    individual bodies occupy discrete positions, so any one bin is often
    dominated by whichever single body happened to sample it -- smoothing
    reveals the underlying mixed-vs-segregated pattern instead of per-bin
    noise) and saturation-boosted (this project's brand colors are muted/
    desaturated by design for bar charts etc.; a straight blend of two
    muted colors reads as muddy gray rather than a legible two-color mix)."""
    shape = next(iter(hs_by_species.values())).shape
    smoothed = {sp: gaussian_filter(hs_by_species[sp], sigma=smooth_sigma) for sp in species_list}
    total = np.zeros(shape)
    for sp in species_list:
        total += smoothed[sp]
    nonzero = total > 0
    blended = np.zeros(shape + (3,))
    for sp in species_list:
        frac = np.zeros(shape)
        frac[nonzero] = smoothed[sp][nonzero] / total[nonzero]
        color = np.array(to_rgb(species_colors[sp]))
        blended += frac[..., None] * color[None, None, :]
    hsv = rgb_to_hsv(np.clip(blended, 0, 1))
    hsv[..., 1] = np.clip(hsv[..., 1] * saturation_boost, 0, 1)
    blended = hsv_to_rgb(hsv)
    alpha = np.zeros(shape)
    raw_total = sum(hs_by_species[sp] for sp in species_list)
    raw_total_smoothed = gaussian_filter(raw_total, sigma=smooth_sigma)
    alpha = np.clip(norm(np.clip(raw_total_smoothed, 1e-9, None)), 0, 1)
    # composite blended color over white background using alpha
    out_rgb = blended * alpha[..., None] + 1.0 * (1 - alpha[..., None])
    return np.clip(out_rgb, 0, 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--system", nargs=3, action="append", required=True,
                     metavar=("LABEL", "PRMTOP", "TRAJ"),
                     help="one row of the figure; repeat for multiple systems")
    ap.add_argument("--species", nargs=2, action="append", required=True,
                     metavar=("NAME", "HEXCOLOR"), help="repeat per species")
    ap.add_argument("--body-mask", nargs=2, action="append", required=True,
                     metavar=("SPECIES", "MASK"),
                     help="AMBER mask for one solute body's COM, tagged with its species; repeat per body")
    ap.add_argument("--anchor-mask", required=True,
                     help="single-molecule mask for periodic-image resolution (does not bias the plot)")
    ap.add_argument("--windows", type=int, default=6)
    ap.add_argument("--dt-ps", type=float, default=10.0, help="ps between saved frames")
    ap.add_argument("--box-half", type=float, default=32.0, help="plot range, +/- Angstrom")
    ap.add_argument("--bins", type=int, default=40)
    ap.add_argument("--title", default="Solute-body COM density over time, by species\n"
                     "(relative to the solutes' own collective center)")
    ap.add_argument("--out", default="com_density_timeseries.png")
    args = ap.parse_args()

    species_colors = {name: color for name, color in args.species}
    species_list = [name for name, _ in args.species]
    species_of_body = [sp for sp, _ in args.body_mask]
    all_masks = [mask for _, mask in args.body_mask]
    n_bodies = len(all_masks)

    all_hs, all_labels, all_frames = [], [], []
    with tempfile.TemporaryDirectory() as workdir:
        for label, prmtop, traj in args.system:
            com_file = extract_com(prmtop, traj, all_masks, args.anchor_mask, workdir)
            frames, bodies = load(com_file, n_bodies)
            hs, win_ranges = species_histograms(frames, bodies, species_of_body, species_list,
                                                 args.windows, args.box_half, args.bins)
            all_hs.append(hs)
            all_labels.append(label)
            all_frames.append((frames, win_ranges))

    global_max = max(sum(hs[sp] for sp in species_list).max()
                      for row_hs in all_hs for hs in row_hs)
    norm = LogNorm(vmin=1, vmax=global_max)

    n_sys = len(all_hs)
    fig_w, fig_h = 2.6 * args.windows, 3.1 * n_sys
    fig, axes = plt.subplots(n_sys, args.windows, figsize=(fig_w, fig_h), squeeze=False)
    top_margin_in = 0.85
    fig.subplots_adjust(top=1 - top_margin_in / fig_h, left=0.07, right=0.88,
                         wspace=0.08, hspace=0.25)

    for row, (hs_list, label, (frames, win_ranges)) in enumerate(zip(all_hs, all_labels, all_frames)):
        for i, (hs, (lo, hi)) in enumerate(zip(hs_list, win_ranges)):
            ax = axes[row][i]
            rgb = composite_rgba(hs, species_colors, species_list, norm)
            ax.imshow(np.transpose(rgb, (1, 0, 2)), origin="lower",
                       extent=[-args.box_half, args.box_half, -args.box_half, args.box_half],
                       aspect="equal")
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

    legend_handles = [Line2D([0], [0], marker="s", linestyle="none", markersize=10,
                              markerfacecolor=species_colors[sp], markeredgecolor="none", label=sp)
                       for sp in species_list]
    fig.legend(handles=legend_handles, loc="upper right", bbox_to_anchor=(0.995, 0.995),
               ncol=len(species_list), fontsize=10, frameon=False)

    fig.suptitle(args.title, fontsize=13, x=0.46)
    fig.text(0.90, 0.5, "color = local species mix; brightness = total density (log scale)",
              rotation=90, va="center", fontsize=8.5)
    fig.savefig(args.out, dpi=160)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
