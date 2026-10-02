#!/usr/bin/env python3
"""
Browse and export color palettes from Sanzo Wada's "A Dictionary of Colour
Combinations" (1933) for use in publication figures.

Data source: the open-sourced digitization of the book by Dain M. Blodorn
Kim, as republished by mattdesl at
https://github.com/mattdesl/dictionary-of-colour-combinations (MIT license).
That dataset provides 159 named colors (hex/RGB/CMYK/Lab), each tagged with
which of the book's 348 combinations (2-4 colors each) it belongs to; the
combinations themselves are reconstructed here from that tagging, following
the exact method documented in that repo's README.

The 348 original combinations are historical curated color pairs/triples/
quadruples -- they are a source of *aesthetically coordinated* palettes, not
inherently colorblind-safe or validated for data visualization. Use
--cvd-preview to get a quick (approximate) sanity check of how a chosen
combination looks under the three common forms of color vision deficiency
before using it as a *categorical/identity* palette in a real figure; this
is not a substitute for a dedicated scientific-palette validator (e.g. the
`dataviz` skill's `validate_palette.js`, which computes acutal CVD Delta-E
separations rather than just rendering a preview).

Usage:
    # list all 348 combinations, optionally filtered
    python3 sanzo_wada_palette.py list [--n-colors 2|3|4] [--search sienna]

    # show one combination's colors in the terminal
    python3 sanzo_wada_palette.py show 42

    # export a combination as a swatch figure + hex codes, ready to paste
    # into a matplotlib script or LaTeX document
    python3 sanzo_wada_palette.py export 42 --out my_palette [--cvd-preview]

    # grab a random combination with a given number of colors
    python3 sanzo_wada_palette.py random --n-colors 3

    # visualize many combinations at once as a browsable grid
    python3 sanzo_wada_palette.py gallery --n-colors 3 --limit 30 --out gallery

On first use the script downloads colors.json from GitHub and caches it
next to this file (sanzo_wada_colors.json); use --refresh to re-download.
"""
import argparse
import json
import os
import random
import sys
import urllib.request

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE_PATH = os.path.join(HERE, "sanzo_wada_colors.json")
SOURCE_URL = ("https://raw.githubusercontent.com/mattdesl/"
              "dictionary-of-colour-combinations/master/colors.json")

# approximate linear-RGB CVD simulation matrices (Brettel/Vienot-style,
# widely used for quick previews -- not a precise/validated transform)
CVD_MATRICES = {
    "protanopia": np.array([[0.567, 0.433, 0.000],
                             [0.558, 0.442, 0.000],
                             [0.000, 0.242, 0.758]]),
    "deuteranopia": np.array([[0.625, 0.375, 0.000],
                               [0.700, 0.300, 0.000],
                               [0.000, 0.300, 0.700]]),
    "tritanopia": np.array([[0.950, 0.050, 0.000],
                             [0.000, 0.433, 0.567],
                             [0.000, 0.475, 0.525]]),
}


def fetch_colors(refresh=False):
    if refresh or not os.path.exists(CACHE_PATH):
        print(f"Downloading color data from {SOURCE_URL} ...", file=sys.stderr)
        with urllib.request.urlopen(SOURCE_URL, timeout=20) as resp:
            data = resp.read()
        with open(CACHE_PATH, "wb") as f:
            f.write(data)
    with open(CACHE_PATH) as f:
        return json.load(f)


def build_combinations(colors):
    """Reconstruct {combination_id: [color, ...]} from each color's own
    'combinations' tag list, per the dataset's documented method."""
    combos = {}
    for color in colors:
        for cid in color["combinations"]:
            combos.setdefault(cid, []).append(color)
    return combos


def hex_to_rgb01(hexcode):
    hexcode = hexcode.lstrip("#")
    return tuple(int(hexcode[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


def rgb01_to_hex(rgb01):
    r, g, b = (int(round(np.clip(c, 0, 1) * 255)) for c in rgb01)
    return f"#{r:02x}{g:02x}{b:02x}"


# --- sRGB <-> CIELAB (D65), used for perceptually-uniform gradient
# interpolation and for farthest-point selection of distinct colors ---
_SRGB2XYZ = np.array([[0.4124564, 0.3575761, 0.1804375],
                       [0.2126729, 0.7151522, 0.0721750],
                       [0.0193339, 0.1191920, 0.9503041]])
_XYZ2SRGB = np.linalg.inv(_SRGB2XYZ)
_WHITE = np.array([95.047, 100.0, 108.883])  # D65 reference white


def rgb01_to_lab(rgb01):
    rgb = np.clip(np.asarray(rgb01, dtype=float), 0, 1)
    lin = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    xyz = (_SRGB2XYZ @ lin) * 100.0
    t = xyz / _WHITE
    delta = 6.0 / 29.0
    ft = np.where(t > delta ** 3, np.cbrt(t), t / (3 * delta ** 2) + 4.0 / 29.0)
    L = 116 * ft[1] - 16
    a = 500 * (ft[0] - ft[1])
    b = 200 * (ft[1] - ft[2])
    return np.array([L, a, b])


def lab_to_rgb01(lab):
    L, a, b = lab
    fy = (L + 16) / 116
    fx = fy + a / 500
    fz = fy - b / 200
    delta = 6.0 / 29.0

    def finv(t):
        return np.where(t > delta, t ** 3, 3 * delta ** 2 * (t - 4.0 / 29.0))

    xyz = np.array([finv(fx), finv(fy), finv(fz)]) * _WHITE / 100.0
    lin = _XYZ2SRGB @ xyz
    lin = np.clip(lin, 0, 1)
    srgb = np.where(lin <= 0.0031308, lin * 12.92, 1.055 * lin ** (1 / 2.4) - 0.055)
    return tuple(np.clip(srgb, 0, 1))


def interpolate_lab(anchor_hexes, n_steps):
    """Piecewise-linear interpolation through a sequence of anchor colors in
    CIELAB space (perceptually much smoother than interpolating raw RGB),
    returning n_steps evenly-spaced hex colors spanning the full sequence."""
    anchors_lab = np.array([rgb01_to_lab(hex_to_rgb01(h)) for h in anchor_hexes])
    t_anchors = np.linspace(0, 1, len(anchors_lab))
    t_samples = np.linspace(0, 1, n_steps)
    lab_samples = np.column_stack([
        np.interp(t_samples, t_anchors, anchors_lab[:, dim]) for dim in range(3)
    ])
    rgb_samples = [lab_to_rgb01(lab) for lab in lab_samples]
    return [rgb01_to_hex(rgb) for rgb in rgb_samples], rgb_samples


def simulate_cvd(rgb01, kind):
    M = CVD_MATRICES[kind]
    return tuple(np.clip(M @ np.array(rgb01), 0, 1))


def readable_text_color(rgb01):
    """Black or white, whichever contrasts better against this background
    (simple relative-luminance threshold)."""
    r, g, b = rgb01
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return "black" if luminance > 0.55 else "white"


def cmd_list(args, colors, combos):
    ids = sorted(combos)
    for cid in ids:
        members = combos[cid]
        if args.n_colors and len(members) != args.n_colors:
            continue
        names = [c["name"] for c in members]
        if args.search and not any(args.search.lower() in n.lower() for n in names):
            continue
        hexes = ",".join(c["hex"] for c in members)
        print(f"{cid:>4}  [{len(members)}]  {hexes}   {', '.join(names)}")


def cmd_show(args, colors, combos):
    if args.combination_id not in combos:
        sys.exit(f"No such combination id: {args.combination_id} (valid range 1-348)")
    for c in combos[args.combination_id]:
        print(f"  {c['hex']}   rgb{tuple(c['rgb'])}   {c['name']}")


def render_swatch(members, out_prefix, cvd_preview=False):
    n = len(members)
    rows = 1 + (3 if cvd_preview else 0)
    fig, axes = plt.subplots(rows, 1, figsize=(1.6 * n + 1.0, 1.7 * rows),
                              squeeze=False)
    axes = axes[:, 0]

    row_defs = [("original", None)]
    if cvd_preview:
        row_defs += [("protanopia", "protanopia"),
                     ("deuteranopia", "deuteranopia"),
                     ("tritanopia", "tritanopia")]

    for ax, (label, kind) in zip(axes, row_defs):
        for i, c in enumerate(members):
            rgb = hex_to_rgb01(c["hex"])
            if kind:
                rgb = simulate_cvd(rgb, kind)
            ax.add_patch(Rectangle((i, 0), 1, 1, facecolor=rgb, edgecolor="white",
                                    linewidth=1.5))
            if label == "original":
                ax.text(i + 0.5, -0.15, c["hex"], ha="center", va="top", fontsize=8)
                ax.text(i + 0.5, -0.32, c["name"], ha="center", va="top", fontsize=7,
                         color="0.3")
        ax.set_xlim(0, n)
        ax.set_ylim(-0.5 if label == "original" else 0, 1)
        ax.axis("off")
        ax.text(-0.15, 0.5, label, ha="right", va="center", fontsize=8, color="0.4")

    fig.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=300)
    fig.savefig(f"{out_prefix}.pdf")
    plt.close(fig)


def render_gallery(filtered, out_prefix, max_colors):
    """filtered: list of (combination_id, members), one row per combination."""
    rows = len(filtered)
    fig_h = max(0.42 * rows + 0.6, 1.5)
    fig_w = 1.25 * max_colors + 1.3
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    for row_i, (cid, members) in enumerate(filtered):
        y = rows - row_i - 1  # first result at the top
        for i, c in enumerate(members):
            rgb = hex_to_rgb01(c["hex"])
            ax.add_patch(Rectangle((i, y), 1, 0.85, facecolor=rgb,
                                    edgecolor="white", linewidth=0.8))
            ax.text(i + 0.5, y + 0.425, c["hex"], ha="center", va="center",
                     fontsize=6, color=readable_text_color(rgb), family="monospace")
        ax.text(-0.15, y + 0.425, f"#{cid}", ha="right", va="center",
                 fontsize=7, color="0.35", family="monospace")

    ax.set_xlim(0, max_colors)
    ax.set_ylim(0, rows)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=300)
    fig.savefig(f"{out_prefix}.pdf")
    plt.close(fig)


def cmd_gallery(args, colors, combos):
    ids = sorted(combos)
    filtered = []
    for cid in ids:
        members = combos[cid]
        if args.n_colors and len(members) != args.n_colors:
            continue
        names = [c["name"] for c in members]
        if args.search and not any(args.search.lower() in n.lower() for n in names):
            continue
        filtered.append((cid, members))

    if not filtered:
        sys.exit("No combinations match the given filters")
    if args.limit:
        filtered = filtered[:args.limit]

    max_colors = max(len(m) for _, m in filtered)
    render_gallery(filtered, args.out, max_colors)
    print(f"Wrote {args.out}.png, {args.out}.pdf  ({len(filtered)} combinations, "
          f"ids: {filtered[0][0]}-{filtered[-1][0]})")


def render_gradient(anchor_hexes, step_hexes, out_prefix):
    """Top: a smooth, continuous colormap bar (256-step Lab interpolation)
    suitable as a publication colorbar. Bottom: the discrete requested
    steps, labeled with hex codes, for use as a categorical/ordinal list."""
    smooth_hexes, _ = interpolate_lab(anchor_hexes, 256)
    smooth_rgb = np.array([hex_to_rgb01(h) for h in smooth_hexes])[None, :, :]

    n = len(step_hexes)
    fig, (ax_bar, ax_steps) = plt.subplots(
        2, 1, figsize=(max(6, 0.9 * n), 2.6), gridspec_kw={"height_ratios": [1, 1.3]})

    ax_bar.imshow(smooth_rgb, aspect="auto", extent=[0, n, 0, 1])
    ax_bar.set_xlim(0, n)
    ax_bar.axis("off")
    ax_bar.set_title("continuous", fontsize=8, color="0.4", loc="left")

    for i, h in enumerate(step_hexes):
        rgb = hex_to_rgb01(h)
        ax_steps.add_patch(Rectangle((i, 0), 1, 1, facecolor=rgb, edgecolor="white",
                                      linewidth=1.5))
        ax_steps.text(i + 0.5, 0.5, h, ha="center", va="center", fontsize=7,
                       color=readable_text_color(rgb), family="monospace")
    ax_steps.set_xlim(0, n)
    ax_steps.set_ylim(0, 1)
    ax_steps.axis("off")
    ax_steps.set_title(f"{n} discrete steps", fontsize=8, color="0.4", loc="left")

    fig.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=300)
    fig.savefig(f"{out_prefix}.pdf")
    plt.close(fig)


def write_gradient_snippets(anchor_hexes, step_hexes, out_prefix):
    with open(f"{out_prefix}_python.txt", "w") as f:
        f.write("# Sanzo Wada gradient -- paste into a matplotlib script\n")
        f.write(f"# Lab-interpolated through: {', '.join(anchor_hexes)}\n\n")
        f.write("from matplotlib.colors import LinearSegmentedColormap\n")
        f.write(f"GRADIENT_STEPS = {step_hexes!r}\n")
        f.write(f"GRADIENT_CMAP = LinearSegmentedColormap.from_list("
                f"'sanzo_wada_gradient', {anchor_hexes!r}, N=256)\n")
        f.write("# GRADIENT_CMAP interpolates in raw RGB; for the exact Lab-space\n")
        f.write("# curve this script used, sample interpolate_lab() directly instead.\n")

    with open(f"{out_prefix}_latex.tex", "w") as f:
        f.write("% Sanzo Wada gradient steps -- paste into your LaTeX preamble\n")
        f.write("\\usepackage{xcolor}\n")
        for i, h in enumerate(step_hexes):
            f.write(f"\\definecolor{{GradStep{i}}}{{HTML}}{{{h.lstrip('#').upper()}}}\n")

    with open(f"{out_prefix}.json", "w") as f:
        json.dump({"anchors": anchor_hexes, "steps": step_hexes}, f, indent=2)


def cmd_gradient(args, colors, combos):
    if args.from_combination:
        if args.from_combination not in combos:
            sys.exit(f"No such combination id: {args.from_combination} (valid range 1-348)")
        anchor_hexes = [c["hex"] for c in combos[args.from_combination]]
    elif args.colors:
        anchor_hexes = [h.strip() for h in args.colors.split(",")]
    else:
        sys.exit("Provide anchor colors via --colors '#hex,#hex,...' or --from-combination <id>")

    if len(anchor_hexes) < 2:
        sys.exit("Need at least 2 anchor colors to make a gradient")

    step_hexes, _ = interpolate_lab(anchor_hexes, args.steps)
    render_gradient(anchor_hexes, step_hexes, args.out)
    write_gradient_snippets(anchor_hexes, step_hexes, args.out)
    print(f"Anchors: {', '.join(anchor_hexes)}")
    print(f"{args.steps} steps: {', '.join(step_hexes)}")
    print(f"Wrote {args.out}.png, {args.out}.pdf, {args.out}_python.txt, "
          f"{args.out}_latex.tex, {args.out}.json")


def farthest_point_select(colors, n):
    """Greedy farthest-point sampling in CIELAB space: start from the color
    most extreme relative to the full set's centroid, then repeatedly add
    whichever remaining color maximizes its minimum distance to everything
    already chosen -- a simple, deterministic way to get n mutually
    well-separated colors out of the 159-color set."""
    labs = np.array([c["lab"] for c in colors])
    centroid = labs.mean(axis=0)
    chosen = [int(np.argmax(np.linalg.norm(labs - centroid, axis=1)))]
    while len(chosen) < n:
        d = np.linalg.norm(labs[:, None, :] - labs[chosen][None, :, :], axis=-1).min(axis=1)
        d[chosen] = -np.inf
        chosen.append(int(np.argmax(d)))
    return [colors[i] for i in chosen]


def min_pairwise_lab_distance(labs):
    labs = np.asarray(labs)
    d = np.linalg.norm(labs[:, None, :] - labs[None, :, :], axis=-1)
    np.fill_diagonal(d, np.inf)
    return float(d.min())


def cmd_distinct(args, colors, combos):
    if args.n < 2:
        sys.exit("--n must be at least 2")
    if args.n > len(colors):
        sys.exit(f"--n cannot exceed the {len(colors)} available colors")

    picked = farthest_point_select(colors, args.n)
    labs = np.array([c["lab"] for c in picked])

    print(f"{args.n} maximally-distinct colors (greedy farthest-point, CIELAB):")
    for c in picked:
        print(f"  {c['hex']}   {c['name']}")
    print(f"Min pairwise CIE76 Delta-E (normal vision): {min_pairwise_lab_distance(labs):.1f}"
          "  (as a rule of thumb, >10 is usually safely distinguishable, >20 very safe)")

    for kind in ("protanopia", "deuteranopia", "tritanopia"):
        sim_labs = [rgb01_to_lab(simulate_cvd(hex_to_rgb01(c["hex"]), kind)) for c in picked]
        print(f"Min pairwise Delta-E under simulated {kind}: "
              f"{min_pairwise_lab_distance(sim_labs):.1f}")

    render_swatch(picked, args.out, cvd_preview=args.cvd_preview)
    write_snippets(picked, args.out)
    print(f"Wrote {args.out}.png, {args.out}.pdf, {args.out}_python.txt, "
          f"{args.out}_latex.tex, {args.out}.json")


def write_snippets(members, out_prefix):
    hexes = [c["hex"] for c in members]
    names = [c["name"] for c in members]

    with open(f"{out_prefix}_python.txt", "w") as f:
        f.write("# Sanzo Wada palette -- paste into a matplotlib script\n")
        f.write("PALETTE = [\n")
        for h, n in zip(hexes, names):
            f.write(f'    "{h}",  # {n}\n')
        f.write("]\n")

    with open(f"{out_prefix}_latex.tex", "w") as f:
        f.write("% Sanzo Wada palette -- paste into your LaTeX preamble\n")
        f.write("\\usepackage{xcolor}\n")
        for i, (h, n) in enumerate(zip(hexes, names)):
            safe_name = "".join(ch for ch in n if ch.isalpha())
            f.write(f"\\definecolor{{{safe_name}}}{{HTML}}{{{h.lstrip('#').upper()}}}"
                    f"  % {n}\n")

    with open(f"{out_prefix}.json", "w") as f:
        json.dump(members, f, indent=2)


def cmd_export(args, colors, combos):
    if args.combination_id not in combos:
        sys.exit(f"No such combination id: {args.combination_id} (valid range 1-348)")
    members = combos[args.combination_id]
    render_swatch(members, args.out, cvd_preview=args.cvd_preview)
    write_snippets(members, args.out)
    print(f"Wrote {args.out}.png, {args.out}.pdf, {args.out}_python.txt, "
          f"{args.out}_latex.tex, {args.out}.json")


def cmd_random(args, colors, combos):
    candidates = [cid for cid, members in combos.items()
                  if not args.n_colors or len(members) == args.n_colors]
    if not candidates:
        sys.exit(f"No combinations with n_colors={args.n_colors}")
    cid = random.choice(candidates)
    print(f"Combination {cid}:")
    cmd_show(argparse.Namespace(combination_id=cid), colors, combos)
    if args.out:
        render_swatch(combos[cid], args.out, cvd_preview=args.cvd_preview)
        write_snippets(combos[cid], args.out)
        print(f"Wrote {args.out}.png, {args.out}.pdf, {args.out}_python.txt, "
              f"{args.out}_latex.tex, {args.out}.json")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--refresh", action="store_true",
                    help="re-download the color data instead of using the cache")
    sub = p.add_subparsers(dest="command", required=True)

    p_list = sub.add_parser("list", help="list combinations")
    p_list.add_argument("--n-colors", type=int, choices=[2, 3, 4])
    p_list.add_argument("--search", help="filter by color name substring")

    p_show = sub.add_parser("show", help="show one combination's colors")
    p_show.add_argument("combination_id", type=int)

    p_export = sub.add_parser("export", help="export a combination as figure + code snippets")
    p_export.add_argument("combination_id", type=int)
    p_export.add_argument("--out", default="sanzo_wada_palette", help="output basename")
    p_export.add_argument("--cvd-preview", action="store_true",
                           help="add protanopia/deuteranopia/tritanopia preview rows")

    p_random = sub.add_parser("random", help="pick a random combination")
    p_random.add_argument("--n-colors", type=int, choices=[2, 3, 4])
    p_random.add_argument("--out", help="also export this random pick (basename)")
    p_random.add_argument("--cvd-preview", action="store_true")

    p_gallery = sub.add_parser(
        "gallery", help="render a grid of many combinations at once for visual browsing")
    p_gallery.add_argument("--n-colors", type=int, choices=[2, 3, 4])
    p_gallery.add_argument("--search", help="filter by color name substring")
    p_gallery.add_argument("--limit", type=int, default=30,
                            help="max combinations to include (default 30)")
    p_gallery.add_argument("--out", default="sanzo_wada_gallery", help="output basename")

    p_gradient = sub.add_parser(
        "gradient", help="build a smooth, perceptually-interpolated gradient/colormap")
    g_src = p_gradient.add_mutually_exclusive_group(required=True)
    g_src.add_argument("--colors", help="comma-separated anchor hex colors, e.g. '#112233,#ffee00'")
    g_src.add_argument("--from-combination", type=int,
                        help="use an existing combination's colors as the anchors, in order")
    p_gradient.add_argument("--steps", type=int, default=8,
                             help="number of discrete steps to also output (default 8)")
    p_gradient.add_argument("--out", default="sanzo_wada_gradient", help="output basename")

    p_distinct = sub.add_parser(
        "distinct", help="pick N maximally-separated colors for a categorical/identity palette")
    p_distinct.add_argument("--n", type=int, required=True,
                             help="how many distinct colors to pick")
    p_distinct.add_argument("--out", default="sanzo_wada_distinct", help="output basename")
    p_distinct.add_argument("--cvd-preview", action="store_true",
                             help="add protanopia/deuteranopia/tritanopia preview rows")

    args = p.parse_args()
    colors = fetch_colors(refresh=args.refresh)
    combos = build_combinations(colors)

    {"list": cmd_list, "show": cmd_show, "export": cmd_export,
     "random": cmd_random, "gallery": cmd_gallery, "gradient": cmd_gradient,
     "distinct": cmd_distinct}[args.command](args, colors, combos)


if __name__ == "__main__":
    main()
