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

    # sequential: one hue, light -> dark (magnitude data, e.g. a heatmap)
    python3 sanzo_wada_palette.py sequential --color "#1c4286" --steps 8

    # diverging: two hues + a neutral gray midpoint (polarity, e.g. +/- data)
    python3 sanzo_wada_palette.py diverging --colors "#cc1236,#00978d" --steps 9

    # qualitative: N maximally-separated colors (categorical/identity data)
    python3 sanzo_wada_palette.py qualitative --n 6

    # test any palette's distinguishability under color vision deficiency
    python3 sanzo_wada_palette.py cvd-test --from-combination 42
    python3 sanzo_wada_palette.py cvd-test --colors "#cc1236,#00978d,#e2b540"

    # render example line/scatter/heatmap plots using one combination's colors
    python3 sanzo_wada_palette.py demo --from-combination 121 --out demo

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
from matplotlib.colors import ListedColormap

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


def sequential_anchors(hue_hex):
    """One hue, light -> dark: build a 3-point Lab path (light tint, the
    hue itself, dark shade) from a single base color, per the standard
    'sequential = one hue' convention for magnitude data."""
    L, a, b = rgb01_to_lab(hex_to_rgb01(hue_hex))
    light = lab_to_rgb01(np.array([95.0, a * 0.12, b * 0.12]))
    dark = lab_to_rgb01(np.array([22.0, a * 1.05, b * 1.05]))
    return [rgb01_to_hex(light), hue_hex, rgb01_to_hex(dark)]


def diverging_anchors(hex1, hex2, midpoint_hex=None):
    """Two hues + a neutral gray midpoint, per the standard 'diverging'
    convention for polarity data (e.g. positive/negative, above/below a
    reference) -- never a hue at the midpoint."""
    mid = midpoint_hex or "#f2f2f0"
    return [hex1, mid, hex2]


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


def cmd_sequential(args, colors, combos):
    if args.from_combination:
        if args.from_combination not in combos:
            sys.exit(f"No such combination id: {args.from_combination} (valid range 1-348)")
        hue_hex = combos[args.from_combination][0]["hex"]
    else:
        hue_hex = args.color.strip()

    anchor_hexes = sequential_anchors(hue_hex)
    step_hexes, _ = interpolate_lab(anchor_hexes, args.steps)
    render_gradient(anchor_hexes, step_hexes, args.out)
    write_gradient_snippets(anchor_hexes, step_hexes, args.out)
    print(f"Sequential ramp from hue {hue_hex} (light -> dark)")
    print(f"{args.steps} steps: {', '.join(step_hexes)}")
    print(f"Wrote {args.out}.png, {args.out}.pdf, {args.out}_python.txt, "
          f"{args.out}_latex.tex, {args.out}.json")


def cmd_diverging(args, colors, combos):
    if args.from_combination:
        if args.from_combination not in combos:
            sys.exit(f"No such combination id: {args.from_combination} (valid range 1-348)")
        members = combos[args.from_combination]
        if len(members) < 2:
            sys.exit("That combination has fewer than 2 colors; need 2 poles for diverging")
        hex1, hex2 = members[0]["hex"], members[-1]["hex"]
    else:
        parts = [h.strip() for h in args.colors.split(",")]
        if len(parts) != 2:
            sys.exit("--colors must give exactly 2 hex colors (the two poles), e.g. '#aa0000,#0055aa'")
        hex1, hex2 = parts

    anchor_hexes = diverging_anchors(hex1, hex2, args.midpoint)
    step_hexes, _ = interpolate_lab(anchor_hexes, args.steps)
    render_gradient(anchor_hexes, step_hexes, args.out)
    write_gradient_snippets(anchor_hexes, step_hexes, args.out)
    print(f"Diverging ramp: {hex1} -> {anchor_hexes[1]} (neutral midpoint) -> {hex2}")
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


def print_cvd_diagnostics(hexes):
    """Print the minimum pairwise CIE76 Delta-E among a set of hex colors,
    both for normal vision and simulated under each CVD type -- a real
    numeric distinguishability check, not just a rendered preview."""
    labs = np.array([rgb01_to_lab(hex_to_rgb01(h)) for h in hexes])
    print(f"Min pairwise Delta-E (normal vision): {min_pairwise_lab_distance(labs):.1f}"
          "  (rule of thumb: >10 usually safely distinguishable, >20 very safe)")
    for kind in ("protanopia", "deuteranopia", "tritanopia"):
        sim_labs = [rgb01_to_lab(simulate_cvd(hex_to_rgb01(h), kind)) for h in hexes]
        verdict = ""
        d = min_pairwise_lab_distance(sim_labs)
        if d < 10:
            verdict = "  <-- WARNING: two or more colors likely collide here"
        print(f"Min pairwise Delta-E under simulated {kind}: {d:.1f}{verdict}")


def cmd_qualitative(args, colors, combos):
    if args.n < 2:
        sys.exit("--n must be at least 2")
    if args.n > len(colors):
        sys.exit(f"--n cannot exceed the {len(colors)} available colors")

    picked = farthest_point_select(colors, args.n)

    print(f"{args.n} maximally-distinct colors (greedy farthest-point, CIELAB):")
    for c in picked:
        print(f"  {c['hex']}   {c['name']}")
    print_cvd_diagnostics([c["hex"] for c in picked])

    render_swatch(picked, args.out, cvd_preview=args.cvd_preview)
    write_snippets(picked, args.out)
    print(f"Wrote {args.out}.png, {args.out}.pdf, {args.out}_python.txt, "
          f"{args.out}_latex.tex, {args.out}.json")


def cmd_cvd_test(args, colors, combos):
    if args.from_combination:
        if args.from_combination not in combos:
            sys.exit(f"No such combination id: {args.from_combination} (valid range 1-348)")
        members = combos[args.from_combination]
    else:
        hexes = [h.strip() for h in args.colors.split(",")]
        members = [{"hex": h, "name": h} for h in hexes]

    if len(members) < 2:
        sys.exit("Need at least 2 colors to test distinguishability")

    print(f"Testing {len(members)} colors: {', '.join(c['hex'] for c in members)}")
    print_cvd_diagnostics([c["hex"] for c in members])

    render_swatch(members, args.out, cvd_preview=True)
    print(f"Wrote {args.out}.png, {args.out}.pdf")


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


def _demo_lines(n):
    rng = np.random.default_rng(42)
    x = np.linspace(0, 10, 200)
    return [(x, np.sin(x + i * 0.7) * np.exp(-0.05 * x) + rng.normal(scale=0.03, size=x.shape))
            for i in range(n)]


def _demo_clusters(n, points_per=40):
    rng = np.random.default_rng(7)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    clusters = []
    for i in range(n):
        cx, cy = 2.2 * np.cos(angles[i]), 2.2 * np.sin(angles[i])
        clusters.append((rng.normal(cx, 0.45, points_per), rng.normal(cy, 0.45, points_per)))
    return clusters


def _demo_field():
    x = np.linspace(-3, 3, 150)
    y = np.linspace(-3, 3, 150)
    X, Y = np.meshgrid(x, y)
    Z = (np.exp(-(X ** 2 + Y ** 2) / 4)
         + 0.5 * np.exp(-((X - 1.6) ** 2 + (Y - 1.0) ** 2) / 1.5)
         - 0.4 * np.exp(-((X + 1.6) ** 2 + (Y + 1.2) ** 2) / 2.0))
    return X, Y, Z


def render_demo(hexes, names, label, out_prefix):
    fig, (ax_line, ax_scatter, ax_heat) = plt.subplots(1, 3, figsize=(13, 4.2))

    for (x, y), h in zip(_demo_lines(len(hexes)), hexes):
        ax_line.plot(x, y, color=h, lw=1.8)
    ax_line.set_title("line plot", fontsize=10)
    ax_line.set_xlabel("x", fontsize=9)
    ax_line.set_ylabel("y", fontsize=9)
    ax_line.tick_params(labelsize=8)

    for (xs, ys), h, name in zip(_demo_clusters(len(hexes)), hexes, names):
        ax_scatter.scatter(xs, ys, color=h, s=18, alpha=0.85, edgecolor="white",
                            linewidth=0.3, label=name)
    ax_scatter.set_title("scatter (categories)", fontsize=10)
    ax_scatter.set_xticks([])
    ax_scatter.set_yticks([])
    ax_scatter.legend(fontsize=6, frameon=False, loc="upper right",
                       handletextpad=0.3, borderaxespad=0.2)

    X, Y, Z = _demo_field()
    smooth_hexes, smooth_rgb = interpolate_lab(hexes, 256)
    cmap = ListedColormap(smooth_rgb)
    im = ax_heat.imshow(Z, cmap=cmap, origin="lower", extent=[-3, 3, -3, 3], aspect="auto")
    fig.colorbar(im, ax=ax_heat, fraction=0.046, pad=0.04)
    ax_heat.set_title("heatmap (palette as colormap)", fontsize=10)
    ax_heat.set_xticks([])
    ax_heat.set_yticks([])

    fig.suptitle(label, fontsize=10, color="0.3")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(f"{out_prefix}.png", dpi=300)
    fig.savefig(f"{out_prefix}.pdf")
    plt.close(fig)


def cmd_demo(args, colors, combos):
    if args.from_combination:
        if args.from_combination not in combos:
            sys.exit(f"No such combination id: {args.from_combination} (valid range 1-348)")
        members = combos[args.from_combination]
        hexes = [c["hex"] for c in members]
        names = [c["name"] for c in members]
        label = f"Sanzo Wada combination #{args.from_combination}: {', '.join(names)}"
    else:
        hexes = [h.strip() for h in args.colors.split(",")]
        names = hexes
        label = f"Custom palette: {', '.join(hexes)}"

    if len(hexes) < 2:
        sys.exit("Need at least 2 colors for the demo plots")

    render_demo(hexes, names, label, args.out)
    print(f"Wrote {args.out}.png, {args.out}.pdf")


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

    p_seq = sub.add_parser(
        "sequential", help="one hue, light->dark: a magnitude colormap (e.g. for a heatmap)")
    s_src = p_seq.add_mutually_exclusive_group(required=True)
    s_src.add_argument("--color", help="base hue, e.g. '#1c4286'")
    s_src.add_argument("--from-combination", type=int,
                        help="use an existing combination's first color as the base hue")
    p_seq.add_argument("--steps", type=int, default=8,
                        help="number of discrete steps to also output (default 8)")
    p_seq.add_argument("--out", default="sanzo_wada_sequential", help="output basename")

    p_div = sub.add_parser(
        "diverging", help="two hues + neutral midpoint: a polarity colormap (e.g. +/- data)")
    d_src = p_div.add_mutually_exclusive_group(required=True)
    d_src.add_argument("--colors", help="exactly 2 comma-separated hex poles, e.g. '#aa0000,#0055aa'")
    d_src.add_argument("--from-combination", type=int,
                        help="use an existing combination's first and last colors as the two poles")
    p_div.add_argument("--midpoint", help="override the neutral midpoint hex (default a light gray)")
    p_div.add_argument("--steps", type=int, default=9,
                        help="number of discrete steps to also output (default 9)")
    p_div.add_argument("--out", default="sanzo_wada_diverging", help="output basename")

    p_qual = sub.add_parser(
        "qualitative", help="pick N maximally-separated colors for a categorical/identity palette")
    p_qual.add_argument("--n", type=int, required=True,
                         help="how many distinct colors to pick")
    p_qual.add_argument("--out", default="sanzo_wada_qualitative", help="output basename")
    p_qual.add_argument("--cvd-preview", action="store_true",
                         help="add protanopia/deuteranopia/tritanopia preview rows")

    p_cvd = sub.add_parser(
        "cvd-test", help="test any palette's distinguishability under color vision deficiency")
    c_src = p_cvd.add_mutually_exclusive_group(required=True)
    c_src.add_argument("--colors", help="comma-separated hex colors to test")
    c_src.add_argument("--from-combination", type=int, help="test an existing combination's colors")
    p_cvd.add_argument("--out", default="sanzo_wada_cvd_test", help="output basename")

    p_demo = sub.add_parser(
        "demo", help="render example line/scatter/heatmap plots using a palette")
    de_src = p_demo.add_mutually_exclusive_group(required=True)
    de_src.add_argument("--from-combination", type=int, help="combination id to use as the palette")
    de_src.add_argument("--colors", help="comma-separated hex colors to use as the palette")
    p_demo.add_argument("--out", default="sanzo_wada_demo", help="output basename")

    args = p.parse_args()
    colors = fetch_colors(refresh=args.refresh)
    combos = build_combinations(colors)

    {"list": cmd_list, "show": cmd_show, "export": cmd_export,
     "random": cmd_random, "gallery": cmd_gallery, "sequential": cmd_sequential,
     "diverging": cmd_diverging, "qualitative": cmd_qualitative,
     "cvd-test": cmd_cvd_test, "demo": cmd_demo}[args.command](args, colors, combos)


if __name__ == "__main__":
    main()
