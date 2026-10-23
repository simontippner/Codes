#!/usr/bin/env python3
"""Overlay several WHAM PMF curves (e.g. an old vs a re-run analysis), no error bands.

Each input is LABEL=FILE, where FILE is a whitespace table with z in column 1 and
the PMF in column 2 (Grossfield pmf.dat / pmf_shifted.dat, '#' comments ignored).
Each curve is shifted to its own minimum unless --no-shift is given.

    python3 compare_pmf.py "old=wham_analysis/pmf_shifted.dat" \
        "new=wham_analysis_2026-10-23/pmf_shifted.dat" --out pmf_old_vs_new
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLORS = ["#9a9a9a", "#21b0fe", "#fe218b", "#2c8f8f", "#e6a100"]
STYLES = ["--", "-", "-", "-.", ":"]

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("curves", nargs="+", help="LABEL=FILE pairs")
p.add_argument("--out", default="pmf_compare", help="output prefix (.png and .pdf)")
p.add_argument("--title", default="WHAM PMF comparison")
p.add_argument("--xlabel", default="z (Å)")
p.add_argument("--no-shift", action="store_true", help="do not shift each curve to its minimum")
a = p.parse_args()

fig, ax = plt.subplots(figsize=(10, 5))
for i, item in enumerate(a.curves):
    label, fn = item.split("=", 1)
    d = np.loadtxt(fn, comments="#", usecols=(0, 1))
    d = d[np.isfinite(d[:, 1])]
    f = d[:, 1] if a.no_shift else d[:, 1] - d[:, 1].min()
    ax.plot(d[:, 0], f, color=COLORS[i % len(COLORS)], ls=STYLES[i % len(STYLES)], lw=2.2, label=label)
ax.set_xlabel(a.xlabel, fontsize=12)
ax.set_ylabel("PMF (kcal/mol)", fontsize=12)
ax.set_title(a.title, fontsize=13, fontweight="bold")
ax.legend(fontsize=11)
ax.grid(True, linestyle="--", alpha=0.4)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{a.out}.{ext}", dpi=200)
print(f"wrote {a.out}.png / .pdf")
