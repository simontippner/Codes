#!/usr/bin/env python3
"""Images of the tutorial (docs/TUTORIAL.md) from the output of a brush_builder run.
Usage:  python3 make_tutorial_images.py <output folder> [<image folder>]
        (default: ../examples/pdmaema/output_quick -> images/)
Chains as sticks (heavy atoms; central chains wine, periodic copies grey), N of protonated units amber, grafted atoms
black, water O light blue, counterions green; orthographic views."""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE.parent / "examples" / "pdmaema" / "output_quick"
IMG = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE / "images"
IMG.mkdir(parents=True, exist_ok=True)
WINE, ROSE, AMBER, GREY, INK = "#8C1D40", "#E76D83", "#E9A23B", "#b9b3ad", "#1a1a1a"
WATER, CL, GRID = "#7da2ee", "#2fa55a", "#d9d4cf"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#555555", "axes.linewidth": 0.7, "savefig.dpi": 130,
                     "savefig.bbox": "tight", "font.family": "DejaVu Sans"})


def read(path):
    rows = [l for l in Path(path).read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]
    X = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in rows])
    res = np.array([l[17:20].strip() for l in rows]); name = np.array([l[12:16].strip() for l in rows])
    return X, res, name


def split_chains(X, res, name, sizes):
    """consecutive polymer blocks of the given atom counts -> (heavy xyz, bonds, protonated N xyz, grafted atom)"""
    out, k = [], 0
    for n in sizes:
        x, r, nm = X[k:k + n], res[k:k + n], name[k:k + n]
        heavy = np.array([not s.startswith("H") for s in nm])
        P = x[heavy]
        bonds = np.array(sorted(cKDTree(P).query_pairs(1.75)))
        pn = x[(nm == "N1") & np.array([len(s) == 3 and s[1] == "P" for s in r])]
        out.append((P, bonds, pn, x[1])); k += n
    return out


def sticks2d(ax, chains, dims, color, lw=0.8, zorder=2):
    segs = [c[0][b][:, :, dims] for c in chains for b in [c[1]] if len(b)]
    ax.add_collection(LineCollection(np.vstack(segs), colors=color, linewidths=lw, zorder=zorder))
    pn = np.vstack([c[2] for c in chains if len(c[2])]) if any(len(c[2]) for c in chains) else np.zeros((0, 3))
    return pn


def sticks3d(ax, chains, color, lw=0.7, alpha=1.0):
    segs = np.vstack([c[0][c[1]] for c in chains if len(c[1])])
    ax.add_collection3d(Line3DCollection(segs, colors=color, linewidths=lw, alpha=alpha))


def box3d(ax, x0, y0, L1, L2, H, color, lw=1.0):
    c = np.array([[x0, y0], [x0 + L1, y0], [x0 + L1, y0 + L2], [x0, y0 + L2], [x0, y0]])
    for z in (0, H):
        ax.plot(c[:, 0], c[:, 1], [z] * 5, color=color, lw=lw)
    for x, y in c[:4]:
        ax.plot([x, x], [y, y], [0, H], color=color, lw=lw)


def view3d(ax, lo, hi):
    ax.set_proj_type("ortho"); ax.view_init(elev=22, azim=-60)
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
    ax.set_box_aspect(np.array(hi) - np.array(lo)); ax.set_axis_off()


chain_files = sorted((OUT / "chains").glob("chain_*.pdb"), key=lambda p: int(p.stem.split("_")[1]))
sizes = [len(read(f)[0]) for f in chain_files]
seqs = json.loads((OUT / "chains" / "sequences.json").read_text())
P = json.loads((OUT / "pack" / "pack.json").read_text())
S = json.loads((OUT / "solvate" / "system.json").read_text())
L, R = P["cell"], P["region"]

# 1 chains: aligned, side by side
fig, axs = plt.subplots(1, len(chain_files), figsize=(1.6 * len(chain_files), 3.4))
for k, (ax, f) in enumerate(zip(np.atleast_1d(axs), chain_files)):
    X, res, name = read(f)
    (c,) = split_chains(X, res, name, [len(X)])
    pn = sticks2d(ax, [c], [0, 2], WINE)
    ax.scatter(pn[:, 0], pn[:, 2], s=14, color=AMBER, edgecolor="white", linewidth=0.4, zorder=3)
    ax.scatter(c[3][0], c[3][2], s=22, color=INK, zorder=4)
    ax.set_aspect("equal"); ax.autoscale(); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    n_p = sum(1 for r in seqs[k] if r[1] == "P")
    ax.set_title(f"chain {k + 1}\n{n_p} of {len(seqs[k])} protonated", fontsize=8)
fig.suptitle("step 2  chains: random sequences, aligned along z (black: grafted atom)", fontsize=9, y=1.02)
fig.savefig(IMG / "02_chains.png"); plt.close(fig)

# 3a smallest box (top and side view)
X, res, name = read(OUT / "pack" / "box.pdb")
ch = split_chains(X, res, name, sizes)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.2), gridspec_kw=dict(width_ratios=[1, 1.2]))
L0 = P["smallest_box"]
for ax, dims, title in ((a1, [0, 1], f"top view: {L0:.1f} x {L0:.1f} A"), (a2, [0, 2], "side view")):
    pn = sticks2d(ax, ch, dims, WINE)
    ax.scatter(pn[:, dims[0]], pn[:, dims[1]], s=8, color=AMBER, zorder=3)
    g = np.array([c[3] for c in ch]); ax.scatter(g[:, dims[0]], g[:, dims[1]], s=18, color=INK, zorder=4)
    if dims == [0, 1]:
        ax.add_patch(plt.Rectangle((0, 0), L0, L0, fill=False, ec=ROSE, lw=1.4))
    ax.set_aspect("equal"); ax.autoscale(); ax.set_title(title, fontsize=9); ax.set_xlabel("x (A)")
    ax.set_ylabel("y (A)" if dims == [0, 1] else "z (A)")
fig.suptitle("step 3a  smallest square box holding all chains upright at the tolerance", fontsize=9, y=1.0)
fig.savefig(IMG / "03a_smallest_box.png"); plt.close(fig)

# 3b periodic packing: central chains and copies in the 3 x 3 region
X, res, name = read(OUT / "pack" / "region.pdb")
blocks = split_chains(X, res, name, [sizes[c] for c in P["block_chain"]])
cen, cop = blocks[:P["n_central"]], blocks[P["n_central"]:]
fig = plt.figure(figsize=(8.0, 3.8))
a1 = fig.add_axes([0.0, 0.05, 0.42, 0.85]); a2 = fig.add_axes([0.45, 0.0, 0.55, 0.95], projection="3d")
for i in range(4):
    a1.plot([i * L, i * L], [0, R], color=GRID, lw=0.8); a1.plot([0, R], [i * L, i * L], color=GRID, lw=0.8)
a1.add_patch(plt.Rectangle((L, L), L, L, fill=False, ec=ROSE, lw=1.6, zorder=5))
sticks2d(a1, cop, [0, 1], GREY, lw=0.6, zorder=1)
pn = sticks2d(a1, cen, [0, 1], WINE, lw=0.8)
g = np.array([c[3] for c in cen]); a1.scatter(g[:, 0], g[:, 1], s=16, color=INK, zorder=6)
a1.set_xlim(-2, R + 2); a1.set_ylim(-2, R + 2); a1.set_aspect("equal"); a1.set_xticks([]); a1.set_yticks([])
a1.set_title(f"top view: region {R:.0f} A, central cell {L:.1f} A", fontsize=9)
sticks3d(a2, cop, GREY, 0.5, 0.8); sticks3d(a2, cen, WINE, 0.8)
box3d(a2, L, L, L, L, P["height"], ROSE, 1.2)
for i in range(4):
    a2.plot([i * L, i * L], [0, R], [0, 0], color=GRID, lw=0.6); a2.plot([0, R], [i * L, i * L], [0, 0], color=GRID, lw=0.6)
view3d(a2, (0, 0, 0), (R, R, P["height"]))
fig.suptitle(f"step 3b  periodic packing: {P['n_central']} chains (wine) and {P['n_copies']} copies (grey) "
             f"in the 3 x 3 region", fontsize=9, y=1.0)
fig.savefig(IMG / "03b_periodic_packing.png"); plt.close(fig)

# 4a solvation of the whole region (side view, slab through the central cell)
X, res, name = read(OUT / "solvate" / "solvated.pdb")
n_sol = sum(sizes[c] for c in P["block_chain"])
blocks = split_chains(X[:n_sol], res[:n_sol], name[:n_sol], [sizes[c] for c in P["block_chain"]])
wo = (res == "WAT") & (name == "O"); cl = np.char.startswith(res, "Cl")
slab = (X[:, 1] > L) & (X[:, 1] < 2 * L)
fig, ax = plt.subplots(figsize=(7.2, 3.6))
ax.scatter(X[wo & slab, 0], X[wo & slab, 2], s=1.2, color=WATER, alpha=0.6, lw=0, zorder=1)
sticks2d(ax, [b for b in blocks[P["n_central"]:]], [0, 2], GREY, lw=0.5, zorder=2)
sticks2d(ax, blocks[:P["n_central"]], [0, 2], WINE, lw=0.8, zorder=3)
ax.scatter(X[cl, 0], X[cl, 2], s=10, color=CL, zorder=4)
ax.axvline(L, color=ROSE, lw=1.4); ax.axvline(2 * L, color=ROSE, lw=1.4)
ax.set_aspect("equal"); ax.set_xlabel("x (A)"); ax.set_ylabel("z (A)")
ax.set_title("step 4  the whole region solvated (water O in a slab through the central cell), "
             "counterions (green) in the central cell", fontsize=9)
fig.savefig(IMG / "04a_solvated_region.png"); plt.close(fig)

# 4b final periodic box
X, res, name = read(OUT / "solvate" / "system.pdb")
box = S["box"]
n_pol = sum(sizes)
ch = split_chains(X[:n_pol], res[:n_pol], name[:n_pol], sizes)
wo = (res == "WAT") & (name == "O"); cl = np.char.startswith(res, "Cl")
fig = plt.figure(figsize=(7.6, 4.2))
a1 = fig.add_axes([0.0, 0.05, 0.36, 0.88]); a2 = fig.add_axes([0.40, 0.0, 0.6, 1.0], projection="3d")
a1.scatter(X[wo, 0], X[wo, 2], s=2, color=WATER, alpha=0.6, lw=0, zorder=1)
sticks2d(a1, ch, [0, 2], WINE, lw=0.9, zorder=2)
a1.scatter(X[cl, 0], X[cl, 2], s=12, color=CL, zorder=3)
a1.add_patch(plt.Rectangle((0, 0), box[0], box[2], fill=False, ec=ROSE, lw=1.4))
a1.set_aspect("equal"); a1.set_xlabel("x (A)"); a1.set_ylabel("z (A)"); a1.set_title("side view", fontsize=9)
sticks3d(a2, ch, WINE, 0.9)
a2.scatter(X[wo, 0], X[wo, 1], X[wo, 2], s=1, color=WATER, alpha=0.35, lw=0)
a2.scatter(X[cl, 0], X[cl, 1], X[cl, 2], s=10, color=CL)
box3d(a2, 0, 0, box[0], box[1], box[2], ROSE, 1.2)
view3d(a2, (0, 0, 0), box)
fig.suptitle(f"step 4  periodic box {box[0]:.1f} x {box[1]:.1f} x {box[2]:.1f} A: {S['n_chains']} chains, "
             f"{S['n_water']} water, {S['n_ions']} Cl-", fontsize=9, y=1.0)
fig.savefig(IMG / "04b_final_box.png"); plt.close(fig)
print("images in", IMG, sorted(p.name for p in IMG.glob("*.png")))
