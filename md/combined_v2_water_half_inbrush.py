#!/usr/bin/env python3
"""
BRUSH-ONLY version of combined_v2_water_half.py: every RDF uses only the umbrella
windows with the solute inside the brush (z < 58 A: windows 0.0-57.5 A = frames
1-23200 of all.nc); inputs in .../traj/analysis/inbrush/.
Combined BODIPY + Mo3S13 (CAT) donor-H solvation figure, water-only solvent.
Two-column version of combined_v2_water.py (BODIPY-only), now including MoS_in_brush
RETRY_half_water alongside BODIPY_in_brush RETRY_half_water (same brush protonation,
50%, so the two systems are directly comparable). No methanol in either system, so
the preferential-solvation map (panel d of the dry combined_v2.py) is dropped, same
reasoning as the single-system water figure.
Row order:
   a  brush - water   : unprot. N | prot. NH+     (per system)
   b  solute - brush  : sites vs unprot. N | prot. NH+   (broken g-axis, per system)
   c  solute - water  : sites                      (per system)
Y-axis ranges are set from each system's own (brush-only) peak values;
lower-g and solute-water ranges are shared between the two solutes so panels
without their own tick labels read off the left-hand scale correctly.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

CS = ["#fe218b", "#fed700", "#21b0fe"]
WCOL = "#1f7fb5"
LW = 1.8
FS_T, FS_AX, FS_P, FS_L, FS_H = 10, 12.5, 11, 9.5, 15

BOD = "/gpfs/data/stippner/MoS_BODIPY/BODIPY_in_brush/BIGGER/bigger_bottom/RETRY_half_water/MD/traj/analysis"
MOS = "/gpfs/data/stippner/MoS_BODIPY/MoS_in_brush/RETRY_half_water/MD/traj/analysis"

SYS = [
    dict(name="PS: BODIPY", A=f"{BOD}/inbrush",
         sites=[("COO", r"$-$COO$^-$"), ("F", r"BF$_2$"), ("I", "diiodo I")],
         GHI=(26, 35), authk=[30], GLO=11.0, GMAX=2.0, GBS=1.5),
    dict(name=r"CAT: Mo$_3$S$_{13}^{2-}$", A=f"{MOS}/inbrush",
         sites=[("Sap", "apical S"), ("Sbr", "bridging S"), ("Stm", "terminal S")],
         GHI=(11, 15), authk=[13], GLO=11.0, GMAX=2.0, GBS=1.5),
]
GBS_TICKS = [0, 0.5, 1.0, 1.5]


def ld(A, subdir, key, sm=3):
    d = np.loadtxt(f"{A}/{subdir}/rdf_{key}.dat", comments=["#", "@"])
    r, y = d[:, 0], d[:, 1]
    if sm:
        y = np.convolve(y, np.ones(sm) / sm, mode="same")
    return r, y


def xstyle(ax, xlab=False):
    ax.set_xlim(0, 15); ax.set_xticks([0, 5, 10, 15])
    ax.grid(True, color="0.9", lw=0.5)
    ax.tick_params(labelsize=FS_T)
    if xlab:
        ax.set_xlabel(r"r (Å)", fontsize=FS_AX)
    else:
        ax.tick_params(labelbottom=False)


def drop_seam(axl):
    for i, a in enumerate(axl):
        lb = a.get_xticklabels()
        if len(lb) >= 4:
            (lb[3] if i % 2 == 0 else lb[0]).set_visible(False)


fig = plt.figure(figsize=(8.3, 9.9))
gs = fig.add_gridspec(6, 5, width_ratios=[1, 1, 0.30, 1, 1],
                      height_ratios=[2.2, 0.42, 1.85, 0.02, 2.2, 0.02],
                      hspace=0.12, wspace=0.08,
                      left=0.145, right=0.905, top=0.930, bottom=0.052)
DC = [0, 1, 3, 4]

fig.text(0.322, 0.975, SYS[0]["name"], ha="center", fontsize=FS_H, fontweight="bold")
fig.text(0.74, 0.975, SYS[1]["name"], ha="center", fontsize=FS_H, fontweight="bold")
fig.text(0.525, 0.948, "water-only solvent (no methanol), 50% protonated brush",
         ha="center", fontsize=FS_AX - 1, color="0.3")

# ============ (a) brush - water ============
a1 = [fig.add_subplot(gs[0, c]) for c in DC]
for si, S in enumerate(SYS):
    for j, ntag in enumerate(["n3", "n4"]):
        ax = a1[2 * si + j]
        r, g = ld(S["A"], "solvent_brush", f"wat-{ntag}")
        ax.plot(r, g, color=WCOL, lw=LW)
        ax.axhline(1, color="0.6", lw=0.9, ls=":")
        ax.set_ylim(0, S["GBS"]); ax.set_yticks(GBS_TICKS)
        xstyle(ax, xlab=False)
        if j == 1:
            ax.tick_params(labelleft=False)
        ax.text(0.5, 0.92, [r"unprot. N", r"prot. NH$^+$"][j], transform=ax.transAxes,
                ha="center", va="top", fontsize=FS_L, fontweight="bold")
a1[0].set_ylabel(r"$\mathbf{g(r)}$" + "\n" + r"brush N/NH$^+$ - water", fontsize=FS_AX - 2)
a1[0].text(0.96, 0.08, r"$\mathbf{H}_2$O", transform=a1[0].transAxes, ha="right",
           va="bottom", fontsize=FS_L, color=WCOL, fontweight="bold")
drop_seam(a1)

# ============ (b) solute - brush  (broken g-axis) ============
bu = [fig.add_subplot(gs[1, c]) for c in DC]
bl = [fig.add_subplot(gs[2, c], sharex=bu[0]) for c in DC]
for si, S in enumerate(SYS):
    for j, ntag in enumerate(["n3", "n4"]):
        u, l = bu[2 * si + j], bl[2 * si + j]
        for k, (st, sl) in enumerate(S["sites"]):
            r, g = ld(S["A"], "solute_brush", f"{ntag}-{st}")
            u.plot(r, g, color=CS[k], lw=LW); l.plot(r, g, color=CS[k], lw=LW)
        l.axhline(1, color="0.6", lw=0.9, ls=":")
        u.set_ylim(*S["GHI"]); u.set_yticks(S["authk"])
        l.set_ylim(0, S["GLO"]); l.set_yticks([0, 5, 10] if S["GLO"] >= 10 else [0, 3, 6] if S["GLO"] >= 6 else [0, 3])
        for a in (u, l):
            a.set_xlim(0, 15); a.set_xticks([0, 5, 10, 15])
            a.grid(True, color="0.9", lw=0.5)
            a.tick_params(labelsize=FS_T, labelbottom=False)
        u.tick_params(bottom=False)
        for sp in (u.spines["bottom"], l.spines["top"]):
            sp.set_visible(True); sp.set_linestyle((0, (3, 2))); sp.set_linewidth(0.9)
        _ang = np.deg2rad(150.0)
        bk = dict(marker=[(-np.cos(_ang), -np.sin(_ang)), (np.cos(_ang), np.sin(_ang))],
                  markersize=9, linestyle="none", color="k", mec="k", mew=0.9, clip_on=False)
        u.plot([0, 1], [0, 0], transform=u.transAxes, **bk)
        l.plot([0, 1], [1, 1], transform=l.transAxes, **bk)
        if j == 1:
            u.tick_params(labelleft=False); l.tick_params(labelleft=False)
        l.tick_params(labelbottom=False)
        u.text(0.5, 0.82, [r"unprot. N", r"prot. NH$^+$"][j], transform=u.transAxes,
               ha="center", va="top", fontsize=FS_L, fontweight="bold")
    (bl[1] if si == 0 else bl[2]).legend(
        handles=[Line2D([], [], color=CS[k], lw=2.4, label=S["sites"][k][1])
                 for k in range(3)],
        loc="upper right", fontsize=FS_L, frameon=True, framealpha=0.95,
        handlelength=1.4, labelspacing=0.3, borderpad=0.4, handletextpad=0.5)
bl[0].set_ylabel(r"$\mathbf{g(r)}$" + "\n" + r"brush N/NH$^+$ - solute", fontsize=FS_AX - 2)
bl[0].yaxis.set_label_coords(-0.20, 0.72)
drop_seam(bl)

# ============ (c) solute - water RDF (one wide axes per system) ============
c1 = [fig.add_subplot(gs[4, 0:2]), fig.add_subplot(gs[4, 3:5])]
for si, S in enumerate(SYS):
    ax = c1[si]
    for k, (st, sl) in enumerate(S["sites"]):
        r, g = ld(S["A"], "solute_solvent", f"wat-{st}")
        ax.plot(r, g, color=CS[k], lw=LW, label=sl)
    r, g = ld(S["A"], "solute_solvent", "wat-COM")
    ax.plot(r, g, color="0.3", lw=LW, ls="--", label="whole solute")
    ax.axhline(1, color="0.6", lw=0.9, ls=":")
    ax.set_ylim(0, S["GMAX"])
    xstyle(ax, xlab=True)
    if si == 1:
        ax.tick_params(labelleft=False)
    ax.legend(loc="upper right", fontsize=FS_L, frameon=True, framealpha=0.95,
              handlelength=1.6, labelspacing=0.3, borderpad=0.4, handletextpad=0.5)
c1[0].set_ylabel(r"$\mathbf{g(r)}$" + "\n" + r"water - solute", fontsize=FS_AX - 2)

# ============ align y-axis labels ============
fig.canvas.draw()
_TGT = 0.086
for _ax in (a1[0], bl[0], c1[0]):
    _p = _ax.get_position()
    _y = _ax.yaxis.label.get_position()[1]
    _ax.yaxis.set_label_coords((_TGT - _p.x0) / _p.width, _y)

LETTERS = [(a1[0], "A"), (bl[0], "B"), (c1[0], "C"),
           (a1[2], "D"), (bl[2], "E"), (c1[1], "F")]
for _ax, _lt in LETTERS:
    _ax.text(0.028, 0.95, _lt, transform=_ax.transAxes, ha="left", va="top",
             fontsize=FS_L, fontweight="bold", zorder=6,
             bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none", alpha=0.7))

fig.savefig("combined_v2_water_half_inbrush.png", dpi=220)
fig.savefig("combined_v2_water_half_inbrush.svg")
fig.savefig("combined_v2_water_half_inbrush.pdf")
print("wrote combined_v2_water_half_inbrush.png / .svg / .pdf")
