#!/usr/bin/env python3
"""
Combined PS + CAT donor-H solvation figure, v2 -- BRUSH-ONLY version.
Same as combined_v2.py, but every RDF uses only the umbrella windows with the
solute inside the brush (z < 58 A: windows 0.0-57.5 A = frames 1-23200 of all.nc);
inputs in .../traj/analysis/inbrush/.
All the RDFs stacked on top, in the order:
   a  brush - solvent   : H2O H & MeOH H  vs  unprot. N | prot. NH+
   b  solute - brush    : sites  vs  unprot. N | prot. NH+   (broken g-axis)
   c  solute - solvent  : sites  vs  H2O H | MeOH
   d  preferential-solvation map (one per system)
Columns = system.  Sized to fill an A4 sheet width with room for a caption below.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
from matplotlib.lines import Line2D

CS = ["#fe218b", "#fed700", "#21b0fe"]
WCOL, MCOL = "#1f7fb5", "#2ca02c"
LW = 1.8
FS_T, FS_AX, FS_P, FS_L, FS_H = 10, 12.5, 11, 9.5, 15

BOD = "/gpfs/data/stippner/MoS_BODIPY/BODIPY_in_brush/BIGGER/bigger_bottom/RETRY_half/MD/traj/analysis"
MOS = "/gpfs/data/stippner/MoS_BODIPY/MoS_in_brush/RETRY_half/MD/traj/analysis"

SYS = [
    dict(name="PS", A=f"{BOD}/inbrush/solute_solvent_H", BR=f"{BOD}/inbrush/solute_brush",
         sites=[("COO", r"$-$COO$^-$"), ("F", r"BF$_2$"), ("I", "diiodo I")],
         NW=19231, NM=8590, GHI=(40, 66), authk=[50, 60]),
    dict(name="CAT", A=f"{MOS}/inbrush/solute_solvent_H", BR=f"{MOS}/inbrush/solute_brush",
         sites=[("Sap", "apical S"), ("Sbr", "bridging S"), ("Stm", "terminal S")],
         NW=19238, NM=8593, GHI=(17, 32), authk=[20, 30]),
]
GMAX, GLO, GBS = 3.5, 10.5, 2.0   # brush-only peaks are higher than in combined_v2
ANN = {3.: (6, -13), 4.: (7, 5), 5.: (8, -13), 7.: (17, 4), 10.: (9, -13)}


def ld(A, prefix, key, sm=3):
    d = np.loadtxt(f"{A}/{prefix}_{key}.dat", comments=["#", "@"])
    r, y = d[:, 0], d[:, 1]
    if sm and prefix == "rdf":
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


fig = plt.figure(figsize=(8.3, 9.6))
gs = fig.add_gridspec(9, 5, width_ratios=[1, 1, 0.30, 1, 1],
                      height_ratios=[2.2, 0.05, 0.42, 1.85, 0.02, 2.2, 0.60, 3.5, 0.02],
                      hspace=0.10, wspace=0.05,
                      left=0.145, right=0.905, top=0.935, bottom=0.052)
DC = [0, 1, 3, 4]

fig.text(0.322, 0.949, "PS: BODIPY", ha="center", fontsize=FS_H, fontweight="bold")
fig.text(0.74, 0.949, r"CAT: Mo$_3$S$_{13}^{2-}$", ha="center", fontsize=FS_H, fontweight="bold")

# ============ (a) brush - solvent ============
a1 = [fig.add_subplot(gs[0, c]) for c in DC]
SA = [("watH", r"$\mathbf{H}_2$O", WCOL), ("meoH", r"MeO$\mathbf{H}$", MCOL)]
for si, S in enumerate(SYS):
    for j, ntag in enumerate(["n3", "n4"]):
        ax = a1[2 * si + j]
        for sp, sl, col in SA:
            r, g = ld(S["A"], "rdf", f"{sp}-{ntag}")
            ax.plot(r, g, color=col, lw=LW)
        ax.axhline(1, color="0.6", lw=0.9, ls=":")
        ax.set_ylim(0, GBS); ax.set_yticks([0, 1, 2])
        xstyle(ax, xlab=False)          # ticks only, no "r (Å)" (block b/c carry it)
        if j == 1:
            ax.tick_params(labelleft=False)
        ax.text(0.5, 0.92, [r"unprot. N", r"prot. NH$^+$"][j], transform=ax.transAxes,
                ha="center", va="top", fontsize=FS_L, fontweight="bold")
a1[0].set_ylabel(r"$\mathbf{g(r)}$" + "\n" + r"brush N/NH$^+$ - solvent",
                 fontsize=FS_AX - 2)
a1[0].legend(handles=[Line2D([], [], color=c, lw=2.2, label=l) for _, l, c in SA],
             loc="lower left", bbox_to_anchor=(0.30, 0.0), fontsize=FS_L - 1,
             frameon=True, framealpha=0.9, handlelength=1.3, labelspacing=0.22,
             borderpad=0.35, handletextpad=0.45)
drop_seam(a1)

# ============ (b) solute - brush  (broken g-axis) ============
bu = [fig.add_subplot(gs[2, c]) for c in DC]
bl = [fig.add_subplot(gs[3, c], sharex=bu[0]) for c in DC]
for si, S in enumerate(SYS):
    for j, ntag in enumerate(["n3", "n4"]):
        u, l = bu[2 * si + j], bl[2 * si + j]
        for k, (st, sl) in enumerate(S["sites"]):
            r, g = ld(S["BR"], "rdf", f"{ntag}-{st}")
            u.plot(r, g, color=CS[k], lw=LW); l.plot(r, g, color=CS[k], lw=LW)
        l.axhline(1, color="0.6", lw=0.9, ls=":")
        u.set_ylim(*S["GHI"]); u.set_yticks(S["authk"])
        l.set_ylim(0, GLO); l.set_yticks([0, 5, 10])
        for a in (u, l):
            a.set_xlim(0, 15); a.set_xticks([0, 5, 10, 15])
            a.grid(True, color="0.9", lw=0.5)
            a.tick_params(labelsize=FS_T, labelbottom=False)
        u.tick_params(bottom=False)
        for sp in (u.spines["bottom"], l.spines["top"]):
            sp.set_visible(True); sp.set_linestyle((0, (3, 2))); sp.set_linewidth(0.9)
        # fixed-size slanted cut marks -> identical & parallel on both panels
        _ang = np.deg2rad(150.0)
        bk = dict(marker=[(-np.cos(_ang), -np.sin(_ang)), (np.cos(_ang), np.sin(_ang))],
                  markersize=9, linestyle="none",
                  color="k", mec="k", mew=0.9, clip_on=False)
        u.plot([0, 1], [0, 0], transform=u.transAxes, **bk)
        l.plot([0, 1], [1, 1], transform=l.transAxes, **bk)
        if j == 1:
            u.tick_params(labelleft=False); l.tick_params(labelleft=False)
        l.tick_params(labelbottom=False)         # no ticks; r(Å) label on row c
        u.text(0.5, 0.82, [r"unprot. N", r"prot. NH$^+$"][j], transform=u.transAxes,
               ha="center", va="top", fontsize=FS_L, fontweight="bold")
bl[0].set_ylabel(r"$\mathbf{g(r)}$" + "\n" + r"brush N/NH$^+$ - solute",
                 fontsize=FS_AX - 2)
bl[0].yaxis.set_label_coords(-0.20, 0.72)
drop_seam(bl)

# ============ (c) solute - solvent RDF ============
c1 = [fig.add_subplot(gs[5, c]) for c in DC]
for si, S in enumerate(SYS):
    for j, tag in enumerate(["wat", "meo"]):
        ax = c1[2 * si + j]
        for k, (st, sl) in enumerate(S["sites"]):
            r, g = ld(S["A"], "rdf", f"{tag}-{st}")
            ax.plot(r, g, color=CS[k], lw=LW)
        ax.axhline(1, color="0.6", lw=0.9, ls=":")
        ax.set_ylim(0, GMAX); ax.set_yticks([0, 1, 2, 3])
        xstyle(ax, xlab=True)
        if j == 1:
            ax.tick_params(labelleft=False)
        pos = (0.05, "left") if j == 1 else (0.96, "right")
        ax.text(pos[0], 0.92, [r"$\mathbf{H}_2$O", r"MeO$\mathbf{H}$"][j],
                transform=ax.transAxes, ha=pos[1], va="top", fontsize=FS_P,
                fontweight="normal")
    # site legend parked in row b: PS -> prot. NH+ panel, CAT -> unprot. N panel
    (bl[1] if si == 0 else bl[2]).legend(
        handles=[Line2D([], [], color=CS[k], lw=2.4, label=S["sites"][k][1])
                 for k in range(3)],
        loc="upper right", fontsize=FS_L + (1 if si == 0 else 0.2),
        frameon=True, framealpha=0.95,
        handlelength=1.4, labelspacing=0.3, borderpad=0.4, handletextpad=0.5)
c1[0].set_ylabel(r"$\mathbf{g(r)}$" + "\n" + r"solvent - solute", fontsize=FS_AX - 2)
drop_seam(c1)

# ============ (d) preferential-solvation map ============
clo, chi = 2e-5, 3e-2
d1 = [fig.add_subplot(gs[7, 0:2]), fig.add_subplot(gs[7, 3:5])]
for si, S in enumerate(SYS):
    ax = d1[si]
    ax.plot([1e-5, 1e-1], [1e-5, 1e-1], color="0.55", lw=1.1, ls=":")
    ax.fill_between([1e-5, 1e-1], [1e-5, 1e-1], 1e-1, color="#21b0fe", alpha=0.05)
    ax.fill_between([1e-5, 1e-1], 1e-5, [1e-5, 1e-1], color="#2ca02c", alpha=0.05)
    NT = (2 * S["NW"], S["NM"])
    last = S["sites"][-1][0]
    for k, (st, sl) in enumerate(S["sites"]):
        rr, cw = ld(S["A"], "cn", f"wat-{st}", sm=0); cw = cw / NT[0]
        cm = ld(S["A"], "cn", f"meo-{st}", sm=0)[1] / NT[1]
        ax.plot(cm, cw, color=CS[k], lw=2.0)
        for rm in (3., 4., 5., 7., 10.):
            q = np.searchsorted(rr, rm)
            ax.plot(cm[q], cw[q], "o", color=CS[k], ms=5, mec="white", mew=0.8)
            if st == last:
                off = ANN[rm]
                if rm == 4.:
                    off = (17, -3) if si == 0 else (-28, 2)
                ax.annotate(f"{rm:.0f} Å", (cm[q], cw[q]), textcoords="offset points",
                            xytext=off, fontsize=FS_L - 1, color="0.35")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(clo, chi); ax.set_ylim(clo, chi)
    ax.tick_params(labelsize=FS_T)
    ax.grid(True, which="major", color="0.9", lw=0.5)
    ax.set_xlabel(r"MeO$\mathbf{H}$ fraction   ($N(r)/N_\mathrm{tot}$)", fontsize=FS_AX)
    ax.text(0.13, 0.95, "water-enriched", transform=ax.transAxes, fontsize=FS_L,
            color="#1f7fb5", va="top", fontweight="bold")
    ax.text(0.97, 0.05, "methanol-enriched", transform=ax.transAxes, fontsize=FS_L,
            color="#2ca02c", ha="right", fontweight="bold")
    if si == 1:
        ax.tick_params(labelleft=False)
d1[0].set_ylabel(r"$\mathbf{H}_2$O fraction   ($N(r)/N_\mathrm{tot}$)",
                 fontsize=FS_AX, labelpad=8)

# ============ align every row's y-axis label to a common vertical line ========
fig.canvas.draw()
_TGT = 0.086   # common y-label x, figure coords
for _ax in (a1[0], bl[0], c1[0], d1[0]):
    _p = _ax.get_position()
    _y = _ax.yaxis.label.get_position()[1]
    _ax.yaxis.set_label_coords((_TGT - _p.x0) / _p.width, _y)

# ============ panel letters ============
fig.canvas.draw()
# every row carries its interaction pair in the y-axis label; letter each panel
#   row-wise, left (PS) then right (CAT): A B / C D / E F / G H
LETTERS = [(a1[0], "A"), (a1[2], "B"), (bl[0], "C"), (bl[2], "D"),
           (c1[0], "E"), (c1[2], "F"), (d1[0], "G"), (d1[1], "H")]
for _ax, _lt in LETTERS:
    _ax.text(0.028, 0.95, _lt, transform=_ax.transAxes, ha="left", va="top",
             fontsize=FS_L, fontweight="bold", zorder=6,
             bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none", alpha=0.7))

fig.savefig("combined_v2_inbrush.png", dpi=220)
fig.savefig("combined_v2_inbrush.svg")
fig.savefig("combined_v2_inbrush.pdf")
print("wrote combined_v2_inbrush.png / .svg / .pdf")
