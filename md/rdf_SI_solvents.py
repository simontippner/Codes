#!/usr/bin/env python3
"""
SI RDF figure: one brush protonation state, both solvents (50:50 H2O:MeOH vs
100% H2O), both solutes -- styled like ../combined_figure/combined_v2_inbrush.py.

    python3 rdf_SI_solvents.py half     -> rdf_SI_solvents_half.{png,svg,pdf}
    python3 rdf_SI_solvents.py full     -> rdf_SI_solvents_full.{png,svg,pdf}

BRUSH-ONLY: every RDF uses only the umbrella windows with the solute inside the
brush (z < 58 A: windows 0.0-57.5 A = frames 1-23200 of all.nc); inputs in each
system's .../traj/analysis/inbrush/.  Water is counted by its H atoms
(:WAT@H1,H2) in both solvents, so mixture and pure-water curves are comparable.
Solid = 50:50 H2O:MeOH, dashed = 100% H2O.
Rows:
   A/B  brush - solvent : H2O H & MeOH H  vs  unprot. N | prot. NH+
   C/D  solute - brush  : sites  vs  unprot. N | prot. NH+   (broken g-axis)
   E/F  solute - solvent: sites  vs  H2O H | MeOH H
The fully protonated brush has no neutral amine, so rows A-D carry a single
NH+ panel per solute there.  The preferential-solvation map of combined_v2_inbrush
is left out (it needs two solvents, so it has no pure-water counterpart).
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

PROT = sys.argv[1] if len(sys.argv) > 1 else "half"
assert PROT in ("half", "full")

CS = ["#fe218b", "#fed700", "#21b0fe"]
WCOL, MCOL = "#1f7fb5", "#2ca02c"
C_PROT = {"half": "#fe218b", "full": "#21b0fe"}[PROT]   # as in the PMF figure
LW = 1.8
FS_T, FS_AX, FS_P, FS_L, FS_H = 10, 12.5, 11, 9.5, 15
LS_MIX, LS_WAT = "-", (0, (3.2, 1.6))

B = "/gpfs/data/stippner/MoS_BODIPY"
BODB = f"{B}/BODIPY_in_brush/BIGGER/bigger_bottom"
MOSB = f"{B}/MoS_in_brush"
RUN = {
    "half": dict(bod_mix=f"{BODB}/RETRY_half/MD/traj/analysis",
                 bod_wat=f"{BODB}/RETRY_half_water/MD/traj/analysis",
                 mos_mix=f"{MOSB}/RETRY_half/MD/traj/analysis",
                 mos_wat=f"{MOSB}/RETRY_half_water/MD/traj/analysis"),
    "full": dict(bod_mix=f"{BODB}/RETRY_full/MD/new/mixture_run/traj/analysis",
                 bod_wat=f"{BODB}/RETRY_full_water/MD/traj/analysis",
                 mos_mix=f"{MOSB}/RETRY_full/new_final/tleap_run/traj/analysis",
                 mos_wat=f"{MOSB}/RETRY_full_water/MD/traj/analysis"),
}[PROT]
# y-ranges: broken solute-brush axis (upper part GHI, ticks authk), lower part 0-GLO
YR = {
    "half": dict(GLO=10.5, LOTK=[0, 5, 10], GMAX=3.5, GBS=2.0,
                 PS=dict(GHI=(28, 66), authk=[40, 60]),
                 CAT=dict(GHI=(11, 32), authk=[20, 30])),
    "full": dict(GLO=7.0, LOTK=[0, 3, 6], GMAX=3.5, GBS=2.0,
                 PS=dict(GHI=(18, 44), authk=[20, 40]),
                 CAT=dict(GHI=(7.0, 12), authk=[8, 10])),
}[PROT]

SYS = [
    dict(name="PS", key="bod", sites=[("COO", r"$-$COO$^-$"), ("F", r"BF$_2$"), ("I", "diiodo I")]),
    dict(name="CAT", key="mos", sites=[("Sap", "apical S"), ("Sbr", "bridging S"), ("Stm", "terminal S")]),
]
NTAGS = ["n3", "n4"] if PROT == "half" else ["n4"]
NLAB = {"n3": "unprot. N", "n4": r"prot. NH$^+$"}


def path(S, solv, sub):
    return f"{RUN[S['key'] + '_' + solv]}/inbrush/{sub}"


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


fig = plt.figure(figsize=(8.3, 7.4))
gs = fig.add_gridspec(6, 5, width_ratios=[1, 1, 0.30, 1, 1],
                      height_ratios=[2.2, 0.05, 0.42, 1.85, 0.02, 2.2],
                      hspace=0.10, wspace=0.05,
                      left=0.145, right=0.905, top=0.845, bottom=0.075)
DC = [0, 1, 3, 4]


def row_axes(r, sharex=None):
    """one panel per N type and solute (half) or one wide NH+ panel (full)"""
    if PROT == "half":
        return [fig.add_subplot(gs[r, c], sharex=sharex[i] if sharex else None)
                for i, c in enumerate(DC)]
    return [fig.add_subplot(gs[r, 0:2], sharex=sharex[0] if sharex else None),
            fig.add_subplot(gs[r, 3:5], sharex=sharex[1] if sharex else None)]


fig.text(0.525, 0.965, {"half": "partially protonated brush",
                        "full": "fully protonated brush"}[PROT],
         ha="center", fontsize=FS_H, fontweight="bold", color=C_PROT)
fig.text(0.322, 0.900, "PS: BODIPY", ha="center", fontsize=FS_H, fontweight="bold")
fig.text(0.74, 0.900, r"CAT: Mo$_3$S$_{13}^{2-}$", ha="center", fontsize=FS_H,
         fontweight="bold")
fig.legend(handles=[Line2D([], [], color="0.25", lw=2.0, ls=LS_MIX,
                           label=r"50:50 H$_2$O:MeOH"),
                    Line2D([], [], color="0.25", lw=2.0, ls=LS_WAT,
                           label=r"100% H$_2$O")],
           loc="center", bbox_to_anchor=(0.525, 0.935), ncol=2, fontsize=FS_L + 0.5,
           frameon=False, handlelength=2.4, columnspacing=2.0)

# ============ brush - solvent ============
a1 = row_axes(0)
SA = [("watH", r"$\mathbf{H}_2$O", WCOL), ("meoH", r"MeO$\mathbf{H}$", MCOL)]
for si, S in enumerate(SYS):
    for j, ntag in enumerate(NTAGS):
        ax = a1[len(NTAGS) * si + j]
        for solv, ls in (("mix", LS_MIX), ("wat", LS_WAT)):
            for sp, sl, col in SA:
                if solv == "wat" and sp == "meoH":
                    continue
                r, g = ld(path(S, solv, "solute_solvent_H"), "rdf", f"{sp}-{ntag}")
                ax.plot(r, g, color=col, lw=LW, ls=ls)
        ax.axhline(1, color="0.6", lw=0.9, ls=":")
        ax.set_ylim(0, YR["GBS"]); ax.set_yticks([0, 1, 2])
        xstyle(ax, xlab=False)
        if j == 1 or (PROT == "full" and si == 1):
            ax.tick_params(labelleft=False)
        ax.text(0.5, 0.92, NLAB[ntag], transform=ax.transAxes,
                ha="center", va="top", fontsize=FS_L, fontweight="bold")
a1[0].set_ylabel(r"$\mathbf{g(r)}$" + "\n" + ("brush N/NH$^+$ - solvent" if PROT == "half"
                                              else "brush NH$^+$ - solvent"),
                 fontsize=FS_AX - 2)
a1[0].legend(handles=[Line2D([], [], color=c, lw=2.2, label=l) for _, l, c in SA],
             loc="lower left", bbox_to_anchor=(0.30 if PROT == "half" else 0.55, 0.0),
             fontsize=FS_L - 1, frameon=True, framealpha=0.9, handlelength=1.3,
             labelspacing=0.22, borderpad=0.35, handletextpad=0.45)
if PROT == "half":
    drop_seam(a1)

# ============ solute - brush  (broken g-axis) ============
bu = row_axes(2)
bl = row_axes(3, sharex=bu)
for si, S in enumerate(SYS):
    for j, ntag in enumerate(NTAGS):
        u, l = bu[len(NTAGS) * si + j], bl[len(NTAGS) * si + j]
        for solv, ls in (("mix", LS_MIX), ("wat", LS_WAT)):
            for k, (st, sl) in enumerate(S["sites"]):
                r, g = ld(path(S, solv, "solute_brush"), "rdf", f"{ntag}-{st}")
                u.plot(r, g, color=CS[k], lw=LW, ls=ls)
                l.plot(r, g, color=CS[k], lw=LW, ls=ls)
        l.axhline(1, color="0.6", lw=0.9, ls=":")
        u.set_ylim(*YR[S["name"]]["GHI"]); u.set_yticks(YR[S["name"]]["authk"])
        l.set_ylim(0, YR["GLO"]); l.set_yticks(YR["LOTK"])
        for a in (u, l):
            a.set_xlim(0, 15); a.set_xticks([0, 5, 10, 15])
            a.grid(True, color="0.9", lw=0.5)
            a.tick_params(labelsize=FS_T, labelbottom=False)
        u.tick_params(bottom=False)
        for sp in (u.spines["bottom"], l.spines["top"]):
            sp.set_visible(True); sp.set_linestyle((0, (3, 2))); sp.set_linewidth(0.9)
        _ang = np.deg2rad(150.0)
        bk = dict(marker=[(-np.cos(_ang), -np.sin(_ang)), (np.cos(_ang), np.sin(_ang))],
                  markersize=9, linestyle="none",
                  color="k", mec="k", mew=0.9, clip_on=False)
        u.plot([0, 1], [0, 0], transform=u.transAxes, **bk)
        l.plot([0, 1], [1, 1], transform=l.transAxes, **bk)
        if j == 1:   # half: right panel of a solute pair shares the left's scale
            u.tick_params(labelleft=False); l.tick_params(labelleft=False)
        u.text(0.62, 0.82, NLAB[ntag], transform=u.transAxes,
               ha="center", va="top", fontsize=FS_L, fontweight="bold")
bl[0].set_ylabel(r"$\mathbf{g(r)}$" + "\n" + ("brush N/NH$^+$ - solute" if PROT == "half"
                                              else "brush NH$^+$ - solute"),
                 fontsize=FS_AX - 2)
bl[0].yaxis.set_label_coords(-0.20, 0.72)
if PROT == "half":
    drop_seam(bl)

# ============ solute - solvent ============
c1 = [fig.add_subplot(gs[5, c]) for c in DC]
for si, S in enumerate(SYS):
    for j, tag in enumerate(["wat", "meo"]):
        ax = c1[2 * si + j]
        for solv, ls in (("mix", LS_MIX), ("wat", LS_WAT)):
            if solv == "wat" and tag == "meo":
                continue
            for k, (st, sl) in enumerate(S["sites"]):
                r, g = ld(path(S, solv, "solute_solvent_H"), "rdf", f"{tag}-{st}")
                ax.plot(r, g, color=CS[k], lw=LW, ls=ls)
        ax.axhline(1, color="0.6", lw=0.9, ls=":")
        ax.set_ylim(0, YR["GMAX"]); ax.set_yticks([0, 1, 2, 3])
        xstyle(ax, xlab=True)
        if j == 1:
            ax.tick_params(labelleft=False)
        if j == 1:   # site legend in the empty upper right of the MeOH panel
            ax.legend(handles=[Line2D([], [], color=CS[k], lw=2.4,
                                      label=S["sites"][k][1]) for k in range(3)],
                      loc="upper right", bbox_to_anchor=(1.0, 0.82),
                      fontsize=FS_L - 0.5, frameon=True,
                      framealpha=0.95, handlelength=1.3, labelspacing=0.25,
                      borderpad=0.35, handletextpad=0.45)
        pos = (0.05, "left") if j == 1 else (0.96, "right")
        ax.text(pos[0], 0.92, [r"$\mathbf{H}_2$O", r"MeO$\mathbf{H}$"][j],
                transform=ax.transAxes, ha=pos[1], va="top", fontsize=FS_P)
c1[0].set_ylabel(r"$\mathbf{g(r)}$" + "\n" + r"solvent - solute", fontsize=FS_AX - 2)
drop_seam(c1)

# ============ common y-label line ============
fig.canvas.draw()
_TGT = 0.086
for _ax in (a1[0], bl[0], c1[0]):
    _p = _ax.get_position()
    _y = _ax.yaxis.label.get_position()[1]
    _ax.yaxis.set_label_coords((_TGT - _p.x0) / _p.width, _y)

# ============ panel letters: row-wise, PS then CAT ============
n = len(NTAGS)
LETTERS = [(a1[0], "A"), (a1[n], "B"), (bl[0], "C"), (bl[n], "D"),
           (c1[0], "E"), (c1[2], "F")]
for _ax, _lt in LETTERS:
    _ax.text(0.028, 0.95, _lt, transform=_ax.transAxes, ha="left", va="top",
             fontsize=FS_L, fontweight="bold", zorder=6,
             bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none", alpha=0.7))
out = f"rdf_SI_solvents_{PROT}"
fig.savefig(f"{out}.png", dpi=220)
fig.savefig(f"{out}.svg")
fig.savefig(f"{out}.pdf")
print(f"wrote {out}.png / .svg / .pdf")
