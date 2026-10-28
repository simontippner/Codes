#!/usr/bin/env python3
"""
SI RDF figure: one brush protonation state, both solvents (50:50 H2O:MeOH and
100% H2O), both solutes -- styled like ../combined_figure/combined_v2_inbrush.py.

    python3 rdf_SI_solvents.py half     -> rdf_SI_solvents_half.{png,svg,pdf}
    python3 rdf_SI_solvents.py full     -> rdf_SI_solvents_full.{png,svg,pdf}

BRUSH-ONLY: every RDF uses only the umbrella windows with the solute inside the
brush (z < 58 A: windows 0.0-57.5 A = frames 1-23200 of all.nc); inputs in each
system's .../traj/analysis/inbrush/.  Water is counted by its H atoms
(:WAT@H1,H2) in both solvents, so mixture and pure-water panels are comparable.
Each interaction block has a mixture row and a pure-water row below it, on the
same y-scale (no overlays):
   A-D  brush - solvent : H2O H & MeOH H  vs  unprot. N | prot. NH+
   E-H  solute - brush  : sites  vs  unprot. N | prot. NH+   (broken g-axis)
   I-L  solute - solvent: sites  vs  H2O H | MeOH H
         (pure water has no MeOH panel; that slot carries the site legend)
The fully protonated brush has no neutral amine, so rows A-H carry a single
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
SOLV = [("mix", r"50:50 H$_2$O:MeOH"), ("wat", r"100% H$_2$O")]

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
# y-ranges (shared by the mixture and water rows): broken solute-brush axis,
# upper part GHI with ticks authk, lower part 0-GLO
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


# grid rows: A mix | sp | A wat | gap | C mix up/lo | sp | C wat up/lo | gap | E mix | sp | E wat
HR = [2.0, 0.10, 2.0, 0.42, 0.42, 1.7, 0.10, 0.42, 1.7, 0.42, 2.0, 0.10, 2.0]
R_A = {"mix": 0, "wat": 2}
R_C = {"mix": (4, 5), "wat": (7, 8)}
R_E = {"mix": 10, "wat": 12}
fig = plt.figure(figsize=(8.3, 12.6))
gs = fig.add_gridspec(len(HR), 5, width_ratios=[1, 1, 0.30, 1, 1], height_ratios=HR,
                      hspace=0.10, wspace=0.05,
                      left=0.145, right=0.885, top=0.915, bottom=0.05)
DC = [0, 1, 3, 4]


def row_axes(r, sharex=None):
    """one panel per N type and solute (half) or one wide NH+ panel (full)"""
    if PROT == "half":
        return [fig.add_subplot(gs[r, c], sharex=sharex[i] if sharex else None)
                for i, c in enumerate(DC)]
    return [fig.add_subplot(gs[r, 0:2], sharex=sharex[0] if sharex else None),
            fig.add_subplot(gs[r, 3:5], sharex=sharex[1] if sharex else None)]


fig.text(0.515, 0.977, {"half": "partially protonated brush",
                        "full": "fully protonated brush"}[PROT],
         ha="center", fontsize=FS_H, fontweight="bold", color=C_PROT)
fig.text(0.31, 0.945, "PS: BODIPY", ha="center", fontsize=FS_H, fontweight="bold")
fig.text(0.72, 0.945, r"CAT: Mo$_3$S$_{13}^{2-}$", ha="center", fontsize=FS_H,
         fontweight="bold")

n = len(NTAGS)
A, BU, BL, C = {}, {}, {}, {}

# ============ brush - solvent ============
SA = [("watH", r"$\mathbf{H}_2$O", WCOL), ("meoH", r"MeO$\mathbf{H}$", MCOL)]
for solv, _ in SOLV:
    A[solv] = row_axes(R_A[solv])
    for si, S in enumerate(SYS):
        for j, ntag in enumerate(NTAGS):
            ax = A[solv][n * si + j]
            for sp, sl, col in SA:
                if solv == "wat" and sp == "meoH":
                    continue
                r, g = ld(path(S, solv, "solute_solvent_H"), "rdf", f"{sp}-{ntag}")
                ax.plot(r, g, color=col, lw=LW)
            ax.axhline(1, color="0.6", lw=0.9, ls=":")
            ax.set_ylim(0, YR["GBS"]); ax.set_yticks([0, 1, 2])
            xstyle(ax, xlab=False)
            if j == 1 or (PROT == "full" and si == 1):
                ax.tick_params(labelleft=False)
            if solv == "mix":
                ax.text(0.5, 0.92, NLAB[ntag], transform=ax.transAxes,
                        ha="center", va="top", fontsize=FS_L, fontweight="bold")
    A[solv][0].set_ylabel(r"$\mathbf{g(r)}$", fontsize=FS_AX - 1)
    if PROT == "half":
        drop_seam(A[solv])
A["mix"][0].legend(handles=[Line2D([], [], color=c, lw=2.2, label=l) for _, l, c in SA],
                   loc="lower left", bbox_to_anchor=(0.30 if PROT == "half" else 0.55, 0.0),
                   fontsize=FS_L - 1, frameon=True, framealpha=0.9, handlelength=1.3,
                   labelspacing=0.22, borderpad=0.35, handletextpad=0.45)

# ============ solute - brush  (broken g-axis) ============
_ang = np.deg2rad(150.0)
BK = dict(marker=[(-np.cos(_ang), -np.sin(_ang)), (np.cos(_ang), np.sin(_ang))],
          markersize=9, linestyle="none", color="k", mec="k", mew=0.9, clip_on=False)
for solv, _ in SOLV:
    BU[solv] = row_axes(R_C[solv][0])
    BL[solv] = row_axes(R_C[solv][1], sharex=BU[solv])
    for si, S in enumerate(SYS):
        for j, ntag in enumerate(NTAGS):
            u, l = BU[solv][n * si + j], BL[solv][n * si + j]
            for k, (st, sl) in enumerate(S["sites"]):
                r, g = ld(path(S, solv, "solute_brush"), "rdf", f"{ntag}-{st}")
                u.plot(r, g, color=CS[k], lw=LW); l.plot(r, g, color=CS[k], lw=LW)
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
            u.plot([0, 1], [0, 0], transform=u.transAxes, **BK)
            l.plot([0, 1], [1, 1], transform=l.transAxes, **BK)
            if j == 1:
                u.tick_params(labelleft=False); l.tick_params(labelleft=False)
            if solv == "mix":
                u.text(0.62, 0.82, NLAB[ntag], transform=u.transAxes,
                       ha="center", va="top", fontsize=FS_L, fontweight="bold")
    BL[solv][0].set_ylabel(r"$\mathbf{g(r)}$", fontsize=FS_AX - 1)
    BL[solv][0].yaxis.set_label_coords(-0.20, 0.62)

# ============ solute - solvent ============
for solv, _ in SOLV:
    C[solv] = [fig.add_subplot(gs[R_E[solv], c]) for c in DC]
    for si, S in enumerate(SYS):
        for j, tag in enumerate(["wat", "meo"]):
            ax = C[solv][2 * si + j]
            if solv == "wat" and tag == "meo":
                # no methanol in pure water: use the slot for the site legend
                ax.axis("off")
                ax.legend(handles=[Line2D([], [], color=CS[k], lw=2.4,
                                          label=S["sites"][k][1]) for k in range(3)],
                          loc="center left", bbox_to_anchor=(0.04, 0.5),
                          fontsize=FS_L + 0.5, frameon=False, title="site",
                          title_fontsize=FS_L + 0.5, handlelength=1.6,
                          labelspacing=0.45, alignment="left")
                continue
            for k, (st, sl) in enumerate(S["sites"]):
                r, g = ld(path(S, solv, "solute_solvent_H"), "rdf", f"{tag}-{st}")
                ax.plot(r, g, color=CS[k], lw=LW)
            ax.axhline(1, color="0.6", lw=0.9, ls=":")
            ax.set_ylim(0, YR["GMAX"]); ax.set_yticks([0, 1, 2, 3])
            xstyle(ax, xlab=(solv == "wat"))
            if j == 1:
                ax.tick_params(labelleft=False)
            pos = (0.05, "left") if j == 1 else (0.96, "right")
            ax.text(pos[0], 0.92, [r"$\mathbf{H}_2$O", r"MeO$\mathbf{H}$"][j],
                    transform=ax.transAxes, ha=pos[1], va="top", fontsize=FS_P)
    C[solv][0].set_ylabel(r"$\mathbf{g(r)}$", fontsize=FS_AX - 1)
drop_seam(C["mix"])

# ============ common y-label line + block and solvent labels ============
fig.canvas.draw()
_TGT = 0.10
for _ax in (A["mix"][0], A["wat"][0], BL["mix"][0], BL["wat"][0], C["mix"][0], C["wat"][0]):
    _p = _ax.get_position()
    _y = _ax.yaxis.label.get_position()[1]
    _ax.yaxis.set_label_coords((_TGT - _p.x0) / _p.width, _y)

BLOCKS = [(A["mix"][0], A["wat"][0], r"brush N/NH$^+$ - solvent" if PROT == "half"
           else r"brush NH$^+$ - solvent"),
          (BU["mix"][0], BL["wat"][0], r"brush N/NH$^+$ - solute" if PROT == "half"
           else r"brush NH$^+$ - solute"),
          (C["mix"][0], C["wat"][0], "solvent - solute")]
for top, bot, lab in BLOCKS:
    y0, y1 = bot.get_position().y0, top.get_position().y1
    fig.text(0.035, (y0 + y1) / 2, lab, rotation=90, ha="center", va="center",
             fontsize=FS_AX - 1, fontweight="bold")

# solvent label on the right of every row (pair, for the broken axis)
xr = A["mix"][-1].get_position().x1 + 0.025
for solv, top, bot in [("mix", A["mix"][-1], A["mix"][-1]), ("wat", A["wat"][-1], A["wat"][-1]),
                       ("mix", BU["mix"][-1], BL["mix"][-1]), ("wat", BU["wat"][-1], BL["wat"][-1]),
                       ("mix", C["mix"][-1], C["mix"][-1]), ("wat", C["wat"][-2], C["wat"][-2])]:
    y0, y1 = bot.get_position().y0, top.get_position().y1
    fig.text(xr, (y0 + y1) / 2, dict(SOLV)[solv], rotation=-90, ha="center",
             va="center", fontsize=FS_L + 0.5, fontweight="bold", color="0.25")

# ============ panel letters: row-wise, PS then CAT ============
ORDER = [(A["mix"], n), (A["wat"], n), (BL["mix"], n), (BL["wat"], n),
         (C["mix"], 2), (C["wat"], 2)]
for i, (row, step) in enumerate(ORDER):
    for _ax, _lt in ((row[0], "ABCDEFGHIJKL"[2 * i]), (row[step], "ABCDEFGHIJKL"[2 * i + 1])):
        _ax.text(0.028, 0.95, _lt, transform=_ax.transAxes, ha="left", va="top",
                 fontsize=FS_L, fontweight="bold", zorder=6,
                 bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none", alpha=0.7))

out = f"rdf_SI_solvents_{PROT}"
fig.savefig(f"{out}.png", dpi=220)
fig.savefig(f"{out}.svg")
fig.savefig(f"{out}.pdf")
print(f"wrote {out}.png / .svg / .pdf")
