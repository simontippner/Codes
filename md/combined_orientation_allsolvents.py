#!/usr/bin/env python3
"""
Solute orientation in the brush: 50:50 H2O:MeOH mixture vs 100% H2O.

Same analysis as ../combined_figure/combined_orientation_v3.py (plane normal vs z,
folded to [0,90]; 0 = flat / face-on, 90 = edge-on), extended to the water runs.
    columns = BODIPY dipyrrin-core plane | BODIPY phenyl-ring plane | Mo3S13 Mo3 plane
    rows    = half mix | half H2O | full mix | full H2O

Writes
  combined_orientation_allsolvents.{png,svg,pdf}  P(theta|z) heat maps, 4 rows
  orientation_mix_vs_water.{png,svg,pdf}          mean angle(z) and S(z), all 4 conditions
  orientation_mix_vs_water_summary.txt            brush/solvent averages with window-level SE

Errors: frames come in umbrella windows of 200 frames (one window = one z). Each
window is reduced to its mean, and the SE of a z-bin or region mean is the spread
of window means / sqrt(n_windows). Inputs: rc.dat + plane-normal files from the
cpptraj vector runs in each system's .../traj/analysis/orientation[/core_phenyl].
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.lines import Line2D

try:
    from scipy.ndimage import gaussian_filter
    smooth = lambda A: gaussian_filter(A, sigma=(0.8, 0.8), mode="nearest")
except Exception:
    smooth = lambda A: A

BOD = "/gpfs/data/stippner/MoS_BODIPY/BODIPY_in_brush/BIGGER/bigger_bottom"
MOS = "/gpfs/data/stippner/MoS_BODIPY/MoS_in_brush"
ICONS = "/gpfs/data/stippner/MoS_BODIPY/combined_figure"
BOD_DIRS = dict(
    half_mix=f"{BOD}/RETRY_half/MD/traj/analysis/orientation/core_phenyl",
    half_wat=f"{BOD}/RETRY_half_water/MD/traj/analysis/orientation/core_phenyl",
    full_mix=f"{BOD}/RETRY_full/MD/new/mixture_run/traj/analysis/orientation/core_phenyl",
    full_wat=f"{BOD}/RETRY_full_water/MD/traj/analysis/orientation/core_phenyl")
MOS_DIRS = dict(
    half_mix=f"{MOS}/RETRY_half/MD/traj/analysis/orientation",
    half_wat=f"{MOS}/RETRY_half_water/MD/traj/analysis/orientation",
    full_mix=f"{MOS}/RETRY_full/new_final/tleap_run/traj/analysis/orientation",
    full_wat=f"{MOS}/RETRY_full_water/MD/traj/analysis/orientation")

COLS = [
    dict(name="BODIPY core", sub=r"BODIPY dipyrrin normal $\angle$ z", vec="core.dat", dirs=BOD_DIRS,
         icon=f"{ICONS}/bodipy_core_annot.png", scale=1.0),
    dict(name="BODIPY phenyl", sub=r"BODIPY phenyl normal $\angle$ z", vec="phenyl.dat", dirs=BOD_DIRS,
         icon=f"{ICONS}/bodipy_phenyl_annot.png", scale=1.0),
    dict(name="Mo3S13 Mo3 plane", sub=r"Mo$_3$ normal $\angle$ z", vec="norm.dat", dirs=MOS_DIRS,
         icon=f"{ICONS}/mos_plane.png", scale=0.85),
]
C_HALF, C_FULL = "#fe218b", "#21b0fe"        # partially / fully protonated (as in the PMF figure)
ROWS = [("half_mix", "partially prot.\n50:50 H$_2$O:MeOH", C_HALF, "-"),
        ("half_wat", "partially prot.\n100% H$_2$O", C_HALF, "--"),
        ("full_mix", "fully prot.\n50:50 H$_2$O:MeOH", C_FULL, "-"),
        ("full_wat", "fully prot.\n100% H$_2$O", C_FULL, "--")]

FRAMES_PER_WIN = 200
BRUSH_EDGE = 58.0
ZMAX = 110.0
zedges = np.arange(0, ZMAX + 5.0, 5.0)
zc = 0.5 * (zedges[:-1] + zedges[1:])
aedges = np.arange(0, 90.001, 5.0)
ISO = np.degrees(1.0)                        # isotropic mean of the folded angle, 57.3 deg


def unit(v):
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def load(D, vecfile):
    z = np.loadtxt(f"{D}/rc.dat")[:, 3]
    n = unit(np.loadtxt(f"{D}/{vecfile}")[:, 1:4])
    cz = np.clip(np.abs(n[:, 2]), 0, 1)
    return z, np.degrees(np.arccos(cz)), cz ** 2


def window_means(z, th, c2):
    nw = len(z) // FRAMES_PER_WIN
    r = lambda a: a[:nw * FRAMES_PER_WIN].reshape(nw, FRAMES_PER_WIN).mean(axis=1)
    return r(z), r(th), r(c2)


def binned(zw, yw):
    """mean and window-level SE of yw per z bin"""
    ib = np.digitize(zw, zedges) - 1
    m = np.full(len(zc), np.nan); se = np.full(len(zc), np.nan)
    for k in range(len(zc)):
        y = yw[ib == k]
        if len(y):
            m[k] = y.mean()
            se[k] = y.std(ddof=1) / np.sqrt(len(y)) if len(y) > 1 else np.nan
    return m, se


def region(zw, yw, lo, hi):
    y = yw[(zw >= lo) & (zw < hi)]
    return y.mean(), y.std(ddof=1) / np.sqrt(len(y)), len(y)


DATA = {}
for ci, C in enumerate(COLS):
    for rk, *_ in ROWS:
        z, th, c2 = load(C["dirs"][rk], C["vec"])
        H, _, _ = np.histogram2d(z, th, bins=[zedges, aedges])
        H = smooth(H)
        with np.errstate(invalid="ignore", divide="ignore"):
            Hn = np.ma.masked_invalid(H / H.sum(axis=1, keepdims=True))
        zw, thw, c2w = window_means(z, th, c2)
        DATA[(ci, rk)] = dict(Hn=Hn, zw=zw, thw=thw, Sw=(3 * c2w - 1) / 2,
                              th_b=binned(zw, thw), S_b=binned(zw, (3 * c2w - 1) / 2), nframes=len(z))

# ---------------------------------------------------------------- heat maps
vmax = max(np.nanmax(d["Hn"].filled(0)) for d in DATA.values())
fig, axes = plt.subplots(4, 3, figsize=(11.4, 14.2), sharex=True, sharey=True,
                         gridspec_kw=dict(hspace=0.12, wspace=0.10))
FS_AX, FS_T, FS_P, FS_ROW = 13, 11, 12, 11.5
for ci, C in enumerate(COLS):
    axes[0, ci].set_title(C["sub"], fontsize=FS_P, pad=8)
    for ri, (rk, rlabel, col, ls) in enumerate(ROWS):
        ax = axes[ri, ci]; d = DATA[(ci, rk)]
        pc = ax.pcolormesh(aedges, zedges, d["Hn"], cmap="magma", shading="flat", vmin=0, vmax=vmax,
                           edgecolors="face", linewidth=0, rasterized=True)  # no cell seams in pdf/svg
        ax.plot(d["th_b"][0], zc, "-", color="cyan", lw=2.2)
        ax.axvline(ISO, ls="--", color="w", lw=1.3)
        ax.axhline(BRUSH_EDGE, ls=":", color="0.85", lw=1.0)
        ax.set_xlim(0, 90); ax.set_ylim(0, ZMAX); ax.set_xticks([0, 30, 60, 90])
        ax.tick_params(labelsize=FS_T)
        if ri == 0:
            ax.text(3, ZMAX - 3, "flat", color="w", fontsize=FS_T, ha="left", va="top", fontweight="bold")
            ax.text(87, ZMAX - 3, "edge-on", color="w", fontsize=FS_T, ha="right", va="top", fontweight="bold")
        if ci == 0:
            ax.set_ylabel("solute z-coordinate (Å)", fontsize=FS_AX)
            ax.text(-0.36, 0.5, rlabel, transform=ax.transAxes, rotation=90, ha="center",
                    va="center", fontsize=FS_ROW, fontweight="bold", color=col)
        if ri == len(ROWS) - 1:
            ax.set_xlabel("angle (deg)", fontsize=FS_AX)

fig.subplots_adjust(left=0.17, right=0.90, top=0.82, bottom=0.05)
ICON_BOTTOM, ICON_H = 0.835, 0.14
for ci, C in enumerate(COLS):
    img = mpimg.imread(C["icon"]); icon_h = ICON_H * C["scale"]
    w_frac = icon_h * fig.get_figheight() * img.shape[1] / img.shape[0] / fig.get_figwidth()
    p = axes[0, ci].get_position(); xc = (p.x0 + p.x1) / 2; yc = ICON_BOTTOM + ICON_H / 2
    a = fig.add_axes([xc - w_frac / 2, yc - icon_h / 2, w_frac, icon_h]); a.imshow(img); a.axis("off")
cax = fig.add_axes([0.915, 0.05, 0.02, 0.77])
cb = fig.colorbar(pc, cax=cax); cb.set_ticks([]); cb.set_label("angle occurrence", fontsize=FS_AX)
for ext in ("png", "svg", "pdf"):
    fig.savefig(f"combined_orientation_allsolvents.{ext}", dpi=220)
plt.close(fig)

# ---------------------------------------------------------------- mean-angle / S(z) comparison
fig, axes = plt.subplots(2, 3, figsize=(12.5, 9.0), sharey=True, gridspec_kw=dict(hspace=0.25, wspace=0.08))
for ci, C in enumerate(COLS):
    for rk, rlabel, col, ls in ROWS:
        d = DATA[(ci, rk)]
        for ri, (key, iso) in enumerate([("th_b", ISO), ("S_b", 0.0)]):
            m, se = d[key]; ax = axes[ri, ci]
            ax.plot(m, zc, color=col, ls=ls, lw=2.0)
            ax.fill_betweenx(zc, m - se, m + se, color=col, alpha=0.15, lw=0)
    for ri, (lab, iso, lim) in enumerate([("mean angle (deg)", ISO, (0, 90)),
                                          (r"$S = \langle 3\cos^2\theta - 1\rangle/2$", 0.0, (-0.5, 1.0))]):
        ax = axes[ri, ci]
        ax.axvline(iso, color="0.5", ls=":", lw=1.2)
        ax.axhspan(0, BRUSH_EDGE, color="0.9", zorder=0)
        ax.set_xlim(*lim); ax.set_ylim(0, ZMAX); ax.set_xlabel(lab, fontsize=12)
        ax.grid(True, ls="--", alpha=0.4)
        if ri == 0:
            ax.set_title(C["sub"], fontsize=12)
    axes[0, 0].set_ylabel("solute z-coordinate (Å)", fontsize=12)
    axes[1, 0].set_ylabel("solute z-coordinate (Å)", fontsize=12)
handles = [Line2D([], [], color=c, ls=l, lw=2, label=lab.replace("\n", ", ")) for _, lab, c, l in ROWS]
fig.legend(handles=handles, loc="upper center", ncol=4, fontsize=10.5, frameon=False, bbox_to_anchor=(0.5, 0.995))
fig.subplots_adjust(left=0.07, right=0.98, top=0.90, bottom=0.07)
for ext in ("png", "svg", "pdf"):
    fig.savefig(f"orientation_mix_vs_water.{ext}", dpi=200)
plt.close(fig)

# ---------------------------------------------------------------- numbers
out = ["Orientation: 50:50 H2O:MeOH vs 100% H2O  (mean ± window-level SE; n = windows)",
       f"brush: z < {BRUSH_EDGE:.0f} Å, solvent: z ≥ {BRUSH_EDGE + 5:.0f} Å; isotropic: <theta> = 57.3, S = 0", ""]
for ci, C in enumerate(COLS):
    out.append(f"== {C['name']}")
    for rk, rlabel, *_ in ROWS:
        d = DATA[(ci, rk)]; row = f"  {rlabel.replace(chr(10), ', ').replace('$_2$', '2'):30s} frames {d['nframes']:6d}"
        for nm, lo, hi in [("brush", 0, BRUSH_EDGE), ("solvent", BRUSH_EDGE + 5, 1e9)]:
            t = region(d["zw"], d["thw"], lo, hi); s = region(d["zw"], d["Sw"], lo, hi)
            row += f" | {nm}: theta {t[0]:5.1f}±{t[1]:3.1f}  S {s[0]:+.3f}±{s[1]:.3f} (n={t[2]})"
        out.append(row)
    for prot in ("half", "full"):
        a, b = DATA[(ci, f"{prot}_mix")], DATA[(ci, f"{prot}_wat")]
        diff = b["th_b"][0] - a["th_b"][0]
        err = np.sqrt(a["th_b"][1] ** 2 + b["th_b"][1] ** 2)
        k = np.nanargmax(np.abs(diff))
        sig = np.abs(diff) > 2 * err
        out.append(f"  {prot}: water - mix mean angle per 5-Å bin: max |Δ| = {diff[k]:+.1f}±{err[k]:.1f} deg at "
                   f"z={zc[k]:.1f} Å; bins with |Δ|>2SE: {', '.join(f'{z:.0f}' for z in zc[sig]) or 'none'}")
    out.append("")
txt = "\n".join(out)
open("orientation_mix_vs_water_summary.txt", "w").write(txt + "\n")
print(txt)
print("wrote combined_orientation_allsolvents.*, orientation_mix_vs_water.*, orientation_mix_vs_water_summary.txt")
