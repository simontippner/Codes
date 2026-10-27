#!/usr/bin/env python3
"""
Combined P(theta | z) orientation heat maps, both brush protonation states and both
solvents in one figure (same style as ../combined_figure/combined_orientation_v3.py):
    columns = BODIPY dipyrrin-core plane | BODIPY phenyl-ring plane
              | [Mo3S13]2- basal (Mo3) plane
    rows    = partially prot. 50:50 H2O:MeOH | partially prot. 100% H2O
              | fully prot. 50:50 H2O:MeOH | fully prot. 100% H2O

theta = angle(plane normal, z), folded to [0,90]
    0  -> plane parallel to substrate (flat / face-on)
    90 -> plane orthogonal to substrate (edge-on / standing up)
Per-panel conditional probability P(theta|z) (each z-bin normalised to sum 1),
magma colormap on one common 0-vmax scale across all panels; cyan = z-binned
mean angle; dashed white = isotropic folded-cosine mean (57.3 deg).
Source: rc.dat + plane-normal vector file in each system's
.../analysis/orientation[/core_phenyl] directory.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

try:
    from scipy.ndimage import gaussian_filter
    smooth = lambda A: gaussian_filter(A, sigma=(0.8, 0.8), mode="nearest")
except Exception:
    smooth = lambda A: A

BOD = "/gpfs/data/stippner/MoS_BODIPY/BODIPY_in_brush/BIGGER/bigger_bottom"
MOS = "/gpfs/data/stippner/MoS_BODIPY/MoS_in_brush"

COLS = [
    dict(group="BODIPY", sub=r"dipyrrin normal $\angle$ z", vec="core.dat",
         half_mix=f"{BOD}/RETRY_half/MD/traj/analysis/orientation/core_phenyl",
         half_wat=f"{BOD}/RETRY_half_water/MD/traj/analysis/orientation/core_phenyl",
         full_mix=f"{BOD}/RETRY_full/MD/new/mixture_run/traj/analysis/orientation/core_phenyl",
         full_wat=f"{BOD}/RETRY_full_water/MD/traj/analysis/orientation/core_phenyl"),
    dict(group="BODIPY", sub=r"phenyl normal $\angle$ z", vec="phenyl.dat",
         half_mix=f"{BOD}/RETRY_half/MD/traj/analysis/orientation/core_phenyl",
         half_wat=f"{BOD}/RETRY_half_water/MD/traj/analysis/orientation/core_phenyl",
         full_mix=f"{BOD}/RETRY_full/MD/new/mixture_run/traj/analysis/orientation/core_phenyl",
         full_wat=f"{BOD}/RETRY_full_water/MD/traj/analysis/orientation/core_phenyl"),
    dict(group=r"[Mo$_3$S$_{13}$]$^{2-}$", sub=r"Mo$_3$ normal $\angle$ z", vec="norm.dat",
         half_mix=f"{MOS}/RETRY_half/MD/traj/analysis/orientation",
         half_wat=f"{MOS}/RETRY_half_water/MD/traj/analysis/orientation",
         full_mix=f"{MOS}/RETRY_full/new_final/tleap_run/traj/analysis/orientation",
         full_wat=f"{MOS}/RETRY_full_water/MD/traj/analysis/orientation"),
]
ROWS = [("half_mix", "partially prot. brush\n50:50 H$_2$O:MeOH"),
        ("half_wat", "partially prot. brush\n100% H$_2$O"),
        ("full_mix", "fully prot. brush\n50:50 H$_2$O:MeOH"),
        ("full_wat", "fully prot. brush\n100% H$_2$O")]

ZMAX = 110.0
zedges = np.arange(0, ZMAX + 5.0, 5.0)
zc = 0.5 * (zedges[:-1] + zedges[1:])
aedges = np.arange(0, 90.001, 5.0)


def unit(v):
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def load_theta_z(D, vecfile):
    rc = np.loadtxt(f"{D}/rc.dat")[:, 1:4]
    vec = np.loadtxt(f"{D}/{vecfile}")[:, 1:4]
    z = rc[:, 2]
    th = np.degrees(np.arccos(np.clip(np.abs(unit(vec)[:, 2]), 0, 1)))
    return z, th


PANELS = {}
for ci, C in enumerate(COLS):
    for rk, _ in ROWS:
        z, th = load_theta_z(C[rk], C["vec"])
        ib = np.digitize(z, zedges) - 1
        H, _, _ = np.histogram2d(z, th, bins=[zedges, aedges])
        H = smooth(H)
        with np.errstate(invalid="ignore", divide="ignore"):
            Hn = np.ma.masked_invalid(H / H.sum(axis=1, keepdims=True))
        mean_th = np.array([th[ib == k].mean() if (ib == k).any() else np.nan
                             for k in range(len(zc))])
        PANELS[(ci, rk)] = (Hn, mean_th)

vmax = max(np.nanmax(Hn.filled(0)) for Hn, _ in PANELS.values())

fig, axes = plt.subplots(4, 3, figsize=(11.4, 14.6), sharex=True, sharey=True,
                         gridspec_kw=dict(hspace=0.12, wspace=0.10))
FS_AX, FS_T, FS_H, FS_P, FS_ROW = 13, 11, 15, 12, 12.5

for ci, C in enumerate(COLS):
    axes[0, ci].set_title(C["sub"], fontsize=FS_P, pad=8)
    for ri, (rk, rlabel) in enumerate(ROWS):
        ax = axes[ri, ci]
        Hn, mean_th = PANELS[(ci, rk)]
        pc = ax.pcolormesh(aedges, zedges, Hn, cmap="magma", shading="flat",
                           vmin=0, vmax=vmax, edgecolors="face", linewidth=0,
                           rasterized=True)  # no cell seams in pdf/svg
        ax.plot(mean_th, zc, "-", color="cyan", lw=2.2)
        ax.axvline(np.degrees(1.0), ls="--", color="w", lw=1.3)
        ax.set_xlim(0, 90); ax.set_ylim(0, ZMAX)
        ax.set_xticks([0, 30, 60, 90])
        ax.tick_params(labelsize=FS_T)
        if ri == 0:
            ax.text(3, ZMAX - 3, "flat", color="w", fontsize=FS_T,
                    ha="left", va="top", fontweight="bold")
            ax.text(87, ZMAX - 3, "edge-on", color="w", fontsize=FS_T,
                    ha="right", va="top", fontweight="bold")
        if ci == 0:
            ax.set_ylabel("solute z-coordinate (Å)", fontsize=FS_AX)
        if ri == len(ROWS) - 1:
            ax.set_xlabel("angle (deg)", fontsize=FS_AX)

for ri, (rk, rlabel) in enumerate(ROWS):
    axes[ri, 0].text(-0.36, 0.5, rlabel, transform=axes[ri, 0].transAxes,
                     rotation=90, ha="center", va="center",
                     fontsize=FS_ROW, fontweight="bold", color="0.25")

fig.subplots_adjust(left=0.15, right=0.90, top=0.82, bottom=0.05)

# structure renders above each column, with that column's fitted plane highlighted
ICON_BOTTOM, ICON_H = 0.835, 0.135


def add_icon(ax_col, fname, scale=1.0):
    img = mpimg.imread(fname)
    icon_h = ICON_H * scale
    h_in = icon_h * fig.get_figheight()
    w_in = h_in * (img.shape[1] / img.shape[0])
    w_frac = w_in / fig.get_figwidth()
    p = ax_col.get_position()
    xc = (p.x0 + p.x1) / 2
    yc = ICON_BOTTOM + ICON_H / 2
    a = fig.add_axes([xc - w_frac / 2, yc - icon_h / 2, w_frac, icon_h])
    a.imshow(img); a.axis("off")


add_icon(axes[0, 0], "/gpfs/data/stippner/MoS_BODIPY/combined_figure/bodipy_core_annot.png")
add_icon(axes[0, 1], "/gpfs/data/stippner/MoS_BODIPY/combined_figure/bodipy_phenyl_annot.png")
add_icon(axes[0, 2], "/gpfs/data/stippner/MoS_BODIPY/combined_figure/mos_plane.png", scale=0.85)

cax = fig.add_axes([0.915, 0.05, 0.02, 0.77])
cb = fig.colorbar(pc, cax=cax)
cb.set_ticks([])
cb.set_label("angle occurrence", fontsize=FS_AX)

for ext in ("png", "svg", "pdf"):
    fig.savefig(f"combined_orientation_v3_both.{ext}", dpi=220)
print("wrote combined_orientation_v3_both.png / .svg / .pdf")
