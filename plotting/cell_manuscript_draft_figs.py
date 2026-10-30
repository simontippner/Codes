"""Draft/placeholder main-text figures for the Cell manuscript.
model.pdf: panel A placeholder (simulation snapshot), panel B BODIPY and
[Mo3S13]2- with their analysed sites. barebrush.pdf: placeholder grid A-E."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from PIL import Image

S = "/gpfs/data/stippner/MoS_BODIPY/combined_figure/"
RED = "#d62728"


def placeholder(ax, letter, text):
    ax.set_axis_off()
    ax.add_patch(FancyBboxPatch((0.02, 0.02), 0.96, 0.96, boxstyle="round,pad=0,rounding_size=0.03",
                                fill=False, ls="--", lw=1.2, ec=RED, transform=ax.transAxes))
    ax.text(0.5, 0.5, text, ha="center", va="center", color=RED, fontsize=9, wrap=True,
            transform=ax.transAxes)
    ax.text(0.04, 0.96, letter, ha="left", va="top", fontsize=13, weight="bold", transform=ax.transAxes)


# ---- model figure ----
fig = plt.figure(figsize=(7.0, 3.4))
axA = fig.add_axes([0.0, 0.0, 0.36, 1.0])
placeholder(axA, "A", "TODO\nsimulation snapshot:\nPDMAEMA brush (16 chains\n× 30 units) with solute,\nCl$^-$, H$_2$O / MeOH")
for i, (img, title, sites) in enumerate([
        ("bodipy_icon.png", "PS: BODIPY", "$-$COO$^-$ anchor · BF$_2$ · iodine"),
        ("mos_icon.png", "CAT: [Mo$_3$S$_{13}$]$^{2-}$", "apical $\\mu_3$-S · bridging $\\mu$-S$_2$ · terminal $\\eta^2$-S$_2$")]):
    ax = fig.add_axes([0.38 + i * 0.31, 0.12, 0.30, 0.76])
    ax.imshow(Image.open(S + img))
    ax.set_axis_off()
    ax.set_title(title, fontsize=10, weight="bold")
    ax.text(0.5, -0.06, sites, ha="center", va="top", fontsize=7.5, transform=ax.transAxes)
fig.text(0.38, 0.96, "B", fontsize=13, weight="bold", va="top")
fig.savefig("main_figs/model.pdf")

# ---- bare brush placeholder ----
fig, axs = plt.subplots(2, 3, figsize=(7.0, 4.4))
labels = ["brush height /\npolymer density profile\n(water)", "brush height /\npolymer density profile\n(H$_2$O:MeOH)",
          "Cl$^-$ condensation\n(NH$^+$–Cl$^-$ g(r),\ncontact ion-pair fraction)", "amine solvation\n(N/NH$^+$–solvent g(r))",
          "local MeOH fraction\nvs. z"]
for ax, L, t in zip(axs.flat, "ABCDE", labels):
    placeholder(ax, L, "TODO\n" + t)
axs.flat[-1].set_axis_off()
fig.tight_layout()
fig.savefig("main_figs/barebrush.pdf")
