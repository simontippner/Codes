import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load(path):
    wl, eps = [], []
    with open(path) as f:
        for row in csv.DictReader(f):
            wl.append(float(row["wavelength_nm"]))
            eps.append(float(row["eps_L_mol-1_cm-1"]))
    return wl, eps

rubpy_wl, rubpy_eps = load("/gpfs/data/stippner/Luca/04_uvvis/rubpy/qm/rubpy3_uvvis_spectrum.csv")
mos_wl, mos_eps = load("/gpfs/data/stippner/Luca/04_uvvis/mos/qm/mo3s13_uvvis_spectrum.csv")
comb_wl, comb_eps = load("rubpy_mos_uvvis_220roots_spectrum.csv")

assert rubpy_wl == mos_wl, "grids must match for a point-by-point sum"
naive_sum = [r + m for r, m in zip(rubpy_eps, mos_eps)]

# Lowest wavelength actually reached by the combined (220-root) calc
lowest_computed = min(
    float(row["wavelength_nm"])
    for row in csv.DictReader(open("rubpy_mos_uvvis_220roots_sticks.csv"))
)

fig, ax = plt.subplots(figsize=(7.2, 4.5))
ax.plot(rubpy_wl, naive_sum, "--", color="#666", lw=1.3, label="rubpy + mos (isolated, summed)")
ax.plot(comb_wl, comb_eps, "-", color="#1e4e8c", lw=1.8, label="rubpy…mos combined (computed, 220 roots)")
ax.axvline(lowest_computed, color="#d6a93a", ls=":", lw=1)
ax.axvspan(200, lowest_computed, color="#d6453a", alpha=0.08)
ax.text((200 + lowest_computed) / 2, ax.get_ylim()[1] if False else None, "", visible=False)
ax.set_xlim(200, 700)
ax.set_xlabel("Wavelength (nm)")
ax.set_ylabel(r"$\varepsilon$ (L mol$^{-1}$ cm$^{-1}$)")
ax.set_title("Combined spectrum (220 roots) vs. naive sum of isolated fragments")
ax.legend(loc="upper right", frameon=False, fontsize=9)
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
plt.tight_layout()
plt.savefig("rubpy_mos_uvvis_220roots_comparison.png", dpi=200)
plt.savefig("rubpy_mos_uvvis_220roots_comparison.pdf")
print("lowest_computed_nm", lowest_computed)
print("saved comparison figure")
