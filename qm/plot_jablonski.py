import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Energies in eV, all from this project's TD-DFT/triplet-TD-DFT results
systems = [
    {
        "name": "rubpy",
        "label": "Ru(bpy)$_3^{2+}$\n(isolated)",
        "s1_ev": 2.825, "s1_nm": 438.8, "s1_f": 0.131,
        "t1_ev": 2.313,
        "s1_char": "$^1$MLCT\n(Ru$\\rightarrow$bpy)",
        "t1_char": "$^3$MLCT\n(Ru$\\rightarrow$bpy)",
    },
    {
        "name": "mos",
        "label": "[Mo$_3$S$_{13}$]$^{2-}$\n(isolated)",
        "s1_ev": 2.680, "s1_nm": 462.7, "s1_f": 0.016,
        "t1_ev": 1.875,
        "s1_char": "LMCT\n(S$\\rightarrow$Mo)",
        "t1_char": "LMCT\n(S$\\rightarrow$Mo)",
    },
    {
        "name": "combined",
        "label": "rubpy…mos\ncombined",
        "s1_ev": 2.766, "s1_nm": 448.3, "s1_f": 0.098,
        "t1_ev": 1.747,
        "s1_char": "$^1$MLCT\n(§3e, mostly\nrubpy-localized)",
        "t1_char": "T$_1$\n(not yet\nNTO-checked)",
    },
]

rubpy_c = "#1e4e8c"
mos_c = "#d6a93a"
combined_c = "#4a9f5c"
colors = [rubpy_c, mos_c, combined_c]

fig, ax = plt.subplots(figsize=(10, 6))

xw = 1.3  # half-width of each level line
for i, sysd in enumerate(systems):
    xc = i * 4.3
    c = colors[i]

    # S0
    ax.hlines(0, xc - xw, xc + xw, color="#333", lw=2.2)
    ax.text(xc, -0.17, "S$_0$", ha="center", fontsize=9.5)

    # S1 (bright singlet)
    ax.hlines(sysd["s1_ev"], xc - xw, xc + xw, color=c, lw=2.2)
    ax.text(xc + xw + 0.05, sysd["s1_ev"], f"S$_1$  {sysd['s1_ev']:.3f} eV\n({sysd['s1_nm']} nm, f={sysd['s1_f']})",
            fontsize=7.8, va="center", color="#333")
    ax.text(xc, sysd["s1_ev"] + 0.08, sysd["s1_char"], fontsize=7.3, ha="center", color=c)

    # T1
    ax.hlines(sysd["t1_ev"], xc - xw, xc + xw, color=c, lw=2.2, linestyles="--")
    ax.text(xc + xw + 0.05, sysd["t1_ev"], f"T$_1$  {sysd['t1_ev']:.3f} eV",
            fontsize=7.8, va="center", color="#333")
    ax.text(xc, sysd["t1_ev"] - 0.20, sysd["t1_char"], fontsize=7.3, ha="center", color=c)

    # absorption arrow S0->S1
    ax.annotate("", xy=(xc - xw * 0.45, sysd["s1_ev"]), xytext=(xc - xw * 0.45, 0),
                arrowprops=dict(arrowstyle="-|>", color=c, lw=1.8))

    # ISC arrow S1->T1 (wavy-ish dashed diagonal)
    ax.annotate("", xy=(xc + xw * 0.45, sysd["t1_ev"]), xytext=(xc + xw * 0.2, sysd["s1_ev"]),
                arrowprops=dict(arrowstyle="-|>", color="#888", lw=1.4, linestyle=(0, (3, 2))))
    ax.text(xc + xw * 0.5, (sysd["s1_ev"] + sysd["t1_ev"]) / 2, "ISC",
            fontsize=7, color="#888", rotation=-60)

    ax.text(xc, -0.55, sysd["label"], ha="center", fontsize=10, fontweight="bold", color=c)

ax.set_xlim(-2.3, 11.8)
ax.set_ylim(-0.75, 3.3)
ax.set_ylabel("Energy (eV)")
ax.set_xticks([])
for spine in ["top", "right", "bottom"]:
    ax.spines[spine].set_visible(False)
ax.set_title("Jablonski diagram: S$_0$ / lowest bright S$_1$ / T$_1$\n"
              "(this project's TD-DFT results; ISC drawn schematically, no rate computed)",
              fontsize=10.5)

plt.tight_layout()
plt.savefig("jablonski_diagram.png", dpi=200)
plt.savefig("jablonski_diagram.pdf")
print("saved jablonski_diagram.png/.pdf")
