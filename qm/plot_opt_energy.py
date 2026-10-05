#!/usr/bin/env python3
"""
Plot the SCF energy at every cycle of an ORCA geometry optimization, read
directly from the (possibly still-running) log file.

Works even while the optimization job is still running, since it only needs
the live-written log file on the shared filesystem -- not the job's
node-local scratch directory, which usually isn't accessible until the job
finishes (same property as make_trj_xyz.py).

Usage:
    python3 plot_opt_energy.py FILENAME [--out energy]

Outputs:
    <out>.csv   cycle, energy_Eh, energy_eV, relative_kcal_mol
    <out>.png / <out>.pdf   absolute SCF energy (Eh) vs cycle
"""
import argparse
import os
import re
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EH_TO_KCALMOL = 627.509474


def parse_energies(log_path):
    """Return list of (cycle_index, energy_Eh, energy_eV), one entry per
    'TOTAL SCF ENERGY' block in the file (one such block is printed after
    every geometry optimization cycle's SCF converges)."""
    row_re = re.compile(
        r"Total Energy\s*:\s*(-?\d+\.\d+)\s*Eh\s*(-?\d+\.\d+)\s*eV")
    energies = []
    with open(log_path) as f:
        for line in f:
            m = row_re.search(line)
            if m:
                energies.append((float(m.group(1)), float(m.group(2))))
    if not energies:
        sys.exit(f"No 'Total Energy : ... Eh ... eV' lines found in {log_path}")
    return [(i, e_eh, e_ev) for i, (e_eh, e_ev) in enumerate(energies)]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("filename", help="ORCA geometry-optimization output/log file")
    p.add_argument("--out", default=None,
                    help="output basename (default: derived from the input filename)")
    args = p.parse_args()
    if args.out is None:
        args.out = re.sub(r"\.(log|out)$", "", os.path.basename(args.filename)) + "_energy"

    rows = parse_energies(args.filename)
    print(f"Parsed {len(rows)} cycles from {args.filename}")

    e0_eh = rows[0][1]
    with open(f"{args.out}.csv", "w") as f:
        f.write("cycle,energy_Eh,energy_eV,relative_kcal_mol\n")
        for cycle, e_eh, e_ev in rows:
            rel_kcal = (e_eh - e0_eh) * EH_TO_KCALMOL
            f.write(f"{cycle},{e_eh:.8f},{e_ev:.4f},{rel_kcal:.4f}\n")

    cycles = [r[0] for r in rows]
    energies_eh = [r[1] for r in rows]
    final_eh = rows[-1][1]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(cycles, energies_eh, color="#1e4e8c", lw=1.5, marker="o", markersize=3)
    ax.axhline(final_eh, color="0.7", lw=0.8, ls=(0, (4, 3)), zorder=0)
    ax.set_xlabel("Optimization cycle")
    ax.set_ylabel("SCF energy (Eh)", fontsize=9)
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    ax.set_title(f"{os.path.basename(args.filename)}  "
                  f"(current: {final_eh:.6f} Eh, cycle {cycles[-1]})", fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{args.out}.png", dpi=300)
    fig.savefig(f"{args.out}.pdf")
    print(f"Wrote {args.out}.csv, {args.out}.png, {args.out}.pdf")


if __name__ == "__main__":
    main()
