#!/usr/bin/env python3
"""
Parse an ORCA TD-DFT output file (rubpy3_uvvis.log) and produce a
Gaussian-broadened UV-vis absorption spectrum from the computed vertical
excitation energies and oscillator strengths.

Usage:
    python3 plot_uvvis.py [--log rubpy3_uvvis.log] [--fwhm 0.4]
                           [--xmin 200] [--xmax 700] [--out spectrum]

Broadening follows the standard Gaussian-convolution formula used to turn
discrete TD-DFT oscillator strengths into a continuous molar-absorptivity-like
spectrum (same convention as used by Multiwfn/GaussView "UV-Vis spectrum"
tools):

    eps(E) = 1.3062974e8 * sum_i [ f_i / FWHM ] * exp(-4 ln2 * ((E-E_i)/FWHM)^2)

with E and FWHM in eV and eps in L mol^-1 cm^-1. This is a standard
convention for comparing relative band shapes/positions to experiment; it is
not a substitute for an explicit vibronic/Franck-Condon simulation.

Outputs:
    <out>_sticks.csv   raw transitions (state, nm, eV, fosc)
    <out>_spectrum.csv broadened spectrum (nm, eV, eps)
    <out>.png          plot: broadened spectrum + stick transitions
"""
import argparse
import re
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EV_PER_CM1 = 1.0 / 8065.54429  # cm^-1 -> eV
NM_PER_EV = 1239.84198         # nm = NM_PER_EV / eV


def parse_transitions(log_path):
    """Return list of (state, energy_cm1, wavelength_nm, fosc) from the
    dipole-LENGTH form absorption-spectrum table ("ABSORPTION SPECTRUM VIA
    TRANSITION ELECTRIC DIPOLE MOMENTS") -- the standard convention for
    reporting oscillator strengths. ORCA also prints a velocity-dipole-form
    table afterwards (gauge-origin independence check); that one is
    deliberately skipped.

    Row format (ORCA 6.x), note the transition label is e.g. "0-1A -> 12-1A"
    (not a bare integer), so the state index is parsed out of that label:
        Transition          Energy     Energy  Wavelength   fosc(D2)   D2   DX  DY  DZ
                             (eV)      (cm-1)    (nm)
        0-1A  ->  1-1A     2.479374   19997.5   500.1      0.001834871 ...
    """
    with open(log_path) as f:
        lines = f.readlines()

    header_idx = [i for i, l in enumerate(lines) if "ELECTRIC DIPOLE MOMENTS" in l.upper()]
    if not header_idx:
        header_idx = [i for i, l in enumerate(lines) if "ABSORPTION SPECTRUM" in l.upper()]
    if not header_idx:
        sys.exit(f"No 'ABSORPTION SPECTRUM' table found in {log_path}")
    start = header_idx[0]

    row_re = re.compile(
        r"^\s*\d+-\S*\s*->\s*(\d+)-\S*\s+"   # transition label, capture upper-state index
        r"(-?\d+\.\d+)\s+"                    # energy, eV
        r"(-?\d+\.\d+)\s+"                    # energy, cm-1
        r"(-?\d+\.\d+)\s+"                    # wavelength, nm
        r"([\d.eE+-]+)"                       # fosc
    )
    transitions = []
    in_table = False
    dash_count = 0
    for line in lines[start:]:
        if set(line.strip()) == {"-"} and line.strip():
            dash_count += 1
            # table is framed by dashed lines above+below the header/units;
            # once we've passed the second dashed line, data rows follow
            if dash_count >= 2:
                in_table = True
            continue
        if not in_table:
            continue
        m = row_re.match(line)
        if m:
            state = int(m.group(1))
            energy_cm1 = float(m.group(3))
            wavelength_nm = float(m.group(4))
            fosc = float(m.group(5))
            transitions.append((state, energy_cm1, wavelength_nm, fosc))
        else:
            if transitions:
                break  # blank/dashed line (or next table) after data -> end

    if not transitions:
        sys.exit(f"Found an ABSORPTION SPECTRUM header but no data rows in {log_path}")
    return transitions


def broaden(transitions, fwhm_ev, xmin_nm, xmax_nm, npoints):
    energies_ev = np.array([EV_PER_CM1 * t[1] for t in transitions])
    fosc = np.array([t[3] for t in transitions])

    x_nm = np.linspace(xmin_nm, xmax_nm, npoints)
    x_ev = NM_PER_EV / x_nm

    prefactor = 1.3062974e8
    eps = np.zeros_like(x_ev)
    for E_i, f_i in zip(energies_ev, fosc):
        eps += prefactor * (f_i / fwhm_ev) * np.exp(
            -4.0 * np.log(2.0) * ((x_ev - E_i) / fwhm_ev) ** 2
        )
    return x_nm, eps


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--log", default="rubpy3_uvvis.log",
                    help="ORCA TD-DFT output/log file (default: rubpy3_uvvis.log)")
    p.add_argument("--fwhm", type=float, default=0.4,
                    help="Gaussian FWHM in eV (default: 0.4)")
    p.add_argument("--xmin", type=float, default=200.0, help="plot xmin, nm")
    p.add_argument("--xmax", type=float, default=700.0, help="plot xmax, nm")
    p.add_argument("--npoints", type=int, default=2000)
    p.add_argument("--out", default="rubpy3_uvvis",
                    help="output basename (default: rubpy3_uvvis)")
    p.add_argument("--exp-ref", type=float, nargs="*", default=[286.0, 450.0],
                    help="experimental reference wavelengths (nm) to mark, "
                         "default: 286 450 (bpy pi-pi*, 1MLCT max)")
    args = p.parse_args()

    transitions = parse_transitions(args.log)
    print(f"Parsed {len(transitions)} transitions from {args.log}")

    with open(f"{args.out}_sticks.csv", "w") as f:
        f.write("state,energy_cm-1,wavelength_nm,fosc\n")
        for state, e_cm1, wl_nm, fosc in transitions:
            f.write(f"{state},{e_cm1:.2f},{wl_nm:.2f},{fosc:.6f}\n")

    x_nm, eps = broaden(transitions, args.fwhm, args.xmin, args.xmax, args.npoints)
    with open(f"{args.out}_spectrum.csv", "w") as f:
        f.write("wavelength_nm,eV,eps_L_mol-1_cm-1\n")
        for wl, e in zip(x_nm, eps):
            f.write(f"{wl:.3f},{NM_PER_EV/wl:.4f},{e:.4f}\n")

    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.plot(x_nm, eps, color="#1e4e8c", lw=1.8, zorder=3)
    ax.fill_between(x_nm, eps, color="#1e4e8c", alpha=0.12, zorder=2)
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel(r"$\varepsilon$ (L mol$^{-1}$ cm$^{-1}$)")
    ax.set_xlim(args.xmin, args.xmax)
    ax.set_ylim(bottom=0)

    ax2 = ax.twinx()
    stick_nm = [t[2] for t in transitions if args.xmin <= t[2] <= args.xmax]
    stick_f = [t[3] for t in transitions if args.xmin <= t[2] <= args.xmax]
    ax2.vlines(stick_nm, 0, stick_f, color="0.5", lw=1.0, zorder=1)
    ax2.set_ylabel("Oscillator strength (stick)", color="0.4")
    ax2.set_ylim(bottom=0)
    ax2.tick_params(axis="y", colors="0.4")

    for ref_nm in args.exp_ref:
        if args.xmin <= ref_nm <= args.xmax:
            ax.axvline(ref_nm, color="#c0392b", ls=(0, (4, 3)), lw=1.0, zorder=4)
            ax.text(ref_nm, ax.get_ylim()[1] * 0.97, f"{ref_nm:.0f} nm (exp.)",
                     rotation=90, ha="right", va="top", fontsize=7.5, color="#c0392b")

    fig.tight_layout()
    fig.savefig(f"{args.out}.png", dpi=300)
    fig.savefig(f"{args.out}.pdf")
    print(f"Wrote {args.out}_sticks.csv, {args.out}_spectrum.csv, "
          f"{args.out}.png, {args.out}.pdf")


if __name__ == "__main__":
    main()
