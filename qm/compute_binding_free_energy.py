#!/usr/bin/env python3
"""Compute a complex's binding free energy and association constant from
three ORCA frequency-job outputs (donor, acceptor, complex), computational
analogue of a Benesi-Hildebrand-derived K_assoc.

dG_bind = G(complex) - G(fragment1) - G(fragment2)
K_assoc = exp(-dG_bind / RT)

Usage:
    python3 compute_binding_free_energy.py FRAGMENT1_FREQ.log FRAGMENT2_FREQ.log COMPLEX_FREQ.log [--temp 298.15]

Each input must be an ORCA %freq job's .log/.out (same level of theory /
same solvent model for all three, or the comparison is not meaningful).
"""
import argparse
import math

from extract_thermo import extract_thermo

HARTREE_TO_KCAL = 627.5094740631
R_KCAL_MOL_K = 1.98720425864e-3  # kcal/(mol*K)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("fragment1")
    ap.add_argument("fragment2")
    ap.add_argument("complex_")
    ap.add_argument("--temp", type=float, default=None,
                     help="Temperature in K (default: read from the complex "
                          "job's output)")
    args = ap.parse_args()

    t1 = extract_thermo(args.fragment1)
    t2 = extract_thermo(args.fragment2)
    tc = extract_thermo(args.complex_)

    g1 = t1["gibbs_free_energy_Eh"]
    g2 = t2["gibbs_free_energy_Eh"]
    gc = tc["gibbs_free_energy_Eh"]

    dG_Eh = gc - g1 - g2
    dG_kcal = dG_Eh * HARTREE_TO_KCAL

    temp = args.temp or tc.get("temperature_K") or 298.15
    K_assoc = math.exp(-dG_kcal / (R_KCAL_MOL_K * temp))

    print(f"Fragment 1 ({args.fragment1}): G = {g1: .8f} Eh")
    print(f"Fragment 2 ({args.fragment2}): G = {g2: .8f} Eh")
    print(f"Complex    ({args.complex_}):  G = {gc: .8f} Eh")
    print()
    print(f"dG_bind = {dG_Eh: .8f} Eh  =  {dG_kcal: .3f} kcal/mol")
    print(f"T       = {temp:.2f} K")
    print(f"K_assoc = exp(-dG_bind / RT) = {K_assoc:.3e}  (dimensionless, "
          f"standard-state-convention caveats apply)")
    print()
    if dG_kcal < 0:
        print("Negative dG_bind: binding is predicted to be favorable "
              "(associates spontaneously at this level of theory).")
    else:
        print("Positive dG_bind: binding is predicted to be unfavorable "
              "at this level of theory.")
    print()
    print("Caveat: single rigid docked geometry, no conformational search, "
          "continuum (not explicit) solvation -- treat this as an "
          "order-of-magnitude estimate, not a number directly comparable "
          "to an experimental Benesi-Hildebrand-fit K_assoc.")


if __name__ == "__main__":
    main()
