#!/usr/bin/env python3
"""Rehm-Weller-style charge-transfer driving force for a pre-formed
(docked/complexed) donor-acceptor pair, from 4 vertical single-point ORCA
jobs at the SAME frozen geometry (no optimization):
  donor neutral (its ground charge state), donor oxidized (+1 e- removed),
  acceptor neutral (its ground charge state), acceptor reduced (+1 e- added).

    IP_vertical(donor)    = E(donor, oxidized) - E(donor, neutral)
    EA_vertical(acceptor) = E(acceptor, neutral) - E(acceptor, reduced)
    dG_CT (no Coulomb)    = IP_vertical(donor) - EA_vertical(acceptor)

This is the cost to move one electron donor->acceptor at infinite
separation. For a complex at a real through-space separation r, the
resulting ion pair's Coulomb attraction lowers this; two bounds are
reported since how much the surrounding solvent screens a THROUGH-SPACE
term between two already-CPCM-embedded fragments is genuinely ambiguous:
  - vacuum (unscreened, e/(4pi eps0 r)): upper-bound stabilization
  - bulk continuum dielectric screened (divide by the solvent's static
    eps_r): lower-bound stabilization, likely an over-correction at
    short/contact-ion-pair separations where the continuum approximation
    is least valid.

Usage:
    python3 compute_ct_driving_force.py DONOR_NEUTRAL.log DONOR_OXIDIZED.log \
        ACCEPTOR_NEUTRAL.log ACCEPTOR_REDUCED.log --distance 6.5 [--eps 36.6]
"""
import argparse
import re

HARTREE_TO_EV = 27.211386245988
HARTREE_TO_KCAL = 627.5094740631
BOHR_PER_ANG = 1 / 0.529177210903


def final_sp_energy(path):
    with open(path) as f:
        text = f.read()
    matches = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text)
    if not matches:
        raise ValueError(f"{path}: no 'FINAL SINGLE POINT ENERGY' found")
    return float(matches[-1])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("donor_neutral")
    ap.add_argument("donor_oxidized")
    ap.add_argument("acceptor_neutral")
    ap.add_argument("acceptor_reduced")
    ap.add_argument("--distance", type=float, required=True,
                     help="donor-acceptor center-to-center distance, Angstrom")
    ap.add_argument("--eps", type=float, default=36.6,
                     help="solvent static dielectric constant (default: acetonitrile)")
    args = ap.parse_args()

    e_d0 = final_sp_energy(args.donor_neutral)
    e_d1 = final_sp_energy(args.donor_oxidized)
    e_a0 = final_sp_energy(args.acceptor_neutral)
    e_a1 = final_sp_energy(args.acceptor_reduced)

    ip = e_d1 - e_d0
    ea = e_a0 - e_a1
    dg_noCoulomb_eV = (ip - ea) * HARTREE_TO_EV

    r_bohr = args.distance * BOHR_PER_ANG
    e_coul_vac_eV = (-1.0 / r_bohr) * HARTREE_TO_EV
    e_coul_screened_eV = e_coul_vac_eV / args.eps

    dg_vac = dg_noCoulomb_eV + e_coul_vac_eV
    dg_screened = dg_noCoulomb_eV + e_coul_screened_eV

    print(f"IP_vertical(donor)    = {ip:.6f} Eh = {ip*HARTREE_TO_EV:.4f} eV")
    print(f"EA_vertical(acceptor) = {ea:.6f} Eh = {ea*HARTREE_TO_EV:.4f} eV")
    print(f"dG_CT (no Coulomb)    = {dg_noCoulomb_eV:.4f} eV "
          f"= {(ip-ea)*HARTREE_TO_KCAL:.2f} kcal/mol")
    print(f"E_coulomb (vacuum, r={args.distance} A)        = {e_coul_vac_eV:.4f} eV")
    print(f"E_coulomb (bulk-screened, eps={args.eps})  = {e_coul_screened_eV:.4f} eV")
    print(f"dG_CT (vacuum-Coulomb bound, upper stabiliz.)  = {dg_vac:.4f} eV")
    print(f"dG_CT (bulk-screened bound, lower stabiliz.)   = {dg_screened:.4f} eV")
    print()
    print("Caveat: both are approximations. Vacuum likely overestimates the "
          "driving force (no screening at all); bulk dielectric likely "
          "underestimates it (continuum screening is least valid at short, "
          "contact-ion-pair separations). True value lies somewhere between.")


if __name__ == "__main__":
    main()
