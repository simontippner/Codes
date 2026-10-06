#!/usr/bin/env python3
"""Estimate a purely-radiative (fluorescence) rate constant and lifetime
from a TD-DFT oscillator strength via the Strickler-Berg relation:

    k_r [s^-1] = 0.668 * nu_avg^2 [cm^-1] * n^2 * f

where nu_avg is the transition wavenumber, n is the solvent refractive
index, and f is the (dimensionless) oscillator strength.

IMPORTANT CAVEATS (read before using this number for anything):
  - This is an UPPER BOUND on the true excited-state lifetime
    (tau = 1/(k_r + k_ISC + k_nr) <= 1/k_r always) -- any real
    non-radiative or intersystem-crossing channel can only shorten the
    true lifetime further. This is NOT a prediction of the observed
    lifetime.
  - Uses the VERTICAL absorption oscillator strength as a stand-in for
    the (geometry-relaxed) emission oscillator strength -- an
    approximation. Strickler-Berg formally wants emission-band data.
  - Only valid for spin-allowed (e.g. singlet-singlet) transitions. A
    spin-forbidden transition (e.g. S0-T1 phosphorescence without
    explicit spin-orbit coupling) has f=0 in a plain TD-DFT calculation
    and this formula gives a meaningless infinite lifetime -- check f
    is not (numerically) zero before trusting the result.

Usage:
    python3 radiative_lifetime_estimate.py WAVELENGTH_NM FOSC [--n 1.344]
"""
import argparse


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("wavelength_nm", type=float)
    ap.add_argument("fosc", type=float)
    ap.add_argument("--n", type=float, default=1.344,
                     help="solvent refractive index (default: acetonitrile, 1.344)")
    args = ap.parse_args()

    if args.fosc <= 0:
        print("f is zero (or negative) -- this transition is spin-forbidden or "
              "otherwise has no usable oscillator strength; Strickler-Berg does "
              "not apply (would give an infinite/meaningless lifetime). A "
              "spin-orbit-coupling calculation is needed for a phosphorescence "
              "rate instead.")
        return

    nu = 1e7 / args.wavelength_nm  # cm^-1
    kr = 0.668 * nu**2 * args.n**2 * args.fosc  # s^-1
    tau = 1 / kr  # s

    print(f"Wavelength: {args.wavelength_nm:.1f} nm  ({nu:.1f} cm^-1)")
    print(f"f: {args.fosc:.4f}   n: {args.n}")
    print(f"k_r = {kr:.3e} s^-1")
    print(f"tau_r (radiative-only, upper bound on true lifetime) = "
          f"{tau*1e9:.2f} ns = {tau*1e6:.4f} us")


if __name__ == "__main__":
    main()
