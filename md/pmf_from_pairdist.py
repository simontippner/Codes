#!/usr/bin/env python3
"""Approximate CAT...PS binding free energy (potential of mean force, PMF)
from an EXISTING unbiased MD trajectory's pairwise center-of-mass distance
distribution -- no new simulation needed, just re-analysis of the
`com_distances.dat` files this project's own aggregation pipeline
(01_generate_cpptraj/gen_mix_cpptraj.py / gen_mix10_cpptraj.py) already
writes.

Method: for N independent CAT...PS pairs sampled every frame, the expected
pair count in a shell [r, r+dr) under a uniform (uncorrelated / "ideal gas")
reference is  n_frames * N_pairs * (4/3 pi (r+dr)^3 - r^3) / V_box . The
ratio of the OBSERVED count to that expectation is g(r) (same quantity as a
site-site radial distribution function, but for whole-body centers of mass
rather than atom pairs). PMF(r) = -RT ln[g(r) / g(r->ref)], shifted so the
plateau at large r (where the pair is far enough apart to behave like two
independent, uncorrelated bodies) sits at zero -- the depth of the resulting
well at the contact-distance minimum is the pairwise association free
energy relative to the "dissociated" reference state, estimated directly
from how often the unbiased trajectory visits close vs. far separations.

This is NOT a replacement for proper umbrella-sampling/WHAM PMF along a
well-defined reaction coordinate, nor for a true single-pair (infinite
dilution) binding free energy -- see the caveats in the module docstring of
the calling report section. It is a legitimate, standard way to extract an
approximate association free energy from an unbiased trajectory that
already happens to sample both bound and unbound configurations, without
running any new simulation.

Usage:
    python3 pmf_from_pairdist.py COM_DISTANCES_DAT BOX_VOLUME_A3 \\
        --pair-regex 'd_(CAT\\d*)_(RUB\\d+)' --r-max 25 --nbins 100 \\
        --out-csv pmf.csv
"""
import argparse
import re
import numpy as np


def compute_pmf(dat_path, box_vol, pair_regex, r_max=25.0, nbins=100, T=300.0):
    with open(dat_path) as f:
        header = f.readline().split()[1:]
    pat = re.compile(pair_regex)
    idx = [i for i, name in enumerate(header) if pat.match(name)]
    if not idx:
        raise ValueError(f"No columns in {dat_path} matched regex {pair_regex!r}")

    data = np.loadtxt(dat_path, comments="#")
    n_frames = data.shape[0]
    dists = data[:, [i + 1 for i in idx]]  # +1: skip leading Frame column
    all_d = dists.ravel()
    all_d = all_d[all_d < r_max]

    n_pairs = len(idx)
    counts, edges = np.histogram(all_d, bins=nbins, range=(0, r_max))
    r_lo, r_hi = edges[:-1], edges[1:]
    shell_vol = (4.0 / 3.0) * np.pi * (r_hi**3 - r_lo**3)
    expected = n_frames * n_pairs * (shell_vol / box_vol)
    g = counts / expected
    r_mid = 0.5 * (r_lo + r_hi)

    RT = 1.987204e-3 * T  # kcal/mol
    with np.errstate(divide="ignore"):
        pmf = -RT * np.log(g)
    plateau_mask = r_mid > 0.7 * r_max
    finite_plateau = pmf[plateau_mask & np.isfinite(pmf)]
    pmf_ref = np.nanmean(finite_plateau) if len(finite_plateau) else np.nan
    pmf_shifted = pmf - pmf_ref

    valid = np.isfinite(pmf_shifted)
    min_idx = np.nanargmin(np.where(valid, pmf_shifted, np.inf)) if valid.any() else None
    return dict(r=r_mid, pmf=pmf_shifted, g=g, counts=counts, n_frames=n_frames,
                n_pairs=n_pairs, min_idx=min_idx, pmf_ref=pmf_ref)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("com_distances_dat")
    ap.add_argument("box_volume_a3", type=float)
    ap.add_argument("--pair-regex", default=r"d_(CAT\d*)_(RUB\d+)$",
                     help="regex matching header columns to include as CAT...PS pair distances")
    ap.add_argument("--r-max", type=float, default=25.0)
    ap.add_argument("--nbins", type=int, default=100)
    ap.add_argument("--temperature", type=float, default=300.0)
    ap.add_argument("--out-csv", default=None)
    args = ap.parse_args()

    res = compute_pmf(args.com_distances_dat, args.box_volume_a3, args.pair_regex,
                       args.r_max, args.nbins, args.temperature)
    if res["min_idx"] is None:
        print("No valid PMF minimum found (check r-max / data).")
        return
    i = res["min_idx"]
    print(f"n_frames={res['n_frames']}  n_pairs={res['n_pairs']}")
    print(f"PMF minimum: {res['pmf'][i]:.2f} kcal/mol at r = {res['r'][i]:.2f} A "
          f"(g(r)={res['g'][i]:.2f}, counts={res['counts'][i]})")

    if args.out_csv:
        import csv
        with open(args.out_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["r_A", "g_r", "pmf_kcal_mol", "counts"])
            for r, g, p, c in zip(res["r"], res["g"], res["pmf"], res["counts"]):
                w.writerow([f"{r:.3f}", f"{g:.4f}", f"{p:.4f}" if np.isfinite(p) else "", c])
        print("wrote", args.out_csv)


if __name__ == "__main__":
    main()
