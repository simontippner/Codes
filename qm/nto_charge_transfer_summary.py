#!/usr/bin/env python3
"""Full per-state charge-transfer summary from an ORCA DoNTO TD-DFT job.

Combines three pieces: (1) the NTO pair weights ("n=") ORCA prints in the
main .log/.out for each requested state, (2) the actual NTO orbitals
written to "<base>.sN.nto" and converted to Molden with
`orca_2mkl <base>.sN -molden` for each state, and (3) a two-fragment atom
split (e.g. donor complex vs. acceptor cluster) -- to report, per state,
how much of the leading hole/particle NTO pair's density sits on each
fragment, i.e. whether the state is fragment-localized or genuine
donor->acceptor charge transfer.

Handles the index remapping between ORCA's log-printed pair labels (which
use the *original* canonical-orbital numbering, e.g. "255a -> 256a") and
the Molden/NTO file's own internal numbering (NTOs re-sorted into a new
pseudo-canonical set, symmetric around the hole/particle boundary) --
confirmed empirically: if nocc is that file's hole/particle boundary index,
rank-k pair (k=1 = dominant) sits at Molden indices (nocc-k+1, nocc+k).

Usage:
    python3 nto_charge_transfer_summary.py LOGFILE MOLDEN_PREFIX \
        --frag1 1-61 --frag2 62-77 \
        [--frag1-name rubpy] [--frag2-name mos] \
        [--states 29,30,41,43,44] [--max-ranks 4] [--csv out.csv]

MOLDEN_PREFIX + "_sN.molden.input" (or "<prefix>sN.molden.input" -- both
tried) must exist for each requested state; generate with:
    cp <base>.sN.nto <base>_sN.gbw && orca_2mkl <base>_sN -molden
"""
import argparse
import csv as csvmod
import re
import sys

from nto_fragment_analysis import parse_molden, fragment_weight, parse_atom_range


def parse_log_states(logfile, wanted_states):
    with open(logfile) as f:
        text = f.read()

    # Split the log on each state header; each chunk runs up to the next
    # such header (or EOF) and contains exactly that state's "n=" list
    # (plus unrelated text after it, which is harmless for a findall).
    headers = list(re.finditer(
        r"NATURAL TRANSITION ORBITALS FOR STATE\s+(\d+)\n-+\n", text))

    results = {}
    for i, m in enumerate(headers):
        state = int(m.group(1))
        if wanted_states and state not in wanted_states:
            continue
        start = m.end()
        end = headers[i + 1].start() if i + 1 < len(headers) else len(text)
        block = text[start:end]
        nvals = [float(x) for x in re.findall(r"n=\s*([\d.]+)", block)]
        if nvals:
            results[state] = nvals
    return results


def find_nocc(molist, n1, tol=2e-3):
    for idx in range(1, len(molist) + 1):
        occ = molist[idx - 1][0]
        if occ is not None and abs(occ - n1) < tol:
            return idx
    raise ValueError(f"could not locate rank-1 NTO pair (n={n1}) in Molden MOs")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("logfile")
    ap.add_argument("molden_prefix")
    ap.add_argument("--frag1", required=True)
    ap.add_argument("--frag2", required=True)
    ap.add_argument("--frag1-name", default="frag1")
    ap.add_argument("--frag2-name", default="frag2")
    ap.add_argument("--states", default="", help="comma-separated; default: all found in log")
    ap.add_argument("--max-ranks", type=int, default=4)
    ap.add_argument("--csv", default=None)
    args = ap.parse_args()

    wanted = {int(x) for x in args.states.split(",") if x} if args.states else None
    frag1 = parse_atom_range(args.frag1)
    frag2 = parse_atom_range(args.frag2)
    n1name, n2name = args.frag1_name, args.frag2_name

    log_states = parse_log_states(args.logfile, wanted)
    if not log_states:
        sys.exit("No matching 'NATURAL TRANSITION ORBITALS FOR STATE' blocks found.")

    csv_rows = []
    for state in sorted(log_states):
        nvals = log_states[state][: args.max_ranks]
        for cand in (f"{args.molden_prefix}_s{state}.molden.input",
                     f"{args.molden_prefix}s{state}.molden.input",
                     f"{args.molden_prefix}.s{state}.molden.input"):
            try:
                atoms, ao_atom, molist = parse_molden(cand)
                molden_path = cand
                break
            except FileNotFoundError:
                continue
        else:
            print(f"State {state}: no molden file found (tried "
                  f"{args.molden_prefix}_s{state}.molden.input etc.) -- skipping")
            continue

        nocc = find_nocc(molist, nvals[0])

        print(f"\n=== State {state}  (Molden: {molden_path}, hole/particle "
              f"boundary = {nocc}) ===")
        print(f"{'rank':>4} {'n':>8} {'hole':>6} {n1name+'%':>8} {n2name+'%':>7}"
              f"   {'part':>6} {n1name+'%':>8} {n2name+'%':>7}")
        cum = {"h1": 0.0, "h2": 0.0, "p1": 0.0, "p2": 0.0}
        total_n = sum(nvals)
        for k, n in enumerate(nvals, start=1):
            hole_idx, part_idx = nocc - k + 1, nocc + k
            h_occ, h_coef = molist[hole_idx - 1]
            p_occ, p_coef = molist[part_idx - 1]
            h1 = fragment_weight(h_coef, ao_atom, frag1) * 100
            h2 = fragment_weight(h_coef, ao_atom, frag2) * 100
            p1 = fragment_weight(p_coef, ao_atom, frag1) * 100
            p2 = fragment_weight(p_coef, ao_atom, frag2) * 100
            print(f"{k:>4} {n:>8.4f} {hole_idx:>6} {h1:>7.1f}% {h2:>6.1f}%"
                  f"   {part_idx:>6} {p1:>7.1f}% {p2:>6.1f}%")
            cum["h1"] += n * h1
            cum["h2"] += n * h2
            cum["p1"] += n * p1
            cum["p2"] += n * p2
            if k == 1:
                rank1 = dict(hole1=h1, hole2=h2, part1=p1, part2=p2)

        avg = {k: v / total_n for k, v in cum.items()}
        print(f"  n-weighted avg (top {len(nvals)} pairs, {total_n*100:.1f}% "
              f"of state): hole {avg['h1']:.1f}% {n1name} / {avg['h2']:.1f}% "
              f"{n2name};  particle {avg['p1']:.1f}% {n1name} / "
              f"{avg['p2']:.1f}% {n2name}")

        csv_rows.append({
            "state": state,
            "rank1_n": nvals[0],
            "rank1_hole_%s" % n1name: round(rank1["hole1"], 1),
            "rank1_hole_%s" % n2name: round(rank1["hole2"], 1),
            "rank1_particle_%s" % n1name: round(rank1["part1"], 1),
            "rank1_particle_%s" % n2name: round(rank1["part2"], 1),
            "weighted_avg_hole_%s" % n1name: round(avg["h1"], 1),
            "weighted_avg_hole_%s" % n2name: round(avg["h2"], 1),
            "weighted_avg_particle_%s" % n1name: round(avg["p1"], 1),
            "weighted_avg_particle_%s" % n2name: round(avg["p2"], 1),
        })

    if args.csv and csv_rows:
        with open(args.csv, "w", newline="") as f:
            w = csvmod.DictWriter(f, fieldnames=list(csv_rows[0].keys()))
            w.writeheader()
            w.writerows(csv_rows)
        print(f"\nWrote {args.csv}")


if __name__ == "__main__":
    main()
