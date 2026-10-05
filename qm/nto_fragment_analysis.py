#!/usr/bin/env python3
"""Fragment (donor/acceptor) localization analysis of ORCA Natural
Transition Orbitals (NTOs), from a Molden file written by `orca_2mkl
<base>.sN -molden` on an ORCA NTO file (`DoNTO true` TD-DFT job output,
"<base>.sN.nto").

For each requested MO/NTO index, reports what fraction of that orbital's
electron density (by raw |c_i|^2 AO-coefficient weight, summed per atom --
an unweighted-Mulliken-style measure; good enough for a qualitative
fragment split, not meant as a rigorous atomic-charge analysis) sits on
each of two user-defined fragments (e.g. a donor complex vs. an acceptor
cluster docked together).

Usage:
    python3 nto_fragment_analysis.py FILE.molden.input MO_INDEX [MO_INDEX ...] \
        --frag1 1-61 --frag2 62-77 [--frag1-name rubpy --frag2-name mos]

MO_INDEX is 1-based. Note: ORCA's own log labels NTO pairs using the
*original* canonical-orbital numbering (e.g. "255a -> 256a"), which is NOT
the same as the index of that orbital inside the separately-written
Molden/NTO file -- see nto_charge_transfer_summary.py, which handles that
mapping automatically and is the more convenient entry point for a full
per-state CT-character summary.
"""
import argparse
import re


SHELL_NAO = {"s": 1, "p": 3, "d": 5, "f": 7, "g": 9}  # spherical (5D/7F/9G)


def parse_atom_range(spec):
    atoms = set()
    for part in spec.split(","):
        if "-" in part:
            lo, hi = part.split("-")
            atoms.update(range(int(lo), int(hi) + 1))
        else:
            atoms.add(int(part))
    return atoms


def parse_molden(path):
    with open(path) as f:
        lines = f.readlines()

    i = next(k for k, l in enumerate(lines) if l.startswith("[Atoms]"))
    natom_atoms = []
    k = i + 1
    while k < len(lines) and not lines[k].startswith("["):
        parts = lines[k].split()
        if len(parts) >= 2:
            natom_atoms.append(parts[0])
        k += 1

    i = next(k for k, l in enumerate(lines) if l.startswith("[GTO]"))
    ao_atom = []
    k = i + 1
    atom_i = 0
    while k < len(lines) and not lines[k].startswith("["):
        line = lines[k].strip()
        if line == "":
            k += 1
            continue
        parts = line.split()
        if len(parts) == 2 and parts[1] == "0":
            atom_i = int(parts[0])
            k += 1
            continue
        m = re.match(r"^([spdfg])\s+(\d+)\s+([\d.]+)", line)
        if m:
            shell_type = m.group(1)
            nprim = int(m.group(2))
            nao = SHELL_NAO[shell_type]
            ao_atom.extend([atom_i] * nao)
            k += 1 + nprim
            continue
        k += 1

    nbasis = len(ao_atom)

    i = next(k for k, l in enumerate(lines) if l.startswith("[MO]"))
    mos = []
    k = i + 1
    while k < len(lines):
        if not lines[k].startswith(" Sym="):
            break
        occup = None
        k += 1
        while lines[k].strip().startswith(("Ene=", "Spin=", "Occup=")):
            if lines[k].strip().startswith("Occup="):
                occup = float(lines[k].split("=")[1])
            k += 1
        coeffs = [0.0] * nbasis
        for _ in range(nbasis):
            idx_str, val_str = lines[k].split()
            coeffs[int(idx_str) - 1] = float(val_str)
            k += 1
        mos.append((occup, coeffs))

    return natom_atoms, ao_atom, mos


def fragment_weight(coeffs, ao_atom, frag_atoms):
    num = sum(c * c for c, a in zip(coeffs, ao_atom) if a in frag_atoms)
    den = sum(c * c for c in coeffs)
    return num / den if den else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("molden_file")
    ap.add_argument("mo_indices", type=int, nargs="+")
    ap.add_argument("--frag1", required=True, help="e.g. 1-61")
    ap.add_argument("--frag2", required=True, help="e.g. 62-77")
    ap.add_argument("--frag1-name", default="frag1")
    ap.add_argument("--frag2-name", default="frag2")
    args = ap.parse_args()

    frag1 = parse_atom_range(args.frag1)
    frag2 = parse_atom_range(args.frag2)

    atoms, ao_atom, mos = parse_molden(args.molden_file)
    print(f"{args.molden_file}: {len(atoms)} atoms, {len(ao_atom)} basis "
          f"functions, {len(mos)} MOs")
    print(f"{'MO':>5} {'occup':>9} {args.frag1_name + ' %':>12} "
          f"{args.frag2_name + ' %':>10}")
    for idx in args.mo_indices:
        occup, coeffs = mos[idx - 1]
        w1 = fragment_weight(coeffs, ao_atom, frag1) * 100
        w2 = fragment_weight(coeffs, ao_atom, frag2) * 100
        print(f"{idx:>5} {occup:>9.4f} {w1:>11.1f}% {w2:>9.1f}%")


if __name__ == "__main__":
    main()
