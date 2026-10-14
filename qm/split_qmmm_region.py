#!/usr/bin/env python3
"""Rebuild cat-alone/ps-alone QM/MM files from an ALREADY-EXTRACTED complex
region (frame_qm.xyz + frame_mm.pc), by reclassifying a subset of the
complex's own QM atoms as additional MM point charges -- rather than
re-deriving a brand-new independent MM shell for each sub-region.

Why this matters: calling make_qmmm_frame.py separately for the cat-alone
and ps-alone regions (each re-deriving its own 18 A shell centered on a
DIFFERENT reference point) gives each sub-region a DIFFERENT, inconsistent
realization of the local counter-ion cloud -- confirmed 2026-10-14 on
1_PS_1_CAT/0_vol/buffer/exp frame 1: cat_mm.pc net charge +2.0 (as
expected, just PS's +2) but ps_mm.pc net charge -3.0 (CAT's -2 PLUS a
stray -1 counter-ion picked up only by the PS-centered shell, not the
CAT-centered one). That inconsistency corrupts dE_int = E(complex) -
E(cat) - E(ps) with spurious +100-200 kcal/mol noise, completely masking
the real ~1-5 kcal/mol CAT<->PS interaction. This script instead takes the
ONE already-extracted complex MM shell and reuses it unchanged for every
region -- only the QM/MM boundary moves -- which is what the subtraction
actually requires to isolate the real interaction energy.

Requires: cpptraj (for the atom-index list, same technique as
make_qmmm_frame.py) and ParmEd (for real-unit charges) on PATH/PYTHONPATH.

Usage:
    python3 split_qmmm_region.py PRMTOP FRAMEDIR \\
        --qm-mask ":1-11" --n-atoms-a 16 --tag-a cat --tag-b ps
(FRAMEDIR must already contain frame_qm.xyz + frame_mm.pc, the complex
region's output from make_qmmm_frame.py / setup_qmmm_ensemble.py.
--n-atoms-a: how many of the LEADING atoms in frame_qm.xyz -- in --qm-mask
order -- belong to region A; the rest become region B.)
"""
import argparse
import os
import subprocess
import tempfile

import parmed


def load_xyz(path):
    with open(path) as f:
        n = int(f.readline())
        f.readline()
        atoms = []
        for _ in range(n):
            parts = f.readline().split()
            atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    return atoms


def mask_atom_indices(prmtop, traj, mask):
    """0-based original prmtop atom indices selected by `mask`, in cpptraj's
    own (ascending) order -- same technique as make_qmmm_frame.py's
    qm_list_out, deliberately not ParmEd's own mask parser. Needs an actual
    trajin (any valid coordinate source, e.g. the frame's own
    snapshot.ncrst) -- `mask ... out` is a per-frame action and produces no
    output at all if run against a bare topology with no frames."""
    with tempfile.TemporaryDirectory() as td:
        listing = os.path.join(td, "atoms.dat")
        script = os.path.join(td, "mask.in")
        with open(script, "w") as f:
            f.write(f"parm {prmtop}\ntrajin {traj}\nmask {mask} out {listing}\n")
        result = subprocess.run(["cpptraj", "-i", script], capture_output=True, text=True)
        if result.returncode != 0 or not os.path.exists(listing):
            raise RuntimeError(f"cpptraj mask selection failed for {mask!r}:\n"
                                f"{result.stdout}\n{result.stderr}")
        idx = []
        with open(listing) as f:
            next(f)
            for line in f:
                parts = line.split()
                if parts:
                    idx.append(int(parts[1]) - 1)
    return idx


def write_region(tag, atoms, charges, qm_range, mm_extra_range, mm_lines, n_mm, out_dir):
    qm_path = os.path.join(out_dir, f"{tag}_qm.xyz")
    sel_atoms = [atoms[i] for i in qm_range]
    with open(qm_path, "w") as f:
        f.write(f"{len(sel_atoms)}\n")
        f.write(f"QM region (reclassified from complex), tag={tag}\n")
        for el, x, y, z in sel_atoms:
            f.write(f"{el:<3s}{x:14.6f}{y:14.6f}{z:14.6f}\n")

    extra_lines = []
    net_extra = 0.0
    for i in mm_extra_range:
        el, x, y, z = atoms[i]
        extra_lines.append(f"{charges[i]:10.6f}{x:14.6f}{y:14.6f}{z:14.6f}\n")
        net_extra += charges[i]

    mm_path = os.path.join(out_dir, f"{tag}_mm.pc")
    with open(mm_path, "w") as f:
        f.write(f"{n_mm + len(extra_lines)}\n")
        f.writelines(mm_lines)
        f.writelines(extra_lines)
    return qm_path, mm_path, net_extra


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prmtop")
    ap.add_argument("framedir")
    ap.add_argument("--qm-mask", required=True,
                     help="combined mask matching frame_qm.xyz's atom order, e.g. ':1-11'")
    ap.add_argument("--n-atoms-a", type=int, required=True,
                     help="number of LEADING atoms (in --qm-mask order) belonging to region A")
    ap.add_argument("--tag-a", required=True)
    ap.add_argument("--tag-b", required=True)
    args = ap.parse_args()

    complex_qm = os.path.join(args.framedir, "frame_qm.xyz")
    complex_mm = os.path.join(args.framedir, "frame_mm.pc")
    atoms = load_xyz(complex_qm)

    snap = os.path.join(args.framedir, "snapshot.ncrst")
    idx = mask_atom_indices(args.prmtop, snap, args.qm_mask)
    if len(idx) != len(atoms):
        raise RuntimeError(f"{args.qm_mask} selects {len(idx)} atoms in {args.prmtop}, "
                            f"but {complex_qm} has {len(atoms)} -- mask/order mismatch")

    full = parmed.load_file(args.prmtop)
    charges = [full.atoms[i].charge for i in idx]

    with open(complex_mm) as f:
        n_mm = int(f.readline())
        mm_lines = [f.readline() for _ in range(n_mm)]

    n = len(atoms)
    na = args.n_atoms_a
    _, _, net_b_in_a = write_region(args.tag_a, atoms, charges, range(0, na), range(na, n),
                                     mm_lines, n_mm, args.framedir)
    _, _, net_a_in_b = write_region(args.tag_b, atoms, charges, range(na, n), range(0, na),
                                     mm_lines, n_mm, args.framedir)
    print(f"wrote {args.tag_a}_qm.xyz/{args.tag_a}_mm.pc (MM += {args.tag_b}, net {net_b_in_a:+.3f} e) "
          f"and {args.tag_b}_qm.xyz/{args.tag_b}_mm.pc (MM += {args.tag_a}, net {net_a_in_b:+.3f} e) "
          f"in {args.framedir}")


if __name__ == "__main__":
    main()
