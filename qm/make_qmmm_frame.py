#!/usr/bin/env python3
"""Extract one MD frame as a QM/MM (electrostatic embedding) snapshot:
writes the QM-region atoms as a plain .xyz and every other atom within a
cutoff shell as an ORCA point-charge file, using the real atomic charges
from the AMBER topology (not a guessed/rescaled value).

No link atoms / QM-MM bond-cutting support: this is for the case where
the QM region is one or more whole residues/molecules with nothing
covalently bonded across the boundary (e.g. a non-covalent host...guest
or donor...acceptor pair sitting in explicit solvent + ions) -- see
05_qmmm_ensemble/qmmm_ensemble_workflow.pdf for why that applies here.

Requires: AmberTools (cpptraj, for the frame extraction/autoimage step)
and ParmEd (for real-unit atomic charges) on PATH/PYTHONPATH -- e.g.
`module load amber`.

Usage:
    python3 make_qmmm_frame.py PRMTOP TRAJ FRAME \
        --qm-mask "^8,12" --shell 18.0 \
        --out-prefix frame146_cat8_rub12

Use "^" (molecule number) for --qm-mask when the QM region is one or more
whole AMBER-sense "molecules" with internal multi-residue topology (e.g. a
metal cluster split into several residues per molecule) -- a residue mask
(":") would only grab one sub-residue, not the whole molecule. Check with
`cpptraj -p PRMTOP --resmask ':*'` (Mnum column) which molecule number
each piece belongs to before picking a mask.

Writes <out-prefix>_qm.xyz and <out-prefix>_mm.pc (ORCA point-charge
format: first line = count, then "charge x y z" per atom, Angstrom).
"""
import argparse
import re
import subprocess
import tempfile
import os
import parmed


def parse_mask_tokens(mask):
    """Split a single-prefix cpptraj mask ("^8,12" or "^8,^12" or
    ":8,16") into (prefix, [token, ...]), so each token can be
    re-combined with the SAME prefix later (needed because cpptraj's
    distance operator does not reliably combine with a raw comma list --
    see run_cpptraj). Accepts the prefix repeated per-token or only once.
    """
    m = re.match(r"^([:^@])", mask)
    if not m:
        raise ValueError(f"mask {mask!r} must start with one of : ^ @")
    prefix = m.group(1)
    body = mask[1:]
    tokens = [t[1:] if t.startswith(prefix) else t for t in body.split(",")]
    return prefix, tokens


def _run(script):
    with tempfile.NamedTemporaryFile("w", suffix=".in", delete=False) as f:
        f.write(script)
        script_path = f.name
    result = subprocess.run(["cpptraj", "-i", script_path],
                             capture_output=True, text=True)
    os.unlink(script_path)
    return result


def run_cpptraj(prmtop, traj, frame, qm_mask, shell, out_pdb, mask_list_out,
                 qm_list_out, anchor_mask, fixed_mask):
    fixed_clause = f"fixed {fixed_mask}" if fixed_mask else ""
    # cpptraj's distance operator does not combine correctly with a
    # comma-separated mask list (":8,16<:18.0" silently selects almost
    # nothing) -- OR together one distance clause per token instead,
    # which is the syntax that actually works.
    prefix, tokens = parse_mask_tokens(qm_mask)
    sel = "|".join(f"{prefix}{t}<:{shell}" for t in tokens)

    # Pass 1: autoimage + get the authoritative kept-atom list from `mask`
    # (confirmed correct). `strip` with a NEGATED distance-based OR mask is
    # separately confirmed BUGGY in this cpptraj version (silently strips
    # almost nothing) -- so don't ask strip to re-evaluate the same
    # expression; use its explicit atom-number output in pass 2 instead.
    # Also list the plain QM-mask atoms here via cpptraj itself -- ParmEd's
    # own mask parser does not support "^" (molecule) masks at all, so it
    # cannot be used for this lookup when qm_mask uses that prefix.
    script1 = f"""parm {prmtop}
trajin {traj} {frame} {frame}
autoimage anchor {anchor_mask} {fixed_clause} origin
mask {sel} out {mask_list_out}
mask {qm_mask} out {qm_list_out}
trajout {out_pdb}.full.rst7 restart
run
"""
    result = _run(script1)
    if result.returncode != 0 or not os.path.exists(mask_list_out):
        print(result.stdout)
        print(result.stderr)
        raise RuntimeError("cpptraj pass 1 (autoimage/mask listing) failed")

    atom_nums = []
    with open(mask_list_out) as f:
        next(f)
        for line in f:
            parts = line.split()
            if parts:
                atom_nums.append(parts[1])
    atom_list = ",".join(atom_nums)

    # Pass 2: strip to the EXPLICIT atom-number list from pass 1 (an @-mask
    # with a plain comma list, not a distance expression -- the form that
    # works), starting from the already-autoimaged restart coordinates.
    full_rst = f"{out_pdb}.full.rst7"
    script2 = f"""parm {prmtop}
trajin {full_rst}
strip !(@{atom_list})
trajout {out_pdb} pdb
run
"""
    result = _run(script2)
    os.unlink(full_rst)
    if result.returncode != 0 or not os.path.exists(out_pdb):
        print(result.stdout)
        print(result.stderr)
        raise RuntimeError("cpptraj pass 2 (strip to explicit atom list) failed")
    return result.stdout


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prmtop")
    ap.add_argument("traj")
    ap.add_argument("frame", type=int, help="1-based frame index in TRAJ")
    ap.add_argument("--qm-mask", required=True,
                     help='cpptraj atom/residue mask for the QM region, e.g. ":8,16"')
    ap.add_argument("--anchor-mask", default=None,
                     help='autoimage anchor (single molecule), default: first residue in --qm-mask')
    ap.add_argument("--fixed-mask", default=None,
                     help='autoimage "fixed" mask (imaged close to anchor), default: rest of --qm-mask')
    ap.add_argument("--shell", type=float, default=18.0,
                     help="MM shell radius around the QM region, Angstrom (default 18.0)")
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    qm_prefix, qm_tokens = parse_mask_tokens(args.qm_mask)
    anchor_mask = args.anchor_mask or f"{qm_prefix}{qm_tokens[0]}"
    fixed_mask = args.fixed_mask
    if fixed_mask is None and len(qm_tokens) > 1:
        fixed_mask = qm_prefix + ",".join(qm_tokens[1:])

    tmp_pdb = f"{args.out_prefix}_kept.pdb"
    mask_list_out = f"{args.out_prefix}_kept_atoms.dat"
    qm_list_out = f"{args.out_prefix}_qm_atoms.dat"
    run_cpptraj(args.prmtop, args.traj, args.frame, args.qm_mask, args.shell, tmp_pdb,
                mask_list_out, qm_list_out, anchor_mask, fixed_mask)

    full = parmed.load_file(args.prmtop)

    # qm_idx: straight from cpptraj's own "mask" listing for the plain
    # qm_mask (ParmEd's mask parser does not support "^" molecule masks at
    # all, so it cannot be used here when qm_mask uses that prefix).
    qm_idx = set()
    with open(qm_list_out) as f:
        next(f)
        for line in f:
            parts = line.split()
            if parts:
                qm_idx.add(int(parts[1]) - 1)

    # Original (1-based -> convert to 0-based) atom indices that survived
    # the strip, straight from cpptraj's own "mask ... out" listing -- this
    # is the SAME imaged frame the stripped PDB's coordinates came from, so
    # the two are guaranteed consistent (no re-deriving the selection).
    kept_idx = []
    with open(mask_list_out) as f:
        next(f)  # header
        for line in f:
            parts = line.split()
            if not parts:
                continue
            kept_idx.append(int(parts[1]) - 1)

    kept = parmed.load_file(tmp_pdb)
    if len(kept.atoms) == 0:
        raise RuntimeError("No atoms survived the cpptraj strip -- check --qm-mask/--shell")
    if len(kept_idx) != len(kept.atoms):
        raise RuntimeError(
            f"Atom count mismatch: stripped PDB has {len(kept.atoms)} atoms, "
            f"mask listing has {len(kept_idx)}. Inspect {tmp_pdb} and "
            f"{mask_list_out} directly.")

    qm_lines = []
    mm_lines = []
    for i, orig_idx in enumerate(kept_idx):
        atom = kept.atoms[i]
        x, y, z = atom.xx, atom.xy, atom.xz
        elem = full.atoms[orig_idx].element_name
        charge = full.atoms[orig_idx].charge
        if orig_idx in qm_idx:
            qm_lines.append(f"{elem:<3s}{x:14.6f}{y:14.6f}{z:14.6f}")
        else:
            mm_lines.append(f"{charge:10.6f}{x:14.6f}{y:14.6f}{z:14.6f}")

    qm_path = f"{args.out_prefix}_qm.xyz"
    with open(qm_path, "w") as f:
        f.write(f"{len(qm_lines)}\n")
        f.write(f"QM region, frame {args.frame}, mask {args.qm_mask}\n")
        f.write("\n".join(qm_lines) + "\n")

    mm_path = f"{args.out_prefix}_mm.pc"
    with open(mm_path, "w") as f:
        f.write(f"{len(mm_lines)}\n")
        f.write("\n".join(mm_lines) + "\n")

    net_mm_charge = sum(full.atoms[i].charge for i in kept_idx if i not in qm_idx)
    print(f"QM region: {len(qm_lines)} atoms -> {qm_path}")
    print(f"MM shell:  {len(mm_lines)} atoms (net charge {net_mm_charge:+.3f} e) -> {mm_path}")
    os.unlink(tmp_pdb)
    os.unlink(mask_list_out)
    os.unlink(qm_list_out)


if __name__ == "__main__":
    main()
