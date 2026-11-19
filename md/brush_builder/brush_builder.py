#!/usr/bin/env python3
"""brush_builder: all-atom models of randomly grafted, atactic polyelectrolyte brushes for Amber, with periodic
boundary conditions imposed during the packing.

Method of S. Tippner, D. Hernandez-Castillo, F. H. Schacher, L. Gonzalez, J. Phys. Chem. B 129, 2105 (2025)
(PDMAEMA brushes); this is a rewrite of the original collection of bash and Python scripts as one tool.

Steps (each a subcommand; "all" runs 1 to 4, "single" is optional):
  1 params   residues for every monomer type (head, middle and tail unit, from a RESP-charged antechamber file) with
             prepgen and GAFF2 parameters with parmchk2
  2 chains   n_chains random sequences with a fixed number of protonated units and random stereochemistry, built
             with tleap and aligned along z (longest axis onto z, grafted unit at the bottom, lowest atom at z = 0)
  3 pack     (a) smallest square box in which PACKMOL places all chains upright at the given tolerance (the minimum
             distance between chains, which sets the grafting density), widened in steps of box_step;
             (b) packing with periodic boundaries in a region of 3 x 3 cells: every chain is placed with its grafted
             atom on the plane z = anchor_z inside the central cell, against all chains placed before and their
             copies, and is then copied into the 8 neighbouring cells (copies kept if any atom lies within half a cell
             of the central cell). After max_attempts failures for one chain the packing restarts with chain 1;
             after max_fail failed packings the region is enlarged by region_step and packing starts again.
  4 solvate  the whole region, copies included, is filled with water (shell_xy beyond the polymer laterally,
             shell_z above it) and counterions are added to the central cell; then the copies and all water outside
             the central cell are removed and tleap writes the Amber topology (GAFF2, TIP3P) of the periodic box.
  5 single   (optional) every chain solvated on its own with counterions, as reference systems

Requirements: Python >= 3.11 with numpy; AmberTools (prepgen, parmchk2, tleap) on PATH; PACKMOL (path in the
config or on PATH).  Usage:  brush_builder.py <config.toml> {params,chains,pack,solvate,all,single}
"""
import argparse
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path

import numpy as np

AMU_PER_A3 = 0.602214            # g/cm3 -> amu/A3
M_WATER = 18.015
WATER_PDB = ("HETATM    1  O   WAT     1       0.000   0.000   0.000\n"
             "HETATM    2  H1  WAT     1       0.957   0.000   0.000\n"
             "HETATM    3  H2  WAT     1      -0.240   0.927   0.000\nEND\n")
MASS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "S": 32.06, "Cl": 35.45, "Na": 22.99}


# ------------------------------------------------------------------ helpers
def log(msg):
    print(f"[brush_builder] {msg}", flush=True)


def run(cmd, cwd, stdin=None):
    """run a program in cwd, return stdout; raise with its output on failure"""
    r = subprocess.run(cmd, cwd=cwd, input=stdin, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"{' '.join(map(str, cmd))} failed in {cwd}:\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
    return r.stdout


def packmol(cfg, inp, cwd):
    """run PACKMOL with input text; True if it reports success"""
    (cwd / "packmol.inp").write_text(inp)
    with open(cwd / "packmol.inp") as f:                        # as packmol < packmol.inp (a pipe is not read)
        r = subprocess.run([cfg["tools"]["packmol"]], cwd=cwd, stdin=f, capture_output=True, text=True)
    (cwd / "packmol.log").write_text(r.stdout + r.stderr)
    return "Success!" in r.stdout


def read_pdb(path):
    """atoms as list of dicts (record, name, resname, chain, resid, xyz) and coordinates array"""
    atoms = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            atoms.append(dict(rec=line[:6], name=line[12:16], res=line[17:20] + " ", chain=line[21],
                              resid=line[22:26], tail=line[54:]))
    X = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                  for l in Path(path).read_text().splitlines() if l.startswith(("ATOM", "HETATM"))])
    return atoms, X


def write_pdb(path, blocks, cryst=None):
    """blocks: list of (atoms, X) written in order, TER after each block; atom serials wrap at 99999"""
    out = [f"CRYST1{cryst[0]:9.3f}{cryst[1]:9.3f}{cryst[2]:9.3f}  90.00  90.00  90.00 P 1           1"] if cryst else []
    n = 0
    for atoms, X in blocks:
        for a, x in zip(atoms, X):
            n += 1
            out.append(f"{a['rec']}{n % 100000:5d} {a['name']} {a['res']}{a['chain']}{a['resid']}    "
                       f"{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}{a.get('tail', '')}".rstrip())
        out.append("TER")
    out.append("END")
    Path(path).write_text("\n".join(out) + "\n")


def element(name):
    n = name.strip()
    return "Cl" if n.upper().startswith("CL") else n[0]


def mass(atoms):
    return sum(MASS.get(element(a["name"]), 12.0) for a in atoms)


def load_config(path):
    cfg = tomllib.loads(Path(path).read_text())
    cfg["_dir"] = Path(path).resolve().parent
    cfg["_out"] = (cfg["_dir"] / cfg.get("output", "output")).resolve()
    cfg.setdefault("tools", {}).setdefault("packmol", shutil.which("packmol") or "packmol")
    return cfg


# ------------------------------------------------------------------ 1 params
def residue_names(key):
    """monomer key (e.g. 'PA') -> head, middle, tail residue names"""
    return {"head": "H" + key, "middle": "D" + key, "tail": "T" + key}


def step_params(cfg):
    out = cfg["_out"] / "params"; out.mkdir(parents=True, exist_ok=True)
    r = cfg["residues"]
    omit = {"head": [r["cap_atom"]],                                # keeps the end group (stand-in for the wall)
            "middle": r["end_group_atoms"] + [r["cap_atom"]],
            "tail": r["end_group_atoms"]}                            # keeps the H cap
    for key, m in cfg["monomers"].items():
        ac = (cfg["_dir"] / m["ac"]).resolve()
        charge = float(re.search(r"CHARGE\s+([-\d.]+)", ac.read_text()).group(1))
        for pos, rn in residue_names(key).items():
            mc = out / f"mainchain.{rn.lower()}"
            mc.write_text(f"HEAD_NAME {r['head_name']}\nTAIL_NAME {r['tail_name']}\n"
                          + "".join(f"OMIT_NAME {a}\n" for a in omit[pos])
                          + f"PRE_HEAD_TYPE {r['pre_head_type']}\nPOST_TAIL_TYPE {r['post_tail_type']}\n"
                          + f"CHARGE {charge:.1f}\n")
            run(["prepgen", "-i", str(ac), "-o", f"{rn.lower()}.prepi", "-f", "prepi", "-m", mc.name,
                 "-rn", rn, "-rf", f"{rn.lower()}.res"], out)
            run(["parmchk2", "-i", f"{rn.lower()}.prepi", "-f", "prepi", "-o", f"{rn.lower()}.frcmod",
                 "-a", "Y", "-s", "gaff2"], out)
    log(f"params: {3 * len(cfg['monomers'])} residues in {out}")


def leap_header(cfg, water=False):
    p = cfg["_out"] / "params"
    lines = ["source leaprc.gaff2"] + (["source leaprc.water.tip3p"] if water else [])
    for key in cfg["monomers"]:
        for rn in residue_names(key).values():
            lines += [f"loadamberprep {p / (rn.lower() + '.prepi')}", f"loadamberparams {p / (rn.lower() + '.frcmod')}"]
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------ 2 chains
def sequence(cfg, rng):
    """random sequence: n_units units, round(alpha * n_units) of them protonated, stereo and positions random"""
    c = cfg["chains"]; n = c["n_units"]; n_p = round(c["alpha"] * n)
    prot = [k for k, m in cfg["monomers"].items() if m["protonated"]]
    neut = [k for k, m in cfg["monomers"].items() if not m["protonated"]]
    is_p = [True] * n_p + [False] * (n - n_p)
    rng.shuffle(is_p)
    keys = [rng.choice(prot if p else neut) for p in is_p]
    names = [("H" if i == 0 else "T" if i == n - 1 else "D") + k for i, k in enumerate(keys)]
    return names


def step_chains(cfg):
    out = cfg["_out"] / "chains"; out.mkdir(parents=True, exist_ok=True)
    c = cfg["chains"]; b = cfg["brush"]
    rng = random.Random(c.get("seed"))
    seqs = []
    for i in range(1, c["n_chains"] + 1):
        d = out / f"chain_{i}"; d.mkdir(exist_ok=True)
        seq = sequence(cfg, rng); seqs.append(seq)
        (d / "tleap.in").write_text(leap_header(cfg) + f"mol = sequence {{ {' '.join(seq)} }}\n"
                                    f"savepdb mol built.pdb\nquit\n")
        run(["tleap", "-f", "tleap.in"], d)
        # align along z: PACKMOL moves molecules only as rigid bodies, so fitting the chain into a narrow, tall box
        # (as in the original scripts) is a rotation of its longest axis onto z; done here directly
        atoms, X = read_pdb(d / "built.pdb")
        C = X - X.mean(0)
        w, v = np.linalg.eigh(np.cov(C.T))
        Rm = v[:, ::-1][:, [1, 2, 0]]                               # longest axis -> z
        if np.linalg.det(Rm) < 0:
            Rm[:, 0] *= -1
        X = C @ Rm
        if X[b["anchor_atom"] - 1, 2] > X[:, 2].mean():             # grafted atom at the bottom
            X[:, 1:] *= -1
        X[:, 2] -= X[:, 2].min()                                    # lowest atom on z = 0
        write_pdb(out / f"chain_{i}.pdb", [(atoms, X)])
        w = np.ptp(X[:, :2], 0).max()
        log(f"chain {i}: {' '.join(seq)} (aligned, {np.ptp(X[:, 2]):.0f} A long, {w:.0f} A wide)")
    (out / "sequences.json").write_text(json.dumps(seqs, indent=1))


# ------------------------------------------------------------------ 3 pack
def chain_files(cfg):
    return [cfg["_out"] / "chains" / f"chain_{i}.pdb" for i in range(1, cfg["chains"]["n_chains"] + 1)]


def upright(path, box, anchor, z_a, height):
    """PACKMOL structure block: chain upright (rotations about x and y fixed), grafted atom on the plane z_a"""
    x0, y0, x1, y1 = box
    return (f"structure {path}\n  number 1\n  inside box {anchor[0]:.3f} {anchor[1]:.3f} 0. {anchor[2]:.3f} "
            f"{anchor[3]:.3f} {height:.3f}\n  constrain_rotation x 0. 0.\n  constrain_rotation y 0. 0.\n"
            f"  atoms {{atom}}\n    inside box {x0:.3f} {y0:.3f} {z_a:.3f} {x1:.3f} {y1:.3f} {z_a:.3f}\n"
            f"  end atoms\nend structure\n")


def smallest_box(cfg, files, d):
    """(a) smallest square box holding all chains upright at the tolerance, widened in steps of box_step"""
    b = cfg["brush"]; tol = b["tolerance"]
    height = max(read_pdb(f)[1][:, 2].max() for f in files) + 4.0
    L = b.get("box_start", math.sqrt(len(files)) * (tol + 4.0))
    while True:
        body = "".join(upright(f, (0, 0, L, L), (0, 0, L, L), b["anchor_z"], height).replace("{atom}", str(b["anchor_atom"]))
                       for f in files)
        if packmol(cfg, f"tolerance {tol}\nfiletype pdb\noutput box.pdb\nnloop 20\nmaxit 20\nadd_amber_ter\n{body}", d):
            log(f"pack (a): all chains fit in {L:.2f} x {L:.2f} A at tolerance {tol} A")
            return L, height
        L += b["box_step"]


def copies(atoms, X, L, R):
    """copies of a chain in the 8 neighbouring cells, kept if any atom lies within half a cell of the central cell"""
    out = []
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            if i == j == 0:
                continue
            Y = X + [i * L, j * L, 0]
            inside = (Y[:, 0] > 0.5 * L) & (Y[:, 0] < 2.5 * L) & (Y[:, 1] > 0.5 * L) & (Y[:, 1] < 2.5 * L)
            if inside.any():
                out.append(Y)
    return out


def step_pack(cfg):
    out = cfg["_out"] / "pack"; out.mkdir(parents=True, exist_ok=True)
    b = cfg["brush"]; files = chain_files(cfg)
    L0, height = smallest_box(cfg, files, out)
    R = math.ceil(3 * L0) if "region_start" not in b else b["region_start"]
    while True:                                            # (b) periodic packing, region enlarged after max_fail
        L = R / 3; fails = 0; placed = None
        log(f"pack (b): region {R:.1f} A, central cell {L:.2f} A")
        while fails < b["max_fail"]:
            central, copies_all, i, attempts = [], [], 0, 0
            while i < len(files):
                d = out / "work"; d.mkdir(exist_ok=True)
                fixed = ""
                if central:
                    write_pdb(d / "fixed.pdb", [(a, X) for a, X in central] + [(a, Y) for a, Y, _ in copies_all])
                    fixed = "structure fixed.pdb\n  number 1\n  fixed 0. 0. 0. 0. 0. 0.\nend structure\n"
                new = upright(files[i], (L, L, 2 * L, 2 * L), (0, 0, R, R), b["anchor_z"], height)
                ok = packmol(cfg, f"tolerance {b['tolerance']}\nfiletype pdb\noutput placed.pdb\nnloop 20\nmaxit 20\n"
                                  f"add_amber_ter\nseed -1\n{fixed}{new.replace('{atom}', str(b['anchor_atom']))}", d)
                if not ok:
                    fails += 1; attempts += 1
                    if fails >= b["max_fail"]:
                        break
                    if attempts >= b["max_attempts"]:        # restart from the first chain
                        central, copies_all, i, attempts = [], [], 0, 0
                    continue
                atoms_new, _ = read_pdb(files[i]); _, Xall = read_pdb(d / "placed.pdb")
                Xn = Xall[-len(atoms_new):]                 # the new chain is the last molecule of the output
                central.append((atoms_new, Xn))
                copies_all += [(atoms_new, Y, i) for Y in copies(atoms_new, Xn, L, R)]
                log(f"  chain {i + 1}/{len(files)} placed ({len(copies_all)} copies so far, {fails} failed packings)")
                i += 1; attempts = 0
            if i == len(files):
                placed = (central, copies_all); break
        if placed:
            break
        R += b["region_step"]
    central, copies_all = placed
    write_pdb(out / "region.pdb", [(a, X) for a, X in central] + [(a, Y) for a, Y, _ in copies_all])
    write_pdb(out / "central.pdb", central)
    (out / "pack.json").write_text(json.dumps({"region": R, "cell": R / 3, "smallest_box": L0, "height": height,
                                               "n_central": len(central), "n_copies": len(copies_all),
                                               "block_chain": list(range(len(central))) + [c for _, _, c in copies_all],
                                               "grafting_density_nm2": len(central) / (R / 3 / 10) ** 2}, indent=1))
    log(f"pack: {len(central)} chains and {len(copies_all)} copies, cell {R / 3:.2f} A, "
        f"grafting density {len(central) / (R / 30) ** 2:.2f} nm^-2")


# ------------------------------------------------------------------ 4 solvate
def ion_pdb(name, d):
    p = d / f"{name.replace('-', 'm').replace('+', 'p')}.pdb"
    p.write_text(f"HETATM    1 {name:<4}{name:>4}     1       0.000   0.000   0.000\nEND\n")
    return p.name


def solvate_box(cfg, d, solute_blocks, region_xy, n_ion, ion_box, top_z):
    """fill region_xy x [0, top_z] with water around the fixed solute and n_ion counterions inside ion_box.
    The first water count is mass based (box volume at the target density minus the solute mass, as solvate.tcl of
    PACKMOL); as the solute excludes more volume than its mass suggests, the count is then corrected until the density
    of the free water layer above the solute is within 3% of the target (at most 3 passes). Returns the PACKMOL
    output path and the number of water molecules."""
    s = cfg["solvent"]
    write_pdb(d / "solute.pdb", solute_blocks)
    (d / "water.pdb").write_text(WATER_PDB)
    ion = ion_pdb(s["counterion"], d)
    x0, y0, x1, y1 = region_xy
    area = (x1 - x0) * (y1 - y0)
    z_solute = max(x[:, 2].max() for _, x in solute_blocks)
    m_solute = sum(mass(a) for a, _ in solute_blocks) + n_ion * MASS.get(element(s["counterion"]), 35.45)
    n_w = max(1, round((s["density"] * AMU_PER_A3 * area * top_z - m_solute) / M_WATER))
    for it in range(3):
        inp = (f"tolerance 2.0\nfiletype pdb\noutput solvated.pdb\nadd_amber_ter\nseed -1\n"
               f"structure solute.pdb\n  number 1\n  fixed 0. 0. 0. 0. 0. 0.\nend structure\n"
               f"structure water.pdb\n  number {n_w}\n  inside box {x0:.3f} {y0:.3f} 0. {x1:.3f} {y1:.3f} {top_z:.3f}\n"
               f"end structure\n")
        if n_ion:
            inp += (f"structure {ion}\n  number {n_ion}\n  inside box {ion_box[0]:.3f} {ion_box[1]:.3f} 1. "
                    f"{ion_box[2]:.3f} {ion_box[3]:.3f} {top_z - 1:.3f}\nend structure\n")
        if not packmol(cfg, inp, d):
            raise RuntimeError(f"PACKMOL could not solvate the system in {d} (see packmol.log)")
        atoms, X = read_pdb(d / "solvated.pdb")
        o = np.array([a["res"].strip() == "WAT" and a["name"].strip() == "O" for a in atoms])
        layer = o & (X[:, 2] > z_solute + 3) & (X[:, 2] < top_z - 3)
        rho = layer.sum() * M_WATER / AMU_PER_A3 / (area * max(top_z - z_solute - 6, 1.0))
        log(f"solvation pass {it + 1}: {n_w} water molecules, {rho:.3f} g/cm3 in the water layer above the solute")
        if abs(rho / s["density"] - 1) < 0.03 or top_z - z_solute < 12:
            break
        n_w = round(n_w * s["density"] / rho)
    return d / "solvated.pdb", n_w


def remove_clashes(blocks, box, cut=2.0):
    """drop water molecules (blocks of residue WAT) with a heavy atom closer than cut to a heavy atom of another
    molecule under periodic boundaries; chains and ions are kept"""
    from scipy.spatial import cKDTree
    X = np.vstack([x for _, x in blocks]); owner = np.concatenate([[k] * len(x) for k, (_, x) in enumerate(blocks)])
    names = [a["name"].strip() for a, _ in blocks for a in a]
    heavy = np.array([not n.startswith("H") for n in names])
    is_w = np.array([a[0]["res"].strip() == "WAT" for a, _ in blocks])
    idx = np.where(heavy)[0]
    t = cKDTree(X[idx] % box, boxsize=box)
    drop = set()
    for i, j in t.query_pairs(cut, output_type="ndarray"):
        bi, bj = owner[idx[i]], owner[idx[j]]
        if bi == bj or bi in drop or bj in drop:
            continue
        if is_w[bj]:
            drop.add(bj)
        elif is_w[bi]:
            drop.add(bi)
    return [b for k, b in enumerate(blocks) if k not in drop], len(drop)


def tleap_box(cfg, d, pdb, box):
    (d / "tleap.in").write_text(leap_header(cfg, water=True) + f"mol = loadpdb {pdb}\n"
                                f"set mol box {{ {box[0]:.3f} {box[1]:.3f} {box[2]:.3f} }}\n"
                                "charge mol\nsavepdb mol system.pdb\nsaveamberparm mol system.prmtop system.inpcrd\nquit\n")
    log_text = run(["tleap", "-f", "tleap.in"], d)
    (d / "leap.out").write_text(log_text)
    q = re.findall(r"Total unperturbed charge:\s*([-\d.]+)", log_text)
    log(f"tleap: system.prmtop written, total charge {q[-1] if q else '?'}")


def charge_of(cfg, atoms_blocks):
    """number of protonated units (net charge) of the chains"""
    prot = {rn for k, m in cfg["monomers"].items() if m["protonated"] for rn in residue_names(k).values()}
    return sum(len({a["resid"] for a in atoms if a["res"].strip() in prot}) for atoms, _ in atoms_blocks)


def step_solvate(cfg):
    out = cfg["_out"] / "solvate"; out.mkdir(parents=True, exist_ok=True)
    s = cfg["solvent"]; P = json.loads((cfg["_out"] / "pack" / "pack.json").read_text())
    L = P["cell"]
    atoms, X = read_pdb(cfg["_out"] / "pack" / "region.pdb")
    sizes = [len(read_pdb(f)[0]) for f in chain_files(cfg)]
    n_c, all_blocks, k = P["n_central"], [], 0
    for c in P["block_chain"]:                                  # central chains first, then the copies
        n = sizes[c]; all_blocks.append((atoms[k:k + n], X[k:k + n])); k += n
    blocks, copy_blocks = all_blocks[:n_c], all_blocks[n_c:]
    q = charge_of(cfg, blocks)
    allX = np.vstack([x for _, x in blocks + copy_blocks])
    lo = allX[:, :2].min(0) - s["shell_xy"]; hi = allX[:, :2].max(0) + s["shell_xy"]
    top = allX[:, 2].max() + s["shell_z"]
    path, n_w = solvate_box(cfg, out, blocks + copy_blocks, (lo[0], lo[1], hi[0], hi[1]), q, (L, L, 2 * L, 2 * L), top)
    # cut: central chains, all ions (they are in the central cell), water molecules whose O lies in the central cell
    satoms, SX = read_pdb(path)
    n_sol = sum(len(a) for a, _ in blocks + copy_blocks)
    keep = [(a, x - [L, L, 0]) for a, x in blocks]
    W_atoms, W_X, I_atoms, I_X = [], [], [], []
    k = n_sol
    while k < len(satoms):
        if satoms[k]["res"].strip() == "WAT":
            o = SX[k]
            if L <= o[0] < 2 * L and L <= o[1] < 2 * L:
                W_atoms += satoms[k:k + 3]; W_X.append(SX[k:k + 3] - [L, L, 0])
            k += 3
        else:
            I_atoms.append(satoms[k]); I_X.append(SX[k] - [L, L, 0]); k += 1
    box = (L, L, top + 2.0)                                     # a gap of one tolerance at the z boundary
    blocks_final = keep + [([a], np.array([x])) for a, x in zip(I_atoms, I_X)] + \
        [(W_atoms[3 * i:3 * i + 3], W_X[i]) for i in range(len(W_X))]
    blocks_final, n_drop = remove_clashes(blocks_final, np.array(box))
    log(f"cut: {n_drop} water molecules removed that clashed across the periodic boundaries")
    W_X = [x for a, x in blocks_final if a[0]["res"].strip() == "WAT"]
    write_pdb(out / "cut.pdb", blocks_final, cryst=box)
    tleap_box(cfg, out, "cut.pdb", box)
    # residue numbers of the grafted (head) units in the topology: chains come first, in order
    first, r = [], 1
    for a_, _ in keep:
        first.append(r); r += len({x["resid"] for x in a_})
    mask = ":" + ",".join(map(str, first))
    log(f"restraint mask of the grafted units (keep them on their plane during MD): {mask}")
    (out / "system.json").write_text(json.dumps({"box": box, "restraint_mask": mask, "n_chains": len(keep),
                                                 "n_water": len(W_X),
                                                 "n_water_removed_clashes": n_drop,
                                                 "n_ions": len(I_X), "n_water_region": n_w,
                                                 "grafting_density_nm2": len(keep) / (L / 10) ** 2}, indent=1))
    log(f"solvate: box {box[0]:.2f} x {box[1]:.2f} x {box[2]:.2f} A, {len(keep)} chains, {len(W_X)} water, "
        f"{len(I_X)} ions")


# ------------------------------------------------------------------ 5 single
def step_single(cfg):
    out = cfg["_out"] / "single"; out.mkdir(parents=True, exist_ok=True)
    shell = cfg["solvent"].get("shell_single", 30.0)
    for i, f in enumerate(chain_files(cfg), 1):
        d = out / f"chain_{i}"; d.mkdir(exist_ok=True)
        atoms, X = read_pdb(f)
        X = X - X.min(0) + shell
        hi = X.max(0) + shell
        q = charge_of(cfg, [(atoms, X)])
        path, _ = solvate_box(cfg, d, [(atoms, X)], (0, 0, hi[0], hi[1]), q, (0, 0, hi[0], hi[1]), hi[2])
        satoms, SX = read_pdb(path)
        n, k, rest = len(atoms), len(atoms), []
        while k < len(satoms):                                     # water molecules (3 atoms) and ions (1 atom)
            m = 3 if satoms[k]["res"].strip() == "WAT" else 1
            rest.append((satoms[k:k + m], SX[k:k + m])); k += m
        box = (hi[0] + 2.0, hi[1] + 2.0, hi[2] + 2.0)
        blocks, n_drop = remove_clashes([(atoms, SX[:n])] + rest, np.array(box))
        write_pdb(d / "single.pdb", blocks, cryst=box)
        tleap_box(cfg, d, "single.pdb", box)


STEPS = {"params": step_params, "chains": step_chains, "pack": step_pack, "solvate": step_solvate,
         "single": step_single}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", help="TOML configuration (see examples/pdmaema/config.toml)")
    ap.add_argument("step", choices=list(STEPS) + ["all"])
    a = ap.parse_args()
    cfg = load_config(a.config)
    for s in (["params", "chains", "pack", "solvate"] if a.step == "all" else [a.step]):
        log(f"=== {s}")
        STEPS[s](cfg)
