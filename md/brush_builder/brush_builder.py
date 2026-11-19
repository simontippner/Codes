#!/usr/bin/env python3
"""brush_builder: all-atom models of randomly grafted, atactic polyelectrolyte brushes for Amber, with periodic
boundary conditions imposed during the packing, in water or solvent mixtures.

Method of S. Tippner, D. Hernandez-Castillo, F. H. Schacher, L. Gonzalez, J. Phys. Chem. B 129, 2105 (2025)
(PDMAEMA brushes); this is a rewrite of the original collection of bash and Python scripts as one tool.

Steps (each a subcommand; "all" runs 1 to 4, "single" is optional):
  1 params   residues for every monomer type (head, middle and tail unit, from a RESP-charged antechamber file) with
             prepgen and GAFF2 parameters with parmchk2
  2 chains   chains_per_box random sequences with a fixed number of protonated units and random stereochemistry,
             built with tleap and aligned along z (longest axis onto z, grafted unit at the bottom)
  3 pack     the box side L follows from grafting_density and chains_per_box (L = sqrt(chains / density)); without
             grafting_density, the densest packing at min_distance is searched as in the paper (smallest box first).
             Packing with periodic boundaries in a region of 3 x 3 cells: every chain is placed with its grafted atom
             on the plane z = anchor_z inside the central cell, against all chains placed before and their copies, and
             is then copied into the 8 neighbouring cells (copies kept if any atom lies within half a cell of the
             central cell). After max_attempts failures for one chain the packing restarts with chain 1; after
             max_fail failed packings min_distance is lowered (fixed density) or the region enlarged (densest mode).
  4 solvate  counterions are placed in the central cell against the chains and their copies and are copied into the
             neighbouring cells as well; then the central cell is filled with the solvent (water or a mixture by
             volume or mole fraction) around the chains, the ions and their copies, with the layer above the brush; the number of molecules follows
             from the density, corrected until the free solvent above the brush has the target density. The cell is
             cut out, molecules clashing across the periodic boundaries removed, and tleap writes the topology.
  5 single   (optional) every chain in the same solvent on its own, with its counterions, as reference systems

Requirements: Python >= 3.11 with numpy and scipy; AmberTools (prepgen, parmchk2, tleap) on PATH; PACKMOL 20.x.
Usage:  brush_builder.py <config.toml> {params,chains,pack,solvate,all,single}
"""
import argparse
import json
import math
import random
import re
import shutil
import subprocess
import tomllib
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

AMU_PER_A3 = 0.602214            # g/cm3 -> amu/A3
MASS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "S": 32.06, "P": 30.974, "F": 18.998, "Cl": 35.45,
        "Na": 22.99, "K": 39.098, "Br": 79.904, "I": 126.90}
WATER_PDB = ("HETATM    1  O   WAT     1       0.000   0.000   0.000\n"
             "HETATM    2  H1  WAT     1       0.957   0.000   0.000\n"
             "HETATM    3  H2  WAT     1      -0.240   0.927   0.000\nEND\n")
WATER_MODELS = {                 # tleap commands for the built-in water models (each also loads its ion parameters)
    "tip3p": ["source leaprc.water.tip3p"],
    "spce": ["source leaprc.water.spce"],
    "opc": ["source leaprc.water.opc"],
    "spcfw": ["source leaprc.water.spce", "WAT = SPF", "loadamberparams frcmod.spcfw", "set default FlexibleWater on"],
}
WATER_DENSITY = {"tip3p": 0.98, "spce": 0.997, "opc": 0.997, "spcfw": 1.01}   # densities of the models at 298 K


# ------------------------------------------------------------------ helpers
def log(msg):
    print(f"[brush_builder] {msg}", flush=True)


def run(cmd, cwd):
    """run a program in cwd, return stdout; raise with its output on failure"""
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"{' '.join(map(str, cmd))} failed in {cwd}:\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
    return r.stdout


def packmol(cfg, inp, cwd):
    """run PACKMOL on the input text (written to packmol.inp and read as packmol < packmol.inp); True on success"""
    (cwd / "packmol.inp").write_text(inp)
    with open(cwd / "packmol.inp") as f:
        r = subprocess.run([cfg["tools"]["packmol"]], cwd=cwd, stdin=f, capture_output=True, text=True)
    (cwd / "packmol.log").write_text(r.stdout + r.stderr)
    return "Success!" in r.stdout


def read_pdb(path):
    """atoms (record, name, residue, chain, residue number, rest of line) and coordinates"""
    lines = [l for l in Path(path).read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]
    atoms = [dict(rec=l[:6], name=l[12:16], res=l[17:20] + " ", chain=l[21], resid=l[22:26], tail=l[54:])
             for l in lines]
    X = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lines]).reshape(-1, 3)
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


ION_ELEMENT = {"CL-": "Cl", "CL": "Cl", "NA+": "Na", "NA": "Na", "K+": "K", "BR-": "Br", "I-": "I"}


def element(name):
    """element from a PDB atom name: monatomic ions by their names, otherwise the first letter"""
    n = name.strip()
    return ION_ELEMENT.get(n.upper(), n[0].upper())


def mass(atoms):
    return sum(MASS.get(element(a["name"]), 12.011) for a in atoms)


def load_config(path):
    cfg = tomllib.loads(Path(path).read_text())
    cfg["_dir"] = Path(path).resolve().parent
    cfg["_out"] = (cfg["_dir"] / cfg.get("output", "output")).resolve()
    cfg.setdefault("tools", {}).setdefault("packmol", shutil.which("packmol") or "packmol")
    b = cfg["brush"]
    b.setdefault("min_distance", b.get("tolerance", 2.0))         # "tolerance" accepted from version 1 configs
    b.setdefault("chains_per_box", cfg.get("chains", {}).get("n_chains"))
    s = cfg.setdefault("solvent", {})
    s.setdefault("water_model", "tip3p"); s.setdefault("mode", "volume"); s.setdefault("margin", 0.0)
    s.setdefault("layer", s.get("shell_z", 30.0)); s.setdefault("counterion", "Cl-")
    s.setdefault("components", [{"name": "water", "fraction": 100}])
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


def leap_header(cfg, solvent=False):
    p = cfg["_out"] / "params"
    lines = ["source leaprc.gaff2"]
    if solvent:
        lines += WATER_MODELS[cfg["solvent"]["water_model"]]
        for c in cfg["solvent"]["components"]:
            for f in c.get("params", []):
                f = (cfg["_dir"] / f).resolve()
                lines.append({".prepi": "loadamberprep", ".prep": "loadamberprep", ".lib": "loadoff",
                              ".off": "loadoff", ".frcmod": "loadamberparams"}[f.suffix] + f" {f}")
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
    return [("H" if i == 0 else "T" if i == n - 1 else "D") + k for i, k in enumerate(keys)]


def step_chains(cfg):
    out = cfg["_out"] / "chains"; out.mkdir(parents=True, exist_ok=True)
    b = cfg["brush"]
    rng = random.Random(cfg["chains"].get("seed"))
    seqs = []
    for i in range(1, b["chains_per_box"] + 1):
        d = out / f"chain_{i}"; d.mkdir(exist_ok=True)
        seq = sequence(cfg, rng); seqs.append(seq)
        (d / "tleap.in").write_text(leap_header(cfg) + f"mol = sequence {{ {' '.join(seq)} }}\n"
                                    "savepdb mol built.pdb\nquit\n")
        run(["tleap", "-f", "tleap.in"], d)
        # align along z: PACKMOL moves molecules only as rigid bodies, so fitting the chain into a narrow, tall box
        # (as in the original scripts) is a rotation of its longest axis onto z; done here directly
        atoms, X = read_pdb(d / "built.pdb")
        C = X - X.mean(0)
        v = np.linalg.eigh(np.cov(C.T))[1]
        Rm = v[:, ::-1][:, [1, 2, 0]]                               # longest axis -> z
        if np.linalg.det(Rm) < 0:
            Rm[:, 0] *= -1
        X = C @ Rm
        if X[b["anchor_atom"] - 1, 2] > X[:, 2].mean():             # grafted atom at the bottom
            X[:, 1:] *= -1
        X[:, 2] -= X[:, 2].min()                                    # lowest atom on z = 0
        write_pdb(out / f"chain_{i}.pdb", [(atoms, X)])
        log(f"chain {i}: {' '.join(seq)} (aligned, {np.ptp(X[:, 2]):.0f} A long)")
    (out / "sequences.json").write_text(json.dumps(seqs, indent=1))


# ------------------------------------------------------------------ 3 pack
def chain_files(cfg):
    return [cfg["_out"] / "chains" / f"chain_{i}.pdb" for i in range(1, cfg["brush"]["chains_per_box"] + 1)]


def upright(path, plane, region, z_a, height, atom):
    """PACKMOL structure block: chain upright (rotations about x and y fixed), grafted atom on the plane z_a"""
    x0, y0, x1, y1 = plane
    return (f"structure {path}\n  number 1\n  inside box {region[0]:.3f} {region[1]:.3f} 0. {region[2]:.3f} "
            f"{region[3]:.3f} {height:.3f}\n  constrain_rotation x 0. 0.\n  constrain_rotation y 0. 0.\n"
            f"  atoms {atom}\n    inside box {x0:.3f} {y0:.3f} {z_a:.3f} {x1:.3f} {y1:.3f} {z_a:.3f}\n"
            f"  end atoms\nend structure\n")


def smallest_box(cfg, files, d, height):
    """densest mode: smallest square box holding all chains upright at min_distance, widened by box_step"""
    b = cfg["brush"]; tol = b["min_distance"]
    L = b.get("box_start", math.sqrt(len(files)) * (tol + 4.0))
    while True:
        body = "".join(upright(f, (0, 0, L, L), (0, 0, L, L), b["anchor_z"], height, b["anchor_atom"]) for f in files)
        if packmol(cfg, f"tolerance {tol}\nfiletype pdb\noutput box.pdb\nnloop 20\nmaxit 20\nadd_amber_ter\n{body}", d):
            log(f"pack (a): all chains fit in {L:.2f} x {L:.2f} A at {tol} A")
            return L
        L += b.get("box_step", 0.5)


def copies(X, L, keep):
    """copies of X in the 8 neighbouring cells of the central cell [L, 2L]^2 that have an atom within `keep` of it"""
    out = []
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            if i == j == 0:
                continue
            Y = X + [i * L, j * L, 0]
            near = (Y[:, 0] > L - keep) & (Y[:, 0] < 2 * L + keep) & (Y[:, 1] > L - keep) & (Y[:, 1] < 2 * L + keep)
            if near.any():
                out.append(Y)
    return out


def periodic_pack(cfg, files, L, tol, height, d):
    """pack all chains with periodic copies at cell size L; returns (central, copies) or None after max_fail failures"""
    b = cfg["brush"]; R = 3 * L; fails = 0
    while fails < b["max_fail"]:
        central, copies_all, i, attempts = [], [], 0, 0
        while i < len(files):
            fixed = ""
            if central:
                write_pdb(d / "fixed.pdb", [(a, X) for a, X in central] + [(a, Y) for a, Y, _ in copies_all])
                fixed = "structure fixed.pdb\n  number 1\n  fixed 0. 0. 0. 0. 0. 0.\nend structure\n"
            new = upright(files[i], (L, L, 2 * L, 2 * L), (0, 0, R, R), b["anchor_z"], height, b["anchor_atom"])
            if not packmol(cfg, f"tolerance {tol}\nfiletype pdb\noutput placed.pdb\nnloop 20\nmaxit 20\n"
                                f"add_amber_ter\nseed -1\n{fixed}{new}", d):
                fails += 1; attempts += 1
                if fails >= b["max_fail"]:
                    return None
                if attempts >= b["max_attempts"]:                  # restart from the first chain
                    central, copies_all, i, attempts = [], [], 0, 0
                continue
            atoms_new, _ = read_pdb(files[i]); Xall = read_pdb(d / "placed.pdb")[1]
            Xn = Xall[-len(atoms_new):]                            # the new chain is the last molecule of the output
            central.append((atoms_new, Xn))
            copies_all += [(atoms_new, Y, i) for Y in copies(Xn, L, 0.5 * L)]
            log(f"  chain {i + 1}/{len(files)} placed ({len(copies_all)} copies, {fails} failed packings)")
            i += 1; attempts = 0
        return central, copies_all
    return None


def step_pack(cfg):
    out = cfg["_out"] / "pack"; out.mkdir(parents=True, exist_ok=True)
    work = out / "work"; work.mkdir(exist_ok=True)
    b = cfg["brush"]; files = chain_files(cfg); n = len(files)
    height = max(read_pdb(f)[1][:, 2].max() for f in files) + 4.0
    tol = b["min_distance"]
    if b.get("grafting_density"):                                 # fixed density: box from chains and density
        L = 10.0 * math.sqrt(n / b["grafting_density"])
        floor, step = b.get("min_distance_floor", 1.5), b.get("min_distance_step", 0.25)
        log(f"pack: {n} chains at {b['grafting_density']} nm^-2 -> box {L:.2f} x {L:.2f} A")
        while True:
            log(f"pack: periodic packing at a minimum distance of {tol:.2f} A")
            placed = periodic_pack(cfg, files, L, tol, height, work)
            if placed:
                break
            tol = round(tol - step, 3)
            if tol < floor - 1e-9:
                raise RuntimeError(
                    f"{n} chains do not pack at {b['grafting_density']} nm^-2 with a minimum distance of {floor} A or "
                    f"more. Lower grafting_density, or run without it (densest mode) to find the highest density.")
    else:                                                         # densest mode, as in the paper
        L = L_small = smallest_box(cfg, files, out, height)
        R = b.get("region_start", math.ceil(3 * L))
        while True:
            log(f"pack: region {R:.1f} A, cell {R / 3:.2f} A")
            placed = periodic_pack(cfg, files, R / 3, tol, height, work)
            if placed:
                L = R / 3; break
            R += b.get("region_step", 1.0)
    central, copies_all = placed
    write_pdb(out / "region.pdb", [(a, X) for a, X in central] + [(a, Y) for a, Y, _ in copies_all])
    write_pdb(out / "central.pdb", central)
    (out / "pack.json").write_text(json.dumps({
        "cell": L, "region": 3 * L, "height": height, "min_distance": tol, "n_central": len(central),
        "n_copies": len(copies_all), "block_chain": list(range(len(central))) + [c for _, _, c in copies_all],
        "grafting_density_nm2": len(central) / (L / 10) ** 2,
        **({} if b.get("grafting_density") else {"smallest_box": L_small})}, indent=1))
    log(f"pack: {len(central)} chains and {len(copies_all)} copies, box {L:.2f} A, "
        f"grafting density {len(central) / (L / 10) ** 2:.3f} nm^-2, minimum distance {tol:.2f} A")


# ------------------------------------------------------------------ 4 solvate
def components(cfg, d):
    """solvent components: template file, atom count, molar mass, pure density and the relative number of molecules
    per composition unit (from vol% via n ~ phi rho / M, or from mole fractions)"""
    s = cfg["solvent"]; comps = []
    for c in s["components"]:
        if c["name"] == "water":
            (d / "water.pdb").write_text(WATER_PDB)
            pdb, rho = d / "water.pdb", c.get("density", WATER_DENSITY[s["water_model"]])
        else:
            pdb, rho = (cfg["_dir"] / c["structure"]).resolve(), c["density"]
        atoms, X = read_pdb(pdb)
        M = mass(atoms)
        comps.append(dict(name=c["name"], pdb=pdb, n_atoms=len(atoms), M=M, rho=rho, f=c["fraction"],
                          res=atoms[0]["res"].strip()))
    tot = sum(c["f"] for c in comps)
    for c in comps:
        c["f"] /= tot
        c["unit"] = c["f"] * c["rho"] / c["M"] if s["mode"] == "volume" else c["f"]
    # density of the mixture: ideal mixing of volumes (or given)
    v = sum(c["unit"] * c["M"] / c["rho"] for c in comps); m = sum(c["unit"] * c["M"] for c in comps)
    s["_rho_mix"] = s.get("density_mix", m / v)
    return comps


def counts(comps, mass_target):
    """molecules per component for a total solvent mass (amu), keeping the composition"""
    m_unit = sum(c["unit"] * c["M"] for c in comps)
    return [max(1, round(c["unit"] * mass_target / m_unit)) for c in comps]


def ion_pdb(name, d):
    p = d / "ion.pdb"
    p.write_text(f"HETATM    1 {name:<4} {name:<3} X   1       0.000   0.000   0.000\nEND\n")
    return p


def fill(cfg, d, fixed_blocks, box, comps, z_free):
    """fill box (x0, y0, x1, y1, top) with the solvent components around the fixed blocks, in two zones:
      brush zone (0 to z_free): solvent mass = target density x zone volume - mass of the fixed molecules in it (the
        rule of PACKMOL's solvate.tcl; for a PDMAEMA brush it gives about 0.7 g/cm3 of water next to 0.3 g/cm3 of
        polymer), started at brush_fill (0.9) of it and lowered by 10% while PACKMOL cannot place it;
      layer (z_free to top): the target density, corrected until the measured density of the layer is within 3%.
    Both zones get the same composition. Returns per-component lists of molecule blocks and the counts."""
    s = cfg["solvent"]; rho = s["_rho_mix"]
    x0, y0, x1, y1, top = box
    write_pdb(d / "fixed.pdb", fixed_blocks)
    area = (x1 - x0) * (y1 - y0)
    Xf = np.vstack([x for _, x in fixed_blocks]); mf = np.array([MASS.get(element(a["name"]), 12.011)
                                                                  for atoms, _ in fixed_blocks for a in atoms])
    inside = (Xf[:, 0] > x0) & (Xf[:, 0] < x1) & (Xf[:, 1] > y0) & (Xf[:, 1] < y1)
    m_brush = mf[inside & (Xf[:, 2] < z_free)].sum(); m_layer = mf[inside & (Xf[:, 2] >= z_free)].sum()
    fill_b = s.get("brush_fill", 0.9)
    n_b = counts(comps, fill_b * max(rho * AMU_PER_A3 * area * z_free - m_brush, 0.0)) if z_free > 0 else [0] * len(comps)
    n_l = counts(comps, rho * AMU_PER_A3 * area * (top - z_free) - m_layer)
    log(f"solvent: brush zone 0 to {z_free:.1f} A: " + ", ".join(f"{k} {c['name']}" for c, k in zip(comps, n_b))
        + f" ({100 * fill_b:.0f}% of the solvate.tcl estimate); layer to {top:.1f} A: "
        + ", ".join(f"{k} {c['name']}" for c, k in zip(comps, n_l)))
    for it in range(8):
        inp = (f"tolerance 2.0\nfiletype pdb\noutput solvated.pdb\nadd_amber_ter\nseed -1\n"
               f"structure fixed.pdb\n  number 1\n  fixed 0. 0. 0. 0. 0. 0.\nend structure\n")
        groups = []
        for c, kb, kl in zip(comps, n_b, n_l):
            groups += [(c, kb, 0.0, z_free), (c, kl, z_free, top)]
        for c, k, z0, z1 in groups:
            if k > 0:
                inp += (f"structure {c['pdb']}\n  number {k}\n  inside box {x0:.3f} {y0:.3f} {z0:.3f} {x1:.3f} {y1:.3f} "
                        f"{z1:.3f}\nend structure\n")
        if not packmol(cfg, inp, d):
            n_b = [int(0.9 * k) for k in n_b]                       # brush zone too full: 10% less, try again
            log(f"solvent pass {it + 1}: PACKMOL could not place it; brush zone lowered to "
                + ", ".join(f"{k} {c['name']}" for c, k in zip(comps, n_b)))
            if it == 7:
                raise RuntimeError(f"PACKMOL could not place the solvent in {d} (see packmol.log)")
            continue
        atoms, X = read_pdb(d / "solvated.pdb")
        k0 = sum(len(a) for a, _ in fixed_blocks)
        mols = {c["name"]: [] for c in comps}; layer_mass = 0.0
        for c, k, z0, z1 in groups:
            for _ in range(max(k, 0)):
                a, x = atoms[k0:k0 + c["n_atoms"]], X[k0:k0 + c["n_atoms"]]
                mols[c["name"]].append((a, x)); k0 += c["n_atoms"]
                if z_free + 3 < x.mean(0)[2] < top - 3:
                    layer_mass += c["M"]
        rho_layer = layer_mass / AMU_PER_A3 / (area * max(top - z_free - 6, 1.0))
        log(f"solvent pass {it + 1}: brush zone " + ", ".join(f"{k} {c['name']}" for c, k in zip(comps, n_b))
            + "; layer " + ", ".join(f"{k} {c['name']}" for c, k in zip(comps, n_l))
            + f"; {rho_layer:.3f} g/cm3 in the layer (target {rho:.3f})")
        if abs(rho_layer / rho - 1) < 0.03 or top - z_free < 12:
            break
        n_l = [max(1, round(k * rho / rho_layer)) for k in n_l]
    n = [a + b for a, b in zip(n_b, n_l)]
    return [mols[c["name"]] for c in comps], n


def remove_clashes(blocks, removable, box, cut=2.0):
    """drop removable molecules (solvent) with a heavy atom closer than cut to a heavy atom of another molecule under
    periodic boundaries; removable is 0 for fixed molecules, otherwise a priority (the abundance of the molecule's
    component), so that of two clashing solvent molecules the one of the more abundant component goes"""
    X = np.vstack([x for _, x in blocks]); owner = np.concatenate([[k] * len(x) for k, (_, x) in enumerate(blocks)])
    heavy = np.array([not a["name"].strip().startswith("H") for b in blocks for a in b[0]])
    idx = np.where(heavy)[0]
    t = cKDTree(X[idx] % box, boxsize=box)
    drop = set()
    for i, j in t.query_pairs(cut, output_type="ndarray"):
        bi, bj = owner[idx[i]], owner[idx[j]]
        if bi == bj or bi in drop or bj in drop:
            continue
        if removable[bi] and removable[bj]:                       # both solvent: drop the more abundant component
            drop.add(bi if removable[bi] >= removable[bj] else bj)
        elif removable[bj]:
            drop.add(bj)
        elif removable[bi]:
            drop.add(bi)
    return [b for k, b in enumerate(blocks) if k not in drop], [r for k, r in enumerate(removable) if k not in drop], len(drop)


def exact_composition(blocks, comps):
    """after the cut the numbers of molecules of a mixture deviate from the composition by chance; remove the surplus
    of the over-represented components, highest molecules first (from the free solvent above the brush)"""
    if len(comps) < 2:
        return blocks, 0
    idx = {c["res"]: [k for k, (a, _) in enumerate(blocks) if a[0]["res"].strip() == c["res"]] for c in comps}
    scale = min(len(idx[c["res"]]) / c["unit"] for c in comps)
    drop = set()
    for c in comps:
        surplus = len(idx[c["res"]]) - int(math.floor(scale * c["unit"] + 1e-9))
        allowed = max(1, int(0.05 * len(idx[c["res"]])))         # at most 5% of a component (at least 1 molecule)
        if surplus > allowed:
            log(f"composition: {c['name']} would need {surplus} removals (> 5%); kept, see the composition in system.json")
            surplus = allowed
        if surplus > 0:
            top_first = sorted(idx[c["res"]], key=lambda k: -blocks[k][1][:, 2].mean())
            drop.update(top_first[:surplus])
    return [b for k, b in enumerate(blocks) if k not in drop], len(drop)


def charge_of(cfg, blocks):
    """number of protonated units (net charge) of the chains"""
    prot = {rn for k, m in cfg["monomers"].items() if m["protonated"] for rn in residue_names(k).values()}
    return sum(len({a["resid"] for a in atoms if a["res"].strip() in prot}) for atoms, _ in blocks)


def place_ions(cfg, d, fixed_blocks, n_ion, box, L=None):
    """counterions inside box (x0, y0, x1, y1, top), 1 A from its faces, against the fixed blocks; with L, the ions are
    placed in the central cell and their periodic copies (within the margin) are returned as well"""
    s = cfg["solvent"]
    if n_ion == 0:
        return [], []
    ion = ion_pdb(s["counterion"], d)
    write_pdb(d / "fixed.pdb", fixed_blocks)
    x0, y0, x1, y1, top = box
    for attempt in range(10):
        inp = (f"tolerance 2.0\nfiletype pdb\noutput ions.pdb\nseed -1\n"
               f"structure fixed.pdb\n  number 1\n  fixed 0. 0. 0. 0. 0. 0.\nend structure\n"
               f"structure {ion}\n  number {n_ion}\n  inside box {x0 + 1:.3f} {y0 + 1:.3f} 1. {x1 - 1:.3f} {y1 - 1:.3f} "
               f"{top - 1:.3f}\nend structure\n")
        if not packmol(cfg, inp, d):
            raise RuntimeError(f"PACKMOL could not place the counterions in {d} (see packmol.log)")
        atoms, X = read_pdb(d / "ions.pdb")
        ions = [([a], x[None]) for a, x in zip(atoms[-n_ion:], X[-n_ion:])]
        if L is None:
            return ions, []
        # ions must also keep their distance to the periodic images of the other ions
        P = np.array([x[0] for _, x in ions])
        t = cKDTree((P - [L, L, 0]) % [L, L, 1e6], boxsize=[L, L, 1e6])
        if not len(t.query_pairs(2.0)):
            break
    images = [(a, y) for a, x in ions for y in copies(x, L, 0.5 * L)]   # same range as the chain copies
    return ions, images


def solvent_report(cfg, comps, n):
    tot = sum(n); out = {}
    vols = [k * c["M"] / c["rho"] for c, k in zip(comps, n)]
    for c, k, v in zip(comps, n, vols):
        out[c["name"]] = {"molecules": k, "mole_fraction": round(k / tot, 4), "volume_percent": round(100 * v / sum(vols), 2)}
    return out


def step_solvate(cfg):
    out = cfg["_out"] / "solvate"; out.mkdir(parents=True, exist_ok=True)
    s = cfg["solvent"]; P = json.loads((cfg["_out"] / "pack" / "pack.json").read_text())
    L, m = P["cell"], s["margin"]
    atoms, X = read_pdb(cfg["_out"] / "pack" / "region.pdb")
    sizes = [len(read_pdb(f)[0]) for f in chain_files(cfg)]
    blocks, k = [], 0
    for c in P["block_chain"]:                                    # central chains first, then the copies
        blocks.append((atoms[k:k + sizes[c]], X[k:k + sizes[c]])); k += sizes[c]
    central, copy_blocks = blocks[:P["n_central"]], blocks[P["n_central"]:]
    z_poly = max(x[:, 2].max() for _, x in blocks)
    layer = max(30.0, 0.5 * z_poly) if s["layer"] == "auto" else float(s["layer"])
    top = z_poly + layer
    comps = components(cfg, out)
    log(f"solvent: " + ", ".join(f"{c['name']} {100 * c['f']:.1f} {'vol%' if s['mode'] == 'volume' else 'mol%'}"
                                 for c in comps) + f", mixture density {s['_rho_mix']:.3f} g/cm3, water model "
        f"{s['water_model']}; layer {layer:.1f} A above the brush (top {top:.1f} A)")
    # counterions: in the central cell against chains and copies, then copied into the margin like the chains
    q = charge_of(cfg, central)
    ions, ion_images = place_ions(cfg, out, central + copy_blocks, q, (L, L, 2 * L, 2 * L, top), L)
    log(f"counterions: {len(ions)} {s['counterion']} in the cell, {len(ion_images)} periodic copies around it")
    # solvent: central cell plus margin
    mols, n = fill(cfg, out, central + copy_blocks + ions + ion_images, (L - m, L - m, 2 * L + m, 2 * L + m, top),
                   comps, z_poly)
    # cut: central chains and ions, solvent molecules whose centre lies in the central cell; shift to the origin
    shift = np.array([L, L, 0.0])
    keep = [(a, x - shift) for a, x in central] + [(a, x - shift) for a, x in ions]
    removable = [0] * len(keep)
    for blocks_c in mols:
        for a, x in blocks_c:
            c = x.mean(0)
            if L <= c[0] < 2 * L and L <= c[1] < 2 * L:
                keep.append((a, x - shift)); removable.append(len(blocks_c))
    box = (L, L, top + 2.0)                                       # a gap of one tolerance at the z boundary
    keep, removable, n_drop = remove_clashes(keep, removable, np.array(box))
    log(f"cut: central cell kept, {n_drop} solvent molecules removed that clashed across the periodic boundaries")
    keep, n_fix = exact_composition(keep, comps)
    if n_fix:
        log(f"composition: {n_fix} surplus solvent molecules removed from the top of the box to restore the ratio")
    final_n = [sum(1 for a, _ in keep if a[0]["res"].strip() == c["res"]) for c in comps]
    write_pdb(out / "cut.pdb", keep, cryst=box)
    tleap_box(cfg, out, "cut.pdb", box)
    first, r = [], 1                                               # residue numbers of the grafted (head) units
    for a_, _ in central:
        first.append(r); r += len({x["resid"] for x in a_})
    mask = ":" + ",".join(map(str, first))
    rep = solvent_report(cfg, comps, final_n)
    (out / "system.json").write_text(json.dumps({
        "box": box, "restraint_mask": mask, "n_chains": len(central), "n_ions": len(ions),
        "grafting_density_nm2": len(central) / (L / 10) ** 2, "solvent": rep,
        "solvent_removed_clashes": n_drop, "water_model": s["water_model"]}, indent=1))
    log(f"restraint mask of the grafted units: {mask}")
    log(f"solvate: box {box[0]:.2f} x {box[1]:.2f} x {box[2]:.2f} A, {len(central)} chains, {len(ions)} ions, "
        + ", ".join(f"{v['molecules']} {k} ({100 * v['mole_fraction']:.1f} mol%, {v['volume_percent']:.1f} vol%)"
                    for k, v in rep.items()))


def tleap_box(cfg, d, pdb, box):
    (d / "tleap.in").write_text(leap_header(cfg, solvent=True) + f"mol = loadpdb {pdb}\n"
                                f"set mol box {{ {box[0]:.3f} {box[1]:.3f} {box[2]:.3f} }}\n"
                                "charge mol\nsavepdb mol system.pdb\nsaveamberparm mol system.prmtop system.inpcrd\nquit\n")
    out = run(["tleap", "-f", "tleap.in"], d)
    (d / "leap.out").write_text(out)
    errors = re.findall(r"Errors = (\d+)", out)
    if not errors or int(errors[-1]) > 0 or not (d / "system.prmtop").exists():
        raise RuntimeError(f"tleap failed in {d} (see leap.out)")
    q = re.findall(r"Total unperturbed charge:\s*([-\d.]+)", out)
    log(f"tleap: system.prmtop written, total charge {q[-1] if q else '?'}")


# ------------------------------------------------------------------ 5 single
def step_single(cfg):
    out = cfg["_out"] / "single"; out.mkdir(parents=True, exist_ok=True)
    shell = cfg["solvent"].get("shell_single", 30.0)
    for i, f in enumerate(chain_files(cfg), 1):
        d = out / f"chain_{i}"; d.mkdir(exist_ok=True)
        atoms, X = read_pdb(f)
        X = X - X.min(0) + shell
        hi = X.max(0) + shell
        comps = components(cfg, d)
        ions, _ = place_ions(cfg, d, [(atoms, X)], charge_of(cfg, [(atoms, X)]), (0, 0, hi[0], hi[1], hi[2]))
        mols, n = fill(cfg, d, [(atoms, X)] + ions, (0, 0, hi[0], hi[1], hi[2]), comps, X[:, 2].max())
        blocks = [(atoms, X)] + ions + [b for m in mols for b in m]
        removable = [0] * (1 + len(ions)) + [len(m) for m in mols for _ in m]
        box = (hi[0] + 2.0, hi[1] + 2.0, hi[2] + 2.0)
        blocks, removable, n_drop = remove_clashes(blocks, removable, np.array(box))
        write_pdb(d / "single.pdb", blocks, cryst=box)
        tleap_box(cfg, d, "single.pdb", box)


STEPS = {"params": step_params, "chains": step_chains, "pack": step_pack, "solvate": step_solvate,
         "single": step_single}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", help="TOML configuration (see examples/pdmaema/)")
    ap.add_argument("step", choices=list(STEPS) + ["all"])
    a = ap.parse_args()
    cfg = load_config(a.config)
    for s in (["params", "chains", "pack", "solvate"] if a.step == "all" else [a.step]):
        log(f"=== {s}")
        STEPS[s](cfg)
