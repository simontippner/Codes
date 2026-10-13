#!/usr/bin/env python3
"""Build a point-charge-embedding QM/MM ensemble from an existing AMBER MD
trajectory: equally-spaced snapshots, each turned into a QM-region .xyz +
MM point-charge .pc file (via make_qmmm_frame.py) and a ready-to-submit
ORCA .inp, laid out one subfolder per frame the same way as this group's
established Redoxpotentials QM-ensemble projects (FRAMES/<i>/...).

This is step 4+6+8+9 of 05_qmmm_ensemble/qmmm_ensemble_workflow.pdf: pick N
frames spaced across the trajectory, autoimage/strip to a finite MM shell
around the QM region (make_qmmm_frame.py), and write an ORCA ground-state
single-point input at the same B3LYP-D3(BJ)/def2-TZVP/RIJCOSX/def2-J level
used throughout 04_uvvis, but WITHOUT CPCM -- the MM point charges supply
the solvent electrostatics instead (combining both would double-count the
solvent response).

Frame extraction is a two-stage pipeline:
  Stage A (this script, one single cpptraj pass per --system): load the
  trajectory, autoimage-center on the QM region's first body (the CAT/
  Mo3S13 cluster, "MoS", for this project's default --qm-mask "^1,2"), and
  `trajout ... ncrestart multi` to dump every selected frame as its own
  full-system restart file straight away -- cpptraj numbers these
  sequentially (1, 2, 3, ...) by selection order, which is also the FRAMES
  numbering used below. This is the persisted, inspectable raw snapshot
  (same "don't just use a throwaway tempdir" convention as md/
  com_density_timeseries.py in this repo).
  Stage B (per frame, parallel): make_qmmm_frame.py reads that small
  single-frame snapshot.ncrst (not the original multi-GB trajectory) and
  truncates it to the QM .xyz + MM .pc shell.

Frame selection is deliberately dumb equal-spacing over the WHOLE
trajectory (frame 1, 1+stride, 1+2*stride, ...), not clustering -- this
script only sets up the cheap ground-state-SP ensemble path from the
workflow doc's decision #5/#10. Run the companion RMSD-autocorrelation
check (05_qmmm_ensemble/decorrelation/) first to know how many of the N
frames are actually statistically independent; this script does not
decide that for you, it just extracts what you ask for.

Does NOT submit any ORCA jobs -- it only writes a submit_array.slurm
(SLURM array over 1..N, one array task per frame, reusing
Codes/qm/orca_submit.slurm unmodified) for you to sbatch when ready.

Requires cpptraj + ParmEd on PATH/PYTHONPATH (module load amber) and
make_qmmm_frame.py next to this script.

Usage:
    python3 setup_qmmm_ensemble.py \\
        --system "0% water" PRMTOP1 TRAJ1 \\
        --system "50% water" PRMTOP2 TRAJ2 \\
        --qm-mask "^1,2" --shell 18.0 \\
        --qm-charge 0 --qm-mult 1 \\
        --n-frames 2000 --nprocs 16 \\
        --out-dir QM_ensemble_SPE --workers 16
"""
import argparse
import os
import re
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
MAKE_FRAME = os.path.join(HERE, "make_qmmm_frame.py")

# SlowConv KDIIS SOSCF: convergence stabilizers, not part of the original
# 04_uvvis method -- added after frame 1 of 1_PS_1_CAT/0_vol/buf=exp's SCF
# oscillated for 50 iterations (never settling, energy still swinging ~30
# Eh step to step) and then fell back to ORCA's slow TRAH algorithm.
# Matches the convergence aids the Redoxpotentials reference pipeline
# already used for its own (different) difficult transition-metal SCF.
ORCA_METHOD = "! B3LYP D3BJ def2-TZVP def2/J RIJCOSX TightSCF DefGrid2 SlowConv KDIIS SOSCF"

INP_TEMPLATE = """# Ground-state single point, QM/MM point-charge embedding
# QM region: {qm_mask} (from {label})
# MM region: everything within the shell written to {pc_name}
# Same functional/basis as 04_uvvis (B3LYP-D3(BJ)/def2-TZVP), WITHOUT CPCM:
# the MM point charges supply the solvent electrostatics instead.
{method}
%pal nprocs {nprocs} end
%maxcore {maxcore}
%pointcharges "{pc_name}"

* xyzfile {charge} {mult} {xyz_name}
"""


def n_frames_in_traj(prmtop, traj):
    with tempfile.NamedTemporaryFile("w", suffix=".in", delete=False) as f:
        f.write(f"parm {prmtop}\ntrajin {traj}\nrun\n")
        script = f.name
    result = subprocess.run(["cpptraj", "-i", script], capture_output=True, text=True)
    os.unlink(script)
    import re
    m = re.search(r"Coordinate processing will occur on (\d+) frames", result.stdout)
    if not m:
        raise RuntimeError(f"Could not determine frame count for {traj}:\n{result.stdout}\n{result.stderr}")
    return int(m.group(1))


def anchor_from_qm_mask(qm_mask):
    """First body named in --qm-mask -- the CAT/Mo3S13 ("MoS") cluster for
    this project's convention of listing CAT before PS, e.g. "^1,2" -> "^1".
    """
    m = re.match(r"^([:^@])", qm_mask)
    if not m:
        raise ValueError(f"qm_mask {qm_mask!r} must start with one of : ^ @")
    prefix = m.group(1)
    first_token = qm_mask[1:].split(",")[0]
    return f"{prefix}{first_token}"


def extract_ncrst_batch(prmtop, traj, frame_indices, anchor_mask, framesdir):
    """Stage A: one cpptraj pass -- load the trajectory, autoimage-center
    on anchor_mask (the CAT/MoS cluster), and dump every selected frame as
    its own full-system restart via `trajout ... ncrestart multi`.
    cpptraj numbers the written files sequentially by selection order
    (snapshot.ncrst.1, .2, ... -- NOT by original frame number), which
    conveniently matches the 1-based FRAMES/<i> numbering used throughout
    this pipeline.
    """
    tag_out = os.path.join(framesdir, "snapshot.ncrst")
    lines = [f"parm {prmtop}"]
    for fidx in frame_indices:
        lines.append(f"trajin {traj} {fidx} {fidx}")
    lines.append(f"autoimage anchor {anchor_mask} origin")
    lines.append(f"trajout {tag_out} ncrestart multi")
    lines.append("run")
    script_path = os.path.join(framesdir, "extract_ncrst.in")
    with open(script_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    result = subprocess.run(["cpptraj", "-i", script_path], capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise RuntimeError(f"cpptraj ncrestart extraction failed for {traj}")

    for i in range(1, len(frame_indices) + 1):
        src = f"{tag_out}.{i}"
        framedir = os.path.join(framesdir, str(i))
        os.makedirs(framedir, exist_ok=True)
        dst = os.path.join(framedir, "snapshot.ncrst")
        if not os.path.exists(src):
            raise RuntimeError(f"expected {src} after cpptraj ncrestart extraction, not found")
        os.replace(src, dst)
    return script_path


def extract_one(prmtop, framedir, qm_mask, shell):
    """Stage B: strip the already-centered per-frame snapshot.ncrst (from
    extract_ncrst_batch) down to the QM .xyz + MM .pc shell. Frame index is
    always 1 -- snapshot.ncrst is a single-frame restart, not the original
    multi-thousand-frame trajectory."""
    snap = os.path.join(framedir, "snapshot.ncrst")
    prefix = os.path.join(framedir, "frame")
    result = subprocess.run(
        ["python3", MAKE_FRAME, prmtop, snap, "1",
         "--qm-mask", qm_mask, "--shell", str(shell), "--out-prefix", prefix],
        capture_output=True, text=True)
    if result.returncode != 0:
        return (framedir, False, result.stdout + result.stderr)
    return (framedir, True, result.stdout.strip())


def write_inp(framedir, label, qm_mask, charge, mult, nprocs, maxcore):
    inp_path = os.path.join(framedir, "frame.inp")
    with open(inp_path, "w") as f:
        f.write(INP_TEMPLATE.format(
            qm_mask=qm_mask, label=label, pc_name="frame_mm.pc",
            method=ORCA_METHOD, nprocs=nprocs, maxcore=maxcore,
            charge=charge, mult=mult, xyz_name="frame_qm.xyz"))
    return inp_path


QMMM_IN_TEMPLATE = """qmmm single point via ORCA EXTERN (ifqnt), already-centered snapshot.ncrst
&cntrl
 nstlim = 0, dt = 0.0005,
 ntx = 1, irest = 0, iwrap=0, cut = 999.0,
 ntpr = 1, ntwr = 1, ntwx = 1,
 tempi = 300., temp0 = 300.,
 ntt = 0, gamma_ln=0.,
 ntb = 0, ntp = 0,
 ntc = 2, ntf = 2,
 noshakemask = "{qm_resmask}",
 ifqnt=1,
&end

&qmmm
qmmask='{qm_resmask}',
qmcharge={charge},
spin={mult},
qm_theory='EXTERN',
qmshake=0,
qm_ewald=0,
qm_pme=0,
qmcut=999.,
adjust_q=0,
&end
&orc
use_template=1,
&end
"""

ORC_TPL_TEMPLATE = """{method}
%maxcore {maxcore}

%scf
  MaxIter 5000
end
"""


def write_sander_inputs(framesdir, prmtop, n_frames, regions, maxcore):
    """Write the shared sander qmmm_<tag>.in (one per region: complex,
    and optionally cat-alone/ps-alone for an interaction-energy decomposition,
    dE_int = E(complex) - E(cat) - E(ps), same embedding for all three)
    + ORCA orc_job.tpl (same
    method for every region, so one shared file) once at the FRAMES/
    level, then symlink them -- plus the system topology -- into every
    FRAMES/<i>/ so each frame folder is self-contained for `sander -p
    system.prmtop -c snapshot.ncrst -i qmmm_<tag>.in` without duplicating
    the (several-MB) topology file N times.

    regions: list of dicts {tag, resmask, charge, mult}. The "complex"
    tag additionally gets a bare "qmmm.in" alias symlink for backward
    compatibility with earlier single-region frame folders/jobs.
    """
    orc_tpl = os.path.join(framesdir, "orc_job.tpl")
    with open(orc_tpl, "w") as f:
        f.write(ORC_TPL_TEMPLATE.format(method=ORCA_METHOD, maxcore=maxcore))

    region_files = {}
    for r in regions:
        path = os.path.join(framesdir, f"qmmm_{r['tag']}.in")
        with open(path, "w") as f:
            f.write(QMMM_IN_TEMPLATE.format(qm_resmask=r["resmask"], charge=r["charge"], mult=r["mult"]))
        region_files[r["tag"]] = path

    prmtop_abs = os.path.abspath(prmtop)
    for i in range(1, n_frames + 1):
        framedir = os.path.join(framesdir, str(i))
        links = [("orc_job.tpl", orc_tpl), ("system.prmtop", prmtop_abs)]
        for tag, path in region_files.items():
            links.append((f"qmmm_{tag}.in", path))
        if "complex" in region_files:
            links.append(("qmmm.in", region_files["complex"]))  # backward-compat alias
        for name, target in links:
            link = os.path.join(framedir, name)
            if os.path.islink(link) or os.path.exists(link):
                os.remove(link)
            os.symlink(os.path.abspath(target), link)
    return orc_tpl, region_files


SLURM_ARRAY_TEMPLATE = """#!/bin/bash
#SBATCH -J qmmm_{tag}
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --ntasks-per-node={nprocs}
#SBATCH --array=1-{n_frames}%{concurrency}
#SBATCH --output=FRAMES/%a/slurm-%A_%a.out
#SBATCH --error=FRAMES/%a/slurm-%A_%a.err

cd "FRAMES/$SLURM_ARRAY_TASK_ID"
# orca_submit.slurm does its own `cd $SLURM_SUBMIT_DIR` internally (SLURM
# sets that to the ORIGINAL sbatch directory, i.e. one level up from here)
# -- override it to the frame dir we just cd'd into, or that internal cd
# undoes this one and orca_submit.slurm looks for frame.inp/*.xyz in the
# wrong directory (confirmed failure mode, 2026-10-13: "cannot stat
# frame.inp", "Cannot open input file: frame.inp").
export SLURM_SUBMIT_DIR="$(pwd)"
bash {orca_submit} frame.inp
"""


def write_submit_array(outdir, n_frames, nprocs, concurrency, tag):
    path = os.path.join(outdir, "submit_array.slurm")
    orca_submit = os.path.join(HERE, "orca_submit.slurm")
    with open(path, "w") as f:
        f.write(SLURM_ARRAY_TEMPLATE.format(
            tag=tag, nprocs=nprocs, n_frames=n_frames, concurrency=concurrency,
            orca_submit=orca_submit))
    os.chmod(path, 0o755)
    return path


SANDER_ARRAY_TEMPLATE = """#!/bin/bash
#SBATCH -J qmmmS_{tag}_{region}
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-{n_frames}%{concurrency}
#SBATCH --output=FRAMES/%a/slurm-{region}-%A_%a.out
#SBATCH --error=FRAMES/%a/slurm-{region}-%A_%a.err

cd "FRAMES/$SLURM_ARRAY_TASK_ID"
bash {run_sander_orca} {region}
"""


def write_submit_sander_array(outdir, n_frames, concurrency, tag, region):
    """SLURM array for one region (complex/cat/ps) of the sander+ORCA-EXTERN
    path, one array task per frame via Codes/qm/run_sander_orca.sh <region>.
    Submit a subset with `sbatch --array=1-K%C submit_sander_array_<region>.slurm`
    (SLURM CLI flags override the file's own #SBATCH directives) to test
    before committing the full array.
    """
    path = os.path.join(outdir, f"submit_sander_array_{region}.slurm")
    run_sander_orca = os.path.join(HERE, "run_sander_orca.sh")
    with open(path, "w") as f:
        f.write(SANDER_ARRAY_TEMPLATE.format(
            tag=tag, region=region, n_frames=n_frames, concurrency=concurrency,
            run_sander_orca=run_sander_orca))
    os.chmod(path, 0o755)
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--system", nargs=3, action="append", required=True,
                     metavar=("LABEL", "PRMTOP", "TRAJ"))
    ap.add_argument("--qm-mask", default="^1,2",
                     help="cpptraj molecule-number mask for make_qmmm_frame.py (shell-truncation path)")
    ap.add_argument("--qm-resmask", default=None,
                     help="standard AMBER residue mask for sander's &qmmm qmmask/noshakemask "
                          "(sander does not understand '^' molecule masks) -- required to write "
                          "qmmm.in/orc_job.tpl, e.g. \":1-11\". This is the COMPLEX (CAT+PS) region.")
    ap.add_argument("--cat-resmask", default=None,
                     help="standard AMBER residue mask for the CAT-alone region (interaction-energy "
                          "decomposition), e.g. \":1-10\" -- PS becomes classical/MM in this region's "
                          "qmmm_cat.in. Requires --qm-resmask; --cat-charge/--cat-mult to match.")
    ap.add_argument("--cat-charge", type=int, default=-2)
    ap.add_argument("--cat-mult", type=int, default=1)
    ap.add_argument("--ps-resmask", default=None,
                     help="standard AMBER residue mask for the PS-alone region, e.g. \":11\" -- CAT "
                          "becomes classical/MM in this region's qmmm_ps.in. Requires --qm-resmask; "
                          "--ps-charge/--ps-mult to match.")
    ap.add_argument("--ps-charge", type=int, default=2)
    ap.add_argument("--ps-mult", type=int, default=1)
    ap.add_argument("--shell", type=float, default=18.0)
    ap.add_argument("--qm-charge", type=int, default=0)
    ap.add_argument("--qm-mult", type=int, default=1)
    ap.add_argument("--n-frames", type=int, default=2000)
    ap.add_argument("--nprocs", type=int, default=16, help="ORCA %%pal nprocs per job")
    ap.add_argument("--maxcore", type=int, default=4000, help="ORCA %%maxcore (MB) per core")
    ap.add_argument("--concurrency", type=int, default=100,
                     help="SLURM array %% throttle (max concurrently running tasks)")
    ap.add_argument("--out-dir", default="QM_ensemble_SPE",
                     help="subfolder created inside each system's MD dir")
    ap.add_argument("--workers", type=int, default=16,
                     help="parallel worker processes for the (cheap, CPU-only) frame extraction")
    ap.add_argument("--sander-only", action="store_true",
                     help="skip Stage A/B (assume FRAMES/<i>/ already populated from a previous run) "
                          "and only (re)write qmmm.in/orc_job.tpl + symlinks -- requires --qm-resmask")
    args = ap.parse_args()

    if args.qm_resmask is None and args.sander_only:
        ap.error("--sander-only requires --qm-resmask")

    def build_regions():
        regions = []
        if args.qm_resmask:
            regions.append({"tag": "complex", "resmask": args.qm_resmask,
                             "charge": args.qm_charge, "mult": args.qm_mult})
        if args.cat_resmask:
            regions.append({"tag": "cat", "resmask": args.cat_resmask,
                             "charge": args.cat_charge, "mult": args.cat_mult})
        if args.ps_resmask:
            regions.append({"tag": "ps", "resmask": args.ps_resmask,
                             "charge": args.ps_charge, "mult": args.ps_mult})
        return regions

    for label, prmtop, traj in args.system:
        mddir = os.path.dirname(prmtop)
        outdir = os.path.join(mddir, args.out_dir)
        framesdir = os.path.join(outdir, "FRAMES")
        os.makedirs(framesdir, exist_ok=True)

        if args.sander_only:
            existing = [int(d) for d in os.listdir(framesdir) if d.isdigit()]
            n_frames_here = max(existing) if existing else 0
            regions = build_regions()
            print(f"=== {label}: --sander-only, found {n_frames_here} existing FRAMES/<i>/, "
                  f"regions={[r['tag'] for r in regions]} ===")
            write_sander_inputs(framesdir, prmtop, n_frames_here, regions, args.maxcore)
            tag = re.sub(r"[^A-Za-z0-9]+", "_", label).strip("_") or "sys"
            for r in regions:
                sarray_path = write_submit_sander_array(outdir, n_frames_here, args.concurrency,
                                                          tag, r["tag"])
                print(f"  wrote {sarray_path} (NOT submitted)")
            print(f"  wrote orc_job.tpl + qmmm_<region>.in, symlinked into all {n_frames_here} frame dirs")
            continue

        n_total = n_frames_in_traj(prmtop, traj)
        stride = max(1, n_total // args.n_frames)
        frame_indices = [1 + k * stride for k in range(args.n_frames) if 1 + k * stride <= n_total]

        print(f"=== {label}: {traj} has {n_total} frames, extracting {len(frame_indices)} "
              f"(stride {stride}) into {framesdir} ===")

        manifest_path = os.path.join(outdir, "manifest.txt")
        with open(manifest_path, "w") as mf:
            mf.write(f"label: {label}\nprmtop: {prmtop}\ntraj: {traj}\n"
                     f"n_total_frames: {n_total}\nn_extracted: {len(frame_indices)}\n"
                     f"stride: {stride}\nqm_mask: {args.qm_mask}\nshell: {args.shell}\n"
                     f"qm_charge: {args.qm_charge}\nqm_mult: {args.qm_mult}\n"
                     f"orca_method: {ORCA_METHOD}\n"
                     f"i -> original_frame_index\n")
            for i, fidx in enumerate(frame_indices, start=1):
                mf.write(f"{i}\t{fidx}\n")

        anchor_mask = anchor_from_qm_mask(args.qm_mask)
        print(f"  Stage A: cpptraj ncrestart extraction, anchor={anchor_mask} (CAT/MoS) ...")
        extract_ncrst_batch(prmtop, traj, frame_indices, anchor_mask, framesdir)
        print(f"  Stage A done: {len(frame_indices)} snapshot.ncrst written under {framesdir}/<i>/")

        print("  Stage B: QM .xyz + MM .pc per frame ...")
        jobs = []
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            for i in range(1, len(frame_indices) + 1):
                framedir = os.path.join(framesdir, str(i))
                jobs.append(pool.submit(extract_one, prmtop, framedir, args.qm_mask, args.shell))
            n_ok, n_fail = 0, 0
            for fut in as_completed(jobs):
                framedir, ok, msg = fut.result()
                if ok:
                    n_ok += 1
                else:
                    n_fail += 1
                    print(f"  FAILED ({framedir}):\n{msg}")
        print(f"  Stage B done: {n_ok} ok, {n_fail} failed")

        for i in range(1, len(frame_indices) + 1):
            framedir = os.path.join(framesdir, str(i))
            if os.path.exists(os.path.join(framedir, "frame_qm.xyz")):
                write_inp(framedir, label, args.qm_mask, args.qm_charge, args.qm_mult,
                          args.nprocs, args.maxcore)

        tag = re.sub(r"[^A-Za-z0-9]+", "_", label).strip("_") or "sys"
        regions = build_regions()
        if regions:
            write_sander_inputs(framesdir, prmtop, len(frame_indices), regions, args.maxcore)
            for r in regions:
                sarray_path = write_submit_sander_array(outdir, len(frame_indices), args.concurrency,
                                                          tag, r["tag"])
                print(f"  wrote {sarray_path} (NOT submitted)")
            print(f"  wrote orc_job.tpl + qmmm_<region>.in ({[r['tag'] for r in regions]}), "
                  f"symlinked into all {len(frame_indices)} frame dirs")
        array_path = write_submit_array(outdir, len(frame_indices), args.nprocs,
                                         args.concurrency, tag=tag)
        print(f"  wrote {array_path} (NOT submitted -- sbatch it yourself when ready)")


if __name__ == "__main__":
    main()
