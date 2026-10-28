#!/usr/bin/env python3
"""Per-window adaptive equilibration trim for umbrella-sampling dist_XXX.X.dat
files (format: step x: v y: v z: v |r|), replacing a flat global trim
percentage with a trim derived from each window's own trajectory.

Rationale (established empirically + confirmed against the production mdin,
03_prod.sh): most windows settle onto their restrained target almost
instantly; a handful of windows (specifically those whose target is far from
where a fixed, window-independent 100 ps pre-equilibration stage happened to
leave the solute) take substantially longer -- up to several thousand frames
-- to actually arrive at the target. A single global trim fraction (10%, then
5%) either wastes good data in the former, or under-covers the latter. This
script instead finds, per window, the end of the INITIAL excursion block
(the real relaxation transient from frame 0) and trims exactly that, with a
small safety margin -- while explicitly NOT letting a rare, isolated,
late-trajectory outlier (a genuine large thermal fluctuation well after
settling) inflate the trim.

Algorithm per window:
  1. bad = frames where |z - target| > THRESH (5.0 A, ~15-20 sigma for a
     well-behaved window with std ~0.24-0.28 A -- i.e. essentially never
     triggered by ordinary sampling noise).
  2. If bad is empty: trim = FLOOR (a small fixed minimal burn-in, 500
     frames / 1%, standard conservative practice even for windows that show
     no detectable equilibration problem).
  3. Otherwise: walk the bad-frame list from the start; extend the "settling
     index" through consecutive bad frames, stopping the extension the
     moment a gap of >= GAP (1000 frames = 100 ps) of clean frames separates
     it from the next bad frame. This isolates the initial relaxation
     transient from any later, unrelated rare excursion.
  4. trim = settling_index + 1 + MARGIN (100 frames safety padding), floored
     at FLOOR.

Usage: adaptive_equil_trim.py DIST_DIR OUT_DIR
  Writes OUT_DIR/ts/win_*.dat, OUT_DIR/metadata.dat, OUT_DIR/trim_log.tsv
"""
import sys
import os
import re
import glob

THRESH = 5.0
GAP = 1000
MARGIN = 100
FLOOR = 500
K_GROSSFIELD = 10.0


def adaptive_trim(vals, target):
    bad = [i for i, z in enumerate(vals) if abs(z - target) > THRESH]
    if not bad:
        return FLOOR, 0, 0
    settle = bad[0]
    prev = bad[0]
    for idx in bad[1:]:
        if idx - prev > GAP:
            break
        settle = idx
        prev = idx
    trim = max(settle + 1 + MARGIN, FLOOR)
    return trim, len(bad), settle


def main():
    dist_dir, out_dir = sys.argv[1], sys.argv[2]
    ts_dir = os.path.join(out_dir, "ts")
    os.makedirs(ts_dir, exist_ok=True)

    files = sorted(glob.glob(os.path.join(dist_dir, "dist_*.dat")))
    if not files:
        sys.exit(f"no dist_*.dat found in {dist_dir}")

    meta_lines = []
    log_lines = ["window\tn_frames\ttrim\ttrim_pct\tn_bad_total\tsettle_idx"]
    for f in files:
        m = re.match(r"dist_([0-9.]+)\.dat", os.path.basename(f))
        c = m.group(1)
        target = float(c)
        vals = []
        with open(f) as fh:
            lines = fh.readlines()
        steps = []
        for line in lines:
            parts = line.split()
            if len(parts) < 8:
                continue
            steps.append(parts[0])
            vals.append(float(parts[6]))
        n = len(vals)
        trim, n_bad, settle = adaptive_trim(vals, target)
        trim = min(trim, n - 1)  # never trim the whole window

        with open(os.path.join(ts_dir, f"win_{c}.dat"), "w") as out:
            for i in range(trim, n):
                out.write(f"{steps[i]}\t{vals[i]}\n")
        meta_lines.append(f"ts/win_{c}.dat\t{target:.4f}\t{K_GROSSFIELD}\n")
        log_lines.append(f"{c}\t{n}\t{trim}\t{trim/n:.4f}\t{n_bad}\t{settle}")

    with open(os.path.join(out_dir, "metadata.dat"), "w") as mf:
        mf.writelines(meta_lines)
    with open(os.path.join(out_dir, "trim_log.tsv"), "w") as lf:
        lf.write("\n".join(log_lines) + "\n")

    trims = [float(l.split("\t")[3]) for l in log_lines[1:]]
    print(f"{len(files)} windows, trim range {min(trims)*100:.2f}%-{max(trims)*100:.2f}%, "
          f"mean {sum(trims)/len(trims)*100:.2f}%")


if __name__ == "__main__":
    main()
