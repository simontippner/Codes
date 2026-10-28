#!/usr/bin/env python3
"""Scan a directory of dist_XXX.X.dat files (format: step x: v y: v z: v |r|)
for windows whose z column strays far (>5 A) from the window's own restraint
center, and report whether the anomalous frames are a contiguous block at the
start (equilibration lag, benign/fixable by trimming) or scattered throughout
(would indicate something else, e.g. imaging or genuine instability)."""
import sys
import glob
import os
import re

DIST_DIR = sys.argv[1]
THRESH = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0

files = sorted(glob.glob(os.path.join(DIST_DIR, "dist_*.dat")))
if not files:
    print(f"NO dist_*.dat FOUND in {DIST_DIR}")
    sys.exit(0)

n_anom_windows = 0
worst = []
for f in files:
    m = re.match(r"dist_([0-9.]+)\.dat", os.path.basename(f))
    center = float(m.group(1))
    vals = []
    with open(f) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 8:
                continue
            try:
                z = float(parts[6])
            except ValueError:
                continue
            vals.append(z)
    n = len(vals)
    if n == 0:
        print(f"  EMPTY: {os.path.basename(f)}")
        continue
    bad = [i for i, z in enumerate(vals) if abs(z - center) > THRESH]
    if not bad:
        continue
    n_anom_windows += 1
    contiguous_from_0 = (bad == list(range(bad[0], bad[-1] + 1))) and bad[0] == 0
    frac = len(bad) / n
    worst.append((os.path.basename(f), center, len(bad), n, frac, bad[0], bad[-1], contiguous_from_0))

print(f"{DIST_DIR}")
print(f"  {len(files)} windows scanned, threshold |z-center|>{THRESH} A")
print(f"  {n_anom_windows} windows with >=1 anomalous frame")
if worst:
    worst.sort(key=lambda w: -w[4])
    print(f"  worst {min(5,len(worst))} by fraction:")
    for name, center, nbad, ntot, frac, b0, b1, contig0 in worst[:5]:
        tag = "EQUIL-LAG(contig,frame0)" if contig0 else "SCATTERED/NON-CONTIGUOUS"
        print(f"    {name}: center={center:.1f} frac={frac:.4f} ({nbad}/{ntot}) "
              f"span=[{b0},{b1}] -> {tag}")
    n_scattered = sum(1 for w in worst if not w[7])
    if n_scattered:
        print(f"  *** {n_scattered} window(s) NOT simple equilibration lag -- needs manual look ***")
else:
    print("  clean -- no anomalies above threshold")
print()
