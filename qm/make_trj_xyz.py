#!/usr/bin/env python3
"""
Extract every "CARTESIAN COORDINATES (ANGSTROEM)" block from a (possibly
still-running) ORCA optimization log and write them out as a standard
multi-frame xyz trajectory, readable directly in VMD (or any xyz-trajectory
viewer), without needing access to the job's scratch directory.
"""
import sys

LOG = sys.argv[1] if len(sys.argv) > 1 else "rubpy_mos_opt.log"
OUT = sys.argv[2] if len(sys.argv) > 2 else "rubpy_mos_trj_live.xyz"

with open(LOG) as f:
    lines = f.readlines()

frames = []
i = 0
while i < len(lines):
    if "CARTESIAN COORDINATES (ANGSTROEM)" in lines[i]:
        j = i + 2  # skip header line + dashed underline
        atoms = []
        while j < len(lines) and lines[j].strip():
            parts = lines[j].split()
            if len(parts) == 4:
                atoms.append(lines[j])
                j += 1
            else:
                break
        if atoms:
            frames.append(atoms)
        i = j
    else:
        i += 1

with open(OUT, "w") as f:
    for k, atoms in enumerate(frames):
        f.write(f"{len(atoms)}\n")
        f.write(f"cycle {k}\n")
        f.writelines(atoms)

print(f"Wrote {len(frames)} frames ({len(frames[0]) if frames else 0} atoms each) to {OUT}")
