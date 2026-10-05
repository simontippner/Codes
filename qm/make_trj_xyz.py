#!/usr/bin/env python3
"""
Extract every "CARTESIAN COORDINATES (ANGSTROEM)" block from a (possibly
still-running) ORCA optimization log and write them out as a standard
multi-frame xyz trajectory, readable directly in VMD (or any xyz-trajectory
viewer), without needing access to the job's scratch directory.

Usage:
    python3 make_trj_xyz.py FILENAME [--out trajectory.xyz]
"""
import argparse
import os
import re


def extract_frames(log_path):
    with open(log_path) as f:
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
    return frames


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("filename", help="ORCA geometry-optimization output/log file")
    p.add_argument("--out", default=None,
                    help="output xyz path (default: derived from the input filename)")
    args = p.parse_args()
    if args.out is None:
        args.out = re.sub(r"\.(log|out)$", "", os.path.basename(args.filename)) + "_trj.xyz"

    frames = extract_frames(args.filename)
    if not frames:
        raise SystemExit(f"No 'CARTESIAN COORDINATES (ANGSTROEM)' blocks found in {args.filename}")

    with open(args.out, "w") as f:
        for k, atoms in enumerate(frames):
            f.write(f"{len(atoms)}\n")
            f.write(f"cycle {k}\n")
            f.writelines(atoms)

    print(f"Wrote {len(frames)} frames ({len(frames[0])} atoms each) to {args.out}")


if __name__ == "__main__":
    main()
