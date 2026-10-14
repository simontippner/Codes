#!/usr/bin/env python3
"""Patch missing element symbols in a sander/qm2_extern_orc_module-written
ORCA geometry file (inpfile.xyz), in place.

AmberTools' qm2_extern_orc_module atomic-number -> element-symbol lookup
omits Mo (Z=42) and Ru (Z=44) -- confirmed on both amber24 and amber25,
2026-10-13/14 -- so any Mo/Ru atom in the QM region gets written with a
blank label instead of its symbol, and ORCA then rejects the whole file
("Wrong syntax in xyz coordinates"). This project's QM region is always
either the Mo3S13 cluster (CAT, exactly 3 Mo atoms), the Ru(bpy)3 complex
(PS, exactly 1 Ru atom), or both together -- so the blank count alone
tells you which: 3 blanks -> all Mo (cat-alone region), 1 blank -> Ru
(ps-alone region), 4 blanks -> first 3 are Mo, last 1 is Ru (complex
region, matching this project's established atom-ordering convention:
CAT's atoms always precede PS's in the combined QM mask).

Usage (meant to be called by Codes/qm/orca_wrapper/orca before the real
ORCA binary runs, so the patch happens transparently inside sander's own
QM/MM loop):
    python3 patch_mo_ru_labels.py inpfile.xyz
"""
import sys


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "inpfile.xyz"
    with open(path) as f:
        lines = f.readlines()

    # sander writes the missing element as literal NULL bytes (\x00\x00),
    # not empty space -- confirmed 2026-10-14. \x00 is not whitespace, so a
    # plain .split() treats it as a spurious leading token and shifts every
    # coordinate field out of alignment (silently dropping z and writing
    # the null bytes themselves into the "patched" line) unless stripped
    # first.
    lines = [l.replace("\x00", "") for l in lines]

    blanks = [i for i, l in enumerate(lines) if i >= 2 and l.strip() and not l.split()[0][0].isalpha()]
    if not blanks:
        return  # nothing to patch

    if len(blanks) == 3:
        labels = ["Mo"] * 3
    elif len(blanks) == 1:
        labels = ["Ru"]
    elif len(blanks) == 4:
        labels = ["Mo", "Mo", "Mo", "Ru"]
    else:
        raise RuntimeError(f"{path}: {len(blanks)} blank-labeled atom lines, expected 1, 3, or 4 "
                            f"(this patch is specific to the Mo3S13/Ru(bpy)3 system's known "
                            f"Mo/Ru-only blanking pattern) -- refusing to guess, fix manually")

    for i, label in zip(blanks, labels):
        coords = lines[i].split()
        # explicit space separators -- these coordinates are 18-19 chars
        # long (full double precision), so a fixed field width (e.g. 14)
        # adds NO padding at all once exceeded and the fields silently run
        # together with no separator; confirmed 2026-10-14.
        lines[i] = f"{label:<3s} {coords[0]:>20s} {coords[1]:>20s} {coords[2]:>20s}\n"

    with open(path, "w") as f:
        f.writelines(lines)
    print(f"patched {len(blanks)} blank element label(s) in {path}: {labels}", file=sys.stderr)


if __name__ == "__main__":
    main()
