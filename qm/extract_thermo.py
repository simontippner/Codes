#!/usr/bin/env python3
"""Extract thermochemistry (ORCA %freq output) from an ORCA .log/.out file.

Usage:
    python3 extract_thermo.py FILENAME [--json]

Parses the "Final Gibbs free energy", "Total Enthalpy", "Electronic energy",
"Zero point energy" and "Final entropy term" lines from an ORCA frequency
job's output and prints them (Eh and kcal/mol where ORCA gives both).

With --json, prints a single JSON object instead (handy for piping into
another script, e.g. compute_binding_free_energy.py).
"""
import argparse
import json
import re
import sys

PATTERNS = {
    "electronic_energy_Eh": r"Electronic energy\s*\.\.\.\s*(-?\d+\.\d+)\s*Eh",
    "zpe_Eh": r"Zero point energy\s*\.\.\.\s*(-?\d+\.\d+)\s*Eh",
    "total_enthalpy_Eh": r"Total Enthalpy\s*\.\.\.\s*(-?\d+\.\d+)\s*Eh",
    "entropy_term_Eh": r"Final entropy term\s*\.\.\.\s*(-?\d+\.\d+)\s*Eh",
    "gibbs_free_energy_Eh": r"Final Gibbs free energy\s*\.\.\.\s*(-?\d+\.\d+)\s*Eh",
    "temperature_K": r"Temperature\s*\.\.\.\s*(-?\d+\.\d+)\s*K",
}


def extract_thermo(path):
    with open(path) as f:
        text = f.read()

    result = {}
    for key, pattern in PATTERNS.items():
        matches = re.findall(pattern, text)
        if not matches:
            continue
        # Several of these lines also appear earlier in the SCF/property
        # sections for other purposes; the thermochemistry block's own
        # values are the LAST occurrence in the file.
        result[key] = float(matches[-1])

    if "gibbs_free_energy_Eh" not in result:
        raise ValueError(
            f"{path}: no 'Final Gibbs free energy' found -- is this an ORCA "
            f"frequency (%freq) job output?"
        )
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("filename")
    ap.add_argument("--json", action="store_true", help="print as JSON")
    args = ap.parse_args()

    thermo = extract_thermo(args.filename)

    if args.json:
        print(json.dumps(thermo, indent=2))
        return

    print(f"File: {args.filename}")
    for key, val in thermo.items():
        label = key.replace("_", " ").replace("Eh", "").strip()
        print(f"  {label:28s} {val: .8f}")


if __name__ == "__main__":
    main()
