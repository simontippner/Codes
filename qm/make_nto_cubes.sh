#!/bin/bash
# Generate Gaussian cube files for specific MOs/NTOs from an ORCA .gbw
# (or an ORCA NTO file renamed to .gbw, e.g. `cp base.sN.nto base_sN.gbw`)
# using orca_plot's scripted interactive mode.
#
# Usage: make_nto_cubes.sh FILE.gbw MO1 [MO2 MO3 ...] [--res 60]
#
# Writes FILE.moXXXa.cube for each requested (alpha) MO index XXX in the
# current directory. Requires ORCA on PATH (source the same env as
# orca_submit.slurm, or run this on a node where `module load orca6` works).

set -e
gbw="$1"; shift
res=60
mos=()
while [ $# -gt 0 ]; do
  case "$1" in
    --res) res="$2"; shift 2 ;;
    *) mos+=("$1"); shift ;;
  esac
done

if [ -z "$gbw" ] || [ "${#mos[@]}" -eq 0 ]; then
  echo "Usage: $0 FILE.gbw MO1 [MO2 ...] [--res N]" >&2
  exit 1
fi

seq_input=""
for mo in "${mos[@]}"; do
  seq_input+="2\n${mo}\n4\n${res}\n11\n"
done
seq_input+="12\n"

printf "$seq_input" | orca_plot "$gbw" -i
echo "Done. Cube files: $(echo "$gbw" | sed 's/\.gbw$//').mo*a.cube"
