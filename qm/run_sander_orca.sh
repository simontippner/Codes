#!/bin/bash
# Run one sander+ORCA-EXTERN QM/MM single point, inside a FRAMES/<i>/<region>/
# folder (e.g. FRAMES/1/cat/, FRAMES/1/ps/, FRAMES/1/ps_cat/) containing
# system.prmtop, snapshot.ncrst, orc_job.tpl (symlinks shared across the
# whole ensemble) and qmmm.in (specific to that region's qmmask/charge/mult).
#
# Usage: run_sander_orca.sh REGION   (e.g. "cat", "ps", "ps_cat")
#   cd's into the REGION subfolder (relative to $PWD, i.e. the frame's own
#   directory), runs sander there, names the output REGION.out/.rst/.nc.
#   The three regions of an interaction-energy decomposition
#   (dE_int = E(ps_cat) - E(cat) - E(ps), same embedding for all three --
#   see Codes/qm/split_qmmm_region.py for why a SHARED embedding matters)
#   each live in their own subfolder so their ORCA scratch files never
#   collide.
#
# Submitted as a SLURM array task that cd's into FRAMES/$SLURM_ARRAY_TASK_ID
# first, so $PWD here is one frame's own directory (one level above REGION).
#
# ORCA path is hardcoded (not `module load orca6`) for the same reason as
# Codes/qm/orca_submit.slurm: at least one cluster node has been observed
# to be missing the environment-modules installation entirely, which
# silently leaves `orca` off PATH.
set -uo pipefail

REGION="${1:?usage: run_sander_orca.sh REGION (e.g. cat, ps, ps_cat)}"

MPILOC=/usr/remote/mpi/gnu/openmpi-4.1.8
ORCADIR=/usr/license/orca6/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2
# sander (amber25) depends on libfftw3.so.3 (direct link) and libimf.so/
# libsvml.so/etc. (Intel runtime, dlopen'd by MKL at call time, so it
# doesn't show up as a direct `ldd` dependency) -- both resolved via the
# login node's own /lib64 + default search, but NOT present on at least
# one compute node ("indium", confirmed via "error while loading shared
# libraries" for each in turn). /usr/license is a shared filesystem
# visible from every node, so point at copies living there instead of
# relying on each compute node's own default library set.
FFTW_LIB=/usr/license/amber/amber17/lib
INTEL_LIB=/usr/license/intel/2024.1/lib   # matches the MKL 2024.1 sander links against
# orca_wrapper/orca MUST come first: it patches sander's ORCA-EXTERN
# geometry file (inpfile.xyz) for the missing Mo/Ru element symbols
# (AmberTools' qm2_extern_orc_module omits them, confirmed amber24+amber25)
# before handing off to the real binary -- see
# Codes/qm/patch_mo_ru_labels.py for the fix and why it's needed.
WRAPPER="$(dirname "$(readlink -f "$0")")/orca_wrapper"
export PATH="$WRAPPER:$ORCADIR:$MPILOC/bin:/usr/license/amber25/bin:$PATH"
export LD_LIBRARY_PATH="$ORCADIR:$MPILOC/lib:/usr/license/amber25/lib:$FFTW_LIB:$INTEL_LIB:$LD_LIBRARY_PATH"
export PYTHONPATH="/usr/license/amber25/lib/python3.12/site-packages:${PYTHONPATH:-}"
export AMBERHOME=/usr/license/amber25
export OMPI_MCA_btl="^openib"

if [ ! -x "$ORCADIR/orca" ]; then
  echo "FATAL: $ORCADIR/orca not found/executable on $HOSTNAME" >&2
  echo "$(pwd)/$REGION : ORCA binary not found on $HOSTNAME ($(date))" >> I_CRASHED
  exit 1
fi
if ! command -v sander > /dev/null; then
  echo "FATAL: sander not found on PATH on $HOSTNAME" >&2
  echo "$(pwd)/$REGION : sander not found on $HOSTNAME ($(date))" >> I_CRASHED
  exit 1
fi

if [ ! -d "$REGION" ]; then
  echo "FATAL: subfolder $REGION not found in $(pwd)" >&2
  echo "$(pwd)/$REGION : subfolder not found ($(date))" >> I_CRASHED
  exit 1
fi
cd "$REGION"

echo "Entering job in $(pwd) on $HOSTNAME, region=$REGION"
date

TOP=system.prmtop
RST=snapshot.ncrst
INPUT="qmmm.in"
NAME="${REGION}"

if [ ! -f "$INPUT" ]; then
  echo "FATAL: $INPUT not found in $(pwd)" >&2
  echo "$(pwd) : $INPUT not found ($(date))" >> I_CRASHED
  exit 1
fi

export TMPDIR=/public/$(hostname -s)/scratch/tmp/$LOGNAME.$$
mkdir -p "$TMPDIR"
workdir="$TMPDIR/sander_orca_${REGION}/"
mkdir -p "$workdir"
ln -sf "$workdir" "SCRATCH_${REGION}"
cp "$TOP" "$RST" "$INPUT" orc_job.tpl "$workdir"
cd "$workdir"

sander -O -i "$INPUT" -o "$NAME.out" -p "$TOP" -c "$RST" -r "$NAME.rst" -x "$NAME.nc"
status=$?

cd - > /dev/null
# Copy back only genuinely NEW output (never the staged-input basenames,
# which are symlinks here pointing at files SHARED across every frame of
# this system). An earlier version did a blanket `cp "$workdir"/* .`,
# which writes THROUGH those symlinks into the shared target -- harmless
# in isolation, but under a SLURM array with many frames running
# concurrently, many tasks race-wrote the SAME shared system.prmtop
# through their own symlink copies of it and corrupted it (confirmed
# incident, 2026-10-13, 1_PS_1_CAT/0_vol/buffer/exp -- recovered from the
# untouched top-level tleap.in output, but must never happen again).
for f in "$workdir"/*; do
  base=$(basename "$f")
  case "$base" in
    system.prmtop|snapshot.ncrst|orc_job.tpl|qmmm.in) continue ;;
  esac
  cp "$f" .
done
rm -rf "$workdir" "$TMPDIR"

if [ $status -ne 0 ]; then
  echo "$(pwd) : sander exited $status on $HOSTNAME ($(date))" >> I_CRASHED
fi
echo "End time: $(date)"
exit $status
