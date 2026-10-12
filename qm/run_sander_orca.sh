#!/bin/bash
# Run one sander+ORCA-EXTERN QM/MM single point, inside a FRAMES/<i>/
# folder written by setup_qmmm_ensemble.py (needs system.prmtop,
# snapshot.ncrst, qmmm.in, orc_job.tpl present -- all but snapshot.ncrst
# are symlinks shared across the whole ensemble).
#
# Submitted as a SLURM array task (see the submit_sander_array.slurm that
# setup_qmmm_ensemble.py writes next to each system's FRAMES/) that cd's
# into FRAMES/$SLURM_ARRAY_TASK_ID first, so $PWD here is one frame's own
# directory.
#
# ORCA path is hardcoded (not `module load orca6`) for the same reason as
# Codes/qm/orca_submit.slurm: at least one cluster node has been observed
# to be missing the environment-modules installation entirely, which
# silently leaves `orca` off PATH.
set -uo pipefail

MPILOC=/usr/remote/mpi/gnu/openmpi-4.1.8
ORCADIR=/usr/license/orca6/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2
export PATH="$ORCADIR:$MPILOC/bin:/usr/license/amber25/bin:$PATH"
export LD_LIBRARY_PATH="$ORCADIR:$MPILOC/lib:/usr/license/amber25/lib:$LD_LIBRARY_PATH"
export PYTHONPATH="/usr/license/amber25/lib/python3.12/site-packages:${PYTHONPATH:-}"
export AMBERHOME=/usr/license/amber25
export OMPI_MCA_btl="^openib"

if [ ! -x "$ORCADIR/orca" ]; then
  echo "FATAL: $ORCADIR/orca not found/executable on $HOSTNAME" >&2
  echo "$(pwd) : ORCA binary not found on $HOSTNAME ($(date))" >> I_CRASHED
  exit 1
fi
if ! command -v sander > /dev/null; then
  echo "FATAL: sander not found on PATH on $HOSTNAME" >&2
  echo "$(pwd) : sander not found on $HOSTNAME ($(date))" >> I_CRASHED
  exit 1
fi

echo "Entering job in $(pwd) on $HOSTNAME"
date

TOP=system.prmtop
RST=snapshot.ncrst
NAME=qmmm

export TMPDIR=/public/$(hostname -s)/scratch/tmp/$LOGNAME.$$
mkdir -p "$TMPDIR"
workdir="$TMPDIR/sander_orca/"
mkdir -p "$workdir"
ln -sf "$workdir" SCRATCH
cp "$TOP" "$RST" qmmm.in orc_job.tpl "$workdir"
cd "$workdir"

sander -O -i qmmm.in -o "$NAME.out" -p "$TOP" -c "$RST" -r "$NAME.rst" -x "$NAME.nc"
status=$?

cd - > /dev/null
cp "$workdir"/* .
rm -rf "$workdir" "$TMPDIR"

if [ $status -ne 0 ]; then
  echo "$(pwd) : sander exited $status on $HOSTNAME ($(date))" >> I_CRASHED
fi
echo "End time: $(date)"
exit $status
