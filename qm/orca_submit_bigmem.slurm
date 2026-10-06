#!/bin/bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=200G
# Large-NRoots TD-DFT jobs (Davidson expansion space = NRoots*MaxDim trial
# vectors, each a dense 254975-length double vector here) need much more
# memory than the cluster's default --mem-per-cpu allocation provides
# (4096 MB/cpu * 16 cpus = 64 GiB default cap). The 220-root singlet job
# was OOM-killed at exactly that 64 GiB default on its first attempt
# (job 8978388, 2026-10-05) despite ORCA's own %maxcore being set well
# under that nominal figure -- actual peak usage has overhead beyond what
# %maxcore accounts for. 200G is a generous margin on a 500G node.

# Set up ORCA/OpenMPI directly via PATH/LD_LIBRARY_PATH rather than relying
# on `module load` -- at least one cluster node has been observed to be
# missing the environment-modules installation entirely ("couldn't read
# file .../modulecmd.tcl"), which silently leaves `orca` off PATH and the
# job tries to execute its own input file as a shell command. Hardcoding
# the same paths the orca6/ompi modulefiles set makes this work regardless
# of whether Modules itself is intact on the node that gets assigned.
MPILOC=/usr/remote/mpi/gnu/openmpi-4.1.8
ORCADIR=/usr/license/orca6/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2
export PATH="$ORCADIR:$MPILOC/bin:$PATH"
export LD_LIBRARY_PATH="$ORCADIR:$MPILOC/lib:$LD_LIBRARY_PATH"
export OBIN="$ORCADIR/orca"

if [ ! -x "$OBIN" ]; then
  echo "FATAL: $OBIN not found/executable on $HOSTNAME -- aborting before" \
       "attempting to run (would otherwise silently try to execute the" \
       "input file as a shell command)." >&2
  echo "----" >> I_CRASHED
  echo "$1 : ORCA binary not found at $OBIN on $HOSTNAME" >> I_CRASHED
  echo "$(date)" >> I_CRASHED
  echo "----" >> I_CRASHED
  exit 1
fi

# the orca6 modulefile force-enables the InfiniBand (openib) BTL, which has
# been observed to segfault during MPI_Init on at least one cluster node --
# a node-level driver fault, not an ORCA/input problem. These are single-node
# jobs (all ranks communicate via shared memory anyway), so IB is never
# actually needed; excluding it avoids the crash with no performance cost.
export OMPI_MCA_btl="^openib"

input=$1
base=${input%.*}

export TMPDIR=/public/`hostname -s`/scratch/tmp/$LOGNAME.$$
mkdir $TMPDIR

if [ "x$SLURM_SUBMIT_DIR" == "x" ];
then
  SLURM_SUBMIT_DIR=$(pwd)
fi
cd $SLURM_SUBMIT_DIR

workdir=$TMPDIR/orca6/
mkdir -p $workdir
cp *.xyz $workdir
cp $input $workdir
# Also stage any .gbw present (e.g. for MOREAD / TD-DFT "Restart" jobs that
# reuse a previously converged wavefunction) -- harmless no-op if none exist.
cp *.gbw $workdir 2>/dev/null
# Also stage any .pc present (external point-charge files for %pointcharges
# QM/MM jobs) -- harmless no-op if none exist. Missing this causes a
# confusing crash ("ORCA finished by error termination in Startup" right
# after "Reading point charge file") rather than a clear file-not-found
# message, since the input's own directory has the file but the actual
# run happens in this separate scratch workdir.
cp *.pc $workdir 2>/dev/null
ln -sf $workdir SCRATCH

cd $workdir
echo $input
pwd
echo $HOSTNAME
date

$OBIN $input > $SLURM_SUBMIT_DIR/${base}.log

if [ $? == 0 ];
then
  cp * $SLURM_SUBMIT_DIR
  cd $SLURM_SUBMIT_DIR
  rm -r $workdir
  rm -r $TMPDIR
else
  cd $SLURM_SUBMIT_DIR
  echo "----" >> I_CRASHED
  echo ${base}.inp >> I_CRASHED
  echo $(date) >> I_CRASHED
  echo $HOSTNAME >> I_CRASHED
  echo $workdir >> I_CRASHED
  echo "----" >> I_CRASHED
fi
