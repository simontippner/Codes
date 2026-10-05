#!/bin/bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --ntasks-per-node=16

module use /usr/license/modulefiles/
module load orca6
export OBIN=`which orca`

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
