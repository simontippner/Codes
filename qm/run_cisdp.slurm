#!/bin/bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=60G
# TD-DFT (unrelaxed) difference density of one state with orca_plot option 6.
# usage: sbatch run_cisdp.slurm <dir> <base> <vector> <grid_intervals>
# TDA false: singlet k -> vector 2k-1; with Triplets true the triplets follow the singlets.
O=/usr/license/orca6/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2
export LD_LIBRARY_PATH=$O:/usr/remote/mpi/gnu/openmpi-4.1.8/lib:$LD_LIBRARY_PATH
cd $1
printf "4\n$4\n6\ny\n$3\n12\n" | $O/orca_plot $2.gbw -i > orca_plot_v$3.log 2>&1
grep -E "tgain|Output filename" orca_plot_v$3.log
