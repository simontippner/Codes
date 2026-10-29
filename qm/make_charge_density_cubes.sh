#!/bin/bash
## Electron-density and ESP cubes on ONE common grid with Multiwfn, for
## charge-density-difference (drho = rho_AB - rho_A - rho_B) and
## ESP-mapped-density plots of a two-fragment ORCA system.
##
## Usage (in the folder holding the .gbw files):
##   make_charge_density_cubes.sh grid  complex.gbw [spacing_bohr]   -> complex_dens.cub (defines the grid)
##   make_charge_density_cubes.sh dens  frag.gbw out.cub             -> density on the complex_dens.cub grid
##   make_charge_density_cubes.sh esp   complex.gbw out.cub          -> total ESP on the complex_dens.cub grid
##
## Fragment .gbw files must come from single points AT THE COMPLEX GEOMETRY.
## ECP cores are restored by Multiwfn's built-in EDF, consistently for complex
## and fragments, so they cancel in drho. ESP is slow (~2M points): run it
## through sbatch (16 threads), not on the login node.
set -e
ORCADIR=/usr/license/orca6/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2
export LD_LIBRARY_PATH="$ORCADIR:/usr/remote/mpi/gnu/openmpi-4.1.8/lib:$LD_LIBRARY_PATH"
MWFN=/usr/remote/multiwfn/Multiwfn  # the /usr/local/bin link exists only on login nodes
# compute nodes lack the X/Motif libs Multiwfn links against; copies from the login node
export LD_LIBRARY_PATH="/gpfs/data/stippner/lib_multiwfn_compute:$LD_LIBRARY_PATH"
ulimit -s unlimited
export OMP_STACKSIZE=500M
NT=${SLURM_CPUS_ON_NODE:-4}

mode=$1; gbw=$2; base=${gbw%.gbw}
[ -f settings.ini ] || cp /usr/remote/multiwfn/settings.ini .
sed -i "s/^\( *nthreads= *\)[0-9]*/\1$NT/" settings.ini

[ -f $base.molden ] || { $ORCADIR/orca_2mkl $base -molden > /dev/null; mv $base.molden.input $base.molden; }

case $mode in
  grid) sp=${3:-0.25}
        printf "5\n1\n4\n$sp\n2\n0\nq\n" | $MWFN $base.molden > mwfn_${base}_dens.log
        mv density.cub complex_dens.cub ;;
  dens) printf "5\n1\n8\ncomplex_dens.cub\n2\n0\nq\n" | $MWFN $base.molden > mwfn_${base}_dens.log
        mv density.cub $3 ;;
  esp)  printf "5\n12\n8\ncomplex_dens.cub\n2\n0\nq\n" | $MWFN $base.molden > mwfn_${base}_esp.log
        mv totesp.cub $3 ;;
  *) echo "mode must be grid|dens|esp" >&2; exit 1 ;;
esac
grep -A1 "Summing up all value" mwfn_${base}_*.log | tail -1
