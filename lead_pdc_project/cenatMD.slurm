#!/bin/bash
#SBATCH --job-name=cenatMD
#SBATCH --output=cenatMD.o%j
#SBATCH --partition=nu-wide
#SBATCH --ntasks=7
#SBATCH --nodes=7
#SBATCH --time=00:20:00

cd $SLURM_SUBMIT_DIR

module load mpich/3.2.1-gcc-7.2.0
module load gcc/7.2.0

###############Execution commands####################################
export OMP_NUM_THREADS=128
mpiexec -np 16 ./cenatMD 200 7000 0 0 


######################################################################
