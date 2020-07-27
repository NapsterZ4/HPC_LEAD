#!/bin/bash
#SBATCH --job-name=cenatMD
#SBATCH --output=cenatMD.o%j
#SBATCH --partition=nu
#SBATCH --ntasks=8
#SBATCH --time=00:05:00

cd $SLURM_SUBMIT_DIR

module load mpich/3.2.1-gcc-7.2.0
module load gcc/7.2.0

###############Execution commands####################################

mpiexec -np 8 ./cenatMD 100 100 1 1 


######################################################################
