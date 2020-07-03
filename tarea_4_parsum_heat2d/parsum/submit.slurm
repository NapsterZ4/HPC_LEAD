#!/bin/bash
#SBATCH --job-name=parsum_mpi
#SBATCH --output=parsum_mpi.o%j
#SBATCH --partition=nu
#SBATCH --ntasks=64
#SBATCH --time=00:05:00

cd $SLURM_SUBMIT_DIR

module load gcc/7.2.0
module load mpich/3.2.1-gcc-7.2.0

mpiexec -np 64 ./parsum