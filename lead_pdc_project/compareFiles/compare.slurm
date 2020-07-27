#!/bin/bash
#SBATCH --job-name=compare
#SBATCH --output=result_compare.txt
#SBATCH --partition=nu
#SBATCH --ntasks=1
#SBATCH --time=00:05:00

module load gcc/7.2.0

cd $SLURM_SUBMIT_DIR

################################################### Execution commands ##########################################################

./compare test_file.txt cenatMD_100.txt 0.2

################################################################################################################################
