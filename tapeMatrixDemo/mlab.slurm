#!/bin/bash
#SBATCH --job-name=matlab
#SBATCH --output=mout.log
#SBATCH --error=merr.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=01:00:00

module load matlab/2022b
export INPUT_MATRIX_FILE="testData/matrix100x60.txt"
matlab -nodisplay -nosplash -r "run('msolve.m'); exit;" > matlab_run.log 2>&1