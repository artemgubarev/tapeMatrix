#!/bin/bash
#SBATCH --job-name=tapeMatrix
#SBATCH --output=out.log
#SBATCH --error=err.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
module load mpi/impi/2019.4.243
export INPUT_MATRIX_FILE="testData/matrix2000.txt"
mpicxx -o tapeMatrix tapeMatrixDemo.cpp
srun ./tapeMatrix