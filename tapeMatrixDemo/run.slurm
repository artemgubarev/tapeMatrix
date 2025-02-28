#!/bin/bash
#SBATCH --job-name=tapeMatrix
#SBATCH --output=out.log
#SBATCH --error=err.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00

module load mpi/impi/2019.4.243

export INPUT_MATRIX_FILE="testData/matrix100.txt"
export SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

# 0 - serial
# 1 - pthreads
# 2 - openMP
# 3 - MPI
# 4 - python MPI
export MODE = 1

if [ "$MODE" -eq 0 ]; then
    g++ -std=c++11 -o tapeMatrix tapeMatrixDemo.cpp
    srun ./tapeMatrix
elif [ "$MODE" -eq 1 ]; then
    g++ -std=c++11 -pthread -o tapeMatrix tapeMatrixDemo.cpp
    srun ./tapeMatrix
elif [ "$MODE" -eq 2 ]; then
    g++ -std=c++11 -fopenmp -o tapeMatrix tapeMatrixDemo.cpp
    srun ./tapeMatrix
elif [ "$MODE" -eq 3 ]; then
    mpicxx -std=c++11 -o tapeMatrix tapeMatrixDemo.cpp
    mpirun -n $SLURM_NTASKS ./tapeMatrix
else
    echo "Error: Unknown MODE value. Supported values are 0, 1, 2, 3, 4."
    exit 1
fi