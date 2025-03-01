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
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# 0 - C & serial
# 1 - C & pthreads
# 2 - C & openMP
# 3 - C & MPI
# 4 - python & MPI
export MODE=1

if [ "$MODE" -eq 4 ]; then
    mpiexec -n $SLURM_NTASKS tapeMatrixPy.py
elif
    mpicxx -std=c++11 -fopenmp -pthread -o tapeMatrix tapeMatrixDemo.cpp
    if [ "$MODE" -eq 3 ]; then
        mpirun -n $SLURM_NTASKS ./tapeMatrix
    elif [ "$MODE" -eq 0 ] || [ "$MODE" -eq 1 ] || [ "$MODE" -eq 2 ]; then
        srun ./tapeMatrix
    else
        echo "Error: Unknown MODE value. Supported values are 0, 1, 2, 3, 4."
        exit 1
    fi
fi