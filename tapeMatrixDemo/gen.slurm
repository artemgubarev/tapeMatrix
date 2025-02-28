#!/bin/bash
#SBATCH --job-name=gen_data
#SBATCH --output=gen_out.log
#SBATCH --error=gen_err.log
#SBATCH --time=00:10:00

N=100
B=60

module load python/3.8
python3 gen.py $N $B