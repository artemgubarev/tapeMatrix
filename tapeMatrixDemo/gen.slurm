#!/bin/bash
#SBATCH --job-name=gen_data
#SBATCH --output=gen_out.log
#SBATCH --error=gen_err.log
#SBATCH --time=00:10:00

N=2000
B=1200

module load python/3.8
python3 generator.py $N $B