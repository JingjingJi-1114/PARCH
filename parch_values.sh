#!/bin/sh
#SBATCH --job-name=avg-pv
#SBATCH -o ./slurm-outputs/avg-pv-%j.out
#SBATCH --mail-user jji110@syr.edu
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
srun  python ave_pv_5.py $1