#!/bin/sh
#SBATCH --job-name=annealing
#SBATCH -o ./slurm-outputs/ann-%j.out
#SBATCH --mail-user jji110@syr.edu
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
srun python pp_refill_annealing_2.py $1 $2