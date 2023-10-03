#!/bin/sh
#SBATCH --job-name=water-analysis
#SBATCH -o /home/jji110/parch_platform/slurm-outputs/wtanalysis-%j.out
#SBATCH --mail-user jji110@syr.edu
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
srun python water_analysis_3.py $1