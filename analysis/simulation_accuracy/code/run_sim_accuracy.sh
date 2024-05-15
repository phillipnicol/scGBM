#!/bin/bash
#SBATCH --job-name=scGBM_sim_accuracy    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=philnicol740@gmail.com     # Where to send mail
#SBATCH --mem=64gb                    # Job memory request               # Time limit hrs:min:sec
#SBATCH --output=../logs/sim_accuracy.log   # Standard output and error log
#SBATCH -c 16

module load gcc/9.2.0
module load R/4.3.2b
Rscript sim_accuracy.R
