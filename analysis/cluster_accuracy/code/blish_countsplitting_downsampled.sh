#!/bin/bash
#SBATCH --job-name=gbm_clusteracc_blish_countsplit_downsampled    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=philnicol740@gmail.com     # Where to send mail
#SBATCH --mem=256gb                    # Job memory request               # Time limit hrs:min:sec
#SBATCH --output=../logs/blish_downsampled.log   # Standard output and error log
#SBATCH -c 32

module load shared R/4.3.2
module load gcc/9.2.0
Rscript blish_countsplitting_downsampled.R
