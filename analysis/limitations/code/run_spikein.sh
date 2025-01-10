#!/bin/bash
#SBATCH --job-name=scGBM_spikein   # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=philnicol740@gmail.com     # Where to send mail
#SBATCH --mem=32gb                    # Job memory request               # Time limit hrs:min:sec
#SBATCH --output=../logs/spikein.log   # Standard output and error log
#SBATCH -c 4

module load shared R/4.3.2
module load gcc/9.2.0
Rscript spike_in.R
