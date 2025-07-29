#!/bin/bash
#SBATCH --job-name=gbm_projaccuracy_zhengmix    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=philnicol740@gmail.com     # Where to send mail
#SBATCH --mem=128gb                    # Job memory request               # Time limit hrs:min:sec
#SBATCH --output=../logs/zhengmix.log   # Standard output and error log
#SBATCH -c 32

module load R/4.3.2
module load gcc/9.2.0
Rscript zhengmix.R
