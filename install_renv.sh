#!/bin/bash
#SBATCH --job-name=renv_restore      # Job name
#SBATCH --output=renv_restore.log    # Output log file
#SBATCH --error=renv_restore.err     # Error log file
#SBATCH --time=08:00:00              # Time limit (hh:mm:ss)
#SBATCH --ntasks=1                   # Number of tasks
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --mem=4GB                    # Memory per node

# Load necessary modules (adjust this based on your HPC environment)
module load shared R/4.3.2   # Replace with the version of R available on your cluster
#module unload gcc/9.2.0
module load gcc/13.2.0


export PKG_CONFIG_PATH=/cm/shared/apps/curl/lib/pkgconfig:$PKG_CONFIG_PATH
export CFLAGS="-I/cm/shared/apps/curl/include"
export LDFLAGS="-L/cm/shared/apps/curl/lib"

export CC=/path/to/gcc-13.2.0/bin/gcc
export CXX=/path/to/gcc-13.2.0/bin/g++

# Navigate to the directory containing the R project with renv

# Run renv::restore() in R
Rscript -e "renv::rebuild('RcppAnnoy')"
Rscript -e "renv::restore(prompt=FALSE)"
