#!/bin/bash

#SBATCH --account=rubin:developers      # Account name
#SBATCH --partition=milano
#
#SBATCH --job-name=weather3.5
#SBATCH --output=output-%j.txt
#SBATCH --error=output-%j.txt
#
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=20g
#
#SBATCH --time=2-14:10:00
#SBATCH --chdir=/sdf/data/rubin/shared/fbs_sims/sims_featureScheduler_runs3.5/weather

## Set up the evironment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate rubin
export OPENBLAS_NUM_THREADS=1

cat weather.sh | parallel -j 8
