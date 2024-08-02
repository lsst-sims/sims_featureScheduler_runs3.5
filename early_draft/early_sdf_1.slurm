#!/bin/bash

#SBATCH --account=rubin:developers      # Account name
#SBATCH --partition=milano
#
#SBATCH --job-name=one_y_early_draft3.5
#SBATCH --output=output-%j.txt
#SBATCH --error=output-%j.txt
#
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10g
#
#SBATCH --time=0-14:10:00
#SBATCH --chdir= /sdf/group/rubin/u/yoachim/git_repos/sims_featureScheduler_runs3.5/early_draft

## Set up the evironment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate rubin
export OPENBLAS_NUM_THREADS=1

python early_draft_too.py --no_too --survey_length 700
# python early_draft_too.py 