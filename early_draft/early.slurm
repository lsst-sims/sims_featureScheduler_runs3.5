#!/bin/bash

## Job Name

#SBATCH --job-name=early3.5

## Allocation Definition

## On mox and ikt, the account and partition options should be the same.
#SBATCH --account=astro
#SBATCH --partition=compute-bigmem

## Resources

## Nodes

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

## Walltime (hours:min:sec) Do not specify a walltime substantially more than your job needs.

#SBATCH --time=24:00:00

## Memory per node. It is important to specify the memory since the default memory is very small.

## For mox, --mem may be more than 100G depending on the memory of your nodes.

## For ikt, --mem may be 58G or more depending on the memory of your nodes.

## See above section on "Specifying memory" for choices for --mem.

#SBATCH --mem=58G

## Specify the working directory for this job

#SBATCH --chdir=/gscratch/dirac/yoachim/sims_featureScheduler_runs3.5/early_draft

##turn on e-mail notification

#SBATCH --mail-type=ALL

#SBATCH --mail-user=yoachim@uw.edu

## export all your environment variables to the batch job session

#SBATCH --export=all

## Set up the evironment
source /gscratch/astro/yoachim/miniconda3/etc/profile.d/conda.sh
conda activate /gscratch/astro/yoachim/miniconda3/envs/rubin/envs/rubin2

export OPENBLAS_NUM_THREADS=1

cd /gscratch/dirac/yoachim/sims_featureScheduler_runs3.5/early_draft

## run all the baseline commands in parallel
module load gnu/parallel/20210422

# python early_draft_too.py --no_too
python early_draft_too.py
# 
# rm maf.sh
# ls *10yrs.db | xargs -I'{}' echo "glance_dir --db '{}'" > maf.sh
# ls *10yrs.db | xargs -I'{}' echo "scimaf_dir --db '{}'" >> maf.sh
# ls *10yrs.db | xargs -I'{}' echo "ddf_dir --db '{}'" >> maf.sh
# ls *10yrs.db | xargs -I'{}' echo "metadata_dir --db '{}'" >> maf.sh
# 
# generate_ss 
# cat ss_script.sh >> maf.sh
# 
# cat maf.sh | parallel -j 8
