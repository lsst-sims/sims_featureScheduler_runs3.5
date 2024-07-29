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
#SBATCH --ntasks-per-node=8

## Walltime (hours:min:sec) Do not specify a walltime substantially more than your job needs.

#SBATCH --time=36:00:00

## Memory per node. It is important to specify the memory since the default memory is very small.

## For mox, --mem may be more than 100G depending on the memory of your nodes.

## For ikt, --mem may be 58G or more depending on the memory of your nodes.

## See above section on "Specifying memory" for choices for --mem.

#SBATCH --mem=58G

## Specify the working directory for this job

#SBATCH --chdir=/mmfs1/home/yoachim/runs35/early_draft

##turn on e-mail notification

#SBATCH --mail-type=ALL

#SBATCH --mail-user=yoachim@uw.edu

## export all your environment variables to the batch job session

#SBATCH --export=all

## Set up the evironment
# source ~/anaconda3/etc/profile.d/conda.sh
conda activate rubin
export OPENBLAS_NUM_THREADS=1

cd /mmfs1/home/yoachim/runs35/early_draft

## run all the baseline commands in parallel
module load parallel-20170722

python early_draft.py --no_too

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
