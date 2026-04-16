#!/bin/bash
#SBATCH -p long,himem,gpu
#SBATCH --mem=200G
#SBATCH -c 12

source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate cactus2.6.9

#export PATH=$PATH:/public/software/conda_envs/cactus2.6.9/bin
export TOIL_SLURM_ARGS="-p blade,gpu,himem,hugemem --nice=5000"
export TMPDIR=./tmp/

cactus-hal2maf ./js goby24.hal goby24tba.maf.gz \
	--batchCores 12 --batchMemory 50G --batchCount 12 --batchParallelTaf 12 \
	--refGenome Tridentiger_barbatus --chunkSize 1000000 --filterGapCausingDupes --noAncestors --dupeMode single \
	--maxCores 12 --defaultCores 1 --batchSystem slurm --binariesMode local \
	--retryCount=3 --maxPreemptibleServiceJobs 32 --maxServiceJobs 32 --maxNodes 5 --maxJobs 32 --maxLocalJobs 32 \
       	--workDir ./tmp --realTimeLogging \
#	--consMemory 250G

# --configFile config-slurm.xml



