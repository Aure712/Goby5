#!/bin/bash
#SBATCH -p long,himem,gpu,hugemem
#SBATCH --mem=220G
#SBATCH -c 12

source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate cactus2.6.9

#export PATH=$PATH:/public/software/conda_envs/cactus2.6.9/bin
export TOIL_SLURM_ARGS="--nice=5000"
export TMPDIR=./tmp/

cactus --maxCores 12 --defaultCores 1 --batchSystem slurm --binariesMode local ./js ./cactus_input.txt ./goby24.hal \
        --retryCount=3 --maxPreemptibleServiceJobs 32 --maxServiceJobs 32 --maxNodes 5 --maxJobs 32 --maxLocalJobs 32 \
        --workDir ./tmp --realTimeLogging --maxMemory 250G 
#	--restart
#cactus --maxCores 95 --defaultCores 1 --batchSystem slurm --binariesMode local ./js ./cactus_input.txt ./Cyprinodontiformes.hal \
#	--workDir ./tmp --realTimeLogging
# --configFile config-slurm.xml

#cactus-hal2maf ./js Cyprinodontiformes.hal Cyprinodontiformes.maf.gz --refGenome Aplocheilus_lineatus --chunkSize 1000000 --onlyOrthologs --noAncestors --unique --noDupes
