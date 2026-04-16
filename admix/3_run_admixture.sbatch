#!/bin/bash
#SBATCH -p long,blade,himem,hugemem,gpu
#SBATCH -c 10
#SBATCH --mem 200G

source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate ADMIXTURE

Name=Tra

#for K in 2 3 4 5 6 7 8; do
for K in 2 3 4 5 6; do
	admixture -j30 --cv $Name.bed $K | tee log${K}.out
done
