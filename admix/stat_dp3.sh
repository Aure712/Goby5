#!/usr/bin/bash
#SBATCH -p gpu,himem,hugemem,long,blade
#SBATCH -c 10
#SBATCH --mem=220G

source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate sy_py
module load bcftools

bash analyze_depth_distribution4.sh all_merge.g.vcf.gz 


