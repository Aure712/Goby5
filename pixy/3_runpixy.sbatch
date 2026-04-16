#!/bin/bash
#SBATCH -p long,blade,himem
#SBATCH -c 15
#SBATCH --mem 150G

source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate pixy
conda activate --stack beagle

pixy --stats pi fst dxy --vcf tratba_merge.g.filter.vcf.gz --populations tratba_samples.txt --window_size 10000 --n_cores 15 --bypass_invariant_check 'yes' --output_folder tratba_results 
