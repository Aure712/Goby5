#!/usr/bin/bash
#SBATCH -p long,himem,blade
#SBATCH -c 20
#SBATCH --mem=20G
source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate cd-hit
cd-hit-est -i comb_tmi.fa -o tba_comb_tmi_cd.fa -c 0.95 -T 20 -M 20000 
