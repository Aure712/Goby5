#!/usr/bin/bash
#SBATCH -p blade,himem,hugemem,long
#SBATCH -c 40
#SBATCH --mem=180G

source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate /public4/software/conda_env/peaks2utr

peaks2utr --max-distance 3000 --override-utr -p 40 -o Tba_clean.utr3k.gff3 /fast3/group_crf/home/g21shaoy23/goby5/annotations/Tba_clean.gff3 ../bam/merged_sorted.bam







