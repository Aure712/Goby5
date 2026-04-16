#!/bin/bash
#SBATCH -p long,blade,himem,hugemem,gpu
#SBATCH -c 30
#SBATCH --mem 200G

source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate vg

samples=( T1 T2 T3 T4 T5 T6 T7 T8 T9 T10 )

for sm in "${samples[@]}"; do
  vg map -f /fast3/group_crf/home/g21shaoy23/goby5/vg/fq/T_unknow/${sm}.paired_1.fq.gz /fast3/group_crf/home/g21shaoy23/goby5/vg/fq/T_unknow/${sm}.paired_2.fq.gz -x ../tra.xg -g ../tra.gcsa -t 30 --surject-to bam > ${sm}.bam
done
wait;


