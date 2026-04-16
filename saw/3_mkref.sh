#!/usr/bin/bash
#SBATCH -p blade
#SBATCH -c 32

saw=/public4/software/saw/saw-8.1.1/bin/saw
ref=/fast3/group_crf/home/g21shaoy23/goby5/annotations/Tba_sort.fa
gff=/fast3/group_crf/home/g21shaoy23/goby5/peaks2utr/utr3k/Tba_clean.utr3k.gff3

$saw makeRef --threads-num $SLURM_CPUS_PER_TASK --mode STAR --fasta $ref --gtf out.gtf --genome ref
