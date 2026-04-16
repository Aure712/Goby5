#!/usr/bin/bash
#SBATCH -p gpu,blade,himem,hugemem
#SBATCH -c 24

module load bcftools
module load samtools

bcftools concat -O z -o allpops.withoutgroup.norm.vcf.gz -f norm.chromosome.list.txt
tabix allpops.withoutgroup.norm.vcf.gz
