#!/bin/bash
#SBATCH -p long,blade,himem,hugemem,gpu
#SBATCH -c 30
#SBATCH --mem 200G

module load samtools

Pop=allpops

snpeffjar=/public4/group_crf/home/g21shaoy23/snpeff2/snpEff/snpEff.jar
REF=Tridentiger_barbatus_3.0

java -Xmx48g -jar $snpeffjar $REF /public4/group_crf/shareddata/g21shaoy23_data2/population/snpeff/gvcf_tba/snpeff/allpops.withoutgroup.norm.vcf.gz | bgzip -c > $Pop.merged.snpeff.vcf.gz
