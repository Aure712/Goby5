#!/bin/bash
#SBATCH -p long,blade,himem
#SBATCH -c 10
#SBATCH --mem 100G

module load bcftools

#bcftools view -O z  -v "snps"  -e 'SOR>2 || QD<5 || MQRankSum<-12.5 || INFO/DP>250 || INFO/DP<46  || ExcessHet>30' joincalled.genotyped.g.vcf.gz > joincalled.genotyped.g.filtered.snps.vcf.gz

bcftools view -O z -e 'INFO/DP>286 || INFO/DP<128 || INFO/DP="."' /public4/group_crf/home/g21shaoy23/CEGA2/tbatra/all.samples.vcf.gz > all.samples.forpixy.g.vcf.gz

bcftools reheader -s ./header_change.txt -o all.samples.forpixy.newid.g.vcf.gz all.samples.forpixy.g.vcf.gz