#!/bin/bash
#SBATCH -p long,blade,himem,gpu
#SBATCH -c 5
#SBATCH --mem 100G

source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate bgzip
conda activate --stack beagle
module load vcftools
module load bcftools

#bcftools view -O z  -v "snps"  -e 'SOR>2 || QD<5 || MQRankSum<-12.5 || INFO/DP>462 || INFO/DP<66  || ExcessHet>30' joincalled.genotyped.g.vcf.gz > joincalled.genotyped.g.filtered.snps.vcf.gz

##filter
vcftools --gzvcf all.samples.forpixy.newid.g.vcf.gz \
 	--remove-indels \
	--max-missing 0.8 \
	--recode \
  --recode-INFO-all \
  --stdout > tratba_merge.g.filter.vcf;
#--maf 0.05 \

bgzip -f tratba_merge.g.filter.vcf;
tabix tratba_merge.g.filter.vcf.gz;

#bcftools reheader -s header_change.txt -o joincalled.genotyped.forpixy.newid.g.vcf.gz joincalled.genotyped.forpixy.g.vcf.gz