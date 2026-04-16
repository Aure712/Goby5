#!/bin/bash
#SBATCH -p long,blade,himem
#SBATCH -c 10
#SBATCH --mem 100G

source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate sy_py
module load vcftools

bcftools view -O z -o all.filtered.snp.vcf.gz -i 'count(REF)==1 && count(ALT)==1 && QUAL>=60 && TYPE=="snp" && INFO/DP < 476 && FILTER!="LowQual" && F_PASS(GT!="mis")>0.95' ../admix/all_merge.g.vcf.gz

vcftools --gzvcf all.filtered.snp.vcf.gz \
         --minDP 3 \
         --max-missing 0.8 \
         --recode \
         --recode-INFO-all \
         --out all.filtered.snp.vcf

tabix all.filtered.snp.vcf.recode.vcf

bcftools reheader -s ./header_change.txt -o all.filtered.snp.newid.vcf all.filtered.snp.vcf.recode.vcf

bcftools index all.filtered.snp.newid.vcf
