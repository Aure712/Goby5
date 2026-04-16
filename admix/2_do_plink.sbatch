#!/bin/bash
#SBATCH -p long,blade,himem
#SBATCH -c 10
#SBATCH --mem 100G

source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate plink2
Name=Tra

# Generate the input file in plink format
plink2 --vcf all.filtered.snp.vcf.recode.vcf --make-bed --out $Name --allow-extra-chr --const-fid 0

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
awk '{$1="0";print $0}' $Name.bim > $Name.bim.tmp 
mv $Name.bim.tmp $Name.bim


