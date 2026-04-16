#!/usr/bin/bash
#SBATCH -p gpu,himem,hugemem,blade
#SBATCH -c 10

module load bcftools
module load samtools


#bcftools view -O z -M 2 -m 2 -i 'TYPE=="snp" || TYPE=="indel"' allpops.withoutgroup.newid.g.vcf.gz > allpops.withoutgroup.vcf.gz 
tabix allpops.withoutgroup.vcf.gz

