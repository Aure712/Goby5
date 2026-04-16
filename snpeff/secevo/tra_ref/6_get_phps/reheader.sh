module load bcftools
module load samtools
bcftools reheader -s header.txt ../allpops.merged.snpeff.vcf.gz -o allpops.merged.snpeff.renameoutgroup.vcf.gz
tabix allpops.merged.snpeff.renameoutgroup.vcf.gz
