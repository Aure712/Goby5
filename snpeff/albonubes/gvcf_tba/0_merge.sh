bcftools merge -o allpops.withoutgroup.g.vcf.gz -O z --threads 32 *.gz
tabix allpops.withoutgroup.g.vcf.gz
bcftools reheader -s header_change.txt -o allpops.withoutgroup.newid.g.vcf.gz allpops.withoutgroup.g.vcf.gz
tabix allpops.withoutgroup.newid.g.vcf.gz