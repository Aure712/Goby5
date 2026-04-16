bcftools view -H /data2/projects/yshao/population/GATK/bwa/tra/refs/R-for.bam.genotyped.g.vcf.gz | awk '{split($10,arr,":");if(arr[3]>0){print $1"\t"$2"\t"$2}}' > Rfor.mask.bed
bedtools merge -i Rfor.mask.bed > pairwise.simp.ref.Rhinogobius_formosanus.txt

bcftools view -H /data2/projects/yshao/population/GATK/bwa/tra/refs/Tbifasciatus-muscle2.bam.genotyped.g.vcf.gz | awk '{split($10,arr,":");if(arr[3]>0){print $1"\t"$2"\t"$2}}' > Tbi.mask.bed
bedtools merge -i Tbi.mask.bed > pairwise.simp.ref.Tridentiger_bifasciatus.txt

