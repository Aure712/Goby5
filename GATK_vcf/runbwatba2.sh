fa=Tba_sort.fa
samples=(Tba1 Tba2 Tba3 Tba4 )

for sm in "${samples[@]}"; do
  ( bwa mem -t 10 -R "@RG\tID:Tba\tPL:illumina\tSM:${sm}" $fa /data2/projects/yshao/population/trimmed_illumina/Tba/${sm}.paired_1.fq.gz /data2/projects/yshao/population/trimmed_illumina/Tba/${sm}.paired_2.fq.gz | samtools view -Sb > ${sm}.bam ) > ${sm}bwa.log 2>&1 &
done
wait;