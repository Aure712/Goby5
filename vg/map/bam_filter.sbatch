#!/bin/bash
#SBATCH -p long,blade,himem,hugemem,gpu
#SBATCH -c 10
#SBATCH --mem 100G

module load samtools

samples=( T4 T5 T6 T7 T8 T10 ) #T1 T2 T3 T9 

for sm in "${samples[@]}"; do
  samtools sort -@ 10 -o ${sm}.sorted.bam ${sm}.bam
  samtools index ${sm}.sorted.bam
  samtools view -@ 10 -b ${sm}.sorted.bam TraScf_{1..11} -o ${sm}.filtered.bam
  samtools view -H ${sm}.filtered.bam > ${sm}_header.sam
  grep -E '^@HD|^@SQ.*SN:TraScf_(1[0-1]|[1-9])\b' ${sm}_header.sam > ${sm}_new_header.sam
  cat ${sm}_new_header.sam <(samtools view ${sm}.filtered.bam) > ${sm}.comb.sam
  samtools view -@ 10 -b ${sm}.comb.sam -o ${sm}.final.bam
done
wait;


