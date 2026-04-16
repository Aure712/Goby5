#!/usr/bin/bash
#SBATCH -p long,himem,blade
#SBATCH -c 16

updated_gff=update.gff3
origin_gff=predict.gff3
fasta=Tba_sort.fa
sp=orig
module load gffread
 
grep -vP "^#" $origin_gff | \
sed 's/FUN_/'$sp'/g' | \
sed 's/novel_/'$sp'_/g' \
> orig.gff3

cat $updated_gff orig.gff3 > combined.gff

gffread -g $fasta --t-adopt --keep-genes -V -F -M -d dup.info -K -w out.mrna.fa -x out.cds.fa -y out.prot.fa -S -o merged.gff combined.gff
