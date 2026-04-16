source genomedef.sh
module load gffread

tooldir=/public/apps/miniconda3/envs/funannotate/opt/pasa-2.5.2/misc_utilities/
gff=out.gff

cat $gff | grep -v '^#' |  awk -F'\t' 'BEGIN{OFS=FS} {if ($3=="CDS")  {$3="exon"; print $0; $3="CDS"; print $0} else {print $0} }' > $gff.fixed.gff
gffread -g $GENOME -V -N -M -d dup.info -K -Q -x out.cds.fa -y out.prot.fa  -o $gff.dedup.gff $gff.fixed.gff

$tooldir/gff3_file_to_proteins.pl $gff.dedup.gff $GENOME CDS > out.cds.fa
