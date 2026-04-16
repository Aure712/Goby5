#!/usr/bin/bash
#SBATCH -p long,himem,blade
#SBATCH -c 32
#SBATCH --mem=40G

RNA_Evd=../tba_comb_tmi_cd.fa
REF=../Tba_sort.fa
CPU=32
SP=Tba

module load minimap2
module load stringtie
module load samtools

minimap2 -ax splice:hq -t $CPU -u f $REF $RNA_Evd | samtools sort --reference $REF -@ 12 -o ${SP}_aln.bam - > aln.log 2>&1
stringtie -p $CPU -a 5 -E 50 -M 1 -f 0.001 -m 33 -j 0.001 -s 0.001 -c 0.001 -L -t -u -o ${SP}_alnout.gtf ${SP}_aln.bam > stringtie.log 2>&1
#perl gtf_genome_to_cdna_fasta.pl ${SP}_alnout.gtf $REF > ${SP}_aln_trans.fa

