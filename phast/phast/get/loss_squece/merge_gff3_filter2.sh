source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate mixnmatch
module load python3

length=50

#python3 extract_gene_noUTR.py  /fast3/group_crf/home/g21shaoy23/goby5/annotations/cds_phase_fix2/Tra.prechange_longest_fixed.gff3
#cat Tra.prechange_longest_fixed.gff3.gene_noUTR.gff | cut -f1,4,5,9 | cut -f1 -d";" | awk '{print $1, $2, $3, $5}' | sed -e 's/ /\t/g' | sed -e 's/\"//g' | sed -e 's/\t*$//' > gene_noUTR.bed
python3 get_gene_bed.py

sort -k1,1 -k4,4n Traref_validated_TraScf_all.gff > all_Traref.acc.sort.gff
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, "+", $7, $8}' all_Traref.acc.sort.gff > all_Traref.acc.sort.add.gff
bedtools merge -i all_Traref.acc.sort.add.gff > all_Traref.acc.sort.add.gff.bed
awk 'BEGIN {OFS="\t"} {print $1, "PHAST", "4SP", $2, $3, "000", "id", "x"}' all_Traref.acc.sort.add.gff.bed > all_Traref.acc.sort.add.gff.bed.gff

bedtools subtract -a all_Traref.acc.sort.add.gff.bed -b gene_noUTR.bed > Traref.acc.rmgene.bed
#bedtools intersect -a /fast3/group_crf/home/g21shaoy23/goby5/annotations/cds_phase_fix2/Tba.prechange_longest_fixed.gff3 -b all_Tbaref.acc.sort.add.gff.bed > Tbaref.acc.gene.gff3

python3 re_filter_loss_sequence2.py all_Traref.acc.sort.add.gff.bed.gff Traref.acc.rmgene.bed Traref.acc.rmgene.gff $length


#cat NOT_APLref.loss18.gff | awk '{print $1, $4, $5}' | sort -k1,1V -k2,2n -k3,3n | sed -e 's/ /\t/g' | sed -e 's/\"//g' | sed -e 's/\t*$//' > NOT_APLref.loss18.bed
#cat NFZ_APLref.loss18.gff | awk '{print $1, $4, $5}' | sort -k1,1V -k2,2n -k3,3n | sed -e 's/ /\t/g' | sed -e 's/\"//g' | sed -e 's/\t*$//' > NFZ_APLref.loss18.bed
#bedtools intersect -a NFZ_APLref.loss18.bed -b NOT_APLref.loss18.bed > APLref.loss18.bed
#bedtools subtract -a APLref.loss18.bed -b gene_noUTR.bed > APLref.loss18.rmgene.bed
#python re_filter_loss_sequence.py NFZ_APLref.loss18.gff APLref.loss18.rmgene.bed ../APLref.loss18.rmgene.gff $length
