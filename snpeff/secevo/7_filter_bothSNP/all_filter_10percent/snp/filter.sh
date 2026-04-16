#NONSYN_Tba=/public4/group_crf/home/g21shaoy23/snpeff2_reftba/6_get_phps/6_NONSYN_chrom/refTba.NONSYN.polarized.out.txt
#NONSYN_Tra=/public4/group_crf/home/g21shaoy23/snpeff2/6_get_phps/6_NONSYN_chrom/refTra.NONSYN.polarized.out.txt
NONSYN_Tba=../refTba.dedup.NONSYN.polarized.out.txt
NONSYN_Tra=../refTra.dedup.NONSYN.polarized.out.txt
sPop=allpops

#step1: stats AF with statsAF.R

#step2:
#awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$10,$11,$14}' /public4/group_crf/shareddata/g21shaoy23_data2/population/genetic_load/liftover/$sPop.cds.bothDP.bed > temp.$sPop.cds.bothDP.bed
#bedtools intersect -a temp.$sPop.cds.bothDP.bed -b $sPop.refTra.NONSYN.withAF.bed -wa -wb > $sPop.cds.bothDP.TraAF.bed
#awk -v OFS='\t' '{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12,$13,$14}' $sPop.cds.bothDP.TraAF.bed > temp.$sPop.cds.bothDP.TraAF.bed
#bedtools intersect -a temp.$sPop.cds.bothDP.TraAF.bed -b $sPop.refTba.NONSYN.withAF.bed -wa -wb> $sPop.cds.bothDP.bothAF.bed

#step3:filter DP and SNP with filter.R

#step4:
awk -F"\t" -v OFS='\t' '{ start=$2-1;  printf("%s",$1"\t"start"\t"$2"\t");$1="";$2="";print}' $NONSYN_Tba > $sPop.refTba.NONSYN.polarized.out.bed
awk -F"\t" -v OFS='\t' '{ start=$2-1;  printf("%s",$1"\t"start"\t"$2"\t");$1="";$2="";print}' $NONSYN_Tra > $sPop.refTra.NONSYN.polarized.out.bed
bedtools intersect -a $sPop.refTba.NONSYN.polarized.out.bed -b $sPop.refTba.NONSYN.filtered.bed > $sPop.refTba.NONSYN.polarized.out.filtered.bed
bedtools intersect -a $sPop.refTra.NONSYN.polarized.out.bed -b $sPop.refTra.NONSYN.filtered.bed > $sPop.refTra.NONSYN.polarized.out.filtered.bed

#step5:
bedtools intersect -a $sPop.refTra.NONSYN.polarized.out.filtered.bed -b /public4/group_crf/shareddata/g21shaoy23_data2/population/snpeff/gvcf/snpeff/5_pairwise/pairwise_intersect_refTra.bed > $sPop.refTra.NONSYN.polarized.out.filtered.final.bed
bedtools intersect -a $sPop.refTba.NONSYN.polarized.out.filtered.bed -b /public4/group_crf/shareddata/g21shaoy23_data2/population/snpeff/gvcf_tba/snpeff/5_pairwise/pairwise_intersect_refTba.bed > $sPop.refTba.NONSYN.polarized.out.filtered.final.bed
