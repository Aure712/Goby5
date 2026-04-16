#SV_Tba=/public4/group_crf/home/g21shaoy23/snpeff2_reftba/6_get_phps/6_SV_chrom_nopairwise/refTba.SV.w.indelmodifier.polarized.out.txt
#SV_Tra=/public4/group_crf/home/g21shaoy23/snpeff2/6_get_phps/6_SV_chrom_nopairwise/refTra.SV.w.indelmodifier.polarized.out.txt
SV_Tra=../refTra.dedup.SV.w.indelmodifier.polarized.out.txt
SV_Tba=../refTba.dedup.SV.w.indelmodifier.polarized.out.txt
sPop=allpops

#step1: stats AF with statsAF.R

#step2:
#awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$10,$11,$14}' /public4/group_crf/shareddata/g21shaoy23_data2/population/genetic_load/liftover/$sPop.cds.bothDP.bed > temp.$sPop.cds.bothDP.bed
#bedtools intersect -a temp.$sPop.cds.bothDP.bed -b $sPop.refTra.SV.withAF.bed -wa -wb > $sPop.cds.bothDP.TraAF.bed
#awk -v OFS='\t' '{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12,$13,$14}' $sPop.cds.bothDP.TraAF.bed > temp.$sPop.cds.bothDP.TraAF.bed
#bedtools intersect -a temp.$sPop.cds.bothDP.TraAF.bed -b $sPop.refTba.SV.withAF.bed -wa -wb> $sPop.cds.bothDP.bothAF.bed

#step3:filter DP and SNP with filter.R

#step4:
awk -F"\t" -v OFS='\t' '{ start=$2-1;  printf("%s",$1"\t"start"\t"$2"\t");$1="";$2="";print}' $SV_Tra > $sPop.refTra.SV.polarized.out.bed
awk -F"\t" -v OFS='\t' '{ start=$2-1;  printf("%s",$1"\t"start"\t"$2"\t");$1="";$2="";print}' $SV_Tba > $sPop.refTba.SV.polarized.out.bed
bedtools intersect -a $sPop.refTra.SV.polarized.out.bed -b $sPop.refTra.SV.filtered.bed > $sPop.refTra.SV.polarized.out.filtered.bed
bedtools intersect -a $sPop.refTba.SV.polarized.out.bed -b $sPop.refTba.SV.filtered.bed > $sPop.refTba.SV.polarized.out.filtered.bed

#step5:
bedtools intersect -a $sPop.refTra.SV.polarized.out.filtered.bed -b /public4/group_crf/shareddata/g21shaoy23_data2/population/snpeff/gvcf/snpeff/5_pairwise/pairwise_intersect_refTra.bed > $sPop.refTra.SV.polarized.out.filtered.final.bed
bedtools intersect -a $sPop.refTba.SV.polarized.out.filtered.bed -b /public4/group_crf/shareddata/g21shaoy23_data2/population/snpeff/gvcf_tba/snpeff/5_pairwise/pairwise_intersect_refTba.bed > $sPop.refTba.SV.polarized.out.filtered.final.bed