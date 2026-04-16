#step1:
cp allpops_on_MOP2MHK.withOriginName.cds.depth.bed temp.allpops_on_MOP2MHK.withOriginName.cds.depth.bed
sed -i "s/^chr/mhkscf_/" temp.allpops_on_MOP2MHK.withOriginName.cds.depth.bed
awk -v OFS="\t" '{print $1,$2,$3,$5,$7-1,$7}' temp.allpops_on_MOP2MHK.withOriginName.cds.depth.bed > allpops_on_MOP2MHK.renamed.cds.bed

MHK=/data2/projects/zwang/macropodus_compare/genetic_load/all_MHK_MOP/filter_DP/refMHK/allpops.cds.depth.dedup.bed
MOP=/data2/projects/zwang/macropodus_compare/genetic_load/all_MHK_MOP/filter_DP/refMOP/allpops.cds.depth.dedup.bed

bedtools intersect -a allpops_on_MOP2MHK.renamed.cds.bed -b $MHK -wa -wb > allpops.cds.singleDP_MHK.bed

awk -v OFS="\t" '{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10}' allpops.cds.singleDP_MHK.bed > temp.allpops.cds.singleDP_MHK.bed

bedtools intersect -a temp.allpops.cds.singleDP_MHK.bed -b $MOP -wa -wb > allpops.cds.bothDP.bed

#step2:plot hist and find 2.5%~97.5% of difference of DP
