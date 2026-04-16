source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate mixnmatch

#bedtools intersect -a tragene.bed -b ismc.bed > ismc_tra.bed
#bedtools intersect -a Tra_ismc_all_gene.bed -b ismc.bed > ismc_all_tra.bed

grep HIGH indel.tra.bed > indel.tra.high.bed
grep HIGH snp.tra.bed > snp.tra.high.bed

grep synonymous_variant indel.tra.bed > indel.tra.syn.bed
grep synonymous_variant snp.tra.bed > snp.tra.syn.bed

#cat indel.tra.high.bed snp.tra.high.bed | sort -k1,1 -k2,2n | bedtools merge > tra.high.bed
#cat allpops.refTra.NONSYN.polarized.out.filtered.final.high.bed genomeload_indel-snp.bed > allpops.refTra.all.high.bed
bedtools subtract -a indel.tra.syn.bed -b snp.tra.syn.bed > indel-snp_syn.bed
cat snp.tra.syn.bed indel-snp_syn.bed > allpops.refTra.all.syn.bed

bedtools subtract -a indel.tra.high.bed -b snp.tra.high.bed > indel-snp_high.bed
cat snp.tra.high.bed indel-snp_high.bed > allpops.refTra.all.high.bed