genome1=Tba.fa
genome2=Tra.fa
#paftools=/data/software/minimap2-2.20/misc/paftools.js

minimap2 -t 20 -cx asm20 --cs $genome1 $genome2 | sort -k6,6 -k8,8n --parallel=4 > pairwise.aln.paf

cut -f1,3,4,5,6,8,9,12 pairwise.aln.paf | awk '$8>=60 {print}' > pairwise.simp.txt

#$paftools call -q 60 -f $genome1 -L10000 -l1000 pairwise.aln.paf > out.vcf 2>out.vcf.info
