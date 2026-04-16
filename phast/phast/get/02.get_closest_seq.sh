source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate mixnmatch

mkdir -p 100kb
mkdir -p 50kb
mkdir -p 20kb
mkdir -p 10kb
mkdir -p 5kb
mkdir -p 1kb
mkdir -p 500bp

for gff3 in phast_results/*.gff; do
	base=`basename $gff3`
	sed -e 's/+/-/g' $gff3 | sort -k1,1 -k4,4n > phast_results/$base.-gff3
	sort -k1,1 -k4,4n $gff3 > phast_results/$base.+gff3
	# -s: same strandedness; -iu: ignore upstream; -D a: signed distance based on A
	sort -k1,1 -k4,4n Tra.prechange_longest_fixed.gff3.+firstCDS.gff3 | bedtools closest -s -a phast_results/$base.+gff3 -b - -D a -iu > phast_results/$base.closest_gene.+gff3
	sort -k1,1 -k4,4n Tra.prechange_longest_fixed.gff3.-firstCDS.gff3 | bedtools closest -s -a phast_results/$base.-gff3 -b - -D a -iu > phast_results/$base.closest_gene.-gff3
	cat phast_results/$base.closest_gene.*gff3 | sort -k1,1V -k4,4n > phast_results/$base.closest_gene.merged.gff3
	python3 filter_nearest_CDS.py phast_results/$base.closest_gene.merged.gff3 100000 100kb/$base.closest_gene.gff3.filtered.gff3
	python3 filter_nearest_CDS.py phast_results/$base.closest_gene.merged.gff3 50000 50kb/$base.closest_gene.gff3.filtered.gff3
	python3 filter_nearest_CDS.py phast_results/$base.closest_gene.merged.gff3 20000 20kb/$base.closest_gene.gff3.filtered.gff3
	python3 filter_nearest_CDS.py phast_results/$base.closest_gene.merged.gff3 10000 10kb/$base.closest_gene.gff3.filtered.gff3
	python3 filter_nearest_CDS.py phast_results/$base.closest_gene.merged.gff3 5000 5kb/$base.closest_gene.gff3.filtered.gff3
	python3 filter_nearest_CDS.py phast_results/$base.closest_gene.merged.gff3 1000 1kb/$base.closest_gene.gff3.filtered.gff3
	python3 filter_nearest_CDS.py phast_results/$base.closest_gene.merged.gff3 500 500bp/$base.closest_gene.gff3.filtered.gff3
	sort -k1,1V -k4,4n 100kb/$base.closest_gene.gff3.filtered.gff3 > 100kb/$base.closest_gene.final.gff3
	sort -k1,1V -k4,4n 50kb/$base.closest_gene.gff3.filtered.gff3 > 50kb/$base.closest_gene.final.gff3
	sort -k1,1V -k4,4n 20kb/$base.closest_gene.gff3.filtered.gff3 > 20kb/$base.closest_gene.final.gff3
	sort -k1,1V -k4,4n 10kb/$base.closest_gene.gff3.filtered.gff3 > 10kb/$base.closest_gene.final.gff3
	sort -k1,1V -k4,4n 5kb/$base.closest_gene.gff3.filtered.gff3 > 5kb/$base.closest_gene.final.gff3
	sort -k1,1V -k4,4n 1kb/$base.closest_gene.gff3.filtered.gff3 > 1kb/$base.closest_gene.final.gff3
	sort -k1,1V -k4,4n 500bp/$base.closest_gene.gff3.filtered.gff3 > 500bp/$base.closest_gene.final.gff3

done
