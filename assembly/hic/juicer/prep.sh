genome=/data/projects/yshao/hic_raw/ref/Tra/Tra_ref.fasta

mkdir -p references
cp $genome references/ref.fa
bwa index references/ref.fa
fastahack -i references/ref.fa

mkdir -p restriction_sites
cut -f1,2 references/ref.fa.fai > restriction_sites/ref.chrom.sizes

cd restriction_sites
python /data/software/juicer/misc/generate_site_positions.py DpnII ref ../references/ref.fa


