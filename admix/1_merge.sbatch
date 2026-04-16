#!/bin/bash
#SBATCH -p long,blade,himem
#SBATCH -c 30
#SBATCH --mem 200G

module load bcftools
source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate beagle #include tabix

out=all_merge.g.vcf.gz

#bcftools merge -o all.samples.vcf.gz -O  z --threads 30 /public3/group_crf/home/g21shaoy23/goby5/vg/vcf/tra/*.genotyped.g.vcf.gz /public3/group_crf/home/g21shaoy23/goby5/vg/vcf/tba3/Tba2*.genotyped.g.vcf.gz /public3/group_crf/home/g21shaoy23/goby5/vg/vcf/tbi1/*.genotyped.g.vcf.gz /public3/group_crf/home/g21shaoy23/goby5/vg/vcf/tbi2/*.genotyped.g.vcf.gz

#bcftools merge -o all.samples.vcf.gz -O  z --threads 30 /public3/group_crf/home/g21shaoy23/goby5/vg/vcf/tba1/*.genotyped.g.vcf.gz /public3/group_crf/home/g21shaoy23/goby5/vg/vcf/tba2/*.genotyped.g.vcf.gz /public3/group_crf/home/g21shaoy23/goby5/vg/vcf/tba3/Tba1*.genotyped.g.vcf.gz /public3/group_crf/home/g21shaoy23/goby5/vg/vcf/tba3/Tba3*.genotyped.g.vcf.gz /public3/group_crf/home/g21shaoy23/goby5/vg/vcf/tba3/Tba4*.genotyped.g.vcf.gz /public3/group_crf/home/g21shaoy23/goby5/vg/vcf/tbi1/*.genotyped.g.vcf.gz /public3/group_crf/home/g21shaoy23/goby5/vg/vcf/tbi2/*.genotyped.g.vcf.gz


bcftools merge  -o $out -O  z --threads 30 /public3/group_crf/home/g21shaoy23/goby5/vg/vcf/tra/*.genotyped.g.vcf.gz  /public3/group_crf/home/g21shaoy23/goby5/vg/vcf/tba1/*.genotyped.g.vcf.gz /public3/group_crf/home/g21shaoy23/goby5/vg/vcf/tba2/*.genotyped.g.vcf.gz /public3/group_crf/home/g21shaoy23/goby5/vg/vcf/tba3/*.genotyped.g.vcf.gz /public3/group_crf/home/g21shaoy23/goby5/vg/vcf/tbi1/*.genotyped.g.vcf.gz /public3/group_crf/home/g21shaoy23/goby5/vg/vcf/tbi2/*.genotyped.g.vcf.gz
tabix $out
