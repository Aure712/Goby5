#!/bin/bash
#SBATCH -p long,blade,himem,hugemem,gpu
#SBATCH -c 50
#SBATCH --mem 220G

#do gatk.sh for each sample.bam file
module load bwa
module load java
module load samtools
for sSample in /public3/group_crf/home/g21shaoy23/goby5/vg/vcf/tba1/*.bam; do
	#echo $sSample
        sBase=`basename $sSample`
        sName=${sBase/.bam/}
	#echo $sName
        cd /public3/group_crf/home/g21shaoy23/goby5/vg/vcf/gatktba1/$sName
        source /public3/group_crf/home/g21shaoy23/goby5/vg/vcf/gatktba1/gatk.sh $sSample &
done;
wait
