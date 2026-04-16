#!/bin/bash
#do gatk.sh for each sample.bam file
for sSample in /data2/projects/yshao/population/GATK/bwa/tba3/*.bam; do
	#echo $sSample
        sBase=`basename $sSample`
        sName=${sBase/.bam/}
	#echo $sName
        cd /data2/projects/yshao/population/GATK/GATK/tba3/$sName
        source /data2/projects/yshao/population/GATK/GATK/tba3/gatk.sh $sSample &
done;
wait
