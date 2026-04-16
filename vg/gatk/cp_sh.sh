#copy gatk.sh for each sample.bam and put it in the directory named by sample's name
for sR1 in /public3/group_crf/home/g21shaoy23/goby5/vg/vcf/tba1/*.bam; do
	sDir=`dirname $sR1`
	sBase=`basename $sR1`
	sSample=${sBase/.bam/}
	#echo $sSample
  mkdir $sSample
done
