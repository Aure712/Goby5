vcfcaller=/data/software/msmc-tools/vcfAllSiteParser.py

CONVERT=/data/software/msmc-tools/generate_multihetsep.py
BCFTOOLS="bcftools"
sVCFFolder=/data2/projects/yshao/population/msmc/vcf/tba/
sSrcFolder="/data2/projects/yshao/population/msmc/called_separately/tba/formsmc2_in/" #first run vcf2msmc.sh in that folder to produce files before running this script
sVCFSuffix=.bam.genotyped.g.vcf.gz
sChrSuffix="TbaScf_"
sVCFFolder=`realpath $sVCFFolder`

sSampleList=sample.txt

arrSamples=( `cut -d" " -f1 $sSampleList` )

arrChr=( `seq 1 22` )

sOutDir=formsmc2_in/

for nChr in "${arrChr[@]}"; do
         sChr=$sChrSuffix$nChr
         sChrDir=$sOutDir/$sChr
	mkdir -p $sChrDir
	sMaskCmd="";
	sVCFList="";
	for nSample in "${!arrSamples[@]}"; do
		sSample=${arrSamples[$nSample]};
        	sMaskCmd="$sMaskCmd --mask \"$sSrcFolder/$sSample/$sChr/out_mask.bed.gz\""
		sVCFList="$sVCFList \"$sSrcFolder/$sSample/$sChr/out.vcf.gz\""
	done

	#echo $sSample
	#echo $sMaskCmd
	#echo $sVCFList
	sCMD="$CONVERT $sMaskCmd $sVCFList > $sChrDir/formsmc2.multihetsep.txt 2> $sChrDir/formsmc2.multihetsep.log && touch $sChrDir/done.txt"
	echo $sCMD
	( if [ ! -e $sChrDir/done.txt ]; then eval $sCMD; fi; ) 2>$sChrDir/log.txt &
done

wait
