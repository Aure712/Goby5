MSMC2=/data/software/msmc2-2.1.3/build/release/msmc2

sSampleList="sample.txt"
nSamples=`cat $sSampleList | wc -l `
sIn=formsmc2_in
sOut=msmc2ret
recOverMu=1

mkdir -p $sOut

nSamples=$(( nSamples *2 ))
arrStarts=( `seq 0 2 $(( nSamples-1 ))` ) #haplotype starts
arrEnds=( `seq 1 2 $(( nSamples-1 ))` ) #haplotype ends

sI=""
arrHaplotypes=();

for i in "${!arrStarts[@]}";do
	arrHaplotypes[$i]=${arrStarts[i]}"-"${arrEnds[i]}
done

sI=$( IFS=$','; echo "${arrHaplotypes[*]}" )

echo $sI

sFiles=`realpath ${sIn}/*/formsmc2.multihetsep.txt ` #do not include sex chromosome
sCmd="$MSMC2 -I $sI -o $sOut/Tba -r $recOverMu -i 50 -t 24 $sFiles > /dev/null"
echo $sCmd
eval $sCmd > /dev/null 2>&1 &

wait

