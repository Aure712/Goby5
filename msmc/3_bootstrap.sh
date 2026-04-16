bootstrap=/data/software/msmc-tools/multihetsep_bootstrap.py
nReps=30
nChunkSize=5000000


arrPops=( Tba )

sOutDir=bootstrapped

sSrcDir=formsmc2_in

mkdir -p $sOutDir

for sPop in "${arrPops[@]}"; do
	sFiles="$sSrcDir/*/formsmc2.multihetsep.txt";
	sBSFolder="$sOutDir/";
	mkdir -p $sBSFolder
	$bootstrap -n $nReps -s $nChunkSize --nr_chromosomes 1 --chunks_per_chromosome 100 --seed 33333 "$sBSFolder" $sFiles > $sBSFolder/bs.log 2>&1 &
done

wait

