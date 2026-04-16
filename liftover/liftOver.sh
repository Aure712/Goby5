liftOver=/data/software/LiftOver/liftOver

#cp ./chain/tba2tra.chain ./renamed.tba2tra.chain 
#sed -i "s/TbaScf_/chr/g" renamed.tba2tra.chain
#sed -i "s/TraScf_/chr/g" renamed.tba2tra.chain

sF=./bed/all.gvcf.depth2.bed
awk -v OFS="\t" '{print $0,$1,$2,$3}' $sF > temp.all.depth2.bed
sed -i "s/^TraScf_/chr/" temp.all.depth2.bed

$liftOver temp.all.depth2.bed renamed.tba2tra.chain all_tba2tra.bed unmapped.txt
