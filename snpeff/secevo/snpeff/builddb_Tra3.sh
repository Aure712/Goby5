ID=Tridentiger_radiatus_3.0
REF=/fast3/group_crf/home/g21shaoy23/goby5/annotations/Tra_sort.fa
ANN=/fast3/group_crf/home/g21shaoy23/goby5/annotations/Tra_clean_fix.gff3
#cds=/fast3/group_crf/home/g21shaoy23/goby5/annotations/cds_phase_fix/Tra_longest.cds.fa
#prot=/fast3/group_crf/home/g21shaoy23/goby5/annotations/cds_phase_fix/Tra_longest.prot.fa

mkdir -p ./data/${ID}

mkdir -p ./data/genomes

ln -sf `realpath $REF` ./data/genomes/${ID}.fa
cat $ANN | gzip -c >  ./data/${ID}/genes.gff.gz
#ln -sf `realpath $cds` ./data/${ID}/cds.fa
#ln -sf `realpath $prot` ./data/${ID}/protein.fa
java -jar snpEff.jar build -gff3 -noCheckProtein -noCheckCds -v Tridentiger_radiatus_3.0
