ID=Tridentiger_barbatus_3.0
REF=/fast3/group_crf/home/g21shaoy23/goby5/annotations/Tba_sort.fa
ANN=/fast3/group_crf/home/g21shaoy23/goby5/annotations/Tba_clean_fix.gff3
#cds=/data2/projects/yshao/Goby5RNA/annotations/longest_iso/Tba_longest.cds.fa
#prot=/data2/projects/yshao/Goby5RNA/annotations/longest_iso/Tba_longest.prot.fa

mkdir -p ./data/${ID}

mkdir -p ./data/genomes

ln -sf `realpath $REF` ./data/genomes/${ID}.fa
cat $ANN | gzip -c >  ./data/${ID}/genes.gff.gz
#ln -sf `realpath $cds` ./data/${ID}/cds.fa
#ln -sf `realpath $prot` ./data/${ID}/protein.fa
java -jar snpEff.jar build -gff3 -noCheckProtein -noCheckCds -v Tridentiger_barbatus_3.0
