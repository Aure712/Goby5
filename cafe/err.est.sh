CAFE=/data/software/CAFE5/bin/cafe5
CPU=80

$CAFE -c $CPU -t tree.tre -i orthogroups.genecount.txt -p -e > err.est.log 2>&1 
