
CAFE=/data/software/CAFE5/bin/cafe5
CPU=16
k=4

for i in `seq 1 20`; do
	mkdir -p runresults_k$k/$i
$CAFE -c $CPU -t tree.tre -i orthogroups.genecount.txt -p -eresults/Base_error_model.txt -k $k -o  runresults_k$k/$i >  runresults_k$k/$i/est.log 2>&1 
done
