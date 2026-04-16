for i in run*/*/Gamma_results.txt;do
	k=`dirname $i`
	echo $k >> stats.txt
	cat $i | head -1 >> stats.txt
done
