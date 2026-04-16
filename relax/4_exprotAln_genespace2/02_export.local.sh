#!/bin/bash

#SET THIS TO THE CORRECT NUMBER
TOTALPARTS=24
END=$(( $TOTALPARTS - 1 ))

mkdir -p output
mkdir -p logs

for i in $( seq 0 $END ); do 

php 02_exportAln_multiref_withorigseq.php -N ${TOTALPARTS} -f ${i} > logs/export.${i}.log 2>&1 &

done
wait
exit

