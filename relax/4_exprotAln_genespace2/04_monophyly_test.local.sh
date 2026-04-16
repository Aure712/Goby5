#!/bin/bash

#SET THIS TO THE CORRECT NUMBER
TOTALPARTS=32
END=$(( $TOTALPARTS - 1 ))
mkdir -p genetrees_improved
mkdir -p logs

for i in $( seq 0 $END ); do

php monophyly_test.php -N ${TOTALPARTS} -f ${i} > logs/monophylytest.$i.of.${TOTALPARTS}.log 2>&1  &

done
wait
exit

