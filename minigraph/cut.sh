awk 'BEGIN {count=0; output=1} 
     /^>/ {count++; if (count>525) output=0} 
     output {print}' tra.gfa.fa > tra.gfa.cut.fa
