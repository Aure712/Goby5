zcat joincalled.genotyped.g.vcf.gz | awk -v OFS="\t" '{if(index($0,"#")==0){if(index($8,"DP")){split($8,arr,";");for(i=1; i<=length(arr); i++){split(arr[i],arr2,"=");if(arr2[1]=="DP"){print $1,$2-1,$2,arr2[2]}}}else{print $1,$2-1,$2,0}}}' > all.gvcf.depth2.bed

