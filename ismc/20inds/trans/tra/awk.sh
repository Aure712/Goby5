awk 'BEGIN{OFS="\t"} { $8=""; print $0 }' output2.bed | sed 's/\t\t/\t/g' > new_output2.bed
