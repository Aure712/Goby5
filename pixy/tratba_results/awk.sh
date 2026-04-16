awk -v OFS="\t" 'NR > 1 {print $3, $4-1, $5, $6}' pixy_dxy.txt > Tba.Tra.dxy.bed
awk -v OFS="\t" 'NR > 1 {print $3, $4-1, $5, $6}' pixy_fst.txt > Tba.Tra.fst.bed
awk -v OFS="\t" 'NR > 1 {print $3, $4-1, $5, $6}' pixy_pi.txt > Tba.Tra.pi.bed