awk 'BEGIN{OFS="\t"}
NR==1{next}  # 跳过表头
$1=="Tba"{print $2, $3-1, $4, $5 > "Tba.refTra.pi.bed"}
$1=="Tra"{print $2, $3-1, $4, $5 > "Tra.refTra.pi.bed"}
' pixy_pi.txt
