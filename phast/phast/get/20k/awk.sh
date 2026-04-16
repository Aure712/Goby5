awk -F'\t' '{match($(NF-1), /Parent=([^;]+)/, a); if(a[1]!="") print a[1]}' Traref.acc.rmgene.gff.closest_gene.final.gff3 \
| sort -u > parent_unique.txt
