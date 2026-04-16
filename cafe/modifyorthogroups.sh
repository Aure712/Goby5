#grep -v -f family.remove.id.txt orthogroups.genecount.txt > orthogroups.genecount.mod.txt 
for i in `cat family.remove.id.txt`;do sed -i "/$i/d" orthogroups.genecount.txt;done
