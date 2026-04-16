infile = "synorthos.txt"
infile_gff = "Tba_merged.fixed.removeTrna.gff3"
out_rm_gene = "soloGenes_tba.txt"
out_clean_gff = "Tba_clean.gff3"
out_solo_gene_gff = "Tba_solo.gff3"
sp_num1 = 6
sp_num2 = 13
table = []
table2 = []
table3 = []
table4 = []
table5 = []
gff = []
gff_solo = []
gff_clean = []

with open(infile) as f:
    for i in f:
        i = i.strip()
        if all(element == 'NA' for index, element in enumerate(i.split('\t')) if index not in [0,1,2,sp_num1,sp_num2]):
            table.append(i)
for i in table:
    flag = i.split("\t")
    if flag[sp_num1] == 'NA' and flag[sp_num2] == 'NA':
        pass
    else:
        table2.extend([flag[sp_num1]] + [i for i in flag[sp_num2].split(";")])
table3 = [i.split("-")[0] for index, i in enumerate(table2) if i not in table2[:index]]
table3.remove("NA")

with open(infile_gff) as f:
    for i in f:
        gff.append(i)

for i in gff:
    if any(flag in i for flag in table3) and ";Parent" in i:
        tag = i.split("ID=")[1]
        tag = tag.split("-")[0]
        table4.append(tag)

table4.extend(table3)
table5 = [i for index, i in enumerate(table4) if i not in table4[:index]]

for i in gff:
    if any(flag in i for flag in table5):
        gff_solo.append(i)
    else:
        gff_clean.append(i)

with open(out_rm_gene,'w') as o:
    for i in table5:
        o.write(i + "\n")

with open(out_solo_gene_gff,'w') as o:
    for i in gff_solo:
        o.write(i)

with open(out_clean_gff,'w') as o:
    for i in gff_clean:
        o.write(i)