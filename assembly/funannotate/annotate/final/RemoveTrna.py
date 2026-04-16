infile = "Gma_merged.fixed.gff3"
outfile = "Gma_merged.fixed.removeTrna.gff3"
gtf = []
trna = {}

with open(infile) as f:
    for i in f:
        gtf.append(i)
        if "tRNA" in i:
            flag = i.split("=")[1]
            flag = flag.split("-")[0]
            star = i.split("\t")[3]
            end = i.split("\t")[4]
            trna[flag] = star + "\t" + end

with open(outfile,"w") as o:
    for i in gtf[3:]:
        i2 = i.split("=")[1]
        i2 = i2.split("-")[0]
        star = i.split("\t")[3]
        end = i.split("\t")[4]
        flag2 = star + "\t" + end
        if i2 in trna.keys() and flag2 in trna.values():
            pass
        else:
            o.write(i)