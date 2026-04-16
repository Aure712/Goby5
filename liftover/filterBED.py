inbed = "all_tba2tra.bed"
outbed = "all_tba2tra.filter.bed"
bed = []

with open(inbed) as f:
    for i in f:
        i2 = i.split("\t")
        flag1 = i2[0].replace("chr","")
        flag2 = i2[4].replace("TraScf_","")
        if int(flag1) > 22 or int(flag2) > 11:
            continue
        else:
            bed.append(i)

with open(outbed,"w") as o:
    for i in bed:
        o.write(i)