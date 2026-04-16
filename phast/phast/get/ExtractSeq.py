group = []
group_clean = []
family = []
trinity = "/fast3/group_crf/home/g21shaoy23/goby5/annotations/longest_iso/Tba_longest.prot1.fa"

def ExtractSeq(in_file,gene):
    a = []
    b = []
    c = ""
    flag = 0
    with open(in_file) as in_fasta:
        for e1 in in_fasta:
            if gene in e1:
                flag = 1
            if flag == 1:
                a.append(e1)
    for e2 in a:
        if flag <= 2:
            b.append(e2)
        if ">" in e2:
            flag += 1
        if flag == 3:
            break
    for e3 in b[:-1]:
        c += e3
    return c

with open("20k.txt") as in_group:
    for i1 in in_group:
        i1 = i1.strip("\n")
        group.append(i1)
[group_clean.append(i2) for i2 in group if i2 != '']

with open("20k.txt"+".fa","w") as out:
        for gene in group_clean:
            gene = ExtractSeq(trinity,gene)
            out.write(gene)