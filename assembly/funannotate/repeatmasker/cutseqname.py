file = "scf.softmasked.fa"
outfile = "scf.softmasked_cut.fa"
fa = []
with open(file) as f:
    for i in f:
        if ">" in i:
            if "np1212" in i:
                if "fragment" in i:
                    i = i.replace("_np1212___fragment_","f")
                    if "___debris" in i:
                        i = i.replace("___debris","d")
                        fa.append(i)
                    else:
                        fa.append(i)
                else:
                    i = i.replace("_np1212","")
                    fa.append(i)
            else:
                fa.append(i)
        else:
            fa.append(i)
with open(outfile,'w') as o:
    for i1 in fa:
        o.write(i1)