from sys import argv

script, gff = argv

inp_gff = open(gff, 'r')
out_name = gff.split('/')[-1] + '.gene_noUTR.gff'
out_f = open(out_name, 'w')


output = False
for line in inp_gff:
    if not line.startswith('TraScf'):
        pass
    elif not line.startswith('#'):
        line_split = line.split('\t')
        feature = line_split[2]
        start = line_split[3]
        end = line_split[4]
        if 'RNA' in feature:
            output = False
            rna_start = start
            rna_end = end
            out_f.write(line)
        elif feature == 'exon':
            if not output:
                output = True
            elif end == rna_end or start == rna_start:
                pass
            elif output:
                out_f.write(line)
        else:
            out_f.write(line)
inp_gff.close()
out_f.close()