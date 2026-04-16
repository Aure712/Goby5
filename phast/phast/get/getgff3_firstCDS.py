from sys import argv

script, file_gff = argv

inp_gff = open(file_gff, 'r')
out_f = open(file_gff + ".firstCDS.gff", 'w')

out_flag = False
out_list = []
for line in inp_gff.readlines():
    if not line.startswith('TbaScf'):
        pass
    elif line.startswith('TbaScf_23'):
        if out_flag:
            out_flag = False
            if out_strand == '+':
                out_f.write(out_list[0])
            else:
                out_f.write(out_list[-1])
        break
    elif line.startswith('TbaScf'):
        line_split = line.split('\t')
        chrom = line_split[0]
        feature = line_split[2]
        strand = line_split[6]
        start = line_split[3]
        end = line_split[4]
        attributes = line_split[8]
        attributes_split = attributes.strip().split(';')
        if 'RNA' in feature:
            if out_flag:
                out_flag = False
                if out_strand == '+':
                    out_f.write(out_list[0])
                else:
                    out_f.write(out_list[-1])
            rna_id = attributes_split[0].strip('ID=')
            out_flag = True
            out_strand = strand
            out_rna = rna_id
            rna_start = start
            rna_end = end
            out_list = []
        elif feature == 'CDS' and out_flag:
            out_list.append(line)
        else:
            pass

out_f.close()