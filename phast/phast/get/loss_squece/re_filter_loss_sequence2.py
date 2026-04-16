from sys import argv

script, origin_loss_gff, loss_bed, out_name, length = argv
length = int(length)

dic_loss_sequence = {}

with open(origin_loss_gff, 'r') as inp_gff:
    for line in inp_gff:
        if line.startswith('#'):
            continue
        line_split = line.strip().split('\t')
        chrom = line_split[0]
        feature = line_split[2]
        start = int(line_split[3])
        end = int(line_split[4])
        score = line_split[5]
        loss_id = line_split[6]

        if chrom not in dic_loss_sequence:
            dic_loss_sequence[chrom] = {}
        dic_loss_sequence[chrom][start] = [end, feature, score, loss_id]

with open(loss_bed, 'r') as inp_bed, open(out_name, 'w') as out_f:
    for line in inp_bed:
        if line.startswith('#'):
            continue
        line_split = line.strip().split('\t')
        chrom = line_split[0]
        bed_start = int(line_split[1])
        bed_end = int(line_split[2])

        if chrom not in dic_loss_sequence:
            continue

        for check_start in dic_loss_sequence[chrom]:
            check_end, feature, score, loss_id = dic_loss_sequence[chrom][check_start]
            if check_start <= bed_start <= check_end:
                check_length = bed_end - bed_start
                if check_length >= length:
                    gff_start = bed_start + 1  # BED 0-based -> GFF 1-based
                    out_f.write(
                        chrom + '\t' + 'phast' + '\t' + feature + '\t' +
                        str(gff_start) + '\t' + str(bed_end) + '\t' +
                        score + '\t' + '+' + '\t' + '.' + '\t' + loss_id + '\n'
                    )
