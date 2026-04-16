import re

input_file = "Tra.prechange_longest_fixed.1.gff3"
output_file = "Tra.prechange_longest_fixed.gff3.firstCDS.gff3"

cds_records = {}

chromosome_pattern = re.compile(r'^[1-9]$|^1[0-1]$')

with open(input_file, 'r') as infile:
    for line in infile:
        if line.startswith("#"):
            continue
        fields = line.strip().split('\t')
        chromosome = fields[0]
        feature_type = fields[2]
        strand = fields[6]
        attributes = fields[8]
        if not chromosome_pattern.match(chromosome.split("_")[1]):
            continue

        if feature_type == "CDS":
            match = re.search(r'Parent=([^;]+)', attributes)
            if match:
                gene_id = match.group(1)
                if gene_id not in cds_records:
                    cds_records[gene_id] = []
                cds_records[gene_id].append((fields, strand))

result_records = {}

for gene_id, records in cds_records.items():
    records.sort(key=lambda x: int(x[0][3]))
    if records[0][1] == '+':
        result_records[gene_id] = records[0][0]
    else:
        result_records[gene_id] = records[-1][0]

with open(output_file, 'w') as outfile:
    for record in result_records.values():
        outfile.write("\t".join(record) + "\n")