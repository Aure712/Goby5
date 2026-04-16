from sys import argv

script, file_gff, length, out = argv

inp_gff = open(file_gff, 'r')
out_f = open(out, 'w')

compare_flag = False
start_last = ''
end_last = ''
distance_last = int(length)
line_last = ''
output = {}
for line in inp_gff.readlines():
    if not line.startswith('TraScf'):
        pass
    elif line.strip().split('\t')[-1] == '-1':
        pass
    elif line.startswith('TraScf'):
        line_split = line.split('\t')
        strand = line_split[6]
        start = line_split[3]
        end = line_split[4]
        distance = int(line_split[18])
        if start == start_last and end == end_last:
            if distance <= distance_last:
                if distance <= int(length):
                    output[line] = ''
                if distance_last <= int(length):
                    del output[line_last]
        else:
            if distance <= int(length):
                output[line] = ''
        start_last = start
        end_last = end
        distance_last = distance
        line_last = line
inp_gff.close()


for line in output:
    out_f.write(line)
out_f.close()