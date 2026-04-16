infa = "Tba0.fa"
outfa = "Tba.fa"
seq = {}
sp = "Tba"

class InvalidFastaError(ValueError):
    pass

def ReadFasta(fasta):
    sequences = dict()
    active_seq = None
    seq_list = list()
    with open(fasta) as f:
        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == '>':
                if active_seq is not None:
                    sequences[active_seq] = ''.join(seq_list)
                    seq_list = list()
                active_seq = line[1:].split()[0]
                if active_seq in sequences.keys():
                    raise InvalidFastaError('Fasta {0} contains multiple contigs named {1}'.format(fasta, active_seq))
                sequences[active_seq] = ''
            elif active_seq is not None:
                seq_list.append(line)
            else:
                raise InvalidFastaError('Fasta {0} does not begin with a contig name'.format(fasta))
    sequences[active_seq] = ''.join(seq_list)
    return sequences

seq = ReadFasta(infa)
with open(outfa,"w") as o:
    for k,v in seq.items():
        if "Chromosome" in k:
            k = k.replace("Chromosome",sp+"Scf_")
            o.write(">" + k + "\n" + v + "\n")
        else:
            o.write(">" + k + "\n" + v + "\n")
