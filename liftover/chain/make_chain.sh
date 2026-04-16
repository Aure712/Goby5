#tba coordinate to tra coordinate, tba is reference genome and tra is query genome in chain file

query=/data2/projects/yshao/Goby5Compare/synteny/tba_tra2/Tba_sort.fa
ref=/data2/projects/yshao/Goby5Compare/synteny/tba_tra2/Tra_sort.fa
transanno=/data/software/transanno/transanno-x86_64-unknown-linux-musl-v0.4.4/transanno

minimap2 -cx asm5 --cs $query $ref > tba2tra.paf

$transanno minimap2chain tba2tra.paf --output tba2tra.chain
