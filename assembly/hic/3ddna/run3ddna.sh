Contigs=/public3/group_crf/home/g21shaoy23/HiC/Juicer/Tradiatus/references/ref.fa
MND=/public3/group_crf/home/g21shaoy23/HiC/Juicer/Tradiatus/aligned/merged_nodups.txt
MAPQ=30
GAPSIZE=1000

module load java
/public3/group_crf/software/3d-dna-201008/run-asm-pipeline.sh -m haploid -i 100 -r 3 -q $MAPQ --sort-output -g $GAPSIZE $Contigs $MND > 3ddna.log 2>&1