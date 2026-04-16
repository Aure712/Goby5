#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 30
#SBATCH --mem=50G
source /data/software/funannotate/env.sh

GENOME=./scf.fa
SPECIES="Tridentiger barbatus"
STRAIN="PanYu2021"

funannotate sort -i $GENOME -o Tba_sort.fa -b TbaScf