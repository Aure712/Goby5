#!/usr/bin/bash
#SBATCH -c 30
#SBATCH -p long
#SBATCH --mem=50G

ASSEMBLY=scf.fa # assembly to mask
REPBASELIB=`realpath ./all_te.fa` #repeat library, with killifish repeats
THREADS=30

#ln -sf $ASSEMBLY scf.fa


singularity run docker://dfam/tetools:1.3.1 RepeatMasker -gff -html -pa $THREADS -lib $REPBASELIB  scf.fa

