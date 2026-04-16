#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 12
#SBATCH --mem=50G

ASSEMBLY=scf.fa # the assembly to look at
SP=Tba
REPMODEL="singularity run docker://dfam/tetools:1.3.1 " #command of repeatmodeler
THREADS=12

#ln -sf $ASSEMBLY scf.fa

$REPMODEL BuildDatabase -name $SP -engine ncbi $ASSEMBLY

$REPMODEL RepeatModeler -LTRStruct -engine ncbi -pa $THREADS -database $SP
