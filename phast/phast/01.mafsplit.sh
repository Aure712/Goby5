#!/bin/bash
#SBATCH -p blade
#SBATCH --mem=48G
#SBATCH -c 1

mkdir -p maf
/public3/group_crf/software/UCSC/mafSplit _.bed maf/ ../cactus/goby24tra.maf -byTarget -useFullSequenceName
#/public3/group_crf/software/UCSC/mafSplit _.bed maf_ggi_loss/ ../cactus2/goby4tba.maf -byTarget -useFullSequenceName
#/public3/group_crf/software/UCSC/mafSplit _.bed maf_tbi_loss/ ../cactus2/goby4tba.maf -byTarget -useFullSequenceName