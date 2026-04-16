#!/usr/bin/bash
#SBATCH -p long,himem,blade,himem,hugemem
#SBATCH -c 2
#SBATCH --mem=120G

#source /public/apps/miniconda3/etc/profile.d/conda.sh
#conda activate EDTA
module load apptainer
#calcDivergenceFromAlign=/public4/group_crf/home/g21shaoy23/repeatmask/tba/calcDivergenceFromAlign.pl
#createRepeatLandscape=/public/software/conda_envs/EDTA/share/RepeatMasker/util/createRepeatLandscape.pl
repeatmasker_align=/public4/group_crf/home/g21shaoy23/repeatmask/tba/Tba_review_s.fasta.align
divsum=/public4/group_crf/home/g21shaoy23/repeatmask/tba/Tba.genome.divsum
genomeSize=886424161
#awk '{sum += $2} END {print sum}' Tba_sort.fa.fai
#perl $calcDivergenceFromAlign -s $divsum $repeatmasker_align
#perl $createRepeatLandscape -div $divsum  -g 886424161  > replandscape.html # -g genomesize to calculate the percentage of genome covered by repeats
#singularity exec /public4/software/tetools/tetools.sif calcDivergenceFromAlign.pl -s $divsum $repeatmasker_align
#singularity exec /public4/software/tetools/tetools.sif createRepeatLandscape.pl -div $divsum -g $genomeSize | sed 's|www.google.com/jsapi|www.gstatic.com/charts/loader.js|' > tba.landscape.html
perl createRepeatLandscape2.2.pl -div $divsum -g $genomeSize | sed 's|www.google.com/jsapi|www.gstatic.com/charts/loader.js|' > tba.landscape.html