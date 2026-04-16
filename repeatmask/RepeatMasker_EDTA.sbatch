#!/usr/bin/bash
#SBATCH -p long,himem,blade,himem,hugemem,gpu
#SBATCH -c 32
#SBATCH --mem=120G

source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate EDTA

CPU=32
genome=/public4/group_crf/home/g21shaoy23/repeatmask/ref/Tba_review_s.fasta
out_dir=/public4/group_crf/home/g21shaoy23/repeatmask/tba/
TE_lib=/public4/group_crf/home/g21shaoy23/repeatmask/ref/all_te.fa


RepeatMasker -pa $CPU -s -lib $TE_lib -no_is -cutoff 225 -nolow -norna -div 40 -dir $out_dir -a -inv -xsmall $genome
# -s slow mode; 
# -no_is skips bacterial insertion element check
# -cutoff [number] sets cutoff score for masking repeats when using -lib; 
# -norna not mask rna gene; 
# -nolow  does not mask low complexity DNA or simple repeats
# -div [number]  masks only those repeats that are less than [number] percent diverged from the consensus sequence
# -a -inv shows the alignments in a .align output file;
# -xsmall returns repetitive regions in lowercase (rest capitals) rather than masked

