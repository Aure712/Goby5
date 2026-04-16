#!/usr/bin/bash
#SBATCH -p gpu,blade,long,himem,hugemem
#SBATCH -c 16

source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate gpu

ulimit -n 9999

saw=/public4/software/saw/saw-8.1.1/bin/saw
ref=/fast3/group_crf/home/g21shaoy23/goby5/annotations/Tba_sort.fa
gff=/fast3/group_crf/home/g21shaoy23/goby5/peaks2utr/utr3k/Tba_clean.utr3k.gff3
inputfastq=/fast3/group_crf/home/g21shaoy23/goby5/steromics/rawdata/fq/


#$saw makeRef --threads-num $SLURM_CPUS_PER_TASK --mode STAR --fasta $ref --gtf $gff --genome ref
$saw count --threads-num $SLURM_CPUS_PER_TASK \
    --id Tridentiger_barbel \
    --sn D03254E512 \
    --omics transcriptomics \
    --kit-version "Stereo-seq T FF V1.2" \
    --sequencing-type "PE100_50+100" \
    --chip-mask /fast3/group_crf/shareddata/2024-11-13_Tridentiger_barbel_stereoseq_mask/D03254E512.barcodeToPos.h5 \
    --organism Tridentiger_barbatus \
    --tissue barbel \
    --fastqs $inputfastq \
    --reference `pwd`/ref/ \
    --image-tar /fast3/group_crf/home/g21shaoy23/goby5/steromics5/ssdna/D03254E512_SC_20241230_094943_4.1.0.tar.gz
