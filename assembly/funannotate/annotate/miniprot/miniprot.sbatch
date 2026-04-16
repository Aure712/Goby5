#!/usr/bin/bash
#SBATCH -p long,gpu,himem,hugemem,blade
#SBATCH -c 64
#SBATCH --mem=100G

hostname

CPU=$SLURM_CPUS_PER_TASK

module load miniprot

miniprot -I -t $CPU --gff ref.mpi /fast3/group_crf/home/cuirf/goby_flye/Oxyeleotris_marmorata/annotation/prot_evidences/all.prot.fa  > out.gff

