#!/bin/bash
#SBATCH -p himem,hugemem
#SBATCH --mem=500G
#SBATCH -c 1
#SBATCH --array=1-11%11

module load R/4.3.0


#Rscript rphastcon_final5_allsp.R $SLURM_ARRAY_TASK_ID all cache=1 rebuild=1 compat=1 

Rscript rphastcon_final5_bar2_scalesame2_exportmsa2.R $SLURM_ARRAY_TASK_ID arrfreshs cache=1 rebuild=1 compat=1