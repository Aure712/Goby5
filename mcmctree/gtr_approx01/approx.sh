#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 1
#SBATCH --mem=10G

export PATH=/data/software/paml4.9j/src:$PATH

mcmctree mcmctree.ctl > screenlog.log 2>&1
