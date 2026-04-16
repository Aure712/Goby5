#!/usr/bin/bash
#SBATCH -p blade,long,himem,hugemem
#SBATCH -c 10
#SBATCH --mem=180G

#source /public/apps/miniconda3/etc/profile.d/conda.sh
#conda activate /public4/software/conda_env/jupyter/

set -euo pipefail
export MPLBACKEND=Agg

/public4/software/conda_env/stereopy/bin/python3.8 -u st_pipeline2.py

