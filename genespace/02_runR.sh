source /opt/miniconda3/etc/profile.d/conda.sh

export PATH=$PATH:/data/software/MCScanX:/opt/miniconda3/bin/
#conda activate orthofinder
Rscript run.R > run.log 2>&1
