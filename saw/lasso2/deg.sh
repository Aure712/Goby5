#!/usr/bin/bash
#SBATCH -p long,himem,hugemem,blade
#SBATCH -c 16

source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate gpu

ulimit -n 9999

saw=/public4/software/saw/saw-8.1.1/bin/saw


$saw reanalyze diffExp --count-data /fast3/group_crf/home/g21shaoy23/goby5/steromics5/Tridentiger_barbel/ \
	     --diffexp-geojson D03254E512_20251215191126.diffexp.labelvslabel.geojson \
	     --output ./out
