#!/usr/bin/bash
#SBATCH -p himem,long,blade
#SBATCH -c 32
#SBATCH --mem=50G
source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate funannotate1.8.15

GENOME=Tba_sort.fa
SPECIES="Tridentiger barbatus"
STRAIN="PanYu2021_comb_tmi"
CPU=32
STRINGTIE=./stringtie/Tba_alnout.gtf
PROT=../all.prot.fa

funannotate update -i trained_3  --cpus $CPU --pasa_db mysql --alt_transcripts 0.3
