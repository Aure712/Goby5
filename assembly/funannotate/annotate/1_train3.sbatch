#!/usr/bin/bash
#SBATCH -p long,himem,blade
#SBATCH -c 24
#SBATCH --mem=50G
source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate funannotate1.8.15

GENOME=Tba_sort.fa
TRINITY_TRANSCRIPT=tba_comb_tmi_cd.fa
SPECIES="Tridentiger barbatus"
STRAIN="PanYu2021_comb_tmi"
CPU=24

funannotate train -i $GENOME -o trained_3  \
--left left.fq.gz \
--right right.fq.gz  \
--no_trimmomatic \
--no_normalize_reads \
--pasa_min_avg_per_id 72 \
--jaccard_clip \
--memory 100G \
--trinity "$TRINITY_TRANSCRIPT" \
--species "$SPECIES" \
--pasa_db mysql \
--max_intronlen 200000 \
--aligners minimap2 minimap2 \
--strain "$STRAIN" --cpus $CPU
