#!/usr/bin/bash
#SBATCH -p himem,long,blade
#SBATCH -c 24
#SBATCH --mem=50G
source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate funannotate1.8.15

GENOME=Tba_sort.fa
TRINITY_TRANSCRIPT=tba_comb_tmi_cd.fa
SPECIES="Tridentiger barbatus"
STRAIN="PanYu2021_comb_tmi"
CPU=24
STRINGTIE=./stringtie/Tbi_alnout.gtf
PROT=../all.prot.fa


 funannotate predict --force -i $GENOME \
            -o trained_3 -s "$SPECIES" --strain $STRAIN --cpus $CPU \
--max_intronlen 200000 \
--weights augustus:10 hiq:25 pasa:20 snap:0 genemark:0 glimmerhmm:0 proteins:25 transcripts:2 \
--p2g_pident 65 \
--optimize_augustus \
--busco_db actinopterygii \
--busco_seed_species zebrafish \
--organism other \
--protein_evidence $PROT \
--stringtie $STRINGTIE \
--keep_no_stops \
--min_protlen 10 \
--repeat_filter none
--repeats2evm \
#--no-evm-partitions \
# StartWeights = {'augustus': 1, 'hiq': 2, 'genemark': 1, 'pasa': 6,'codingquarry': 0, 'snap': 1, 'glimmerhmm': 1, 'proteins': 1, 'transcripts': 1}
