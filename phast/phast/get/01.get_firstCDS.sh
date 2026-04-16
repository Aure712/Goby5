source /public/apps/miniconda3/etc/profile.d/conda.sh
module load python3
conda activate bedtools

bedtools sort -i Tra.prechange_longest_fixed.gff3 > Tra.prechange_longest_fixed.sort.gff3
sed '/^$/d' Tra.prechange_longest_fixed.gff3 > Tra.prechange_longest_fixed.1.gff3
python3 getFirstCDS.py

grep '+' Tra.prechange_longest_fixed.gff3.firstCDS.gff3 > Tra.prechange_longest_fixed.gff3.+firstCDS.gff3
grep -v '+' Tra.prechange_longest_fixed.gff3.firstCDS.gff3 > Tra.prechange_longest_fixed.gff3.-firstCDS.gff3
