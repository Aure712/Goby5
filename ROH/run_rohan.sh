#!/usr/bin/bash
#SBATCH -p blade,long,himem
#SBATCH -o slurmlog/slurm-%a.out
#SBATCH -e slurmlog/slurm-%a.err
#SBATCH -a 0-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8


module load clusterbasics
module load rohan
ref="/fast3/group_crf/home/g21shaoy23/GATK/bwa/tba/Tba_sort.fa"
sBam_path="/public4/group_crf/home/g21shaoy23/ROH/refbam/"
nJobCount=0
for i in $sBam_path/T-bar.filtered.bam; do
	if (( ${SLURM_ARRAY_TASK_ID} == $nJobCount  )); then
	sF=$(basename $i)
	sStem=${sF/.filtered.bam/}
	echo $sStem
	( mkdir -p $sStem && cd $sStem ; \
	rohan --rohmu 2e-5 -t 8 --size 50000 --step 100 -o $sStem $ref $i > run.log 2>&1 ; cd .. )
	fi
	nJobCount=$(( nJobCount + 1 ))
done;


