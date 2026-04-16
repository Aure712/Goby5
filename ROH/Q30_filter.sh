#!/usr/bin/bash
#SBATCH -p blade,long,himem
#SBATCH -o slurmlog/slurm-%a.out
#SBATCH -e slurmlog/slurm-%a.err
#SBATCH -a 0-2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8


module load samtools

sPath="/public4/group_crf/home/g21shaoy23/ROH/refbam"
nJobCount=0
for i in $sPath/*.dedup.bam; do
	if (( ${SLURM_ARRAY_TASK_ID} == $nJobCount  )); then
	sFile=$(basename $i)
	sStem=${sFile%%\.*}
	echo $sStem
	samtools view -q 30 -o $sStem.filtered.bam -O bam $i
	samtools index $sStem.filtered.bam
	fi
	nJobCount=$(( nJobCount + 1 ))
	sleep 10
done;
