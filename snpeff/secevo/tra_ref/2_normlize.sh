#!/usr/bin/bash
#SBATCH -p gpu,himem,hugemem,blade
#SBATCH -o slurmlog/slurm-%a.out
#SBATCH -e slurmlog/slurm-%a.err
#SBATCH -a 0-524
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load bcftools
module load samtools

REF=Tridentiger_radiatus
genome=/fast3/group_crf/home/g21shaoy23/goby5/annotations/Tra_sort.fa

#bcftools view -O z -M 2 -m 2 -i 'TYPE=="snp" || TYPE=="indel"' allpops.withoutgroup.newid.g.vcf.gz > allpops.withoutgroup.vcf.gz 

nJobCount=0
for i in $(bcftools view -h allpops.withoutgroup.vcf.gz | grep "^##contig=<ID=" | cut -f3 -d'=' | cut -f1 -d','); do
	if (( ${SLURM_ARRAY_TASK_ID} == $nJobCount  )); then
#		echo $i
		( bcftools norm -w 99999999 -O z -f $genome -r $i allpops.withoutgroup.vcf.gz > allpops.withoutgroup.norm.$i.vcf.gz
		tabix allpops.withoutgroup.norm.$i.vcf.gz )
	fi
	nJobCount=$(( nJobCount + 1 ))
done;

#bcftools merge -o allpops.withoutgroup.norm.vcf.gz -O z allpops.withoutgroup.norm.*.vcf.gz
#tabix allpops.withoutgroup.norm.vcf.gz

