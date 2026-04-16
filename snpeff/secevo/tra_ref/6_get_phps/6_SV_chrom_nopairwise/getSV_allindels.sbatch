#!/usr/bin/bash
#SBATCH -p gpu,himem,hugemem,blade,long
#SBATCH -o slurmlog/slurm-%a.out
#SBATCH -e slurmlog/slurm-%a.err
#SBATCH -a 0-524
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load bcftools
module load samtools
module load clusterbasics
module load hyphy

REF=Tridentiger_radiatus

nJobCount=0
for i in $(bcftools view -h ../allpops.merged.snpeff.renameoutgroup.vcf.gz | grep "^##contig=<ID=" | cut -f3 -d'=' | cut -f1 -d','); do
	if (( ${SLURM_ARRAY_TASK_ID} == $nJobCount  )); then
#		echo $i
		( bcftools view -r $i ../allpops.merged.snpeff.renameoutgroup.vcf.gz -o allpops.snpeff.$i.vcf.gz
		tabix allpops.snpeff.$i.vcf.gz
		SnpEffVcf=allpops.snpeff.$i.vcf.gz
		Out=SV.w.indelmodifier.$i.polarized.out.txt
       		php getSV_allindels_nomask.php "$SnpEffVcf" "$Out")
	fi
	nJobCount=$(( nJobCount + 1 ))
done;


