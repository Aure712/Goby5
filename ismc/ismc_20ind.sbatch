#!/bin/bash
#SBATCH -p himem,hugemem,long,blade
#SBATCH -o slurmlog/slurm-%a.out
#SBATCH -e slurmlog/slurm-%a.err
#SBATCH -a 1-22
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
sAllTabs=alltabs.tab
gVCF=/fast3/group_crf/home/g21shaoy23/GATK/bwa/joined_call_tba/call/joincalled.genotyped.g.vcf.newid.gz
nChr=${SLURM_ARRAY_TASK_ID} 
sOutDir=20inds/${nChr}

sChr=TbaScf_${nChr}

sAllTabs=`realpath $sAllTabs`

mkdir -p $sOutDir
ln -sf `realpath ismc_opt.bpp` $sOutDir/
ln -sf `realpath mapper_opt_100k.bpp` $sOutDir/
ln -sf `realpath mapper_opt_10k.bpp` $sOutDir/

cd $sOutDir

module load bcftools

 bcftools view -M2 -r ${sChr} $gVCF | gzip -c > chr.vcf.gz
 grep -P "$sChr\t" $sAllTabs > in.tab

source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate iSMC
ismc params=ismc_opt.bpp
ismc_mapper params=mapper_opt_10k.bpp
ismc_mapper params=mapper_opt_100k.bpp
#rm chr.vcf.gz

exit



