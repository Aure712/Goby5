#!/bin/bash
#SBATCH -p long
#SBATCH -o slurmlogT/slurm-%a.out
#SBATCH -e slurmlogT/slurm-%a.err
#SBATCH -a 0-399
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load hyphy/2.5.1
module load clusterbasics

#SET THIS TO THE CORRECT NUMBER
TreeFolder="./for_relax/Relax_T_only"
OutDIR=`basename $TreeFolder`;


mkdir -p $OutDIR

 php bustedmh.php -N ${SLURM_ARRAY_TASK_COUNT} -f ${SLURM_ARRAY_TASK_ID} -T $TreeFolder -o $OutDIR

exit



