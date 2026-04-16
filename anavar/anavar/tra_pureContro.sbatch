#!/bin/bash
#SBATCH -p blade,long,himem,hugemem,gpu,short
#SBATCH -o slurmlog/tra_out_slurm-%a.out
#SBATCH -e slurmlog/tra_out_slurm-%a.err
#SBATCH -a 0-99
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load anavar
#sPop=ORTDRY.excludehybrid
PREFIX=tra_out_
PREFIX2=tra_out_sametheta_
PREFIX3=tra_out_sametheta_nopolerr_

#$ANAVAR ${PREFIX}control.txt ${PREFIX}out.rep${REP}.txt ${PREFIX}log.rep${REP}.txt $(( 333333 + REP ))
cd ./tra_out
anavar ../${PREFIX}control.txt ${PREFIX}out.rep${SLURM_ARRAY_TASK_ID}.txt ${PREFIX}log.rep${SLURM_ARRAY_TASK_ID}.txt $(( 333333 + ${SLURM_ARRAY_TASK_ID} ))
anavar ../${PREFIX2}control.txt ${PREFIX2}out.rep${SLURM_ARRAY_TASK_ID}.txt ${PREFIX2}log.rep${SLURM_ARRAY_TASK_ID}.txt $(( 333333 + ${SLURM_ARRAY_TASK_ID} ))
anavar ../${PREFIX3}control.txt ${PREFIX3}out.rep${SLURM_ARRAY_TASK_ID}.txt ${PREFIX3}log.rep${SLURM_ARRAY_TASK_ID}.txt $(( 333333 + ${SLURM_ARRAY_TASK_ID} ))