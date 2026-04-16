#!/bin/bash
#SBATCH -p blade
#SBATCH -o slurmlog/RelaxNothos_slurm-%a.out
#SBATCH -e slurmlog/RelaxNothos_slurm-%a.err
#SBATCH -a 0-199
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

#SET THIS TO THE CORRECT NUMBER
TreeFolder="/data2/projects/yshao/relax/4_exprotAln_genespace2/Relax_Tra/"
OutDIR=`basename $TreeFolder`;

TOTALPARTS=32

mkdir -p $OutDIR
mkdir -p logs

for i in `seq 0 31`; do
php hyphyrelax.php -N ${TOTALPARTS} -f ${i} -T $TreeFolder -o $OutDIR > logs/$i.log 2>&1 &

done
wait
exit

