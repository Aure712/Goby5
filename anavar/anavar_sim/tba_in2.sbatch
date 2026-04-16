#!/bin/bash
#SBATCH -p blade,long,himem,hugemem,gpu
#SBATCH -o slurmlog/tba_in_slurm-%a.out
#SBATCH -e slurmlog/tba_in_slurm-%a.err
#SBATCH -a 0-99
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load anavar
PREFIX=tba_in_

cd ./tba_in

# 获取所有匹配的control文件并存储为数组
control_files=($(ls ../${PREFIX}control*.txt 2>/dev/null))

# 检查是否找到文件
if [ ${#control_files[@]} -eq 0 ]; then
    echo "错误：未找到匹配的control文件 ../${PREFIX}control*.txt"
    exit 1
fi

# 对每个control文件运行anavar
for control_file in "${control_files[@]}"; do
    echo "正在处理文件: $control_file"
    anavar "$control_file" ${PREFIX}out.rep${SLURM_ARRAY_TASK_ID}.$(basename "$control_file" .txt).txt ${PREFIX}log.rep${SLURM_ARRAY_TASK_ID}.$(basename "$control_file" .txt).txt $(( 333333 + ${SLURM_ARRAY_TASK_ID} ))
done
