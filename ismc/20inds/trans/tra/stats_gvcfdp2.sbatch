#!/bin/bash
#SBATCH -p long,blade,himem,gpu
#SBATCH -c 5
#SBATCH --mem 200G

module load bcftools
input=/fast3/group_crf/home/g21shaoy23/GATK/bwa/joined_call_tra/call/joincalled.genotyped.g.vcf.newid.gz

# 获取样本名称列表并存储到数组
samples=($(bcftools query -l $input))

bcftools query -f '[%DP\t]\n' $input | \
awk -v samples_str="${samples[*]}" '
BEGIN {
    # 将样本名称字符串转换为数组
    split(samples_str, samples, " ")
}
{
    for (i=1; i<=NF; i++) {
        if ($i != ".") {
            sum[i] += $i
            count[i]++
        }
    }
}
END {
    for (i=1; i<=length(samples); i++) {
        avg = (count[i] > 0) ? sum[i]/count[i] : 0
        printf "%s: Average DP = %.2f\n", samples[i], avg
    }
}'