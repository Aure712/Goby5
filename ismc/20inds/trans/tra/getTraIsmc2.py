import pandas as pd
import pyranges as pr

# 读取 a.bed 文件
bed_columns = ['Chromosome', 'Start', 'End', 'Score', 'Ref_Chr', 'Other1', 'Target_Coord']
bed_df = pd.read_csv('../all_tba2tra.filter.bed', sep='\t', header=None, names=bed_columns)

# 读取 b.txt 文件
txt_columns = ['Chromosome', 'Start', 'End'] + [f'col_{i}' for i in range(4, 18)]
txt_df = pd.read_csv('all.out.rho.100kb.bedgraph', sep='\t', header=None, names=txt_columns)

# 计算第4、5、7、8列的平均值
txt_df['average'] = txt_df[['col_5', 'col_6', 'col_8', 'col_9', 'col_10']].mean(axis=1)

# **步骤 1：转换 a.bed 染色体格式** (TraScf_1 → chr1)
bed_df['Chromosome'] = bed_df['Ref_Chr'].str.replace("TraScf_", "chr", regex=False)

# **步骤 2：创建 pyranges 结构**
# 1. 修正 bed_pr 结构
bed_df['Start'] = bed_df['Target_Coord'].astype(int)
bed_df['End'] = bed_df['Start'] + 1  # 确保 Start < End

# 2. 构造 pyranges 结构
bed_pr = pr.PyRanges(bed_df[['Chromosome', 'Start', 'End']])

# b.txt 作为目标区间，使用计算的平均值列
txt_pr = pr.PyRanges(txt_df[['Chromosome', 'Start', 'End', 'average']])

# **步骤 3：执行高效区间匹配**
result = bed_pr.join(txt_pr, how="left")  # left join 保留所有 bed_pr 记录

# **步骤 4：合并匹配结果**
matched_df = result.df[['Chromosome', 'Start', 'End', 'average']].rename(columns={'Start': 'Target_Coord'})

# **步骤 5：合并回 a.bed**
final_df = bed_df.merge(matched_df, on=['Chromosome', 'Target_Coord'], how='left')

# **步骤 6：填充未匹配的值**
final_df['average'].fillna('NA', inplace=True)

# **步骤 7：输出**
final_df.to_csv('output2.bed', sep='\t', index=False, header=False)
