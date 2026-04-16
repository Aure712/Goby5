import pandas as pd
import pyranges as pr

# 读取 bed 文件
bed_df = pd.read_csv('all_tba2tra.filter.bed', sep='\t', header=None, 
                    names=['Chromosome', 'Start', 'End', 'deep', 'other_genome', 'num', 'position'])

# 读取 b.txt 文件
txt_df = pd.read_csv('all.out.rho.100kb.bedgraph', sep='\t', header=None,
                    names=['Chromosome', 'Start', 'End'] + [f'col_{i}' for i in range(4, 18)])

# 计算指定列的平均值
cols_to_avg = ['col_4', 'col_6', 'col_8', 'col_9', 'col_10', 'col_11', 'col_12', 'col_13', 'col_15']
txt_df['mean_value'] = txt_df[cols_to_avg].mean(axis=1)
txt_df = txt_df[['Chromosome', 'Start', 'End', 'mean_value']]  # 保留关键列

# 确保数据类型
bed_df[['Start', 'End']] = bed_df[['Start', 'End']].astype(int)
txt_df[['Start', 'End']] = txt_df[['Start', 'End']].astype(int)

# 转换为 PyRanges
bed_pr = pr.PyRanges(bed_df)
txt_pr = pr.PyRanges(txt_df)

# 执行区间join
joined_pr = bed_pr.join(txt_pr, suffix='_txt')
joined_df = joined_pr.df

# 调试：打印实际列名
print("当前所有列:", joined_df.columns.tolist())

# 根据实际列名调整（示例使用无后缀的情况）
joined_df = joined_df[['Chromosome', 'Start', 'End', 'deep', 'other_genome', 
                      'num', 'position', 'mean_value']]  # 注意这里没有_txt后缀

joined_df.to_csv('output2.bed', sep='\t', index=False, header=False)
