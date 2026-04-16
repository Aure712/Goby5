import os
import pandas as pd

# 获取当前目录下的所有子目录
base_dir = os.getcwd()
all_data = []

# 遍历所有子目录
for folder in os.listdir(base_dir):
    folder_path = os.path.join(base_dir, folder)
    
    # 检查是否是文件夹，且文件夹名为数字
    if os.path.isdir(folder_path) and folder.isdigit():
        # 构造文件路径
        file_path = os.path.join(folder_path, 'out.rho.100kb.bedgraph')
        
        # 如果文件存在，读取文件
        if os.path.exists(file_path):
            # 读取文件，跳过前两行
            df = pd.read_csv(file_path, sep='\t', header=None, skiprows=2)
            
            # 替换第一列的"chr1"为"chr"加上文件夹的数字
            df[0] = df[0].replace('chr1', f'chr{folder}', regex=False)
            
            # 将当前文件的数据添加到列表中
            all_data.append(df)

# 合并所有数据
merged_data = pd.concat(all_data, ignore_index=True)

# 输出合并后的结果
merged_data.to_csv('all.out.rho.100kb.bedgraph', sep='\t', header=False, index=False)
