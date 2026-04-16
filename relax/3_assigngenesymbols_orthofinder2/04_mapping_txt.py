# 用于存储标题行的映射关系的字典
mapping_dict = {}

# 打开 allproteins.fa 文件
with open('allproteins.fa', 'r') as fa_file:
    for line in fa_file:
        if line.startswith(">"):
            # 提取标题行中的第二列和第三列
            header = line.strip()
            columns = header.split('|')
            if len(columns) >= 3:
                second_column = columns[1]
                third_column = columns[2]
                # 去掉标题中的 "-T1"
                third_column = third_column.replace('-T1', '')
                # 将映射关系添加到字典中
                mapping_dict[header] = (second_column, third_column)

# 打印映射关系
# 将结果写入到文本文件
with open('mapping_results.txt', 'w') as output_file:
    for header, (second_column, third_column) in mapping_dict.items():
        output_file.write(f"{second_column}\t{third_column}\n")

print("映射关系已保存到 mapping_results.txt 文件中。")

