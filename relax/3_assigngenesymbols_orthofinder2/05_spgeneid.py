# 从 mapping_results.txt 中读取映射关系
mapping_dict = {}
with open('mapping_results.txt', 'r') as mapping_file:
    for line in mapping_file:
        line = line.strip()
        if line:
            parts = line.split("\t")
            if len(parts) == 2:
                first_column = parts[0]
                second_column = parts[1]
                mapping_dict[second_column] = first_column

# 打开 Rform.longest_isoform.genesymbol.gff3 文件进行处理
with open('Rform.longest_isoform.genesymbol.gff3', 'r') as gff3_file:
    with open('Rform.longest_isoform.genesymbol.spgeneid.gff3', 'w') as output_file:
        for line in gff3_file:
            if not line.startswith("#"):
                fields = line.strip().split('\t')
                if fields[2] in {"mRNA", "gene"}:
                    attributes = fields[8].split(';')
                    for i in range(len(attributes)):
                        attr = attributes[i].strip()
                        if attr.startswith("ID="):
                            ID = attr[3:]
                            # 查找部分匹配的映射关系
                            spgeneid = None
                            for key, value in mapping_dict.items():
                                if key in ID or ID in key:
                                    spgeneid = value
                                    break
                            if spgeneid:
                                # 添加 ";spgeneid=" 和 spgeneid 标识符
                                attributes.append(f"spgeneid={spgeneid}")
                                # 更新第九列的信息
                                fields[8] = ';'.join(attributes)
                    # 写入处理后的行到输出文件
                    output_file.write('\t'.join(fields) + '\n')
                else:
                    # 对于不是 "mRNA" 或 "gene" 的行，直接写入到输出文件
                    output_file.write(line)

print("处理完成，结果已保存到 output.gff3 文件中。")

