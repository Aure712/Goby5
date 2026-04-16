#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FASTA文件序列名清理脚本
功能：去除序列名中除了'_', '-', ':', '.'外的其他特殊符号
"""

import re
import os


def clean_header(header_line):
    """
    清理FASTA序列名，保留字母、数字和指定的特殊符号
    
    参数:
        header_line: 以'>'开头的序列名行
    
    返回:
        清理后的序列名行
    """
    if not header_line.startswith('>'):
        return header_line
    
    # 保留'>'符号
    prefix = '>'
    # 获取'>'后面的内容
    content = header_line[1:].strip()
    
    # 保留字母、数字、以及'_', '-', ':', '.'
    # 使用正则表达式，替换所有不符合条件的字符为空字符串
    cleaned_content = re.sub(r'[^a-zA-Z0-9_\-:.\s]', '', content)
    
    return prefix + cleaned_content + '\n'


def process_fasta(input_file, output_file=None):
    """
    处理FASTA文件，清理所有序列名
    
    参数:
        input_file: 输入FASTA文件路径
        output_file: 输出文件路径，默认为输入文件名_cleaned.fasta
    """
    # 检查输入文件是否存在
    if not os.path.exists(input_file):
        print(f"错误：文件 '{input_file}' 不存在！")
        return
    
    # 如果没有指定输出文件，则自动生成输出文件名
    if output_file is None:
        base_name = os.path.splitext(input_file)[0]
        output_file = f"{base_name}_cleaned.fasta"
    
    # 统计信息
    total_sequences = 0
    cleaned_count = 0
    
    print(f"正在处理文件: {input_file}")
    print(f"输出文件: {output_file}")
    print("-" * 50)
    
    # 读取并处理文件
    with open(input_file, 'r', encoding='utf-8') as infile, \
         open(output_file, 'w', encoding='utf-8') as outfile:
        
        for line in infile:
            if line.startswith('>'):
                total_sequences += 1
                original_line = line
                cleaned_line = clean_header(line)
                
                # 检查是否有修改
                if original_line != cleaned_line:
                    cleaned_count += 1
                    print(f"序列 {total_sequences}:")
                    print(f"  原始: {original_line.strip()}")
                    print(f"  清理: {cleaned_line.strip()}")
                
                outfile.write(cleaned_line)
            else:
                # 序列数据行直接写入
                outfile.write(line)
    
    # 输出统计信息
    print("-" * 50)
    print(f"处理完成！")
    print(f"总序列数: {total_sequences}")
    print(f"修改的序列名数量: {cleaned_count}")
    print(f"输出文件已保存: {output_file}")


if __name__ == "__main__":
    # ============================================
    # 在这里指定输入文件路径
    # ============================================
    input_fasta = "/public4/group_crf/home/g21shaoy23/cactus/0_species/Rhinogobius_similis/Rhinogobius_similis2.fa"  # 修改为你的输入文件路径
    output_fasta = "/public4/group_crf/home/g21shaoy23/cactus/0_species/Rhinogobius_similis/Rhinogobius_similis3.fa"  # 可选：指定输出文件路径，默认为 input_cleaned.fasta
    
    # 处理FASTA文件
    process_fasta(input_fasta, output_fasta)
