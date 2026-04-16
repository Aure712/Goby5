#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
批量修改TBA配置文件脚本
用于将sampling results数据批量替换到control文件中
"""

import os
import re


def read_file(filepath):
    """读取文件内容"""
    with open(filepath, 'r', encoding='utf-8') as f:
        return f.read()


def read_sampling_results(filepath):
    """读取sampling results文件，返回行列表"""
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    # 去除每行末尾的换行符，保留数据内容
    return [line.strip() for line in lines if line.strip()]


def modify_control_content(control_content, neu_data, sel_data):
    """
    修改control文件内容
    替换neu_sfs和sel_sfs后面的数据
    """
    # 使用正则表达式替换neu_sfs后面的内容
    # 匹配模式：neu_sfs: 后面跟着任意字符直到换行
    neu_pattern = r'(neu_sfs:\s*)[^\n]+'
    control_content = re.sub(neu_pattern, r'\1 ' + neu_data, control_content)

    # 使用正则表达式替换sel_sfs后面的内容
    # 匹配模式：sel_sfs: 后面跟着任意字符直到换行
    sel_pattern = r'(sel_sfs:\s*)[^\n]+'
    control_content = re.sub(sel_pattern, r'\1 ' + sel_data, control_content)

    return control_content


def save_modified_file(content, output_filepath):
    """保存修改后的文件"""
    with open(output_filepath, 'w', encoding='utf-8') as f:
        f.write(content)
    print(f"已保存文件: {output_filepath}")


def main():
    # ===== 配置部分 - 在这里指定输入文件路径 =====
    control_file = "tba_out_control.txt"  # 控制文件
    neu_sampling_file = "tba_pureControl_neu_m_sampling_results.txt"  # neu采样结果文件
    sel_sampling_file = "tba_pureControl_sel_m_sampling_results.txt"  # sel采样结果文件

    # 输出目录（可选，默认为当前目录）
    output_dir = "."  # 可以修改为其他目录路径
    # ==========================================

    # 检查文件是否存在
    for filepath in [control_file, neu_sampling_file, sel_sampling_file]:
        if not os.path.exists(filepath):
            print(f"错误: 文件 {filepath} 不存在！")
            return

    # 读取文件内容
    print("正在读取文件...")
    control_content = read_file(control_file)
    neu_data_lines = read_sampling_results(neu_sampling_file)
    sel_data_lines = read_sampling_results(sel_sampling_file)

    # 检查数据行数是否一致
    if len(neu_data_lines) != len(sel_data_lines):
        print(f"警告: neu文件有{len(neu_data_lines)}行，sel文件有{len(sel_data_lines)}行，行数不一致！")
        print("将使用较少的行数进行配对。")

    # 获取配对数量
    pair_count = min(len(neu_data_lines), len(sel_data_lines))
    print(f"共有 {pair_count} 对数据需要处理")

    # 创建输出目录（如果不存在）
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 处理每一对数据
    print("\n开始处理配对数据...")
    for i in range(pair_count):
        # 获取当前配对的数据
        neu_data = neu_data_lines[i]
        sel_data = sel_data_lines[i]

        # 修改control内容
        modified_content = modify_control_content(control_content, neu_data, sel_data)

        # 生成输出文件名（序号从1开始）
        output_filename = f"tba_out_control_{i + 1}.txt"
        output_filepath = os.path.join(output_dir, output_filename)

        # 保存修改后的文件
        save_modified_file(modified_content, output_filepath)

        # 显示进度
        if (i + 1) % 10 == 0:
            print(f"已处理 {i + 1}/{pair_count} 个文件...")

    print(f"\n处理完成！共生成 {pair_count} 个文件。")
    print(f"文件保存在目录: {os.path.abspath(output_dir)}")


if __name__ == "__main__":
    main()