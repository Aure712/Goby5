#!/usr/bin/env python3
# coding: utf-8

import re

INPUT_PATH = "pairwise_tbi_tra.simp.txt"
OUTPUT_PATH = "pairwise_tbi_tra.simp.sorted.txt"


def natural_sort_key(text):

    """
    将字符串转换为自然排序键
    例如: 'scaf2' 排在 'scaf10' 前面
    """

    def atoi(s):
        return int(s) if s.isdigit() else s

    return [atoi(c) for c in re.split(r'(\d+)', text)]


def main():
    lines = []

    with open(INPUT_PATH, "r") as in_fh:
        for line in in_fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                continue

            # 提取用于排序的字段
            qry_scaf = fields[4]
            qry_start = int(fields[5])
            qry_end = int(fields[6])

            # 保存：(排序键, 原始行)
            lines.append(((natural_sort_key(qry_scaf), qry_start, qry_end), line))

    # 排序：先按 qry_scaf 自然排序，再按 qry_start 数值，最后按 qry_end 数值
    lines.sort(key=lambda x: x[0])

    # 输出排序后的结果
    with open(OUTPUT_PATH, "w") as out_fh:
        for _, original_line in lines:
            out_fh.write(original_line + "\n")

    print(f"排序完成，共 {len(lines)} 行，已写入 {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
