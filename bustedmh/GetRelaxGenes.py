#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
from typing import List, Dict, Set

# ===== 配置（可按需修改） =====
INPUT_ORTHO = "genespace.orthogroups.txt"         # genespace 文件（制表符或逗号分隔均可自动识别）
INPUT_RELAX = "sumT_p0.05.txt"         # RelaxTra 文件（通常是制表符 TSV）
TARGET_COLUMNS = ["Tridentiger_barbatus", "Tridentiger_radiatus"]  # 需要提取的物种列名
# ===========================


def detect_delimiter(path: str) -> str:
    """简单探测分隔符：优先选包含更多的那个（tab 或 comma）。"""
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        first = f.readline()
    if not first:
        return "\t"
    tabs = first.count("\t")
    commas = first.count(",")
    if tabs == 0 and commas == 0:
        return "\t"
    return "\t" if tabs >= commas else ","


def read_first_column_values(path: str) -> Set[str]:
    """读取文件第一列（跳过表头），返回集合。"""
    delim = detect_delimiter(path)
    values: Set[str] = set()
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f, delimiter=delim)
        _header = next(reader, None)  # 跳过表头
        for row in reader:
            if not row:
                continue
            v = row[0].strip()
            if v:
                values.add(v)
    return values


def extract_last_pipe_tokens(cell: str) -> List[str]:
    """
    从单元格中提取最后一个“|”后的内容。
    - 忽略 NA/空值
    - 若同一单元格含有用 ';' 分隔的多个条目，则全部处理
    例如: 'Tridentiger_radiatus|040076|FUN_044172-T1' -> ['FUN_044172-T1']
    """
    if not cell:
        return []
    cell = cell.strip()
    if not cell or cell.upper() == "NA":
        return []

    parts = [p.strip() for p in cell.split(";")] if ";" in cell else [cell]
    out: List[str] = []
    for part in parts:
        if not part or part.upper() == "NA":
            continue
        token = part.split("|")[-1].strip()
        if token and token.upper() != "NA":
            out.append(token)
    return out


def main() -> None:
    # 1) 读取 RelaxTra 的第一列（Group id）
    relax_ids = read_first_column_values(INPUT_RELAX)
    if not relax_ids:
        print("警告：在 RelaxTra 文件第一列未读到任何 ID（可能文件为空或分隔符识别不一致）。")

    # 2) 在 genespace 文件中匹配第一列（OrthoID），并从目标物种列提取最后一个“|”之后的字符串
    delim = detect_delimiter(INPUT_ORTHO)
    results: Dict[str, List[str]] = {col: [] for col in TARGET_COLUMNS}

    with open(INPUT_ORTHO, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f, delimiter=delim)
        header = next(reader, None)
        if not header:
            print("错误：genespace 文件为空或缺少表头。")
            return

        # 映射列名 -> 索引
        name_to_idx = {name.strip(): i for i, name in enumerate(header)}

        # 确定 OrthoID 列索引（优先使用名为 OrthoID 的列，否则用第一列）
        ortho_idx = 0
        for i, name in enumerate(header):
            if name.strip().lower() == "orthoid":
                ortho_idx = i
                break

        missing = [c for c in TARGET_COLUMNS if c not in name_to_idx]
        if missing:
            print(f"警告：genespace 文件中找不到如下列，将跳过：{', '.join(missing)}")

        for row in reader:
            if not row or len(row) <= ortho_idx:
                continue
            gid = row[ortho_idx].strip()
            if gid not in relax_ids:
                continue

            for col in TARGET_COLUMNS:
                idx = name_to_idx.get(col)
                if idx is None or idx >= len(row):
                    continue
                cell = row[idx]
                tokens = extract_last_pipe_tokens(cell)
                if tokens:
                    results[col].extend(tokens)

    # 3) 分列输出到对应文件名（每行一个）
    for col in TARGET_COLUMNS:
        outfile = f"{col}.txt"
        with open(outfile, "w", encoding="utf-8") as wf:
            for item in results[col]:
                wf.write(item + "\n")
        print(f"写入完成：{outfile}（{len(results[col])} 条）")


if __name__ == "__main__":
    main()
