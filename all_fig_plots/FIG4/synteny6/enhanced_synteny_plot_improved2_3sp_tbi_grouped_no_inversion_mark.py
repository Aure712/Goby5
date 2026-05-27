#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
三物种增强版生物信息学共线性图绘制脚本
物种顺序：Tra（上）- Tbi（中）- Tba（下）

主要变更：
1. 新增 Tbi 物种，位于 Tra 和 Tba 之间；
2. 共线性输入文件拆分为两个：
   - Tra-Tbi：文件中 chr1 = Tbi, chr2 = Tra
   - Tba-Tbi：文件中 chr1 = Tbi, chr2 = Tba
3. 新增 tbi_fai 输入；
4. chr_order 增加一列用于记录 Tbi 染色体顺序；
5. 可分别设置 Tbi 的颜色和反转；
6. 端粒信号仍保持原逻辑：只绘制 Tra 与 Tba 的端粒；
7. 断点箭头仍保持原逻辑：在 Tra 端粒层附近标注；
8. 新增 Tbi 分组映射：按 Tra 染色体顺序分组对齐，使 Tbi 与 Tra/Tba 横向布局一致。
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import seaborn as sns
import matplotlib.patheffects as path_effects
from matplotlib.collections import LineCollection
import matplotlib.transforms as transforms

# ==================== 路径与基础配置 ====================
DATA_DIR = "."
OUTPUT_DIR = "."

FILES = {
    'tra_fai': 'Tra_sort.fa.fai',
    'tbi_fai': 'Tbi_sort.fa.fai',
    'tba_fai': 'Tba_sort.fa.fai',

    # 注意：
    # 1) synteny_tbi_tra 文件中 chr1=Tbi, chr2=Tra
    # 2) synteny_tbi_tba 文件中 chr1=Tbi, chr2=Tba
    'synteny_tbi_tra': 'pairwise_tbi_tra.simp.txt',
    'synteny_tbi_tba': 'pairwise_tbi_tba.simp.txt',

    'tra_telomere': 'tidk_default_tra_onlychr.bed',
    'tba_telomere': 'tidk_default_tba_onlychr.bed',
    'chr_order': 'chr.order.txt',
    'breakpoints': 'lines.txt'
}

FIGURE_SIZE = (30, 16)
DPI = 1200

# ==================== chr_order 列名配置 ====================
# 如果你的 chr_order 文件列名不同，只改这里即可
CHR_ORDER_COLUMNS = {
    'tra': 'Tra_chr',      # Tra 染色体顺序列
    'tbi': 'Tbi_chr',      # Tbi 染色体顺序列（新增）
    'tba': 'Tba_chr'   # Tba 染色体顺序列
}

# ==================== 染色体绘制数量限制（None 表示不限制） ====================
TRA_CHR_LIMIT = 11
TBI_CHR_LIMIT = 22
TBA_CHR_LIMIT = 22

# ==================== 染色体边框配置 ====================
CHR_BORDER_CONFIG = {
    'color': '#000000',
    'linewidth': 4.0
}

# ==================== 染色体间隔配置 ====================
TRA_CHR_SPACING = 20000000
TBI_CHR_SPACING = 22450000
TBA_CHR_SPACING = 20000000

# Tbi 分组内间隔（同组内）
TBI_INTRA_GROUP_SPACING = 5000000

# Tbi 分组映射：按 Tra 染色体顺序分组对齐（新增，可选）
# 说明：
# - key 为 Tra 染色体名称；
# - value 为需要排在该 Tra 染色体下方附近的一组 Tbi 染色体；
# - 若不想启用 Tbi 分组对齐，可将 TBI_GROUPS_BY_TRA = {}。
TBI_GROUPS_BY_TRA = {
    'TraScf_1': ['TbiScf_14', 'TbiScf_15', 'TbiScf_13'],
    'TraScf_2': ['TbiScf_4', 'TbiScf_2'],
    'TraScf_3': ['TbiScf_16', 'TbiScf_1'],
    'TraScf_4': ['TbiScf_9', 'TbiScf_3'],
    'TraScf_5': ['TbiScf_7', 'TbiScf_12'],
    'TraScf_6': ['TbiScf_19', 'TbiScf_8'],
    'TraScf_7': ['TbiScf_6', 'TbiScf_21'],
    'TraScf_8': ['TbiScf_20', 'TbiScf_10'],
    'TraScf_9': ['TbiScf_18', 'TbiScf_17'],
    'TraScf_10': ['TbiScf_22', 'TbiScf_5'],
    'TraScf_11': ['TbiScf_11']
}

# Tba 分组内间隔（同组内）
TBA_INTRA_GROUP_SPACING = 5000000

# Tba 分组映射：按 Tra 染色体顺序分组对齐（保留原逻辑，可选）
# 若不想启用 Tba 分组对齐，可将 TBA_GROUPS_BY_TRA = {}。
TBA_GROUPS_BY_TRA = {
    'TraScf_1': ['TbaScf_16', 'TbaScf_15', 'TbaScf_12'],
    'TraScf_2': ['TbaScf_5', 'TbaScf_2'],
    'TraScf_3': ['TbaScf_17', 'TbaScf_1'],
    'TraScf_4': ['TbaScf_4', 'TbaScf_3'],
    'TraScf_5': ['TbaScf_11', 'TbaScf_13'],
    'TraScf_6': ['TbaScf_19', 'TbaScf_9'],
    'TraScf_7': ['TbaScf_6', 'TbaScf_20'],
    'TraScf_8': ['TbaScf_21', 'TbaScf_8'],
    'TraScf_9': ['TbaScf_14', 'TbaScf_18'],
    'TraScf_10': ['TbaScf_22', 'TbaScf_7'],
    'TraScf_11': ['TbaScf_10']
}

MIN_SYNTENY_SIZE = 10000
CHR_HEIGHT = 0.15
CHR_Y_SPACING = 1.2
TELOMERE_LINE_WIDTH = 2.0
TELOMERE_LINE_ALPHA = 0.85
TELOMERE_DOT_SIZE = 40

# ==================== 染色体反转配置 ====================
REVERSED_CHROMOSOMES = {
    'tra': [],
    'tbi': ['TbiScf_14', 'TbiScf_13', 'TbiScf_1', 'TbiScf_7', 'TbiScf_21', 'TbiScf_18',
            'TbiScf_22', 'TbiScf_5', 'TbiScf_11'],  # 新增：可在这里设置 Tbi 反转染色体
    'tba': ['TbaScf_16', 'TbaScf_2', 'TbaScf_4', 'TbaScf_3',
            'TbaScf_11', 'TbaScf_13', 'TbaScf_9', 'TbaScf_6', 'TbaScf_18']
}

# ==================== 染色体显示名映射（可选） ====================
CHR_DISPLAY_NAME = {
    # TRA
    'TraScf_1': 'Tra1',
    'TraScf_2': 'Tra2',
    'TraScf_3': 'Tra3',
    'TraScf_4': 'Tra4',
    'TraScf_5': 'Tra5',
    'TraScf_6': 'Tra6',
    'TraScf_7': 'Tra7',
    'TraScf_8': 'Tra8',
    'TraScf_9': 'Tra9',
    'TraScf_10': 'Tra10',
    'TraScf_11': 'Tra11',

    # TBI（按需启用/补充）
    'TbiScf_1': 'Tbi1',
    'TbiScf_2': 'Tbi2',
    'TbiScf_3': 'Tbi3',
    'TbiScf_4': 'Tbi4',
    'TbiScf_5': 'Tbi5',
    'TbiScf_6': 'Tbi6',
    'TbiScf_7': 'Tbi7',
    'TbiScf_8': 'Tbi8',
    'TbiScf_9': 'Tbi9',
    'TbiScf_10': 'Tbi10',
    'TbiScf_11': 'Tbi11',
    'TbiScf_12': 'Tbi12',
    'TbiScf_13': 'Tbi13',
    'TbiScf_14': 'Tbi14',
    'TbiScf_15': 'Tbi15',
    'TbiScf_16': 'Tbi16',
    'TbiScf_17': 'Tbi17',
    'TbiScf_18': 'Tbi18',
    'TbiScf_19': 'Tbi19',
    'TbiScf_20': 'Tbi20',
    'TbiScf_21': 'Tbi21',
    'TbiScf_22': 'Tbi22',

    # TBA
    'TbaScf_1': 'Tba1',
    'TbaScf_2': 'Tba2',
    'TbaScf_3': 'Tba3',
    'TbaScf_4': 'Tba4',
    'TbaScf_5': 'Tba5',
    'TbaScf_6': 'Tba6',
    'TbaScf_7': 'Tba7',
    'TbaScf_8': 'Tba8',
    'TbaScf_9': 'Tba9',
    'TbaScf_10': 'Tba10',
    'TbaScf_11': 'Tba11',
    'TbaScf_12': 'Tba12',
    'TbaScf_13': 'Tba13',
    'TbaScf_14': 'Tba14',
    'TbaScf_15': 'Tba15',
    'TbaScf_16': 'Tba16',
    'TbaScf_17': 'Tba17',
    'TbaScf_18': 'Tba18',
    'TbaScf_19': 'Tba19',
    'TbaScf_20': 'Tba20',
    'TbaScf_21': 'Tba21',
    'TbaScf_22': 'Tba22'
}

def get_display_name(chr_name):
    return CHR_DISPLAY_NAME.get(chr_name, chr_name)

# 反转标记配置：在名称后追加“*”
REVERSE_LABEL_MARK = {
    'show': True,
    'char': '*',
    'color': 'red',
    'fontsize': 20,
    'scale': 1.0,
    'gap_pts': 5.0,
    'legend_markersize': 14
}

# ==================== 染色体标签位置与字号配置 ====================
LABEL_CONFIG = {
    'bottom_offset': 0.35,
    'top_offset': 0.40,
    'rotation': 45,
    'fontsize': 18,
    'fontweight': 'bold'
}

LABEL_SIZE_PER_CHR = {}
STAR_SIZE_PER_CHR = {}

# ==================== 断点箭头（Tra 端粒层附近的标注）配置 ====================
ARROW_MARK_CONFIG = {
    'color': '#FF00A6',
    'apply_to': 'tra',         # 保持原逻辑：默认只在 Tra 端粒层上画
    'mutation_scale': 40.0,
    'linewidth': 6.0,
    'offset': 0.15,
    'zorder': 11,
    'draw_inside_chr': False,
    'inside_color': '#000000',
    'inside_width': 2.5,

    'show_labels': True,
    'label_color': '#000000',
    'label_fontsize': 24,
    'label_offset': 0.15,
    'label_start': 1
}

# ==================== Tra 分段区域“无 inversion”标记配置 ====================
# 用法：
# - 在 segments 中按 Tra 染色体名称指定需要标记的分段编号；
# - 分段编号为 1-based，按该 Tra 染色体原始坐标从左到右计算：
#   例如有 2 个断点时会形成 3 段，对应编号 1、2、3；
# - 若某条染色体所有 Tra 分段都需要标记，可写成 'all'；
# - 默认不预设具体分段，避免误标；需要时取消下面示例注释并修改即可。
NO_INVERSION_MARK_CONFIG = {
    'show': True,
    'segments': {
        'TraScf_1': [2],
        'TraScf_5': [2],
        'TraScf_8': [2],
        'TraScf_10': [2]
    },
    'marker': 'D',
    'size': 420,
    'facecolor': '#FFFFFF',
    'edgecolor': '#000000',
    'linewidth': 2.0,
    'text': 'NI',
    'text_color': '#000000',
    'fontsize': 11,
    'fontweight': 'bold',
    'y_offset': 0.12,
    'zorder': 13,
    'legend_markersize': 11
}

# ==================== 端粒垂直线配置参数 ====================
TELOMERE_CONFIG = {
    'mode': 'log_quantile',
    'gamma': 2.5,
    'min_height': 0.0008,
    'max_height': 0.30,
    'base_height': 0.02,
    'log_scale_factor': 2.0,
    'show_dots': False
}

COLORS = {
    'tra_chr_base': '#CCCCCC',
    'tbi_chr': sns.color_palette("Set2", 8) + sns.color_palette("Pastel1", 4),
    'tba_chr': sns.color_palette("tab20", 20) + sns.color_palette("Set3", 2),
    'synteny_alpha': 0.3,
    'telomere_tra': '#FF6B35',
    'telomere_tba': '#4A90E2',
    'breakpoint': '#9B59B6',
    'breakpoint_line': '#95A5A6',
    'background': '#FAFAFA',
    'edge_color': '#2C3E50'
}

# ==================== 图表文字配置 ====================
TEXT_CONFIG = {
    'species_label_tra': 'Tra',
    'species_label_tbi': 'Tbi',
    'species_label_tba': 'Tba',
    'species_label_fontsize': 28,
    'title_main': 'Chromosome Synteny Analysis',
    'title_sub': 'Tra - Tbi - Tba with Telomere Distribution and Breakpoint-based Segmentation',
    'title_fontsize': 20,
    'legend_breakpoint': 'Breakpoints',
    'legend_no_inversion': 'NI: no inversion near breakpoint',
    'legend_tra_telomere': 'Tra telomeres',
    'legend_tba_telomere': 'Tba telomeres',
    'legend_tbi_chr': 'Tbi Chr',
    'legend_tba_chr': 'Tba Chr',
    'legend_more_tbi_chrs': '... and {} more Tbi chrs',
    'legend_more_tba_chrs': '... and {} more Tba chrs',
    'legend_synteny_tra_tbi': 'Tra-Tbi synteny (n={})',
    'legend_synteny_tbi_tba': 'Tbi-Tba synteny (n={})',
    'legend_reversed': 'Reversed chromosomes (*)',
    'legend_fontsize': 12,
    'output_filename_base': 'enhanced_synteny_plot_tra_tbi_tba',
    'stats_title': '三物种优化版生物信息学共线性分析结果',
}

# ==================== 自定义染色体区域颜色配置（可选） ====================
# Tra 分段颜色（保留原逻辑，可按断点分段后指定）
CUSTOM_TRA_COLORS = {
    'TraScf_1': ['#aec7e8', '#ffbb78', '#17becf'],
    'TraScf_2': ['#f7b6d2', '#bcbd22'],
    'TraScf_3': ['#c49c94', '#ff7f0e'],
    'TraScf_4': ['#1f77b4', '#ff9896'],
    'TraScf_5': ['#393b79', '#c7c7c7'],
    'TraScf_6': ['#8c564b', '#2ca02c'],
    'TraScf_7': ['#98df8a', '#c5b0d5'],
    'TraScf_8': ['#637939', '#e377c2'],
    'TraScf_9': ['#dbdb8d', '#9467bd'],
    'TraScf_10': ['#9edae5', '#d62728'],
    'TraScf_11': ['#7f7f7f']
}

# Tbi 染色体颜色（新增，可自行配置）
CUSTOM_TBI_COLORS = {
    'TbiScf_14': ['#aec7e8'],
    'TbiScf_15': ['#ffbb78'],
    'TbiScf_13': ['#17becf'],
    'TbiScf_4':  ['#f7b6d2'],
    'TbiScf_2':  ['#bcbd22'],
    'TbiScf_16': ['#c49c94'],
    'TbiScf_1':  ['#ff7f0e'],
    'TbiScf_9':  ['#1f77b4'],
    'TbiScf_3':  ['#ff9896'],
    'TbiScf_7': ['#393b79'],
    'TbiScf_12': ['#c7c7c7'],
    'TbiScf_19':  ['#8c564b'],
    'TbiScf_8': ['#2ca02c'],
    'TbiScf_6':  ['#98df8a'],
    'TbiScf_21': ['#c5b0d5'],
    'TbiScf_20': ['#637939'],
    'TbiScf_10':  ['#e377c2'],
    'TbiScf_18': ['#dbdb8d'],
    'TbiScf_17': ['#9467bd'],
    'TbiScf_22': ['#9edae5'],
    'TbiScf_5':  ['#d62728'],
    'TbiScf_11': ['#7f7f7f']
}

# Tba 染色体颜色（保留原逻辑）
CUSTOM_TBA_COLORS = {
    'TbaScf_16': ['#aec7e8'],
    'TbaScf_15': ['#ffbb78'],
    'TbaScf_12': ['#17becf'],
    'TbaScf_5':  ['#f7b6d2'],
    'TbaScf_2':  ['#bcbd22'],
    'TbaScf_17': ['#c49c94'],
    'TbaScf_1':  ['#ff7f0e'],
    'TbaScf_4':  ['#1f77b4'],
    'TbaScf_3':  ['#ff9896'],
    'TbaScf_11': ['#393b79'],
    'TbaScf_13': ['#c7c7c7'],
    'TbaScf_9':  ['#2ca02c'],
    'TbaScf_19': ['#8c564b'],
    'TbaScf_6':  ['#98df8a'],
    'TbaScf_20': ['#c5b0d5'],
    'TbaScf_21': ['#637939'],
    'TbaScf_8':  ['#e377c2'],
    'TbaScf_14': ['#dbdb8d'],
    'TbaScf_18': ['#9467bd'],
    'TbaScf_22': ['#9edae5'],
    'TbaScf_7':  ['#d62728'],
    'TbaScf_10': ['#7f7f7f']
}

plt.rcParams.update({
    'font.size': 13,
    'font.family': ['Arial', 'DejaVu Sans', 'sans-serif'],
    'axes.linewidth': 1.8,
    'lines.linewidth': 2.5,
    'patch.linewidth': 1.5,
    'figure.facecolor': 'white',
    'axes.facecolor': COLORS['background'],
    'savefig.facecolor': 'white',
    'savefig.edgecolor': 'none'
})

# ==================== 数据加载与预处理 ====================

def load_all_data():
    print("正在加载数据...")

    tra_chr = pd.read_csv(f"{DATA_DIR}/{FILES['tra_fai']}", sep='\t', header=None,
                          names=['chr', 'length', 'offset', 'line_bases', 'line_width'])
    tbi_chr = pd.read_csv(f"{DATA_DIR}/{FILES['tbi_fai']}", sep='\t', header=None,
                          names=['chr', 'length', 'offset', 'line_bases', 'line_width'])
    tba_chr = pd.read_csv(f"{DATA_DIR}/{FILES['tba_fai']}", sep='\t', header=None,
                          names=['chr', 'length', 'offset', 'line_bases', 'line_width'])

    # synteny_tbi_tra: chr1=Tbi, chr2=Tra
    synteny_tbi_tra = pd.read_csv(f"{DATA_DIR}/{FILES['synteny_tbi_tra']}", sep='\t', header=None,
                                  names=['chr1', 'start1', 'end1', 'strand', 'chr2', 'start2', 'end2', 'score'])

    # synteny_tbi_tba: chr1=Tbi, chr2=Tba
    synteny_tbi_tba = pd.read_csv(f"{DATA_DIR}/{FILES['synteny_tbi_tba']}", sep='\t', header=None,
                                  names=['chr1', 'start1', 'end1', 'strand', 'chr2', 'start2', 'end2', 'score'])

    tra_telomere = pd.read_csv(f"{DATA_DIR}/{FILES['tra_telomere']}", sep='\t', header=None,
                               names=['chr', 'start', 'end', 'count'])
    tba_telomere = pd.read_csv(f"{DATA_DIR}/{FILES['tba_telomere']}", sep='\t', header=None,
                               names=['chr', 'start', 'end', 'count'])

    chr_order = pd.read_csv(f"{DATA_DIR}/{FILES['chr_order']}", sep='\t')
    breakpoints = pd.read_csv(f"{DATA_DIR}/{FILES['breakpoints']}", sep='\t', header=None,
                              names=['chr', 'pos', 'start_range', 'end_range'])

    print("数据加载完成！")
    return (tra_chr, tbi_chr, tba_chr,
            synteny_tbi_tra, synteny_tbi_tba,
            tra_telomere, tba_telomere,
            chr_order, breakpoints)


def filter_synteny_by_size(synteny, min_size=MIN_SYNTENY_SIZE, label="synteny"):
    print(f"正在过滤 {label} 共线性数据，最小窗口大小: {min_size:,} bp")
    window_size1 = synteny['end1'] - synteny['start1']
    window_size2 = synteny['end2'] - synteny['start2']
    min_window_size = np.minimum(window_size1, window_size2)
    original_count = len(synteny)
    filtered_synteny = synteny[min_window_size >= min_size].copy()
    filtered_synteny['min_window_size'] = min_window_size[min_window_size >= min_size]
    filtered_count = len(filtered_synteny)

    print(f"{label} 原始共线性区域: {original_count:,}")
    print(f"{label} 过滤后共线性区域: {filtered_count:,}")
    if original_count > 0:
        print(f"{label} 过滤掉: {original_count - filtered_count:,} ({(original_count - filtered_count) / original_count * 100:.1f}%)")
    return filtered_synteny


def _get_ordered_list_from_column(df, col_name):
    if col_name not in df.columns:
        raise KeyError(f"chr_order 文件中未找到列：{col_name}")
    return df[~df[col_name].isna() & (df[col_name] != '')][col_name].astype(str).tolist()


def _build_tbi_order(tra_ordered, chr_order, synteny_tbi_tra, synteny_tbi_tba):
    """按 Tra 染色体顺序构建 TBI 排列；若未配置则按 chr_order 的 TBI 列顺序。"""
    if TBI_GROUPS_BY_TRA:
        ordered = []
        for tra_chr_name in tra_ordered:
            ordered.extend(TBI_GROUPS_BY_TRA.get(tra_chr_name, []))

        fallback = _get_ordered_list_from_column(chr_order, CHR_ORDER_COLUMNS['tbi'])
        # 两个共线性文件中 chr1 均为 Tbi，因此都加入兜底列表，避免遗漏。
        extra = fallback + list(synteny_tbi_tra['chr1'].unique()) + list(synteny_tbi_tba['chr1'].unique())
        for name in extra:
            if name not in ordered:
                ordered.append(name)
        return ordered

    return _get_ordered_list_from_column(chr_order, CHR_ORDER_COLUMNS['tbi'])


def _build_tba_order(tra_ordered, chr_order, synteny_tbi_tba):
    """保留原有 TBA 分组逻辑；若未配置则按 chr_order 的 TBA 列顺序。"""
    if TBA_GROUPS_BY_TRA:
        ordered = []
        for tra_chr_name in tra_ordered:
            ordered.extend(TBA_GROUPS_BY_TRA.get(tra_chr_name, []))

        fallback = _get_ordered_list_from_column(chr_order, CHR_ORDER_COLUMNS['tba'])
        extra = fallback + list(synteny_tbi_tba['chr2'].unique())  # chr2=Tba
        for name in extra:
            if name not in ordered:
                ordered.append(name)
        return ordered

    return _get_ordered_list_from_column(chr_order, CHR_ORDER_COLUMNS['tba'])


def prepare_chromosome_data(tra_chr, tbi_chr, tba_chr, chr_order, synteny_tbi_tra, synteny_tbi_tba):
    print("准备三物种染色体数据...")

    tra_ordered = _get_ordered_list_from_column(chr_order, CHR_ORDER_COLUMNS['tra'])
    tbi_ordered = _build_tbi_order(tra_ordered, chr_order, synteny_tbi_tra, synteny_tbi_tba)
    tba_ordered = _build_tba_order(tra_ordered, chr_order, synteny_tbi_tba)

    if TRA_CHR_LIMIT is not None:
        tra_ordered = tra_ordered[:TRA_CHR_LIMIT]
    if TBI_CHR_LIMIT is not None:
        tbi_ordered = tbi_ordered[:TBI_CHR_LIMIT]
    if TBA_CHR_LIMIT is not None:
        tba_ordered = tba_ordered[:TBA_CHR_LIMIT]

    # 从共线性文件中补充可能出现但不在顺序表前 N 个中的染色体
    tra_keep = list(dict.fromkeys(tra_ordered + list(synteny_tbi_tra['chr2'].unique())))   # chr2=Tra
    tbi_keep = list(dict.fromkeys(tbi_ordered + list(synteny_tbi_tra['chr1'].unique()) + list(synteny_tbi_tba['chr1'].unique())))  # chr1=Tbi
    tba_keep = list(dict.fromkeys(tba_ordered + list(synteny_tbi_tba['chr2'].unique())))   # chr2=Tba

    tra_filtered = tra_chr[tra_chr['chr'].isin(tra_keep)].copy()
    tbi_filtered = tbi_chr[tbi_chr['chr'].isin(tbi_keep)].copy()
    tba_filtered = tba_chr[tba_chr['chr'].isin(tba_keep)].copy()

    tra_order_dict = {chr_name: i for i, chr_name in enumerate(tra_ordered)}
    tbi_order_dict = {chr_name: i for i, chr_name in enumerate(tbi_ordered)}
    tba_order_dict = {chr_name: i for i, chr_name in enumerate(tba_ordered)}

    tra_filtered['order'] = tra_filtered['chr'].map(tra_order_dict)
    tbi_filtered['order'] = tbi_filtered['chr'].map(tbi_order_dict)
    tba_filtered['order'] = tba_filtered['chr'].map(tba_order_dict)

    tra_filtered = tra_filtered.dropna(subset=['order']).sort_values('order')
    tbi_filtered = tbi_filtered.dropna(subset=['order']).sort_values('order')
    tba_filtered = tba_filtered.dropna(subset=['order']).sort_values('order')

    def calc_positions_vectorized(df, spacing):
        df = df.copy()
        lengths = df['length'].values
        cumulative_with_spacing = np.concatenate([[0], np.cumsum(lengths[:-1] + spacing)])
        df['x_start'] = cumulative_with_spacing
        df['x_end'] = cumulative_with_spacing + lengths
        df['x_center'] = (df['x_start'] + df['x_end']) / 2
        return df

    def calc_positions_grouped(df, tra_ordered_list, group_map, group_spacing, intra_spacing):
        lengths = dict(zip(df['chr'], df['length']))
        rows = []
        x_cursor = 0.0
        used = set()

        for tra_name in tra_ordered_list:
            tba_list = group_map.get(tra_name, [])
            group_has = False
            for chr_name in tba_list:
                if chr_name not in lengths:
                    continue
                length = lengths[chr_name]
                rows.append({
                    'chr': chr_name,
                    'length': length,
                    'x_start': x_cursor,
                    'x_end': x_cursor + length,
                    'x_center': x_cursor + length / 2
                })
                used.add(chr_name)
                x_cursor = x_cursor + length + intra_spacing
                group_has = True
            if group_has:
                x_cursor = x_cursor - intra_spacing + group_spacing

        for chr_name in df['chr'].tolist():
            if chr_name in used:
                continue
            length = lengths[chr_name]
            rows.append({
                'chr': chr_name,
                'length': length,
                'x_start': x_cursor,
                'x_end': x_cursor + length,
                'x_center': x_cursor + length / 2
            })
            x_cursor = x_cursor + length + group_spacing

        return pd.DataFrame(rows)

    tra_pos = calc_positions_vectorized(tra_filtered, TRA_CHR_SPACING)

    if TBI_GROUPS_BY_TRA:
        tbi_pos = calc_positions_grouped(
            tbi_filtered, tra_ordered, TBI_GROUPS_BY_TRA,
            TBI_CHR_SPACING, TBI_INTRA_GROUP_SPACING
        )
        tbi_pos['order'] = tbi_pos['chr'].map(tbi_order_dict)
        tbi_pos = tbi_pos.sort_values('order')
    else:
        tbi_pos = calc_positions_vectorized(tbi_filtered, TBI_CHR_SPACING)

    if TBA_GROUPS_BY_TRA:
        tba_pos = calc_positions_grouped(
            tba_filtered, tra_ordered, TBA_GROUPS_BY_TRA,
            TBA_CHR_SPACING, TBA_INTRA_GROUP_SPACING
        )
        tba_pos['order'] = tba_pos['chr'].map(tba_order_dict)
        tba_pos = tba_pos.sort_values('order')
    else:
        tba_pos = calc_positions_vectorized(tba_filtered, TBA_CHR_SPACING)

    if TRA_CHR_LIMIT is not None:
        tra_pos = tra_pos.head(TRA_CHR_LIMIT)
    if TBI_CHR_LIMIT is not None:
        tbi_pos = tbi_pos.head(TBI_CHR_LIMIT)
    if TBA_CHR_LIMIT is not None:
        tba_pos = tba_pos.head(TBA_CHR_LIMIT)

    tra_pos['is_reversed'] = tra_pos['chr'].isin(REVERSED_CHROMOSOMES.get('tra', []))
    tbi_pos['is_reversed'] = tbi_pos['chr'].isin(REVERSED_CHROMOSOMES.get('tbi', []))
    tba_pos['is_reversed'] = tba_pos['chr'].isin(REVERSED_CHROMOSOMES.get('tba', []))

    if tbi_pos['is_reversed'].sum() > 0:
        print(f"TBI反转: {list(tbi_pos[tbi_pos['is_reversed']]['chr'])}")
    if tba_pos['is_reversed'].sum() > 0:
        print(f"TBA反转: {list(tba_pos[tba_pos['is_reversed']]['chr'])}")
    if tra_pos['is_reversed'].sum() > 0:
        print(f"TRA反转: {list(tra_pos[tra_pos['is_reversed']]['chr'])}")

    print(f"处理了 TRA/TBI/TBA 染色体数: {len(tra_pos)}/{len(tbi_pos)}/{len(tba_pos)}")
    print(f"反转数: TRA={tra_pos['is_reversed'].sum()}, TBI={tbi_pos['is_reversed'].sum()}, TBA={tba_pos['is_reversed'].sum()}")

    return tra_pos, tbi_pos, tba_pos


def determine_tra_segments_optimized(tra_pos, breakpoints, synteny_tbi_tra, tbi_colors, custom_colors=None):
    """
    仍只对 Tra 染色体做断点分段着色。
    但由于 Tra-Tbi 文件中 chr1=Tbi, chr2=Tra，
    所以这里要根据 Tra(chr2) 与 Tbi(chr1) 的共线性关系给 Tra 染色体分段着色。
    """
    print("正在分析 Tra 染色体断点分段和着色...")
    chromosome_segments = {}

    # 这里按 Tra 染色体分组：Tra 在 chr2
    synteny_by_chr = synteny_tbi_tra.groupby('chr2')

    for _, chr_row in tra_pos.iterrows():
        chr_name = chr_row['chr']
        chr_length = chr_row['length']
        chr_x_start = chr_row['x_start']
        is_reversed = chr_row.get('is_reversed', False)

        chr_breakpoints = breakpoints[breakpoints['chr'] == chr_name].sort_values('pos')
        segments = []

        has_custom_colors = custom_colors and chr_name in custom_colors
        custom_color_list = custom_colors.get(chr_name, []) if custom_colors else []

        if len(chr_breakpoints) == 0:
            if has_custom_colors and len(custom_color_list) > 0:
                color = custom_color_list[0]
                main_tbi_chr = None
            elif chr_name in synteny_by_chr.groups:
                chr_synteny = synteny_by_chr.get_group(chr_name)
                tbi_chr_counts = chr_synteny['chr1'].value_counts()  # chr1 = Tbi
                main_tbi_chr = tbi_chr_counts.index[0]
                color = tbi_colors.get(main_tbi_chr, list(tbi_colors.values())[0])
            else:
                main_tbi_chr = None
                color = list(tbi_colors.values())[0] if tbi_colors else COLORS['tbi_chr'][0]

            segments.append({
                'start': 0,
                'end': chr_length,
                'x_start': chr_x_start,
                'x_end': chr_x_start + chr_length,
                'color': color,
                'tbi_chr': main_tbi_chr,
                'is_reversed': is_reversed,
                'segment_index': 1,   # 1-based，按 Tra 原始坐标顺序编号
                'plot_index': 1       # 1-based，按图中绘制顺序编号
            })

        else:
            positions = np.concatenate([[0], chr_breakpoints['pos'].values, [chr_length]])
            chr_synteny = synteny_by_chr.get_group(chr_name) if chr_name in synteny_by_chr.groups else pd.DataFrame()
            segment_indices = range(len(positions) - 1)
            if is_reversed:
                segment_indices = reversed(list(segment_indices))

            for i in segment_indices:
                seg_start = positions[i]
                seg_end = positions[i + 1]
                seg_middle = (seg_start + seg_end) / 2

                if is_reversed:
                    actual_x_start = chr_x_start + chr_length - seg_end
                    actual_x_end = chr_x_start + chr_length - seg_start
                else:
                    actual_x_start = chr_x_start + seg_start
                    actual_x_end = chr_x_start + seg_end

                if has_custom_colors and i < len(custom_color_list):
                    color = custom_color_list[i]
                    main_tbi_chr = None
                else:
                    if not chr_synteny.empty:
                        # Tra 在 chr2，所以要用 start2/end2 判断该段中点落在哪些共线性块上
                        mask = (chr_synteny['start2'] <= seg_middle) & (chr_synteny['end2'] >= seg_middle)
                        seg_synteny = chr_synteny[mask]
                        if len(seg_synteny) > 0:
                            tbi_chr_counts = seg_synteny['chr1'].value_counts()
                            main_tbi_chr = tbi_chr_counts.index[0]
                            color = tbi_colors.get(main_tbi_chr, COLORS['tbi_chr'][i % len(COLORS['tbi_chr'])])
                        else:
                            main_tbi_chr = None
                            color = COLORS['tbi_chr'][i % len(COLORS['tbi_chr'])]
                    else:
                        main_tbi_chr = None
                        color = COLORS['tbi_chr'][i % len(COLORS['tbi_chr'])]

                segments.append({
                    'start': seg_start,
                    'end': seg_end,
                    'x_start': actual_x_start,
                    'x_end': actual_x_end,
                    'color': color,
                    'tbi_chr': main_tbi_chr,
                    'is_reversed': is_reversed,
                    'segment_index': i + 1,             # 1-based，按 Tra 原始坐标顺序编号
                    'plot_index': len(segments) + 1     # 1-based，按图中绘制顺序编号
                })

        chromosome_segments[chr_name] = segments

    return chromosome_segments

# ==================== 文本绘制（支持每染色体字号与星标字号） ====================

def _draw_label_with_optional_star(ax, x_center, y, raw_chr_name, display_chr_name, position, is_reversed):
    rotation = LABEL_CONFIG['rotation']
    fontsize_label = LABEL_SIZE_PER_CHR.get(raw_chr_name, LABEL_CONFIG['fontsize'])
    fontweight = LABEL_CONFIG['fontweight']

    if position == 'bottom':
        va = 'top'
        ha = 'right'
    else:
        va = 'bottom'
        ha = 'left'

    main_text = ax.text(
        x_center, y, display_chr_name,
        rotation=rotation, ha=ha, va=va,
        fontsize=fontsize_label, fontweight=fontweight,
        path_effects=[path_effects.withStroke(linewidth=2.5, foreground='white')],
        zorder=6
    )

    if is_reversed and REVERSE_LABEL_MARK.get('show', True):
        star_char = str(REVERSE_LABEL_MARK.get('char', '*'))
        star_color = REVERSE_LABEL_MARK.get('color', 'red')

        star_fontsize = STAR_SIZE_PER_CHR.get(raw_chr_name, None)
        if star_fontsize is None:
            star_fontsize = REVERSE_LABEL_MARK.get('fontsize', None)
        if star_fontsize is None:
            star_fontsize = fontsize_label * float(REVERSE_LABEL_MARK.get('scale', 1.0))

        fig = ax.figure
        renderer = fig.canvas.get_renderer()
        if renderer is None:
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()

        gap_pts = float(REVERSE_LABEL_MARK.get('gap_pts', 5.0))
        dx_pts = gap_pts

        if ha == 'left':
            bbox = main_text.get_window_extent(renderer=renderer)
            width_px = bbox.width
            dx_pts += width_px * 72.0 / fig.dpi

        offset = transforms.ScaledTranslation(dx_pts / 72.0, 0.0, fig.dpi_scale_trans)
        star_transform = ax.transData + offset

        ax.text(
            x_center, y, star_char,
            rotation=rotation, ha='left', va=va,
            fontsize=star_fontsize, fontweight=fontweight,
            color=star_color,
            path_effects=[path_effects.withStroke(linewidth=2.0, foreground='white')],
            transform=star_transform,
            zorder=6
        )

# ==================== 绘图函数 ====================

def _get_no_inversion_selection(chr_name):
    selected = NO_INVERSION_MARK_CONFIG.get('segments', {}).get(chr_name, [])
    if selected is None:
        return set()

    if isinstance(selected, str):
        if selected.lower() == 'all':
            return 'all'
        selected = [selected]

    try:
        return {int(x) for x in selected}
    except (TypeError, ValueError):
        return set()


def _is_no_inversion_segment_marked(chr_name, segment):
    if not NO_INVERSION_MARK_CONFIG.get('show', True):
        return False

    selected = _get_no_inversion_selection(chr_name)
    if selected == 'all':
        return True
    if not selected:
        return False

    return int(segment.get('segment_index', -1)) in selected


def count_configured_no_inversion_segments(chr_segments):
    count = 0
    for chr_name, segments in chr_segments.items():
        count += sum(_is_no_inversion_segment_marked(chr_name, seg) for seg in segments)
    return count


def draw_no_inversion_segment_markers(ax, segments, y_center, height, chr_name):
    if not NO_INVERSION_MARK_CONFIG.get('show', True):
        return

    marked_segments = [
        seg for seg in segments
        if _is_no_inversion_segment_marked(chr_name, seg)
    ]
    if not marked_segments:
        return

    conf = NO_INVERSION_MARK_CONFIG
    marker_y = y_center + height / 2 + float(conf.get('y_offset', 0.10))

    for segment in marked_segments:
        marker_x = (segment['x_start'] + segment['x_end']) / 2

        ax.scatter(
            [marker_x], [marker_y],
            marker=conf.get('marker', 'D'),
            s=conf.get('size', 360),
            facecolors=conf.get('facecolor', '#FFFFFF'),
            edgecolors=conf.get('edgecolor', '#000000'),
            linewidths=conf.get('linewidth', 1.8),
            zorder=conf.get('zorder', 13),
            clip_on=False
        )

        marker_text = str(conf.get('text', '')).strip()
        if marker_text:
            ax.text(
                marker_x, marker_y, marker_text,
                ha='center', va='center',
                fontsize=conf.get('fontsize', 10),
                fontweight=conf.get('fontweight', 'bold'),
                color=conf.get('text_color', '#000000'),
                zorder=conf.get('zorder', 13) + 1,
                clip_on=False,
                path_effects=[path_effects.withStroke(linewidth=1.2, foreground='white')]
            )

def draw_segmented_chromosome(ax, segments, y_center, height, chr_name, position='bottom', is_reversed=False):
    draw_internal_edges = ARROW_MARK_CONFIG.get('draw_inside_chr', True) is True
    seg_edgecolor = CHR_BORDER_CONFIG.get('color', COLORS['edge_color']) if draw_internal_edges else 'none'
    seg_linewidth = CHR_BORDER_CONFIG.get('linewidth', 1.5) if draw_internal_edges else 0.0

    for segment in segments:
        chr_segment = FancyBboxPatch(
            (segment['x_start'], y_center - height / 2),
            segment['x_end'] - segment['x_start'], height,
            boxstyle="round,pad=0.01",
            facecolor=segment['color'],
            edgecolor=seg_edgecolor,
            linewidth=seg_linewidth,
            alpha=0.85,
            zorder=5
        )
        ax.add_patch(chr_segment)

    total_start = min(seg['x_start'] for seg in segments)
    total_end = max(seg['x_end'] for seg in segments)
    x_center = (total_start + total_end) / 2

    if not draw_internal_edges:
        outline = FancyBboxPatch(
            (total_start, y_center - height / 2),
            total_end - total_start, height,
            boxstyle="round,pad=0.01",
            facecolor='none',
            edgecolor=CHR_BORDER_CONFIG.get('color', COLORS['edge_color']),
            linewidth=CHR_BORDER_CONFIG.get('linewidth', 1.5),
            zorder=6
        )
        ax.add_patch(outline)

    draw_no_inversion_segment_markers(ax, segments, y_center, height, chr_name)

    label_y = y_center - height / 2 - LABEL_CONFIG['bottom_offset'] if position == 'bottom' \
        else y_center + height / 2 + LABEL_CONFIG['top_offset']

    _draw_label_with_optional_star(
        ax, x_center, label_y,
        raw_chr_name=chr_name,
        display_chr_name=get_display_name(chr_name),
        position=position,
        is_reversed=is_reversed
    )


def draw_chromosome(ax, x_start, x_end, y_center, height, color, label, position='bottom', is_reversed=False):
    chromosome = FancyBboxPatch(
        (x_start, y_center - height / 2), x_end - x_start, height,
        boxstyle="round,pad=0.01",
        facecolor=color,
        edgecolor=CHR_BORDER_CONFIG.get('color', COLORS['edge_color']),
        linewidth=CHR_BORDER_CONFIG.get('linewidth', 1.5),
        alpha=0.85,
        zorder=5
    )
    ax.add_patch(chromosome)

    x_center = (x_start + x_end) / 2
    label_y = y_center - height / 2 - LABEL_CONFIG['bottom_offset'] if position == 'bottom' \
        else y_center + height / 2 + LABEL_CONFIG['top_offset']

    _draw_label_with_optional_star(
        ax, x_center, label_y,
        raw_chr_name=label,
        display_chr_name=get_display_name(label),
        position=position,
        is_reversed=is_reversed
    )


def _compute_telomere_heights(counts):
    mode = TELOMERE_CONFIG.get('mode', 'log_quantile')
    min_h = TELOMERE_CONFIG.get('min_height', 0.0008)
    max_h = TELOMERE_CONFIG.get('max_height', 0.35)

    counts = np.asarray(counts, dtype=float)
    counts = np.maximum(1.0, counts)

    if mode == 'legacy':
        base = TELOMERE_CONFIG.get('base_height', 0.02)
        factor = TELOMERE_CONFIG.get('log_scale_factor', 2.0)
        heights = base * (1.0 + np.log10(counts) * factor)
        heights = np.clip(heights, min_h, max_h)
        return heights

    c_log = np.log10(counts)
    cmin = float(np.min(c_log))
    cmax = float(np.max(c_log))
    if cmax > cmin:
        norm = (c_log - cmin) / (cmax - cmin)
    else:
        norm = np.zeros_like(c_log)
    gamma = float(TELOMERE_CONFIG.get('gamma', 3.0))
    norm = np.power(norm, gamma)
    heights = min_h + norm * (max_h - min_h)
    return heights


def draw_telomeres_as_lines_optimized(ax, telomere_data, chr_positions, y_position, color, species_name):
    tel_plot = telomere_data.merge(chr_positions[['chr', 'x_start', 'length', 'is_reversed']],
                                   on='chr', how='inner')
    if tel_plot.empty:
        return 0

    def calc_telomere_x(row):
        if row['is_reversed']:
            tel_start_from_right = row['length'] - row['end']
            tel_end_from_right = row['length'] - row['start']
            return row['x_start'] + (tel_start_from_right + tel_end_from_right) / 2
        else:
            return row['x_start'] + (row['start'] + row['end']) / 2

    tel_centers = tel_plot.apply(calc_telomere_x, axis=1).values
    counts = tel_plot['count'].values

    line_heights = _compute_telomere_heights(counts)
    dot_sizes = np.minimum(TELOMERE_DOT_SIZE * (1 + np.log10(np.maximum(1, counts))), 120)

    if len(counts) > 0:
        print(f"  {species_name} 端粒统计: Count范围 {counts.min():.0f}-{counts.max():.0f}, 高度范围 {line_heights.min():.4f}-{line_heights.max():.4f}")

    segments = [[(xc, y_position - h / 2), (xc, y_position + h / 2)]
                for xc, h in zip(tel_centers, line_heights)]

    lc = LineCollection(
        segments,
        colors=color,
        linewidths=TELOMERE_LINE_WIDTH * 1.5,
        alpha=TELOMERE_LINE_ALPHA,
        zorder=10,
        capstyle='round'
    )
    ax.add_collection(lc)

    if TELOMERE_CONFIG.get('show_dots', False):
        ax.scatter(tel_centers, [y_position] * len(tel_centers),
                   s=dot_sizes, c=color, alpha=TELOMERE_LINE_ALPHA * 0.8,
                   edgecolors='white', linewidths=1.5, zorder=11, marker='o')

    return len(tel_plot)


def draw_breakpoint_arrows(ax, x_positions, telomere_y, species, conf, labels=None):
    if x_positions.size == 0:
        return

    offset = float(conf.get('offset', 0.08))
    ms = float(conf.get('mutation_scale', 22.0))
    lw = float(conf.get('linewidth', 2.4))
    color = conf.get('color', '#FF00A6')
    z = conf.get('zorder', 9)

    if species.lower() == 'tra':
        tail_y = telomere_y + offset
        head_y = telomere_y
        label_y = telomere_y + float(conf.get('label_offset', 0.03))
        label_va = 'bottom'
    else:
        tail_y = telomere_y - offset
        head_y = telomere_y
        label_y = telomere_y - float(conf.get('label_offset', 0.03))
        label_va = 'top'

    for i, x in enumerate(x_positions):
        arrow = FancyArrowPatch(
            (x, tail_y), (x, head_y),
            arrowstyle='-|>',
            mutation_scale=ms,
            linewidth=lw,
            color=color,
            shrinkA=0.0, shrinkB=0.0,
            zorder=z
        )
        ax.add_patch(arrow)

        if conf.get('show_labels', True) and labels is not None:
            ax.text(
                x, label_y, str(labels[i]),
                ha='center', va=label_va,
                fontsize=conf.get('label_fontsize', 10),
                color=conf.get('label_color', color),
                fontweight='bold',
                zorder=z + 1,
                path_effects=[path_effects.withStroke(linewidth=2.0, foreground='white')]
            )


def draw_breakpoint_markers_optimized(ax, breakpoints_plot, tra_y, chr_segments,
                                      tra_telomere_y=None, tba_telomere_y=None, arrow_conf=None):
    if breakpoints_plot.empty:
        return

    if arrow_conf is None:
        arrow_conf = {}

    breakpoints_plot = breakpoints_plot.reset_index(drop=True)
    labels = np.arange(
        arrow_conf.get('label_start', 1),
        arrow_conf.get('label_start', 1) + len(breakpoints_plot)
    )

    adjusted_x_positions = []
    for _, bp_row in breakpoints_plot.iterrows():
        chr_name = bp_row['chr']
        segments = chr_segments.get(chr_name, [])

        if len(segments) > 0 and segments[0].get('is_reversed', False):
            chr_length = max(seg['end'] for seg in segments)
            chr_x_start = min(seg['x_start'] for seg in segments)
            adjusted_x = chr_x_start + chr_length - bp_row['pos']
        else:
            adjusted_x = bp_row['x_pos']

        adjusted_x_positions.append(adjusted_x)

    x_positions = np.array(adjusted_x_positions)

    draw_inside_chr = arrow_conf.get('draw_inside_chr', False) is True
    if draw_inside_chr:
        chr_top = tra_y + CHR_HEIGHT / 2
        chr_bottom = tra_y - CHR_HEIGHT / 2
        segments_line = [[(x, chr_bottom), (x, chr_top)] for x in x_positions]
        lc = LineCollection(
            segments_line,
            colors=arrow_conf.get('inside_color', COLORS['breakpoint']),
            linewidths=arrow_conf.get('inside_width', 2.5),
            alpha=0.85,
            linestyles='solid',
            zorder=8
        )
        ax.add_collection(lc)

    apply_to = arrow_conf.get('apply_to', 'tra').lower()
    if apply_to in ('tra', 'both') and tra_telomere_y is not None:
        draw_breakpoint_arrows(ax, x_positions, tra_telomere_y, 'tra', arrow_conf, labels=labels)
    if apply_to in ('tba', 'both') and tba_telomere_y is not None:
        draw_breakpoint_arrows(ax, x_positions, tba_telomere_y, 'tba', arrow_conf, labels=labels)

# ==================== 共线性连接线 ====================

def draw_synteny_connections_generic(ax, synteny_plot, top_y, bottom_y, color_map,
                                     top_pos, bottom_pos,
                                     top_chr_col='chr_top', bottom_chr_col='chr_bottom',
                                     x_top_col='x_top', x_bottom_col='x_bottom',
                                     color_chr_col='color_chr'):
    """
    通用双物种共线性绘图函数：
    - top_pos: 上层物种位置表
    - bottom_pos: 下层物种位置表
    - synteny_plot 中需包含：
        x_top, x_bottom, chr_top, chr_bottom, min_window_size, color_chr
    """
    if len(synteny_plot) == 0:
        return

    top_reversed = dict(zip(top_pos['chr'], top_pos['is_reversed']))
    bottom_reversed = dict(zip(bottom_pos['chr'], bottom_pos['is_reversed']))

    top_lengths = dict(zip(top_pos['chr'], top_pos['length']))
    bottom_lengths = dict(zip(bottom_pos['chr'], bottom_pos['length']))

    top_starts = dict(zip(top_pos['chr'], top_pos['x_start']))
    bottom_starts = dict(zip(bottom_pos['chr'], bottom_pos['x_start']))

    def adjust_x_for_reverse(x_pos, chr_name, chr_reversed, chr_lengths, chr_starts):
        if chr_name in chr_reversed and chr_reversed[chr_name]:
            chr_length = chr_lengths.get(chr_name, 0)
            chr_start = chr_starts.get(chr_name, 0)
            relative_x = x_pos - chr_start
            return chr_start + chr_length - relative_x
        return x_pos

    synteny_plot = synteny_plot.copy()

    synteny_plot['x_top_adj'] = [
        adjust_x_for_reverse(x, chr_name, top_reversed, top_lengths, top_starts)
        for x, chr_name in zip(synteny_plot[x_top_col], synteny_plot[top_chr_col])
    ]
    synteny_plot['x_bottom_adj'] = [
        adjust_x_for_reverse(x, chr_name, bottom_reversed, bottom_lengths, bottom_starts)
        for x, chr_name in zip(synteny_plot[x_bottom_col], synteny_plot[bottom_chr_col])
    ]

    line_widths = pd.Series(
        np.maximum(0.5, np.minimum(2.5, np.log10(synteny_plot['min_window_size']) * 0.4)),
        index=synteny_plot.index
    )

    top_connection_y = top_y - CHR_HEIGHT / 2
    bottom_connection_y = bottom_y + CHR_HEIGHT / 2

    for chr_key, group in synteny_plot.groupby(color_chr_col):
        color = color_map.get(chr_key, '#808080')
        segments = [[(x_top, top_connection_y), (x_bottom, bottom_connection_y)]
                    for x_top, x_bottom in zip(group['x_top_adj'], group['x_bottom_adj'])]
        widths = line_widths.loc[group.index].values

        lc = LineCollection(
            segments,
            colors=color,
            alpha=COLORS['synteny_alpha'],
            linewidths=widths,
            zorder=2
        )
        ax.add_collection(lc)

# ==================== 通用颜色构建 ====================

def build_species_colors(chr_pos, custom_colors, default_palette):
    colors = {}
    palette_len = len(default_palette)
    for i, chr_name in enumerate(chr_pos['chr']):
        if custom_colors and chr_name in custom_colors:
            val = custom_colors[chr_name]
            colors[chr_name] = val[0] if isinstance(val, list) else val
        else:
            colors[chr_name] = default_palette[i % palette_len]
    return colors

# ==================== 主绘制流程 ====================

def create_enhanced_synteny_plot(tra_pos, tbi_pos, tba_pos,
                                 synteny_tbi_tra, synteny_tbi_tba,
                                 tra_telomere, tba_telomere, breakpoints):
    print("正在创建三物种增强版图形...")

    synteny_tbi_tra_filtered = filter_synteny_by_size(synteny_tbi_tra, MIN_SYNTENY_SIZE, label="Tra-Tbi")
    synteny_tbi_tba_filtered = filter_synteny_by_size(synteny_tbi_tba, MIN_SYNTENY_SIZE, label="Tbi-Tba")

    fig, ax = plt.subplots(figsize=FIGURE_SIZE, facecolor='white')
    ax.set_facecolor(COLORS['background'])

    # 物种纵向位置：Tra 上，Tbi 中，Tba 下
    tra_y = CHR_Y_SPACING * 2
    tbi_y = CHR_Y_SPACING * 1
    tba_y = 0

    # 颜色
    tbi_colors = build_species_colors(tbi_pos, CUSTOM_TBI_COLORS, COLORS['tbi_chr'])
    tba_colors = build_species_colors(tba_pos, CUSTOM_TBA_COLORS, COLORS['tba_chr'])

    # Tra 分段（根据 Tra-Tbi 共线性）
    chr_segments = determine_tra_segments_optimized(
        tra_pos, breakpoints, synteny_tbi_tra_filtered,
        tbi_colors, CUSTOM_TRA_COLORS
    )

    # ---------- Tra-Tbi 共线性 ----------
    # 文件中 chr1=Tbi, chr2=Tra
    synteny_plot_tra_tbi = synteny_tbi_tra_filtered.merge(
        tra_pos[['chr', 'x_start']], left_on='chr2', right_on='chr', how='inner'
    ).merge(
        tbi_pos[['chr', 'x_start']], left_on='chr1', right_on='chr', how='inner',
        suffixes=('_tra', '_tbi')
    )

    # 上层是 Tra，用 start2/end2；下层是 Tbi，用 start1/end1
    synteny_plot_tra_tbi['x_top'] = synteny_plot_tra_tbi['x_start_tra'] + (synteny_plot_tra_tbi['start2'] + synteny_plot_tra_tbi['end2']) / 2
    synteny_plot_tra_tbi['x_bottom'] = synteny_plot_tra_tbi['x_start_tbi'] + (synteny_plot_tra_tbi['start1'] + synteny_plot_tra_tbi['end1']) / 2
    synteny_plot_tra_tbi['chr_top'] = synteny_plot_tra_tbi['chr2']
    synteny_plot_tra_tbi['chr_bottom'] = synteny_plot_tra_tbi['chr1']
    synteny_plot_tra_tbi['color_chr'] = synteny_plot_tra_tbi['chr1']  # 用 Tbi 颜色

    draw_synteny_connections_generic(
        ax, synteny_plot_tra_tbi, tra_y, tbi_y, tbi_colors,
        tra_pos, tbi_pos
    )

    # ---------- Tbi-Tba 共线性 ----------
    # 文件中 chr1=Tbi, chr2=Tba
    synteny_plot_tbi_tba = synteny_tbi_tba_filtered.merge(
        tbi_pos[['chr', 'x_start']], left_on='chr1', right_on='chr', how='inner'
    ).merge(
        tba_pos[['chr', 'x_start']], left_on='chr2', right_on='chr', how='inner',
        suffixes=('_tbi', '_tba')
    )

    # 上层是 Tbi，用 start1/end1；下层是 Tba，用 start2/end2
    synteny_plot_tbi_tba['x_top'] = synteny_plot_tbi_tba['x_start_tbi'] + (synteny_plot_tbi_tba['start1'] + synteny_plot_tbi_tba['end1']) / 2
    synteny_plot_tbi_tba['x_bottom'] = synteny_plot_tbi_tba['x_start_tba'] + (synteny_plot_tbi_tba['start2'] + synteny_plot_tbi_tba['end2']) / 2
    synteny_plot_tbi_tba['chr_top'] = synteny_plot_tbi_tba['chr1']
    synteny_plot_tbi_tba['chr_bottom'] = synteny_plot_tbi_tba['chr2']
    synteny_plot_tbi_tba['color_chr'] = synteny_plot_tbi_tba['chr2']  # 用 Tba 颜色

    draw_synteny_connections_generic(
        ax, synteny_plot_tbi_tba, tbi_y, tba_y, tba_colors,
        tbi_pos, tba_pos
    )

    # ---------- 染色体 ----------
    # Tba（下）
    for _, row in tba_pos.iterrows():
        draw_chromosome(
            ax, row['x_start'], row['x_end'], tba_y, CHR_HEIGHT,
            tba_colors[row['chr']], row['chr'], 'bottom',
            is_reversed=row.get('is_reversed', False)
        )

    # Tbi（中）
    for _, row in tbi_pos.iterrows():
        draw_chromosome(
            ax, row['x_start'], row['x_end'], tbi_y, CHR_HEIGHT,
            tbi_colors[row['chr']], row['chr'], 'bottom',
            is_reversed=row.get('is_reversed', False)
        )

    # Tra（上，分段）
    for _, row in tra_pos.iterrows():
        chr_name = row['chr']
        segments = chr_segments[chr_name]
        is_reversed = row.get('is_reversed', False)
        draw_segmented_chromosome(ax, segments, tra_y, CHR_HEIGHT, chr_name, 'top', is_reversed)

    # ---------- 端粒（保持原逻辑：只画 Tra 和 Tba） ----------
    print("\n绘制端粒（log-quantile 缩放，低端大幅缩短）:")
    print(f"  模式: {TELOMERE_CONFIG.get('mode', 'log_quantile')}, gamma={TELOMERE_CONFIG.get('gamma', 3.0)}")
    print(f"  高度上下限: [{TELOMERE_CONFIG['min_height']}, {TELOMERE_CONFIG['max_height']}]")

    tra_telomere_y = tra_y + CHR_HEIGHT / 2 + 0.18
    tba_telomere_y = tba_y - CHR_HEIGHT / 2 - 0.18

    tra_telomere_count = draw_telomeres_as_lines_optimized(
        ax, tra_telomere, tra_pos, tra_telomere_y, COLORS['telomere_tra'], 'TRA'
    )
    tba_telomere_count = draw_telomeres_as_lines_optimized(
        ax, tba_telomere, tba_pos, tba_telomere_y, COLORS['telomere_tba'], 'TBA'
    )

    # ---------- 断点（保持原逻辑：Tra 层） ----------
    breakpoints_plot = breakpoints.merge(tra_pos[['chr', 'x_start']], on='chr', how='inner')
    if not breakpoints_plot.empty:
        breakpoints_plot['x_pos'] = breakpoints_plot['x_start'] + breakpoints_plot['pos']
        draw_breakpoint_markers_optimized(
            ax, breakpoints_plot, tra_y, chr_segments,
            tra_telomere_y=tra_telomere_y, tba_telomere_y=tba_telomere_y,
            arrow_conf=ARROW_MARK_CONFIG
        )

    # ---------- 物种标签 ----------
    ax.text(-0.03, tra_y, TEXT_CONFIG['species_label_tra'],
            fontsize=TEXT_CONFIG['species_label_fontsize'], fontweight='bold',
            ha='left', va='center', transform=ax.get_yaxis_transform(),
            path_effects=[path_effects.withStroke(linewidth=3, foreground='white')])

    ax.text(-0.03, tbi_y, TEXT_CONFIG['species_label_tbi'],
            fontsize=TEXT_CONFIG['species_label_fontsize'], fontweight='bold',
            ha='left', va='center', transform=ax.get_yaxis_transform(),
            path_effects=[path_effects.withStroke(linewidth=3, foreground='white')])

    ax.text(-0.03, tba_y, TEXT_CONFIG['species_label_tba'],
            fontsize=TEXT_CONFIG['species_label_fontsize'], fontweight='bold',
            ha='left', va='center', transform=ax.get_yaxis_transform(),
            path_effects=[path_effects.withStroke(linewidth=3, foreground='white')])

    # ---------- 轴与标题 ----------
    max_x = max(tra_pos['x_end'].max(), tbi_pos['x_end'].max(), tba_pos['x_end'].max())
    ax.set_xlim(-max_x * 0.02, max_x * 1.02)
    ax.set_ylim(-0.9, tra_y + 1.0)
    ax.set_title(
        f"{TEXT_CONFIG['title_main']}\n{TEXT_CONFIG['title_sub']}",
        fontsize=TEXT_CONFIG['title_fontsize'], fontweight='bold', pad=30
    )
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # ---------- 图例（简化版） ----------
    legend_elements = []

    if not breakpoints_plot.empty:
        legend_elements.append(
            plt.Line2D(
                [0], [0], marker='>', linestyle='None',
                color=ARROW_MARK_CONFIG.get('color', '#FF00A6'),
                markersize=max(8, ARROW_MARK_CONFIG.get('mutation_scale', 22.0) * 0.45),
                label=TEXT_CONFIG['legend_breakpoint']
            )
        )

    no_inversion_count = count_configured_no_inversion_segments(chr_segments)
    if no_inversion_count > 0:
        legend_elements.append(
            plt.Line2D(
                [0], [0],
                marker=NO_INVERSION_MARK_CONFIG.get('marker', 'D'),
                linestyle='None',
                color=NO_INVERSION_MARK_CONFIG.get('edgecolor', '#000000'),
                markerfacecolor=NO_INVERSION_MARK_CONFIG.get('facecolor', '#FFFFFF'),
                markeredgecolor=NO_INVERSION_MARK_CONFIG.get('edgecolor', '#000000'),
                markeredgewidth=NO_INVERSION_MARK_CONFIG.get('linewidth', 1.8),
                markersize=NO_INVERSION_MARK_CONFIG.get('legend_markersize', 11),
                label=f"{TEXT_CONFIG['legend_no_inversion']} (n={no_inversion_count})"
            )
        )

    legend_elements.extend([
        plt.Line2D([0], [0], color=COLORS['telomere_tra'], linewidth=TELOMERE_LINE_WIDTH,
                   label=f"{TEXT_CONFIG['legend_tra_telomere']} (n={tra_telomere_count})"),
        plt.Line2D([0], [0], color=COLORS['telomere_tba'], linewidth=TELOMERE_LINE_WIDTH,
                   label=f"{TEXT_CONFIG['legend_tba_telomere']} (n={tba_telomere_count})"),
    ])

    if len(synteny_plot_tra_tbi) > 0:
        legend_elements.append(
            plt.Line2D([0], [0], color='gray', alpha=COLORS['synteny_alpha'], linewidth=2.5,
                       label=TEXT_CONFIG['legend_synteny_tra_tbi'].format(len(synteny_plot_tra_tbi)))
        )

    if len(synteny_plot_tbi_tba) > 0:
        legend_elements.append(
            plt.Line2D([0], [0], color='dimgray', alpha=COLORS['synteny_alpha'], linewidth=2.5,
                       label=TEXT_CONFIG['legend_synteny_tbi_tba'].format(len(synteny_plot_tbi_tba)))
        )

    total_reversed = tra_pos['is_reversed'].sum() + tbi_pos['is_reversed'].sum() + tba_pos['is_reversed'].sum()
    if total_reversed > 0 and REVERSE_LABEL_MARK.get('show', True):
        legend_elements.append(
            plt.Line2D([0], [0], marker='*', linestyle='None',
                       color=REVERSE_LABEL_MARK.get('color', 'red'),
                       markersize=REVERSE_LABEL_MARK.get('legend_markersize', 14),
                       label=f"{TEXT_CONFIG['legend_reversed']} (n={total_reversed})")
        )

    legend = ax.legend(handles=legend_elements, loc='upper right',
                       bbox_to_anchor=(1.18, 1), frameon=True,
                       fancybox=True, shadow=True, fontsize=TEXT_CONFIG['legend_fontsize'])
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(0.95)
    legend.get_frame().set_edgecolor(CHR_BORDER_CONFIG.get('color', COLORS['edge_color']))

    plt.tight_layout()

    return fig, len(synteny_plot_tra_tbi), len(synteny_plot_tbi_tba), tra_telomere_count, tba_telomere_count


def save_enhanced_results(fig, tra_pos, tbi_pos, tba_pos,
                          synteny_tbi_tra_count, synteny_tbi_tba_count,
                          breakpoints_count,
                          filtered_tbi_tra_count, filtered_tbi_tba_count,
                          tra_telomere_count, tba_telomere_count):
    print("正在保存结果...")

    filename_base = TEXT_CONFIG['output_filename_base']
    output_files = {
        'pdf': f'{OUTPUT_DIR}/{filename_base}.pdf',
        'png': f'{OUTPUT_DIR}/{filename_base}.png',
        'svg': f'{OUTPUT_DIR}/{filename_base}.svg'
    }

    fig.savefig(output_files['pdf'], dpi=DPI, bbox_inches='tight',
                facecolor='white', edgecolor='none', format='pdf')
    fig.savefig(output_files['png'], dpi=DPI, bbox_inches='tight',
                facecolor='white', edgecolor='none', format='png')
    fig.savefig(output_files['svg'], bbox_inches='tight',
                facecolor='white', edgecolor='none', format='svg')

    tra_reversed = list(tra_pos[tra_pos['is_reversed']]['chr'])
    tbi_reversed = list(tbi_pos[tbi_pos['is_reversed']]['chr'])
    tba_reversed = list(tba_pos[tba_pos['is_reversed']]['chr'])

    total_synteny_raw = synteny_tbi_tra_count + synteny_tbi_tba_count
    total_synteny_filtered = filtered_tbi_tra_count + filtered_tbi_tba_count

    stats = f"""{TEXT_CONFIG['stats_title']}
======================================
缩放与标注:
- 端粒高度: 模式={TELOMERE_CONFIG.get('mode')}, gamma={TELOMERE_CONFIG.get('gamma')}, 高度范围=[{TELOMERE_CONFIG['min_height']}, {TELOMERE_CONFIG['max_height']}]
- 断点标注: 箭头 color={ARROW_MARK_CONFIG.get('color')}, mutation_scale={ARROW_MARK_CONFIG.get('mutation_scale')}, linewidth={ARROW_MARK_CONFIG.get('linewidth')}, offset={ARROW_MARK_CONFIG.get('offset')}
- 断点编号: show={ARROW_MARK_CONFIG.get('show_labels')}, fontsize={ARROW_MARK_CONFIG.get('label_fontsize')}, color={ARROW_MARK_CONFIG.get('label_color')}
- Tra无inversion标记: show={NO_INVERSION_MARK_CONFIG.get('show')}, configured_chr={len(NO_INVERSION_MARK_CONFIG.get('segments', {}))}, marker={NO_INVERSION_MARK_CONFIG.get('text')}
- Tra/Tbi/Tba 组间间隔: {TRA_CHR_SPACING}/{TBI_CHR_SPACING}/{TBA_CHR_SPACING}; Tbi组内间隔: {TBI_INTRA_GROUP_SPACING}; Tba组内间隔: {TBA_INTRA_GROUP_SPACING}
- Tbi分组数: {len(TBI_GROUPS_BY_TRA)}；Tba分组数: {len(TBA_GROUPS_BY_TRA)}
- 绘制数量限制: TRA={TRA_CHR_LIMIT}, TBI={TBI_CHR_LIMIT}, TBA={TBA_CHR_LIMIT}
- 标签字号: 默认={LABEL_CONFIG['fontsize']}, per-chr={len(LABEL_SIZE_PER_CHR)}条；星标字号：专用={REVERSE_LABEL_MARK.get('fontsize')}, scale={REVERSE_LABEL_MARK.get('scale')}, per-chr={len(STAR_SIZE_PER_CHR)}条
- 显示名映射: {len(CHR_DISPLAY_NAME)} 条
- chr_order 列配置: TRA={CHR_ORDER_COLUMNS['tra']}, TBI={CHR_ORDER_COLUMNS['tbi']}, TBA={CHR_ORDER_COLUMNS['tba']}

染色体反转:
- TRA反转数量: {len(tra_reversed)}；列表: {', '.join(tra_reversed) if tra_reversed else '无'}
- TBI反转数量: {len(tbi_reversed)}；列表: {', '.join(tbi_reversed) if tbi_reversed else '无'}
- TBA反转数量: {len(tba_reversed)}；列表: {', '.join(tba_reversed) if tba_reversed else '无'}

统计:
- TRA染色体: {len(tra_pos)}；TBI染色体: {len(tbi_pos)}；TBA染色体: {len(tba_pos)}
- Tra-Tbi 共线性原始区域: {synteny_tbi_tra_count:,}；过滤后: {filtered_tbi_tra_count:,}
- Tbi-Tba 共线性原始区域: {synteny_tbi_tba_count:,}；过滤后: {filtered_tbi_tba_count:,}
- 共线性总原始区域: {total_synteny_raw:,}；过滤后总计: {total_synteny_filtered:,}
- 端粒：TRA={tra_telomere_count}，TBA={tba_telomere_count}
- 断点条目: {breakpoints_count}

输出:
- {filename_base}.pdf / .png / .svg
- {filename_base}_results.txt
"""
    with open(f'{OUTPUT_DIR}/{filename_base}_results.txt', 'w', encoding='utf-8') as f:
        f.write(stats)

    print("结果保存完成！")

# ==================== 主程序 ====================

if __name__ == "__main__":
    import time

    try:
        start_time = time.time()
        print("开始三物种优化版生物信息学共线性图分析...")
        print("=" * 80)

        (tra_chr, tbi_chr, tba_chr,
         synteny_tbi_tra, synteny_tbi_tba,
         tra_telomere, tba_telomere,
         chr_order, breakpoints) = load_all_data()

        tra_pos, tbi_pos, tba_pos = prepare_chromosome_data(
            tra_chr, tbi_chr, tba_chr, chr_order,
            synteny_tbi_tra, synteny_tbi_tba
        )

        fig, filtered_tbi_tra_count, filtered_tbi_tba_count, tra_telomere_count, tba_telomere_count = create_enhanced_synteny_plot(
            tra_pos, tbi_pos, tba_pos,
            synteny_tbi_tra, synteny_tbi_tba,
            tra_telomere, tba_telomere, breakpoints
        )

        save_enhanced_results(
            fig, tra_pos, tbi_pos, tba_pos,
            len(synteny_tbi_tra), len(synteny_tbi_tba),
            len(breakpoints),
            filtered_tbi_tra_count, filtered_tbi_tba_count,
            tra_telomere_count, tba_telomere_count
        )

        elapsed_time = time.time() - start_time
        print("=" * 80)
        print(f"✓ 完成，总耗时: {elapsed_time:.2f} 秒")
        print(f"✓ TRA/TBI/TBA 染色体: {len(tra_pos)}/{len(tbi_pos)}/{len(tba_pos)}")
        print(f"✓ Tra-Tbi 高质量共线性区域: {filtered_tbi_tra_count:,}")
        print(f"✓ Tbi-Tba 高质量共线性区域: {filtered_tbi_tba_count:,}")
        print(f"✓ 断点条目: {len(breakpoints)}")
        print(f"✓ 端粒（TRA/TBA）: {tra_telomere_count}/{tba_telomere_count}")

    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()

