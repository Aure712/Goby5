#!/usr/bin/env python3
"""
优化版生物信息学共线性图绘制脚本

新增/变更：
- 断点标注：改为在端粒层附近画箭头指向断点（默认颜色 #FF00A6），支持指定箭头大小与粗细。
- 箭头位置：位于端粒层与染色体之间的空间，避免遮挡（zorder 介于染色体与端粒之间）。
- 文本尺寸：支持为每条染色体设置标签字号，支持为反转星标“*”设置字号或按比例缩放。
- 保留端粒高度的 log-quantile + gamma 缩放，低计数端粒显著缩短。
- 新增：端粒箭头编号显示断点序号（加粗）；支持按 Tra 分组设定 Tba 间隔；染色体边框颜色与粗细可配置；
- 新增：可分别限制 Tra/Tba 绘制的染色体数量。
- 新增：支持染色体显示名映射（例如 TraScf_1 -> Tra1），不影响内部计算。
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import seaborn as sns
import matplotlib.patheffects as path_effects
from matplotlib.collections import LineCollection
import matplotlib.transforms as transforms  # 用于屏幕坐标偏移（叠加“*”）

# ==================== 路径与基础配置 ====================
DATA_DIR = "."
OUTPUT_DIR = "."

FILES = {
    'tra_fai': 'Tra_sort.fa.fai',
    'tba_fai': 'Tba_sort.fa.fai',
    'synteny': 'pairwise.simp.txt',
    'tra_telomere': 'tidk_default_tra_onlychr.bed',
    'tba_telomere': 'tidk_default_tba_onlychr.bed',
    'chr_order': 'tba.chr.order.txt',
    'breakpoints': 'lines.txt'
}

FIGURE_SIZE = (28, 12)
DPI = 1200

# ==================== 染色体绘制数量限制（None 表示不限制） ====================
TRA_CHR_LIMIT = 11
TBA_CHR_LIMIT = 22

# ==================== 染色体边框配置 ====================
CHR_BORDER_CONFIG = {
    'color': '#000000',
    'linewidth': 4.0
}

# ==================== 染色体间隔配置 ====================
# Tra 染色体之间的间隔（主要控制）
TRA_CHR_SPACING = 20000000

# Tba 分组内间隔（同组内）
TBA_INTRA_GROUP_SPACING = 5000000  # 可自行调小或调大

# Tba 分组映射：按 Tra 染色体顺序分组对齐
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
TELOMERE_DOT_SIZE = 40  # 已不使用，保留以防需要

# ==================== 染色体反转配置 ====================
REVERSED_CHROMOSOMES = {
    'tra': [],
    'tba': ['TbaScf_16', 'TbaScf_2', 'TbaScf_4', 'TbaScf_3',
            'TbaScf_11', 'TbaScf_13', 'TbaScf_9', 'TbaScf_6', 'TbaScf_18']
}

# ==================== 染色体显示名映射（可选） ====================
# key: 原始名称；value: 图中显示名称
# 不在该字典中的名称将保持原样显示
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

    # TBA（按需启用/补充）
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
    """返回染色体显示名；未配置映射时返回原名。"""
    return CHR_DISPLAY_NAME.get(chr_name, chr_name)

# 反转标记配置：在名称后追加“*”
REVERSE_LABEL_MARK = {
    'show': True,         # 显示星标
    'char': '*',          # 星标字符
    'color': 'red',       # 星标颜色
    'fontsize': 20,     # 指定星标字号；None 表示随标签字号；也可配合 scale 使用
    'scale': 1.0,         # 星标字号 = (标签字号) * scale（当 fontsize 为 None 时生效）
    'gap_pts': 5.0,       # 名称与星标之间的间隙（点）
    'legend_markersize': 14  # 图例中星标大小
}

# ==================== 染色体标签位置与字号配置 ====================
LABEL_CONFIG = {
    'bottom_offset': 0.35,
    'top_offset': 0.40,
    'rotation': 45,
    'fontsize': 18,       # 默认标签字号
    'fontweight': 'bold'
}

# 为每条染色体单独设置标签字号；不设置的采用 LABEL_CONFIG['fontsize']
LABEL_SIZE_PER_CHR = {}
# 为每条染色体单独设置星标字号；不设置的采用 REVERSE_LABEL_MARK['fontsize'] 或标签字号*scale
STAR_SIZE_PER_CHR = {}

# ==================== 断点箭头（端粒层附近的标注）配置 ====================
ARROW_MARK_CONFIG = {
    'color': '#FF00A6',   # 箭头颜色（鲜艳玫红，未在图中其他地方使用）
    'apply_to': 'tra',    # 'tra' | 'tba' | 'both'
    'mutation_scale': 40.0,  # 箭头头部尺寸（点为单位的缩放）
    'linewidth': 6.0,        # 箭杆粗细
    'offset': 0.15,          # 箭尾到箭头的垂直间距（数据坐标，越大箭头越长）
    'zorder': 11,             # 图层顺序：染色体(z=5) < 箭头(z=9) < 端粒(z=10)
    'draw_inside_chr': False, # 是否保留染色体内部的断点竖线
    'inside_color': '#000000',  # Tra 断点线条颜色（默认黑色）
    'inside_width': 2.5,        # Tra 断点线条粗细

    # 断点编号显示
    'show_labels': True,
    'label_color': '#000000',
    'label_fontsize': 24,
    'label_offset': 0.15,   # 标号与端粒层的距离（数据坐标）
    'label_start': 1
}

# ==================== 端粒垂直线配置参数（对比度增强） ====================
TELOMERE_CONFIG = {
    'mode': 'log_quantile',    # 'log_quantile' 或 'legacy'
    'gamma': 2.5,              # 幂指数，>1 拉大差异
    'min_height': 0.0008,      # 极低下限
    'max_height': 0.30,        # 上限
    'base_height': 0.02,       # legacy 使用
    'log_scale_factor': 2.0,   # legacy 使用
    'show_dots': False
}

COLORS = {
    'tra_chr_base': '#CCCCCC',
    'tba_chr': sns.color_palette("tab20", 20) + sns.color_palette("Set3", 2),
    'synteny_alpha': 0.3,
    'telomere_tra': '#FF6B35',
    'telomere_tba': '#4A90E2',
    'breakpoint': '#9B59B6',    # 兼容保留
    'breakpoint_line': '#95A5A6',
    'background': '#FAFAFA',
    'edge_color': '#2C3E50'
}

# ==================== 图表文字配置 ====================
TEXT_CONFIG = {
    'species_label_tra': 'Tra',
    'species_label_tba': 'Tba',
    'species_label_fontsize': 28,
    'title_main': 'Chromosome Synteny Analysis',
    'title_sub': 'with Telomere Distribution and Breakpoint-based Segmentation',
    'title_fontsize': 20,
    'legend_breakpoint': 'Recombination Breakpoints',
    'legend_tra_telomere': 'Tra Telomeres',
    'legend_tba_telomere': 'Tba Telomeres',
    'legend_tba_chr': 'Tba Chr',
    'legend_more_chrs': '... and {} more Tba chrs',
    'legend_synteny': 'Synteny (≥{:,} bp, n={})',
    'legend_reversed': 'Reversed Chromosomes (*)',
    'legend_fontsize': 12,
    'output_filename_base': 'enhanced_synteny_plot_with_reverse',
    'stats_title': '优化版生物信息学共线性分析结果',
}

# ==================== 自定义染色体区域颜色配置（可选） ====================
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
    tba_chr = pd.read_csv(f"{DATA_DIR}/{FILES['tba_fai']}", sep='\t', header=None,
                          names=['chr', 'length', 'offset', 'line_bases', 'line_width'])
    synteny = pd.read_csv(f"{DATA_DIR}/{FILES['synteny']}", sep='\t', header=None,
                          names=['chr1', 'start1', 'end1', 'strand', 'chr2', 'start2', 'end2', 'score'])
    tra_telomere = pd.read_csv(f"{DATA_DIR}/{FILES['tra_telomere']}", sep='\t', header=None,
                               names=['chr', 'start', 'end', 'count'])
    tba_telomere = pd.read_csv(f"{DATA_DIR}/{FILES['tba_telomere']}", sep='\t', header=None,
                               names=['chr', 'start', 'end', 'count'])
    chr_order = pd.read_csv(f"{DATA_DIR}/{FILES['chr_order']}", sep='\t')
    breakpoints = pd.read_csv(f"{DATA_DIR}/{FILES['breakpoints']}", sep='\t', header=None,
                              names=['chr', 'pos', 'start_range', 'end_range'])
    print("数据加载完成！")
    return tra_chr, tba_chr, synteny, tra_telomere, tba_telomere, chr_order, breakpoints


def filter_synteny_by_size(synteny, min_size=MIN_SYNTENY_SIZE):
    print(f"正在过滤共线性数据，最小窗口大小: {min_size:,} bp")
    window_size1 = synteny['end1'] - synteny['start1']
    window_size2 = synteny['end2'] - synteny['start2']
    min_window_size = np.minimum(window_size1, window_size2)
    original_count = len(synteny)
    filtered_synteny = synteny[min_window_size >= min_size].copy()
    filtered_synteny['min_window_size'] = min_window_size[min_window_size >= min_size]
    filtered_count = len(filtered_synteny)
    print(f"原始共线性区域: {original_count:,}")
    print(f"过滤后共线性区域: {filtered_count:,}")
    print(f"过滤掉: {original_count - filtered_count:,} ({(original_count - filtered_count) / original_count * 100:.1f}%)")
    return filtered_synteny


def _build_tba_order(tra_ordered, chr_order, synteny):
    if TBA_GROUPS_BY_TRA:
        ordered = []
        for tra_chr in tra_ordered:
            ordered.extend(TBA_GROUPS_BY_TRA.get(tra_chr, []))

        fallback = chr_order[~chr_order['chrname'].isna() & (chr_order['chrname'] != '')]['chrname'].tolist()
        extra = fallback + list(synteny['chr2'].unique())
        for name in extra:
            if name not in ordered:
                ordered.append(name)
        return ordered

    return chr_order[~chr_order['chrname'].isna() & (chr_order['chrname'] != '')]['chrname'].tolist()


def prepare_chromosome_data(tra_chr, tba_chr, chr_order, synteny):
    print("准备染色体数据...")
    tra_ordered = chr_order[~chr_order['scf'].isna() & (chr_order['scf'] != '')]['scf'].tolist()
    if TRA_CHR_LIMIT is not None:
        tra_ordered = tra_ordered[:TRA_CHR_LIMIT]

    tba_ordered = _build_tba_order(tra_ordered, chr_order, synteny)

    tba_filtered = tba_chr[tba_chr['chr'].isin(tba_ordered + list(synteny['chr2'].unique()))].copy()
    tra_filtered = tra_chr[tra_chr['chr'].isin(tra_ordered + list(synteny['chr1'].unique()))].copy()

    tba_order_dict = {chr_name: i for i, chr_name in enumerate(tba_ordered)}
    tra_order_dict = {chr_name: i for i, chr_name in enumerate(tra_ordered)}

    tba_filtered['order'] = tba_filtered['chr'].map(tba_order_dict)
    tra_filtered['order'] = tra_filtered['chr'].map(tra_order_dict)

    tba_filtered = tba_filtered.dropna(subset=['order']).sort_values('order')
    tra_filtered = tra_filtered.dropna(subset=['order']).sort_values('order')

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

    if TBA_GROUPS_BY_TRA:
        tba_pos = calc_positions_grouped(
            tba_filtered, tra_ordered, TBA_GROUPS_BY_TRA,
            TRA_CHR_SPACING, TBA_INTRA_GROUP_SPACING
        )
        tba_pos['order'] = tba_pos['chr'].map(tba_order_dict)
        tba_pos = tba_pos.sort_values('order')
    else:
        tba_pos = calc_positions_vectorized(tba_filtered, TRA_CHR_SPACING)

    tra_pos = calc_positions_vectorized(tra_filtered, TRA_CHR_SPACING)

    if TBA_CHR_LIMIT is not None:
        tba_pos = tba_pos.head(TBA_CHR_LIMIT)
    if TRA_CHR_LIMIT is not None:
        tra_pos = tra_pos.head(TRA_CHR_LIMIT)

    tba_pos['is_reversed'] = tba_pos['chr'].isin(REVERSED_CHROMOSOMES.get('tba', []))
    tra_pos['is_reversed'] = tra_pos['chr'].isin(REVERSED_CHROMOSOMES.get('tra', []))

    if tba_pos['is_reversed'].sum() > 0:
        print(f"TBA反转: {list(tba_pos[tba_pos['is_reversed']]['chr'])}")
    if tra_pos['is_reversed'].sum() > 0:
        print(f"TRA反转: {list(tra_pos[tra_pos['is_reversed']]['chr'])}")

    print(f"处理了 {len(tra_pos)} 个TRA 和 {len(tba_pos)} 个TBA 染色体；反转：TRA {tra_pos['is_reversed'].sum()}，TBA {tba_pos['is_reversed'].sum()}")
    return tra_pos, tba_pos


def determine_chromosome_segments_optimized(tra_pos, breakpoints, synteny, tba_colors, custom_colors=None):
    print("正在分析染色体段和着色...")
    chromosome_segments = {}
    synteny_by_chr = synteny.groupby('chr1')

    for idx, chr_row in tra_pos.iterrows():
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
            elif chr_name in synteny_by_chr.groups:
                chr_synteny = synteny_by_chr.get_group(chr_name)
                tba_chr_counts = chr_synteny['chr2'].value_counts()
                main_tba_chr = tba_chr_counts.index[0]
                color = tba_colors.get(main_tba_chr, list(tba_colors.values())[0])
            else:
                main_tba_chr = None
                color = list(tba_colors.values())[0] if tba_colors else COLORS['tba_chr'][0]

            segments.append({
                'start': 0,
                'end': chr_length,
                'x_start': chr_x_start,
                'x_end': chr_x_start + chr_length,
                'color': color,
                'tba_chr': main_tba_chr if not has_custom_colors else None,
                'is_reversed': is_reversed
            })
        else:
            positions = np.concatenate([[0], chr_breakpoints['pos'].values, [chr_length]])
            chr_synteny = synteny_by_chr.get_group(chr_name) if chr_name in synteny_by_chr.groups else pd.DataFrame()
            segment_indices = range(len(positions) - 1)
            if is_reversed:
                segment_indices = reversed(segment_indices)

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
                    main_tba_chr = None
                else:
                    if not chr_synteny.empty:
                        mask = (chr_synteny['start1'] <= seg_middle) & (chr_synteny['end1'] >= seg_middle)
                        seg_synteny = chr_synteny[mask]
                        if len(seg_synteny) > 0:
                            tba_chr_counts = seg_synteny['chr2'].value_counts()
                            main_tba_chr = tba_chr_counts.index[0]
                            color = tba_colors.get(main_tba_chr, COLORS['tba_chr'][i % len(COLORS['tba_chr'])])
                        else:
                            main_tba_chr = None
                            color = COLORS['tba_chr'][i % len(COLORS['tba_chr'])]
                    else:
                        main_tba_chr = None
                        color = COLORS['tba_chr'][i % len(COLORS['tba_chr'])]

                segments.append({
                    'start': seg_start,
                    'end': seg_end,
                    'x_start': actual_x_start,
                    'x_end': actual_x_end,
                    'color': color,
                    'tba_chr': main_tba_chr,
                    'is_reversed': is_reversed
                })

        chromosome_segments[chr_name] = segments

    return chromosome_segments

# ==================== 文本绘制（支持每染色体字号与星标字号） ====================

def _draw_label_with_optional_star(ax, x_center, y, chr_name, position, is_reversed):
    rotation = LABEL_CONFIG['rotation']
    fontsize_label = LABEL_SIZE_PER_CHR.get(chr_name, LABEL_CONFIG['fontsize'])
    fontweight = LABEL_CONFIG['fontweight']

    if position == 'bottom':
        va = 'top'
        ha = 'right'
    else:
        va = 'bottom'
        ha = 'left'

    main_text = ax.text(
        x_center, y, chr_name,
        rotation=rotation, ha=ha, va=va,
        fontsize=fontsize_label, fontweight=fontweight,
        path_effects=[path_effects.withStroke(linewidth=2.5, foreground='white')],
        zorder=6
    )

    if is_reversed and REVERSE_LABEL_MARK.get('show', True):
        star_char = str(REVERSE_LABEL_MARK.get('char', '*'))
        star_color = REVERSE_LABEL_MARK.get('color', 'red')

        star_fontsize = STAR_SIZE_PER_CHR.get(chr_name, None)
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
            dx_pts += width_px * 72.0 / fig.dpi  # 像素->点

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

    total_start = segments[0]['x_start']
    total_end = segments[-1]['x_end']
    x_center = (total_start + total_end) / 2

    # 只有在关闭内部边框时，补一个整条染色体的外框
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

    label_y = y_center - height / 2 - LABEL_CONFIG['bottom_offset'] if position == 'bottom' \
        else y_center + height / 2 + LABEL_CONFIG['top_offset']

    _draw_label_with_optional_star(ax, x_center, label_y, get_display_name(chr_name), position, is_reversed)



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

    _draw_label_with_optional_star(ax, x_center, label_y, get_display_name(label), position, is_reversed)


def _compute_telomere_heights(counts):
    mode = TELOMERE_CONFIG.get('mode', 'log_quantile')
    min_h = TELOMERE_CONFIG.get('min_height', 0.0008)
    max_h = TELOMERE_CONFIG.get('max_height', 0.35)

    counts = np.asarray(counts, dtype=float)
    counts = np.maximum(1.0, counts)  # 避免log(0)

    if mode == 'legacy':
        base = TELOMERE_CONFIG.get('base_height', 0.02)
        factor = TELOMERE_CONFIG.get('log_scale_factor', 2.0)
        heights = base * (1.0 + np.log10(counts) * factor)
        heights = np.clip(heights, min_h, max_h)
        return heights

    # log_quantile 模式
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

    lc = LineCollection(segments, colors=color, linewidths=TELOMERE_LINE_WIDTH * 1.5,
                        alpha=TELOMERE_LINE_ALPHA, zorder=10, capstyle='round')
    ax.add_collection(lc)

    if TELOMERE_CONFIG.get('show_dots', False):
        ax.scatter(tel_centers, [y_position] * len(tel_centers),
                   s=dot_sizes, c=color, alpha=TELOMERE_LINE_ALPHA * 0.8,
                   edgecolors='white', linewidths=1.5, zorder=11, marker='o')

    return len(tel_plot)


def draw_breakpoint_arrows(ax, x_positions, telomere_y, species, conf, labels=None):
    """
    在端粒层附近画箭头指向断点：
    - TRA: 端粒层在染色体上方，箭头自更高处指向 telomere_y（向下）
    - TBA: 端粒层在染色体下方，箭头自更低处指向 telomere_y（向上）
    """
    if x_positions.size == 0:
        return
    offset = float(conf.get('offset', 0.08))
    ms = float(conf.get('mutation_scale', 22.0))
    lw = float(conf.get('linewidth', 2.4))
    color = conf.get('color', '#FF00A6')
    z = conf.get('zorder', 9)

    if species.lower() == 'tra':
        tail_y = telomere_y + offset  # 端粒之上
        head_y = telomere_y
        label_y = telomere_y + float(conf.get('label_offset', 0.03))
        label_va = 'bottom'
    else:
        tail_y = telomere_y - offset  # 端粒之下
        head_y = telomere_y
        label_y = telomere_y - float(conf.get('label_offset', 0.03))
        label_va = 'top'

    for i, x in enumerate(x_positions):
        arrow = FancyArrowPatch(
            (x, tail_y), (x, head_y),
            arrowstyle='-|>',  # 简洁箭头
            mutation_scale=ms,  # 控制箭头头部尺寸
            linewidth=lw,
            color=color,  # 同时作为边和面颜色
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
    """
    绘制断点标记：
    - 染色体内部：竖线（可开关）
    - 端粒层附近：箭头指向断点位置（可设颜色/大小/编号）
    """
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

    # 染色体内部竖线（仅当显式设为 True 才绘制）
    draw_inside_chr = arrow_conf.get('draw_inside_chr', False) is True
    if draw_inside_chr:
        chr_top = tra_y + CHR_HEIGHT / 2
        chr_bottom = tra_y - CHR_HEIGHT / 2
        segments_line = [[(x, chr_bottom), (x, chr_top)] for x in x_positions]
        lc = LineCollection(segments_line,
                            colors=arrow_conf.get('inside_color', COLORS['breakpoint']),
                            linewidths=arrow_conf.get('inside_width', 2.5),
                            alpha=0.85, linestyles='solid', zorder=8)
        ax.add_collection(lc)

    # 端粒层箭头 + 编号
    apply_to = arrow_conf.get('apply_to', 'tra').lower()
    if apply_to in ('tra', 'both') and tra_telomere_y is not None:
        draw_breakpoint_arrows(ax, x_positions, tra_telomere_y, 'tra', arrow_conf, labels=labels)
    if apply_to in ('tba', 'both') and tba_telomere_y is not None:
        draw_breakpoint_arrows(ax, x_positions, tba_telomere_y, 'tba', arrow_conf, labels=labels)

# ==================== 共线性连接线 ====================

def draw_synteny_connections_optimized(ax, synteny_plot, tra_y, tba_y, tba_colors, tra_pos, tba_pos):
    if len(synteny_plot) == 0:
        return

    tra_reversed = dict(zip(tra_pos['chr'], tra_pos['is_reversed']))
    tba_reversed = dict(zip(tba_pos['chr'], tba_pos['is_reversed']))
    tra_lengths = dict(zip(tra_pos['chr'], tra_pos['length']))
    tba_lengths = dict(zip(tba_pos['chr'], tba_pos['length']))
    tra_starts = dict(zip(tra_pos['chr'], tra_pos['x_start']))
    tba_starts = dict(zip(tba_pos['chr'], tba_pos['x_start']))

    def adjust_x_for_reverse(x_pos, chr_name, chr_reversed, chr_lengths, chr_starts):
        if chr_name in chr_reversed and chr_reversed[chr_name]:
            chr_length = chr_lengths.get(chr_name, 0)
            chr_start = chr_starts.get(chr_name, 0)
            relative_x = x_pos - chr_start
            return chr_start + chr_length - relative_x
        return x_pos

    x_tra_adjusted = []
    x_tba_adjusted = []

    for _, row in synteny_plot.iterrows():
        x_tra = adjust_x_for_reverse(row['x_tra'], row['chr1'], tra_reversed, tra_lengths, tra_starts)
        x_tba = adjust_x_for_reverse(row['x_tba'], row['chr2'], tba_reversed, tba_lengths, tba_starts)
        x_tra_adjusted.append(x_tra)
        x_tba_adjusted.append(x_tba)

    synteny_plot = synteny_plot.copy()
    synteny_plot['x_tra_adj'] = x_tra_adjusted
    synteny_plot['x_tba_adj'] = x_tba_adjusted

    line_widths = np.maximum(0.5, np.minimum(2.5, np.log10(synteny_plot['min_window_size']) * 0.4))
    tra_connection_y = tra_y - CHR_HEIGHT / 2
    tba_connection_y = tba_y + CHR_HEIGHT / 2

    for chr2, group in synteny_plot.groupby('chr2'):
        color = tba_colors.get(chr2, '#808080')
        segments = [[(x_tra, tra_connection_y), (x_tba, tba_connection_y)]
                    for x_tra, x_tba in zip(group['x_tra_adj'], group['x_tba_adj'])]
        widths = line_widths[group.index].values
        lc = LineCollection(segments, colors=color, alpha=COLORS['synteny_alpha'],
                            linewidths=widths, zorder=2)
        ax.add_collection(lc)

# ==================== 主绘制流程 ====================

def create_enhanced_synteny_plot(tra_pos, tba_pos, synteny, tra_telomere, tba_telomere, breakpoints):
    print("正在创建增强版图形...")

    synteny_filtered = filter_synteny_by_size(synteny, MIN_SYNTENY_SIZE)

    fig, ax = plt.subplots(figsize=FIGURE_SIZE, facecolor='white')
    ax.set_facecolor(COLORS['background'])

    tra_y = CHR_Y_SPACING
    tba_y = 0

    # 为TBA染色体分配颜色
    tba_colors = {}
    for i, chr_name in enumerate(tba_pos['chr']):
        if CUSTOM_TBA_COLORS and chr_name in CUSTOM_TBA_COLORS:
            tba_colors[chr_name] = CUSTOM_TBA_COLORS[chr_name][0] if isinstance(CUSTOM_TBA_COLORS[chr_name], list) else CUSTOM_TBA_COLORS[chr_name]
        else:
            tba_colors[chr_name] = COLORS['tba_chr'][i % len(COLORS['tba_chr'])]

    # TRA段分析
    chr_segments = determine_chromosome_segments_optimized(tra_pos, breakpoints, synteny_filtered,
                                                           tba_colors, CUSTOM_TRA_COLORS)

    # 共线性连接（底层）
    synteny_plot = synteny_filtered.merge(
        tra_pos[['chr', 'x_start']], left_on='chr1', right_on='chr', how='inner'
    ).merge(
        tba_pos[['chr', 'x_start']], left_on='chr2', right_on='chr', how='inner',
        suffixes=('_tra', '_tba')
    )
    synteny_plot['x_tra'] = synteny_plot['x_start_tra'] + (synteny_plot['start1'] + synteny_plot['end1']) / 2
    synteny_plot['x_tba'] = synteny_plot['x_start_tba'] + (synteny_plot['start2'] + synteny_plot['end2']) / 2

    draw_synteny_connections_optimized(ax, synteny_plot, tra_y, tba_y, tba_colors, tra_pos, tba_pos)

    # 染色体（上层）
    for _, row in tba_pos.iterrows():
        draw_chromosome(ax, row['x_start'], row['x_end'], tba_y, CHR_HEIGHT,
                        tba_colors[row['chr']], row['chr'], 'bottom',
                        is_reversed=row.get('is_reversed', False))

    for _, row in tra_pos.iterrows():
        chr_name = row['chr']
        segments = chr_segments[chr_name]
        is_reversed = row.get('is_reversed', False)
        draw_segmented_chromosome(ax, segments, tra_y, CHR_HEIGHT, chr_name, 'top', is_reversed)

    # 端粒（增强高度差异、不显示圆点）
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

    # 断点（染色体内部可选 + 端粒层箭头+编号）
    breakpoints_plot = breakpoints.merge(tra_pos[['chr', 'x_start']], on='chr', how='inner')
    if not breakpoints_plot.empty:
        breakpoints_plot['x_pos'] = breakpoints_plot['x_start'] + breakpoints_plot['pos']
        draw_breakpoint_markers_optimized(
            ax, breakpoints_plot, tra_y, chr_segments,
            tra_telomere_y=tra_telomere_y, tba_telomere_y=tba_telomere_y,
            arrow_conf=ARROW_MARK_CONFIG
        )

    # 物种标签
    ax.text(-0.03, tra_y, TEXT_CONFIG['species_label_tra'],
            fontsize=TEXT_CONFIG['species_label_fontsize'], fontweight='bold',
            ha='left', va='center', transform=ax.get_yaxis_transform(),
            path_effects=[path_effects.withStroke(linewidth=3, foreground='white')])
    ax.text(-0.03, tba_y, TEXT_CONFIG['species_label_tba'],
            fontsize=TEXT_CONFIG['species_label_fontsize'], fontweight='bold',
            ha='left', va='center', transform=ax.get_yaxis_transform(),
            path_effects=[path_effects.withStroke(linewidth=3, foreground='white')])

    # 轴与标题
    max_x = max(tra_pos['x_end'].max(), tba_pos['x_end'].max())
    ax.set_xlim(-max_x * 0.02, max_x * 1.02)
    ax.set_ylim(-0.8, 2.1)
    ax.set_title(
        f"{TEXT_CONFIG['title_main']}\n{TEXT_CONFIG['title_sub']}",
        fontsize=TEXT_CONFIG['title_fontsize'], fontweight='bold', pad=30
    )
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # 图例
    legend_elements = []
    if not breakpoints_plot.empty:
        legend_elements.append(
            plt.Line2D([0], [0], marker='>', linestyle='None',
                       color=ARROW_MARK_CONFIG.get('color', '#FF00A6'),
                       markersize=max(8, ARROW_MARK_CONFIG.get('mutation_scale', 22.0) * 0.5),
                       label=f"{TEXT_CONFIG['legend_breakpoint']} (arrows)")
        )

    legend_elements.extend([
        plt.Line2D([0], [0], color=COLORS['telomere_tra'], linewidth=TELOMERE_LINE_WIDTH,
                   label=f"{TEXT_CONFIG['legend_tra_telomere']} (n={tra_telomere_count})"),
        plt.Line2D([0], [0], color=COLORS['telomere_tba'], linewidth=TELOMERE_LINE_WIDTH,
                   label=f"{TEXT_CONFIG['legend_tba_telomere']} (n={tba_telomere_count})"),
    ])
    for i in range(min(3, len(tba_pos))):
        chr_name = tba_pos.iloc[i]['chr']
        legend_elements.append(
            plt.Rectangle((0, 0), 1, 1, facecolor=tba_colors[chr_name],
                          edgecolor=CHR_BORDER_CONFIG.get('color', COLORS['edge_color']),
                          label=f"{TEXT_CONFIG['legend_tba_chr']} {get_display_name(chr_name)}")
        )
    if len(tba_pos) > 3:
        legend_elements.append(
            plt.Rectangle((0, 0), 1, 1, facecolor='gray', alpha=0.3,
                          edgecolor=CHR_BORDER_CONFIG.get('color', COLORS['edge_color']),
                          label=TEXT_CONFIG['legend_more_chrs'].format(len(tba_pos) - 3))
        )
    if len(synteny_plot) > 0:
        legend_elements.append(
            plt.Line2D([0], [0], color='gray', alpha=COLORS['synteny_alpha'], linewidth=2.5,
                       label=TEXT_CONFIG['legend_synteny'].format(MIN_SYNTENY_SIZE, len(synteny_plot)))
        )
    total_reversed = tra_pos['is_reversed'].sum() + tba_pos['is_reversed'].sum()
    if total_reversed > 0 and REVERSE_LABEL_MARK.get('show', True):
        legend_elements.append(
            plt.Line2D([0], [0], marker='*', linestyle='None',
                       color=REVERSE_LABEL_MARK.get('color', 'red'),
                       markersize=REVERSE_LABEL_MARK.get('legend_markersize', 14),
                       label=f"{TEXT_CONFIG['legend_reversed']} (n={total_reversed})")
        )

    legend = ax.legend(handles=legend_elements, loc='upper right',
                       bbox_to_anchor=(1.15, 1), frameon=True,
                       fancybox=True, shadow=True, fontsize=TEXT_CONFIG['legend_fontsize'])
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(0.95)
    legend.get_frame().set_edgecolor(CHR_BORDER_CONFIG.get('color', COLORS['edge_color']))

    plt.tight_layout()
    return fig, len(synteny_plot), tra_telomere_count, tba_telomere_count


def save_enhanced_results(fig, tra_pos, tba_pos, synteny_count, breakpoints_count,
                          filtered_synteny_count, tra_telomere_count, tba_telomere_count):
    print("正在保存结果...")
    filename_base = TEXT_CONFIG['output_filename_base']
    output_files = {
        'pdf': f'{OUTPUT_DIR}/{filename_base}.pdf',
        'png': f'{OUTPUT_DIR}/{filename_base}.png',
        'svg': f'{OUTPUT_DIR}/{filename_base}.svg'
    }
    fig.savefig(output_files['pdf'], dpi=DPI, bbox_inches='tight', facecolor='white', edgecolor='none', format='pdf')
    fig.savefig(output_files['png'], dpi=DPI, bbox_inches='tight', facecolor='white', edgecolor='none', format='png')
    fig.savefig(output_files['svg'], bbox_inches='tight', facecolor='white', edgecolor='none', format='svg')

    tra_reversed = list(tra_pos[tra_pos['is_reversed']]['chr'])
    tba_reversed = list(tba_pos[tba_pos['is_reversed']]['chr'])

    stats = f"""{TEXT_CONFIG['stats_title']}
======================================
缩放与标注:
- 端粒高度: 模式={TELOMERE_CONFIG.get('mode')}, gamma={TELOMERE_CONFIG.get('gamma')}, 高度范围=[{TELOMERE_CONFIG['min_height']}, {TELOMERE_CONFIG['max_height']}]
- 断点标注: 箭头 color={ARROW_MARK_CONFIG.get('color')}, mutation_scale={ARROW_MARK_CONFIG.get('mutation_scale')}, linewidth={ARROW_MARK_CONFIG.get('linewidth')}, offset={ARROW_MARK_CONFIG.get('offset')}
- 断点编号: show={ARROW_MARK_CONFIG.get('show_labels')}, fontsize={ARROW_MARK_CONFIG.get('label_fontsize')}, color={ARROW_MARK_CONFIG.get('label_color')}
- Tra间隔: {TRA_CHR_SPACING}; Tba组内间隔: {TBA_INTRA_GROUP_SPACING}; Tba分组数: {len(TBA_GROUPS_BY_TRA)}
- 绘制数量限制: TRA={TRA_CHR_LIMIT}, TBA={TBA_CHR_LIMIT}
- 标签字号: 默认={LABEL_CONFIG['fontsize']}, per-chr={len(LABEL_SIZE_PER_CHR)}条；星标字号：专用={REVERSE_LABEL_MARK.get('fontsize')}, scale={REVERSE_LABEL_MARK.get('scale')}, per-chr={len(STAR_SIZE_PER_CHR)}条
- 显示名映射: {len(CHR_DISPLAY_NAME)} 条

染色体反转:
- TRA反转数量: {len(tra_reversed)}；列表: {', '.join(tra_reversed) if tra_reversed else '无'}
- TBA反转数量: {len(tba_reversed)}；列表: {', '.join(tba_reversed) if tba_reversed else '无'}

统计:
- TRA染色体: {len(tra_pos)}；TBA染色体: {len(tba_pos)}
- 共线性原始区域: {synteny_count:,}；过滤后: {filtered_synteny_count:,}；过滤率: {(synteny_count - filtered_synteny_count) / synteny_count * 100:.1f}%
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
        print("开始优化版生物信息学共线性图分析...")
        print("=" * 70)

        tra_chr, tba_chr, synteny, tra_telomere, tba_telomere, chr_order, breakpoints = load_all_data()
        tra_pos, tba_pos = prepare_chromosome_data(tra_chr, tba_chr, chr_order, synteny)

        fig, filtered_synteny_count, tra_telomere_count, tba_telomere_count = create_enhanced_synteny_plot(
            tra_pos, tba_pos, synteny, tra_telomere, tba_telomere, breakpoints
        )

        save_enhanced_results(fig, tra_pos, tba_pos, len(synteny), len(breakpoints),
                              filtered_synteny_count, tra_telomere_count, tba_telomere_count)

        elapsed_time = time.time() - start_time
        print("=" * 70)
        print(f"✓ 完成，总耗时: {elapsed_time:.2f} 秒")
        print(f"✓ TRA/TBA 染色体: {len(tra_pos)}/{len(tba_pos)}")
        print(f"✓ 共线性高质量区域: {filtered_synteny_count:,}")
        print(f"✓ 断点条目: {len(breakpoints)}")
        print(f"✓ 端粒（TRA/TBA）: {tra_telomere_count}/{tba_telomere_count}")

    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()
