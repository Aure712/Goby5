#!/usr/bin/env python3
"""
按物种独立频率分类的区域统计脚本（仅保留主分析，支持多个 region_bed_file）

相对上一版脚本的新增改动：
1. 支持一次性对多个 region_bed_file 进行统计。
2. 每个 region_bed_file 可单独指定输出前缀 prefix。
3. 每组 region_bed_file 的所有输出文件都会自动带上对应前缀，避免结果互相覆盖。
4. 仍然保留：
   - 仅主分析
   - 原有 NA 筛选与频率分箱逻辑
   - 输出统计区间内的原始详细基因型 BED
   - 统计 bed_file 第 8 列 mRNA
   - 每个物种输出统计区间总长度和基因型值总和（0/1/2 直接求和）

输入约定：
- bed_file: 原始位点级 BED，至少前 8 列存在，个体基因型从第 12 列开始（0-based 第 11 列）
- region_bed_file: 统计区间 BED，至少 3 列(chr, start, end)
- sample_file: 两列，species<TAB>individual
- fai_file: 染色体长度文件（这里仅用于基础检查/打印，可不参与统计）
"""

import os
import re
import sys
from bisect import bisect_right
from collections import defaultdict

import pandas as pd


# =========================
# 读取基础文件
# =========================
def read_sample_info(sample_file):
    """读取样本信息文件, 获取物种-个体的对应关系"""
    species_individuals = defaultdict(list)
    individual_to_species = {}
    individual_list = []

    with open(sample_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                species = parts[0]
                individual = parts[1]
                species_individuals[species].append(individual)
                individual_to_species[individual] = species
                individual_list.append(individual)

    species_individuals = dict(species_individuals)

    print(f"读取到 {len(species_individuals)} 个物种:")
    for species, individuals in species_individuals.items():
        print(f"  {species}: {len(individuals)} 个个体")

    return species_individuals, individual_to_species, individual_list


def read_fai_file(fai_file):
    """读取 .fai 文件获取染色体长度信息"""
    chr_lengths = {}
    with open(fai_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chr_name = parts[0]
                chr_length = int(parts[1])
                chr_lengths[chr_name] = chr_length
    return chr_lengths


def read_region_bed(region_bed_file):
    """读取统计区间 BED，并按染色体合并/整理区间"""
    regions_by_chr = defaultdict(list)

    with open(region_bed_file, 'r') as f:
        for line_no, line in enumerate(f, start=1):
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                print(f"警告: region bed 第 {line_no} 行列数不足 3，跳过")
                continue
            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                print(f"警告: region bed 第 {line_no} 行坐标无法解析，跳过")
                continue

            if end < start:
                start, end = end, start
            if end == start:
                end = start + 1

            regions_by_chr[chrom].append((start, end))

    merged_regions_by_chr = {}
    for chrom, intervals in regions_by_chr.items():
        intervals.sort(key=lambda x: (x[0], x[1]))
        merged = []
        for s, e in intervals:
            if not merged or s > merged[-1][1]:
                merged.append([s, e])
            else:
                merged[-1][1] = max(merged[-1][1], e)
        merged_regions_by_chr[chrom] = [(s, e) for s, e in merged]

    total_regions = sum(len(v) for v in merged_regions_by_chr.values())
    print(f"读取统计区间: {len(merged_regions_by_chr)} 条染色体, 合并后 {total_regions} 个区间")

    return merged_regions_by_chr


def calculate_total_region_length(regions_by_chr):
    """计算所有统计区间（合并后）的总长度"""
    total_length = 0
    for intervals in regions_by_chr.values():
        for start, end in intervals:
            total_length += (end - start)
    return total_length


# =========================
# 基因型处理逻辑（保留原逻辑）
# =========================
def filter_na_by_species_row(row_data, species_individuals, individual_to_species,
                             individual_list, na_threshold, species_filtered_stats):
    """对单行数据进行物种级别的 NA 筛选"""
    filtered_row = list(row_data.copy())

    for species, individuals in species_individuals.items():
        species_data = []
        species_indices = []

        for individual in individuals:
            if individual in individual_list:
                idx = individual_list.index(individual)
                if idx < len(row_data):
                    species_data.append(row_data[idx])
                    species_indices.append(idx)

        if species_data:
            na_count = sum(
                1 for val in species_data
                if str(val).strip().upper() == 'NA' or pd.isna(val) or str(val).strip() == ''
            )
            na_ratio = na_count / len(species_data)

            if na_ratio > na_threshold:
                for idx in species_indices:
                    filtered_row[idx] = -999
                species_filtered_stats[species] += 1

    return filtered_row


def convert_genotype_value(value):
    """将基因型值转换为内部数值编码"""
    if value == -999 or str(value) == '-999':
        return -999

    value_str = str(value).strip()

    if value_str.upper() == 'NA' or value_str == 'nan' or pd.isna(value) or value_str == '':
        return -888

    if '.' in value_str:
        try:
            value_str = str(int(float(value_str)))
        except Exception:
            pass

    if value_str == '0':
        return 0
    if value_str == '1':
        return 1
    if value_str == '2':
        return 2

    print(f"警告: 未知基因型值 '{value}'，按 0 处理")
    return 0


def process_genotype_data(df, species_individuals, individual_to_species,
                          individual_list, na_threshold=0.5):
    """统一处理基因型数据，包括物种级别 NA 筛选和格式转换"""
    print(f"处理基因型数据和 NA 筛选(阈值: {na_threshold * 100}%)...")

    processed_df = df.copy()
    species_filtered_stats = {species: 0 for species in species_individuals.keys()}

    for row_idx in range(len(df)):
        row_genotypes = df.iloc[row_idx, 11:].values
        filtered_genotypes = filter_na_by_species_row(
            row_genotypes, species_individuals, individual_to_species,
            individual_list, na_threshold, species_filtered_stats
        )

        processed_genotypes = [convert_genotype_value(val) for val in filtered_genotypes]

        for i, processed_val in enumerate(processed_genotypes):
            processed_df.iloc[row_idx, 11 + i] = processed_val

    print("NA 筛选统计:")
    total_filtered = sum(species_filtered_stats.values())
    print(f"  总共筛选了 {total_filtered} 个物种-行组合")
    for species, count in species_filtered_stats.items():
        percent = (count / len(df)) * 100 if len(df) > 0 else 0
        print(f"  {species}: {count} 行被标记(占 {percent:.2f}%)")

    return processed_df, species_filtered_stats


# =========================
# 频率区间逻辑（保留原逻辑）
# =========================
def create_custom_intervals(interval_percentages):
    """创建自定义区间；保留旧脚本的分箱逻辑"""
    boundaries = {}
    interval_percentages = sorted(interval_percentages)

    has_zero = 0 in interval_percentages
    has_hundred = 100 in interval_percentages

    if has_zero:
        boundaries['exactly_0'] = (0, 0)
        interval_percentages = [p for p in interval_percentages if p != 0]

    if has_hundred:
        interval_percentages = [p for p in interval_percentages if p != 100]

    if interval_percentages:
        first_pct = interval_percentages[0]
        boundaries[f'above_0_to_{first_pct}'] = (0, first_pct, False, True)

        for i in range(len(interval_percentages) - 1):
            lower_pct = interval_percentages[i]
            upper_pct = interval_percentages[i + 1]
            boundaries[f'{lower_pct}_to_{upper_pct}'] = (lower_pct, upper_pct, False, True)

        last_pct = interval_percentages[-1]
        if has_hundred:
            boundaries[f'{last_pct}_to_below_100'] = (last_pct, 100.0, False, False)
            boundaries[f'{last_pct}_to_100'] = (last_pct, 100.0, False, True)
        else:
            boundaries[f'above_{last_pct}'] = (last_pct, float('inf'), False, True)

    if has_hundred:
        boundaries['exactly_100'] = (100.0, 100.0)

    return boundaries


def frequency_in_category(frequency, category, boundary_info):
    """判断频率是否落入某个区间（完全沿用旧脚本边界规则）"""
    if len(boundary_info) == 2:
        lower_pct, upper_pct = boundary_info
        if category == 'exactly_0':
            return frequency == 0
        if category == 'exactly_100':
            return frequency == 100
        if upper_pct == float('inf'):
            return frequency >= lower_pct
        return lower_pct <= frequency < upper_pct

    lower_pct, upper_pct, include_lower, include_upper = boundary_info

    lower_ok = frequency >= lower_pct if include_lower else frequency > lower_pct
    if upper_pct == float('inf'):
        upper_ok = True
    else:
        upper_ok = frequency <= upper_pct if include_upper else frequency < upper_pct
    return lower_ok and upper_ok


def categorize_by_species_frequency(df, species_individuals, individual_to_species,
                                    individual_list, interval_percentages=None):
    """
    按照每个物种各自的频率将数据分为多个区间。
    逻辑与旧版主分析保持一致。

    返回: {species: {category: dataframe}}
    """
    print("按物种独立频率分类数据...")

    if interval_percentages is not None:
        boundaries = create_custom_intervals(interval_percentages)
    else:
        boundaries = {
            'below_20': (0, 20),
            '20_to_40': (20, 40),
            '40_to_60': (40, 60),
            '60_to_80': (60, 80),
            'above_80': (80, float('inf'))
        }

    species_category_indices = {}

    for species, individuals in species_individuals.items():
        print(f"\n  处理物种: {species}")
        indices = [individual_list.index(ind) for ind in individuals if ind in individual_list]
        if not indices:
            continue

        species_category_indices[species] = {category: [] for category in boundaries.keys()}

        for row_idx in range(len(df)):
            row = df.iloc[row_idx]
            species_data = [row.iloc[11 + idx] for idx in indices]
            valid_data = [val for val in species_data if val != -999]

            if not valid_data:
                continue

            has_na = any(val == -888 for val in valid_data)
            if has_na:
                non_na_values = [val for val in valid_data if val != -888]
                if non_na_values:
                    frequency = (sum(non_na_values) / (len(non_na_values) * 2)) * 100
                else:
                    continue
            else:
                frequency = (sum(valid_data) / (len(valid_data) * 2)) * 100

            for category, boundary_info in boundaries.items():
                if frequency_in_category(frequency, category, boundary_info):
                    species_category_indices[species][category].append(row_idx)

    categorized_data = {}
    for species in species_category_indices.keys():
        categorized_data[species] = {}
        for category, indices in species_category_indices[species].items():
            if indices:
                categorized_data[species][category] = df.iloc[indices].copy()
                print(f"    物种 {species} 区间 {category}: {len(indices)} 行")
            else:
                categorized_data[species][category] = pd.DataFrame()

    return categorized_data


# =========================
# 区间筛选逻辑（新）
# =========================
def build_region_index(regions_by_chr):
    """为每条染色体构建 starts 数组，便于二分判断重叠"""
    region_index = {}
    for chrom, intervals in regions_by_chr.items():
        starts = [s for s, _ in intervals]
        region_index[chrom] = {
            'intervals': intervals,
            'starts': starts,
        }
    return region_index


def row_overlaps_regions(chrom, start, end, region_index):
    """判断一个 BED 记录是否与统计区间重叠（半开区间逻辑：[start, end)）"""
    if chrom not in region_index:
        return False

    intervals = region_index[chrom]['intervals']
    starts = region_index[chrom]['starts']

    pos = bisect_right(starts, end - 1) - 1
    while pos >= 0:
        s, e = intervals[pos]
        if e <= start:
            break
        if s < end and e > start:
            return True
        pos -= 1
    return False



def filter_bed_by_regions(df, region_index):
    """从原始/处理后的 BED 中筛选与统计区间重叠的记录"""
    keep_indices = []

    for row_idx in range(len(df)):
        chrom = df.iloc[row_idx, 0]
        try:
            start = int(float(df.iloc[row_idx, 1]))
        except Exception:
            continue

        try:
            end = int(float(df.iloc[row_idx, 2]))
        except Exception:
            end = start + 1

        if end <= start:
            end = start + 1

        if row_overlaps_regions(chrom, start, end, region_index):
            keep_indices.append(row_idx)

    filtered_df = df.iloc[keep_indices].copy()
    print(f"统计区间筛选后保留 {len(filtered_df)} / {len(df)} 行")
    return filtered_df


# =========================
# mRNA 统计（新）
# =========================
def collect_retained_mrnas(df, mrna_col_index=7):
    """统计保留下来的 mRNA（来自 bed_file 第 8 列）"""
    if df.empty:
        return {
            'total_unique_mrnas': 0,
            'unique_mrna_list': []
        }

    mrna_col = df.iloc[:, mrna_col_index]
    unique_mrnas = set()
    for value in mrna_col:
        if pd.notna(value):
            value_str = str(value).strip()
            if value_str:
                unique_mrnas.add(value_str)

    return {
        'total_unique_mrnas': len(unique_mrnas),
        'unique_mrna_list': sorted(unique_mrnas)
    }


def calculate_species_genotype_value_sum(processed_category_df, species_name, species_individuals, individual_list):
    """
    计算某个物种在某个频率区间内的基因型值总和。
    仅累计该物种对应个体列中的有效基因型值 0/1/2；
    -888(NA) 和 -999(因 NA 阈值被整物种屏蔽) 不参与求和。
    """
    if processed_category_df.empty:
        return 0

    individual_indices = [
        individual_list.index(ind) for ind in species_individuals[species_name]
        if ind in individual_list
    ]
    if not individual_indices:
        return 0

    genotype_sum = 0
    for idx in individual_indices:
        col_values = processed_category_df.iloc[:, 11 + idx]
        for value in col_values:
            try:
                int_value = int(value)
            except Exception:
                continue

            if int_value in (0, 1, 2):
                genotype_sum += int_value

    return genotype_sum


# =========================
# 辅助函数（多 region_bed_file 支持）
# =========================
def sanitize_prefix(prefix):
    """将 prefix 规范化为适合文件名前缀的字符串"""
    prefix = str(prefix).strip()
    if not prefix:
        raise ValueError('prefix 不能为空')

    safe_prefix = re.sub(r'[^A-Za-z0-9._-]+', '_', prefix)
    safe_prefix = safe_prefix.strip('._-')
    if not safe_prefix:
        raise ValueError(f'prefix 规范化后为空: {prefix}')
    return safe_prefix



def make_prefixed_filename(prefix, filename):
    return f'{prefix}_{filename}'



def validate_region_configs(region_bed_configs):
    """检查多个 region 配置是否合法"""
    if not region_bed_configs:
        raise ValueError('region_bed_configs 不能为空')

    normalized_configs = []
    seen_prefixes = set()

    for i, config in enumerate(region_bed_configs, start=1):
        if not isinstance(config, dict):
            raise ValueError(f'第 {i} 个 region 配置不是字典')

        if 'prefix' not in config or 'region_bed_file' not in config:
            raise ValueError(f'第 {i} 个 region 配置必须同时包含 prefix 和 region_bed_file')

        prefix = sanitize_prefix(config['prefix'])
        region_bed_file = str(config['region_bed_file']).strip()
        if not region_bed_file:
            raise ValueError(f'第 {i} 个 region 配置的 region_bed_file 为空')

        if prefix in seen_prefixes:
            raise ValueError(f'prefix 重复: {prefix}')
        seen_prefixes.add(prefix)

        normalized_configs.append({
            'prefix': prefix,
            'region_bed_file': region_bed_file,
        })

    return normalized_configs



def run_single_region_analysis(prefix, region_bed_file, raw_df,
                               species_individuals, individual_to_species,
                               individual_list, na_threshold,
                               interval_percentages):
    """对单个 region_bed_file 运行完整主分析"""
    print(f"\n{'#' * 100}")
    print(f"开始处理 region 配置: prefix={prefix}, region_bed_file={region_bed_file}")
    print(f"{'#' * 100}")

    print(f'\n读取统计区间 BED 文件: {region_bed_file}')
    try:
        regions_by_chr = read_region_bed(region_bed_file)
    except FileNotFoundError:
        print(f'错误: 找不到文件 {region_bed_file}')
        return None

    region_index = build_region_index(regions_by_chr)
    total_region_length = calculate_total_region_length(regions_by_chr)
    print(f'统计区间总长度(合并后): {total_region_length}')

    print('\n基于统计区间筛选原始详细 BED 记录...')
    raw_region_df = filter_bed_by_regions(raw_df, region_index)

    if len(raw_region_df) == 0:
        print(f'prefix={prefix} 对应的统计区间内没有任何记录，跳过该 region_bed_file。')
        return {
            'prefix': prefix,
            'region_bed_file': region_bed_file,
            'total_region_length': total_region_length,
            'raw_region_count': 0,
            'species_filtered_stats': {species: 0 for species in species_individuals.keys()},
            'species_region_and_genotype_summaries': {},
            'all_mrna_stats': {},
        }

    print('\n开始处理统计区间内的基因型数据...')
    processed_region_df, species_filtered_stats = process_genotype_data(
        raw_region_df, species_individuals, individual_to_species, individual_list, na_threshold
    )

    print('\n处理后的数据样本(前 5 行, 前 3 个个体):')
    end_col = min(14, processed_region_df.shape[1])
    print(processed_region_df.iloc[:5, 11:end_col])
    print(f'数据类型示例: {type(processed_region_df.iloc[0, 11])}')

    print(f"\n{'#' * 100}")
    print(f'开始主分析: 按物种独立频率分类（prefix={prefix}）')
    print(f"{'#' * 100}")

    categorized_data = categorize_by_species_frequency(
        processed_region_df, species_individuals, individual_to_species,
        individual_list, interval_percentages
    )

    boundaries = create_custom_intervals(interval_percentages)
    print('\n生成的区间:')
    for category, boundary in boundaries.items():
        print(f'  {category}: {boundary}')

    all_mrna_stats = {}
    species_region_and_genotype_summaries = {}

    for species_name, category_dict in categorized_data.items():
        print(f"\n{'=' * 60}")
        print(f'处理物种: {species_name}')
        print(f"包含 {len(species_individuals[species_name])} 个个体: {', '.join(species_individuals[species_name])}")
        print(f"{'=' * 60}")

        species_mrna_stats = []
        species_total_genotype_value_sum = 0
        per_category_genotype_value_sum = {}

        for category, processed_category_df in category_dict.items():
            if len(processed_category_df) == 0:
                print(f'\n物种 {species_name} 区间 {category} 没有数据, 跳过')
                continue

            print(f'\n处理物种 {species_name} 区间: {category}')

            category_indices = processed_category_df.index
            raw_category_df = raw_region_df.loc[category_indices].copy()

            output_file = make_prefixed_filename(prefix, f'detailed_genotypes_{species_name}_{category}.bed')
            raw_category_df.to_csv(output_file, sep='\t', index=False, header=False)
            print(f'  详细基因型结果保存到: {output_file} ({len(raw_category_df)} 行)')

            genotype_value_sum = calculate_species_genotype_value_sum(
                processed_category_df, species_name, species_individuals, individual_list
            )
            per_category_genotype_value_sum[category] = genotype_value_sum
            species_total_genotype_value_sum += genotype_value_sum
            print(f'  该物种在区间 {category} 的基因型值总和: {genotype_value_sum}')

            mrna_stats = collect_retained_mrnas(raw_category_df, mrna_col_index=7)
            mrna_stats['species'] = species_name
            mrna_stats['category'] = category
            mrna_stats['genotype_value_sum'] = genotype_value_sum
            species_mrna_stats.append(mrna_stats)

            mrna_list_file = make_prefixed_filename(prefix, f'mrna_list_{species_name}_{category}.txt')
            with open(mrna_list_file, 'w') as f:
                for mrna in mrna_stats['unique_mrna_list']:
                    f.write(f'{mrna}\n')

            print(f"  保留的唯一 mRNA 数: {mrna_stats['total_unique_mrnas']}")
            print(f'  mRNA 列表保存到: {mrna_list_file}')

        summary_file = make_prefixed_filename(prefix, f'region_length_and_genotype_sum_{species_name}.txt')
        with open(summary_file, 'w') as f:
            f.write(f'prefix\t{prefix}\n')
            f.write(f'region_bed_file\t{region_bed_file}\n')
            f.write(f'species\t{species_name}\n')
            f.write(f'total_region_length\t{total_region_length}\n')
            f.write(f'total_genotype_value_sum_across_all_frequencies\t{species_total_genotype_value_sum}\n')
            for category, genotype_sum in per_category_genotype_value_sum.items():
                f.write(f'genotype_value_sum_{category}\t{genotype_sum}\n')

        print(f'  统计区间总长度和基因型值总和已保存到: {summary_file}')

        all_mrna_stats[species_name] = species_mrna_stats
        species_region_and_genotype_summaries[species_name] = {
            'total_region_length': total_region_length,
            'total_genotype_value_sum_across_all_frequencies': species_total_genotype_value_sum,
            'per_category_genotype_value_sum': per_category_genotype_value_sum,
        }

    summary_main_file = make_prefixed_filename(prefix, 'mrna_statistics_summary_main_analysis.txt')
    with open(summary_main_file, 'w') as f:
        f.write('主分析 mRNA 统计汇总（仅统计 region_bed_file 区间内记录）\n')
        f.write(f'prefix: {prefix}\n')
        f.write(f'统计区间文件: {region_bed_file}\n')
        f.write(f'NA 筛选阈值: {na_threshold * 100}%\n')
        f.write(f'区间百分比: {interval_percentages}\n')
        f.write(f'物种数量: {len(species_individuals)}\n')
        f.write(f'统计区间总长度: {total_region_length}\n')
        f.write(f'统计区间内记录数: {len(raw_region_df)}\n')
        f.write('说明: 仅保留主分析；频率计算逻辑与旧版一致；mRNA 统计来自输出结果中保留的 bed_file 第 8 列\n')
        f.write('=' * 60 + '\n\n')

        for species_name in species_individuals.keys():
            mrna_stats_list = all_mrna_stats.get(species_name, [])
            species_summary = species_region_and_genotype_summaries.get(
                species_name,
                {
                    'total_region_length': total_region_length,
                    'total_genotype_value_sum_across_all_frequencies': 0,
                    'per_category_genotype_value_sum': {},
                }
            )
            f.write(f'物种: {species_name}\n')
            f.write(f'个体数: {len(species_individuals[species_name])}\n')
            f.write(f"个体列表: {', '.join(species_individuals[species_name])}\n")
            f.write(f"NA 筛选影响: {species_filtered_stats.get(species_name, 0)} 行被标记\n")
            f.write(f"统计区间总长度: {species_summary['total_region_length']}\n")
            f.write(
                f"所有频率下该物种基因型值总和: "
                f"{species_summary['total_genotype_value_sum_across_all_frequencies']}\n"
            )
            f.write('-' * 40 + '\n')

            if not mrna_stats_list:
                f.write('  无保留结果\n\n')
                continue

            for stats in mrna_stats_list:
                f.write(f"  区间: {stats['category']}\n")
                f.write(f"  保留的唯一 mRNA 数: {stats['total_unique_mrnas']}\n")
                f.write(f"  该物种该频率下基因型值总和: {stats['genotype_value_sum']}\n")
                f.write('  ' + '-' * 30 + '\n')
            f.write('\n')

    print(f'\n主汇总文件保存到: {summary_main_file}')

    return {
        'prefix': prefix,
        'region_bed_file': region_bed_file,
        'total_region_length': total_region_length,
        'raw_region_count': len(raw_region_df),
        'species_filtered_stats': species_filtered_stats,
        'species_region_and_genotype_summaries': species_region_and_genotype_summaries,
        'all_mrna_stats': all_mrna_stats,
    }


# =========================
# 主程序
# =========================
def main():
    # -------------------------
    # 文件路径（按需修改）
    # -------------------------
    fai_file = 'Tra_sort.fa.fai'
    bed_file = 'allpops.refTra.all.syn.bed' # 'allpops.refTra.all.high.bed'
    sample_file = 'sample_name.txt'

    # 这里改成“多个 region_bed_file + 自定义前缀”
    # 每个元素都是一个字典：
    #   prefix: 该组结果文件的前缀
    #   region_bed_file: 对应的统计区间 bed 文件
    region_bed_configs = [
        {
            'prefix': 'invert_FSR',
            'region_bed_file': 'palette_regions.invert_ismc.merged.bed',
        },
        {
            'prefix': 'invert_NOFSR',
            'region_bed_file': 'palette_regions.invert_non_ismc.merged.bed',
        },
        {
            'prefix': 'NOinvert_FSR',
            'region_bed_file': 'palette_regions.non_invert_ismc.merged.bed',
        },
        {
            'prefix': 'NOinvert_NOFSR',
            'region_bed_file': 'palette_regions.non_invert_non_ismc.merged.bed',
        }
        # 需要更多 region bed 时，按下面格式继续添加：
        # {
        #     'prefix': 'your_prefix_2',
        #     'region_bed_file': 'your_region_file_2.bed',
        # },
    ]

    # -------------------------
    # 参数设置
    # -------------------------
    na_threshold = 0.5
    interval_percentages = [0, 25, 50, 100]

    try:
        region_bed_configs = validate_region_configs(region_bed_configs)
    except ValueError as e:
        print(f'错误: {e}')
        sys.exit(1)

    print('按物种独立频率分类的区域统计脚本（仅主分析，支持多个 region_bed_file）')
    print(f'NA 筛选阈值: {na_threshold * 100}%')
    print(f'区间设置: {interval_percentages}')
    print('region_bed_file 配置:')
    for config in region_bed_configs:
        print(f"  prefix={config['prefix']} -> {config['region_bed_file']}")
    print('=' * 80)
    print('输出说明:')
    print('1. 仅保留主分析')
    print('2. 不再进行滑动窗口统计')
    print('3. 支持一次性对多个 region_bed_file 做统计')
    print('4. 每个 region_bed_file 的输出结果都带对应 prefix')
    print('5. 输出为各物种-各频率区间内、且落在对应 region_bed_file 中的原始详细基因型 BED 记录')
    print('6. mRNA 统计基于 bed_file 第 8 列中被保留下来的 mRNA')
    print('7. 每个物种额外输出统计区间总长度和该物种基因型值总和(0/1/2 直接求和)')
    print('=' * 80)

    print(f'\n读取样本信息文件: {sample_file}')
    try:
        species_individuals, individual_to_species, individual_list = read_sample_info(sample_file)
    except FileNotFoundError:
        print(f'错误: 找不到文件 {sample_file}')
        sys.exit(1)

    if os.path.exists(fai_file):
        print(f'\n读取染色体长度信息: {fai_file}')
        chr_lengths = read_fai_file(fai_file)
        print(f'发现 {len(chr_lengths)} 条染色体')
    else:
        chr_lengths = {}
        print(f'\n提示: 未找到 {fai_file}，继续运行（本脚本不再依赖滑窗统计）')

    print(f'\n读取位点 BED 文件: {bed_file}')
    try:
        raw_df = pd.read_csv(bed_file, sep='\t', header=None, dtype=str)
    except FileNotFoundError:
        print(f'错误: 找不到文件 {bed_file}')
        sys.exit(1)

    print(f'位点 BED 总行数: {len(raw_df)}')
    print(f'位点 BED 总列数: {len(raw_df.columns)}')

    if raw_df.shape[1] < 12:
        print('错误: 输入 bed_file 列数不足，至少应包含 11 列注释信息 + 个体基因型列')
        sys.exit(1)

    n_individuals_total = raw_df.iloc[:, 11:].shape[1]
    print(f'BED 文件中检测到 {n_individuals_total} 个个体列')
    if n_individuals_total != len(individual_list):
        print(f'错误: BED 文件个体数({n_individuals_total})与样本文件个体数({len(individual_list)})不一致')
        sys.exit(1)

    combined_results = []
    for config in region_bed_configs:
        result = run_single_region_analysis(
            prefix=config['prefix'],
            region_bed_file=config['region_bed_file'],
            raw_df=raw_df,
            species_individuals=species_individuals,
            individual_to_species=individual_to_species,
            individual_list=individual_list,
            na_threshold=na_threshold,
            interval_percentages=interval_percentages,
        )
        if result is not None:
            combined_results.append(result)

    combined_summary_file = 'combined_region_prefix_summary.txt'
    with open(combined_summary_file, 'w') as f:
        f.write('多个 region_bed_file 统计结果总览\n')
        f.write(f'位点文件: {bed_file}\n')
        f.write(f'样本文件: {sample_file}\n')
        f.write(f'NA 筛选阈值: {na_threshold * 100}%\n')
        f.write(f'区间百分比: {interval_percentages}\n')
        f.write('=' * 80 + '\n\n')

        for result in combined_results:
            f.write(f"prefix: {result['prefix']}\n")
            f.write(f"region_bed_file: {result['region_bed_file']}\n")
            f.write(f"total_region_length: {result['total_region_length']}\n")
            f.write(f"raw_region_count: {result['raw_region_count']}\n")
            for species_name in species_individuals.keys():
                species_summary = result['species_region_and_genotype_summaries'].get(
                    species_name,
                    {
                        'total_genotype_value_sum_across_all_frequencies': 0,
                    }
                )
                f.write(
                    f"species={species_name}\t"
                    f"total_genotype_value_sum_across_all_frequencies="
                    f"{species_summary['total_genotype_value_sum_across_all_frequencies']}\n"
                )
            f.write('-' * 60 + '\n')

    print('\n分析完成!')
    print('=' * 80)
    print('本版脚本已实现:')
    print('1. 仅保留第一个主分析')
    print('2. 计数/频率分类逻辑保持不变')
    print('3. 支持一次性对多个 region_bed_file 进行统计')
    print('4. 可通过 prefix 为不同 region_bed_file 的结果指定不同输出前缀')
    print('5. 输出仍为统计区间内的详细基因型 BED 记录')
    print('6. 不再使用 mrna_ref_file，而是统计保留下来的 bed_file 第 8 列 mRNA')
    print('7. 每个物种额外输出统计区间总长度和该物种在所有频率下的基因型值总和(0/1/2 直接求和)')
    print(f'8. 额外总览文件: {combined_summary_file}')
    print('=' * 80)


if __name__ == '__main__':
    main()
