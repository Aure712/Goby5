import pandas as pd
import numpy as np
from collections import defaultdict


def read_bed_file(filename):
    """读取完整的BED文件，包括位置信息和基因型数据"""
    df = pd.read_csv(filename, sep='\t', header=None, low_memory=False)

    # 提取位置信息（前3列：chrom, start, end）
    position_info = df.iloc[:, :3].copy()
    position_info.columns = ['chrom', 'start', 'end']
    position_info['start'] = pd.to_numeric(position_info['start'], errors='coerce')
    position_info['end'] = pd.to_numeric(position_info['end'], errors='coerce')

    # 提取从第12列开始的基因型数据（Python索引11开始）
    genotype_data = df.iloc[:, 11:].copy()
    genotype_data = genotype_data.apply(pd.to_numeric, errors='coerce')

    return position_info, genotype_data


def read_sample_file(filename):
    """读取样本文件，获取物种和个体名信息"""
    sample_df = pd.read_csv(filename, sep='\t', header=None, names=['species', 'individual'])
    return sample_df


def read_regions_file(filename):
    """
    读取区域文件（BED格式）
    使用第1列作为染色体，第2列作为起始位置，第3列作为终止位置
    即：BED的第2、3列 -> Python索引的第1、2列
    """
    regions_df = pd.read_csv(
        filename,
        sep='\t',
        header=None,
        comment='#',
        low_memory=False
    )

    if regions_df.shape[1] < 3:
        raise ValueError(f"区域文件 {filename} 列数不足，BED文件至少应包含3列")

    regions_df = regions_df.iloc[:, [0, 1, 2]].copy()
    regions_df.columns = ['chrom', 'start', 'end']
    regions_df['start'] = pd.to_numeric(regions_df['start'], errors='coerce')
    regions_df['end'] = pd.to_numeric(regions_df['end'], errors='coerce')

    # 去掉无效行
    regions_df = regions_df.dropna(subset=['chrom', 'start', 'end']).copy()
    regions_df['start'] = regions_df['start'].astype(int)
    regions_df['end'] = regions_df['end'].astype(int)

    return regions_df


def build_region_dict(regions_df):
    """将区域按染色体整理成字典，便于快速判断位置是否落在区域内"""
    region_dict = defaultdict(list)

    for row in regions_df.itertuples(index=False):
        chrom = str(row.chrom)
        start = int(row.start)
        end = int(row.end)

        # 若区间写反，则自动纠正
        if end < start:
            start, end = end, start

        region_dict[chrom].append((start, end))

    # 对每条染色体上的区间按start排序
    for chrom in region_dict:
        region_dict[chrom].sort(key=lambda x: x[0])

    return region_dict


def is_in_regions(chrom, pos, region_dict):
    """
    判断某个位置是否在指定区域内

    这里按标准BED语义处理为半开区间：[start, end)
    即满足 start <= pos < end 视为在区域内
    """
    chrom = str(chrom)
    regions = region_dict.get(chrom, [])

    for start, end in regions:
        if pos < start:
            break
        if start <= pos < end:
            return True

    return False


def build_in_region_mask(position_info, region_dict):
    """为所有位点预先计算是否在区域内的布尔掩码"""
    mask = []

    for chrom, pos in zip(position_info['chrom'], position_info['start']):
        if pd.isna(pos):
            mask.append(False)
        else:
            mask.append(is_in_regions(chrom, int(pos), region_dict))

    return np.array(mask, dtype=bool)


def check_other_species_zero_or_na(row_data, current_species, all_species_indices):
    """检查除当前物种外，其他物种是否全为0或NA"""
    for other_species, other_indices in all_species_indices.items():
        if other_species == current_species:
            continue

        other_data = row_data[other_indices]
        for value in other_data:
            if not np.isnan(value) and value != 0:
                return False

    return True


def calculate_allele_frequencies_in_regions(genotype_data, sample_df, in_region_mask):
    """
    仅对区域内位点进行统计：
    1. 原始统计
    2. 物种特异性统计（其他物种均为0或NA）
    3. 当前物种自身非零统计（当前物种等位基因计数 > 0）
    """

    # 创建物种到个体索引的映射
    species_to_indices = defaultdict(list)
    for idx, (species, individual) in enumerate(sample_df.values):
        species_to_indices[species].append(idx)

    # 转为numpy数组，提高处理速度
    all_genotype_values = genotype_data.to_numpy(dtype=float, copy=False)
    n_rows = all_genotype_values.shape[0]

    # 存储每个物种的统计结果（仅区域内）
    species_freq_counts_in = {}
    species_freq_counts_specific_in = {}
    species_freq_counts_nonzero_in = {}

    # 对每个物种进行处理
    for species, indices in species_to_indices.items():
        n_individuals = len(indices)
        max_alleles = n_individuals * 2  # 二倍体最大等位基因数

        # 初始化频率计数列表（索引0 ~ max_alleles）
        freq_counts_in = [0] * (max_alleles + 1)
        freq_counts_specific_in = [0] * (max_alleles + 1)
        freq_counts_nonzero_in = [0] * (max_alleles + 1)

        # 当前物种对应的数据
        species_data = all_genotype_values[:, indices]

        for row_idx in range(n_rows):
            # 只统计区域内
            if not in_region_mask[row_idx]:
                continue

            row_data = species_data[row_idx]

            # 如果NA数量超过个体数一半，则跳过
            na_count = np.sum(np.isnan(row_data))
            if na_count > n_individuals / 2:
                continue

            # 计算当前物种的等位基因总数
            has_na = np.any(np.isnan(row_data))

            if not has_na:
                allele_count = int(np.sum(row_data))
            else:
                valid_data = row_data[~np.isnan(row_data)]
                n_valid = len(valid_data)

                if n_valid > 0:
                    sum_valid = np.sum(valid_data)
                    adjusted_count = sum_valid * (n_individuals / n_valid)
                    allele_count = int(round(adjusted_count))
                else:
                    continue

            if 0 <= allele_count <= max_alleles:
                # 1) 原始统计（仅区域内）
                freq_counts_in[allele_count] += 1

                # 2) 物种特异性统计
                entire_row_data = all_genotype_values[row_idx]
                if check_other_species_zero_or_na(entire_row_data, species, species_to_indices):
                    freq_counts_specific_in[allele_count] += 1

                # 3) 当前物种自身非零统计
                if allele_count > 0:
                    freq_counts_nonzero_in[allele_count] += 1

        # 在列表末尾添加总和
        freq_counts_in.append(sum(freq_counts_in))
        freq_counts_specific_in.append(sum(freq_counts_specific_in))
        freq_counts_nonzero_in.append(sum(freq_counts_nonzero_in))

        species_freq_counts_in[species] = freq_counts_in
        species_freq_counts_specific_in[species] = freq_counts_specific_in
        species_freq_counts_nonzero_in[species] = freq_counts_nonzero_in

    return (
        species_freq_counts_in,
        species_freq_counts_specific_in,
        species_freq_counts_nonzero_in
    )


def write_results_by_region(
    species_freq_counts_in,
    species_freq_counts_specific_in,
    species_freq_counts_nonzero_in,
    output_filename,
    region_filename=None
):
    """将仅区域内的统计结果写入文件"""
    with open(output_filename, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("# 仅统计BED区域内位点的结果\n")
        if region_filename is not None:
            f.write(f"# 区域文件: {region_filename}\n")
        f.write("# 每行最后一个值为该类位点总数\n")
        f.write("=" * 80 + "\n")

        # 原始统计
        f.write("\n## 原始统计结果（区域内）\n")
        for species, freq_counts in species_freq_counts_in.items():
            freq_str = ','.join(map(str, freq_counts))
            f.write(f"{species}_in\t{freq_str}\n")

        # 物种特异性统计
        f.write("\n" + "=" * 80 + "\n")
        f.write("# 物种特异性统计结果（其他物种均为0或NA，且仅区域内）\n")
        f.write("=" * 80 + "\n")
        for species, freq_counts in species_freq_counts_specific_in.items():
            freq_str = ','.join(map(str, freq_counts))
            f.write(f"{species}_specific_in\t{freq_str}\n")

        # 当前物种自身非零统计
        f.write("\n" + "=" * 80 + "\n")
        f.write("# 物种自身非零统计结果（当前物种基因计数>0，且仅区域内）\n")
        f.write("=" * 80 + "\n")
        for species, freq_counts in species_freq_counts_nonzero_in.items():
            freq_str = ','.join(map(str, freq_counts))
            f.write(f"{species}_nonzero_in\t{freq_str}\n")

    print(f"结果已保存到 {output_filename}")


def main():
    # 手动指定输入文件名
    bed_file = "allpops.refTra.all.high.bed"     # 等位基因型计数文件
    sample_file = "sample_name.txt"            # 个体物种记录文件

    # 支持多个regions_file输入，并与多个output_file一一对应
    regions_files = [
        "palette_regions.non_invert_ismc.merged.bed",
        "palette_regions.non_invert_non_ismc.merged.bed",
        # "lines3.bed",
    ]

    output_files = [
        "allele_frequency_results_pureFSR_high.txt",
        "allele_frequency_results_pureControl_high.txt",
        # "allele_frequency_results_lines3.txt",
    ]

    print("开始处理文件...")

    try:
        # 兼容单个字符串输入
        if isinstance(regions_files, str):
            regions_files = [regions_files]
        if isinstance(output_files, str):
            output_files = [output_files]

        if len(regions_files) != len(output_files):
            raise ValueError(
                f"regions_files数量({len(regions_files)})与output_files数量({len(output_files)})不一致"
            )

        # 读取基因型数据
        print(f"读取基因型数据文件: {bed_file}")
        position_info, genotype_data = read_bed_file(bed_file)
        print(f"  - 读取到 {len(genotype_data)} 个等位基因位点")
        print(f"  - 包含 {len(genotype_data.columns)} 个个体")

        # 读取样本信息
        print(f"\n读取样本信息文件: {sample_file}")
        sample_df = read_sample_file(sample_file)
        print(f"  - 读取到 {len(sample_df)} 个个体信息")

        # 检查个体数量是否匹配
        if len(genotype_data.columns) != len(sample_df):
            print(
                f"警告：基因型数据中的个体数({len(genotype_data.columns)})"
                f"与样本文件中的个体数({len(sample_df)})不匹配"
            )

        # 逐个处理多个区域文件
        for idx, (regions_file, output_file) in enumerate(zip(regions_files, output_files), start=1):
            print("\n" + "=" * 80)
            print(f"开始处理第 {idx} 组区域文件")
            print(f"读取区域坐标文件(BED): {regions_file}")

            regions_df = read_regions_file(regions_file)
            print(f"  - 读取到 {len(regions_df)} 个目标区域")

            # 构建区域索引
            region_dict = build_region_dict(regions_df)

            # 统计有多少位点在区域内
            in_region_mask = build_in_region_mask(position_info, region_dict)
            in_region_count = int(np.sum(in_region_mask))

            print(f"  - {in_region_count} 个位点在区域内")
            print(f"  - 仅对这些区域内位点进行统计")

            # 计算等位基因频率（仅区域内）
            print("\n计算等位基因频率（仅区域内，包括原始统计、物种特异性统计和自身非零统计）...")
            (
                species_freq_counts_in,
                species_freq_counts_specific_in,
                species_freq_counts_nonzero_in
            ) = calculate_allele_frequencies_in_regions(
                genotype_data,
                sample_df,
                in_region_mask
            )

            # 输出统计信息
            print("\n统计结果:")
            print("-" * 50)

            print("原始统计（仅区域内）：")
            for species in species_freq_counts_in.keys():
                n_individuals = len([s for s in sample_df['species'] if s == species])
                print(f"  物种 {species}: {n_individuals} 个个体, 最大等位基因数 {n_individuals * 2}")
                print(f"    - 区域内: {species_freq_counts_in[species][-1]} 个等位基因位点")

            print("\n物种特异性统计（其他物种均为0或NA，仅区域内）：")
            for species in species_freq_counts_specific_in.keys():
                print(f"  物种 {species}:")
                print(f"    - 区域内: {species_freq_counts_specific_in[species][-1]} 个物种特异性位点")

            print("\n物种自身非零统计（当前物种基因计数>0，仅区域内）：")
            for species in species_freq_counts_nonzero_in.keys():
                print(f"  物种 {species}:")
                print(f"    - 区域内: {species_freq_counts_nonzero_in[species][-1]} 个非零位点")

            # 写入结果
            print(f"\n保存结果到文件: {output_file}")
            write_results_by_region(
                species_freq_counts_in,
                species_freq_counts_specific_in,
                species_freq_counts_nonzero_in,
                output_file,
                region_filename=regions_file
            )

        print("\n全部处理完成！")

    except FileNotFoundError as e:
        print(f"错误：找不到文件 - {e}")
    except Exception as e:
        print(f"错误：{e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
