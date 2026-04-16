import pandas as pd
import numpy as np
from collections import defaultdict


def read_bed_file(filename):
    """读取完整的BED文件，包括位置信息和基因型数据"""
    # 读取文件，使用tab分隔符
    df = pd.read_csv(filename, sep='\t', header=None, low_memory=False)

    # 提取位置信息（通常是前3列：染色体、起始、终止）
    position_info = df.iloc[:, :3].copy()
    position_info.columns = ['chrom', 'start', 'end']

    # 提取从第12列开始的基因型数据（Python索引从0开始，所以是第11列）
    genotype_data = df.iloc[:, 11:].copy()

    # 将字符串'NA'转换为np.nan，数字保持不变
    for col in genotype_data.columns:
        genotype_data[col] = pd.to_numeric(genotype_data[col], errors='coerce')

    return position_info, genotype_data




def read_sample_file(filename):
    """读取样本文件，获取物种和个体名信息"""
    sample_df = pd.read_csv(filename, sep='\t', header=None, names=['species', 'individual'])
    return sample_df


def read_regions_file(filename):
    """读取区域文件，获取目标区域坐标"""
    regions_df = pd.read_csv(filename, sep='\t', header=None,
                             names=['chrom', 'col2', 'start', 'end'])
    # 只保留需要的列（染色体、起始、终止）
    regions_df = regions_df[['chrom', 'start', 'end']]
    return regions_df


def is_in_regions(chrom, pos, regions_df):
    """判断某个位置是否在指定区域内"""
    # 查找相同染色体的区域
    chrom_regions = regions_df[regions_df['chrom'] == chrom]

    # 检查位置是否在任何一个区域内
    for _, region in chrom_regions.iterrows():
        if region['start'] <= pos <= region['end']:
            return True
    return False


def check_other_species_zero_or_na(row_data, current_species_indices, all_species_indices):
    """检查其他物种是否都是0或NA"""
    for other_species, other_indices in all_species_indices.items():
        if other_species == current_species_indices:
            continue  # 跳过当前物种

        # 获取其他物种的数据
        other_data = row_data[other_indices]

        # 检查每个个体是否为0或NA
        for value in other_data:
            if not np.isnan(value) and value != 0:
                return False

    return True


def calculate_allele_frequencies_by_region(position_info, genotype_data, sample_df, regions_df):
    """按物种和区域（内/外）分别计算等位基因频率，包括物种特异性统计和自身非零统计"""

    # 创建物种到个体索引的映射
    species_to_indices = defaultdict(list)
    for idx, (species, individual) in enumerate(sample_df.values):
        species_to_indices[species].append(idx)

    # 存储每个物种的等位基因频率计数
    # 原始统计（区域内外）
    species_freq_counts_in = {}  # 区域内
    species_freq_counts_out = {}  # 区域外

    # 物种特异性统计（只统计其他物种为0或NA的行）
    species_freq_counts_specific_in = {}  # 区域内
    species_freq_counts_specific_out = {}  # 区域外

    # 新增：物种自身非零统计（只统计当前物种基因计数>0的位点）
    species_freq_counts_nonzero_in = {}  # 区域内
    species_freq_counts_nonzero_out = {}  # 区域外

    # 对每个物种进行处理
    for species, indices in species_to_indices.items():
        n_individuals = len(indices)
        max_alleles = n_individuals * 2  # 二倍体，最大等位基因数

        # 初始化频率计数列表，长度增加1以包含0等位基因计数
        # 原始统计
        freq_counts_in = [0] * (max_alleles + 1)  # 0到max_alleles
        freq_counts_out = [0] * (max_alleles + 1)  # 0到max_alleles

        # 物种特异性统计
        freq_counts_specific_in = [0] * (max_alleles + 1)
        freq_counts_specific_out = [0] * (max_alleles + 1)

        # 新增：自身非零统计
        freq_counts_nonzero_in = [0] * (max_alleles + 1)
        freq_counts_nonzero_out = [0] * (max_alleles + 1)

        # 获取该物种对应的基因型数据
        species_genotypes = genotype_data.iloc[:, indices]

        # 对每一行（每个等位基因位点）进行处理
        for row_idx in range(len(species_genotypes)):
            row_data = species_genotypes.iloc[row_idx].values

            # 获取该位点的位置信息（使用起始位置作为位点位置）
            chrom = position_info.iloc[row_idx]['chrom']
            pos = position_info.iloc[row_idx]['start']

            # 判断是否在区域内
            in_region = is_in_regions(chrom, pos, regions_df)

            # 计算NA的数量
            na_count = np.sum(np.isnan(row_data))

            # 如果NA的数量超过个体数的一半，跳过这个位点
            if na_count > n_individuals / 2:
                continue  # 跳过NA数量超过一半的位点

            # 检查是否有NA值
            has_na = np.any(np.isnan(row_data))

            if not has_na:
                # 没有NA，直接计算等位基因总数
                allele_count = int(np.sum(row_data))
            else:
                # 有NA，需要进行缺失值处理
                valid_data = row_data[~np.isnan(row_data)]
                n_valid = len(valid_data)

                if n_valid > 0:
                    # 计算有效数据的和
                    sum_valid = np.sum(valid_data)
                    # 按比例调整到全部个体
                    adjusted_count = sum_valid * (n_individuals / n_valid)
                    # 四舍五入
                    allele_count = int(round(adjusted_count))
                else:
                    # 如果所有数据都是NA，跳过这一行
                    continue

            # 在对应的频率位置加1（注意：频率0对应索引0）
            if 0 <= allele_count <= max_alleles:
                # 原始统计
                if in_region:
                    freq_counts_in[allele_count] += 1
                else:
                    freq_counts_out[allele_count] += 1

                # 检查是否为物种特异性位点（其他物种都是0或NA）
                # 获取整行的所有数据
                entire_row_data = genotype_data.iloc[row_idx].values
                if check_other_species_zero_or_na(entire_row_data, species, species_to_indices):
                    # 物种特异性统计
                    if in_region:
                        freq_counts_specific_in[allele_count] += 1
                    else:
                        freq_counts_specific_out[allele_count] += 1

                # 新增：自身非零统计（只统计当前物种基因计数>0的位点）
                if allele_count > 0:
                    if in_region:
                        freq_counts_nonzero_in[allele_count] += 1
                    else:
                        freq_counts_nonzero_out[allele_count] += 1

        # 添加总和到列表末尾
        freq_counts_in.append(sum(freq_counts_in))
        freq_counts_out.append(sum(freq_counts_out))
        freq_counts_specific_in.append(sum(freq_counts_specific_in))
        freq_counts_specific_out.append(sum(freq_counts_specific_out))
        freq_counts_nonzero_in.append(sum(freq_counts_nonzero_in))
        freq_counts_nonzero_out.append(sum(freq_counts_nonzero_out))

        species_freq_counts_in[species] = freq_counts_in
        species_freq_counts_out[species] = freq_counts_out
        species_freq_counts_specific_in[species] = freq_counts_specific_in
        species_freq_counts_specific_out[species] = freq_counts_specific_out
        species_freq_counts_nonzero_in[species] = freq_counts_nonzero_in
        species_freq_counts_nonzero_out[species] = freq_counts_nonzero_out

    return (species_freq_counts_in, species_freq_counts_out,
            species_freq_counts_specific_in, species_freq_counts_specific_out,
            species_freq_counts_nonzero_in, species_freq_counts_nonzero_out)


def write_results_by_region(species_freq_counts_in, species_freq_counts_out,
                            species_freq_counts_specific_in, species_freq_counts_specific_out,
                            species_freq_counts_nonzero_in, species_freq_counts_nonzero_out,
                            output_filename):
    """将所有统计结果写入文件"""
    with open(output_filename, 'w') as f:
        # 写入原始统计结果
        f.write("=" * 80 + "\n")
        f.write("# 原始统计结果（所有符合条件的位点）\n")
        f.write("=" * 80 + "\n")

        f.write("\n## 区域内的等位基因频率统计\n")
        for species, freq_counts in species_freq_counts_in.items():
            freq_str = ','.join(map(str, freq_counts))
            f.write(f"{species}_in\t{freq_str}\n")

        f.write("\n## 区域外的等位基因频率统计\n")
        for species, freq_counts in species_freq_counts_out.items():
            freq_str = ','.join(map(str, freq_counts))
            f.write(f"{species}_out\t{freq_str}\n")

        # 写入物种特异性统计结果
        f.write("\n" + "=" * 80 + "\n")
        f.write("# 物种特异性统计结果（其他物种均为0或NA的位点）\n")
        f.write("=" * 80 + "\n")

        f.write("\n## 区域内的物种特异性等位基因频率统计\n")
        for species, freq_counts in species_freq_counts_specific_in.items():
            freq_str = ','.join(map(str, freq_counts))
            f.write(f"{species}_specific_in\t{freq_str}\n")

        f.write("\n## 区域外的物种特异性等位基因频率统计\n")
        for species, freq_counts in species_freq_counts_specific_out.items():
            freq_str = ','.join(map(str, freq_counts))
            f.write(f"{species}_specific_out\t{freq_str}\n")

        # 新增：写入自身非零统计结果
        f.write("\n" + "=" * 80 + "\n")
        f.write("# 物种自身非零统计结果（当前物种基因计数>0的位点）\n")
        f.write("=" * 80 + "\n")

        f.write("\n## 区域内的自身非零等位基因频率统计\n")
        for species, freq_counts in species_freq_counts_nonzero_in.items():
            freq_str = ','.join(map(str, freq_counts))
            f.write(f"{species}_nonzero_in\t{freq_str}\n")

        f.write("\n## 区域外的自身非零等位基因频率统计\n")
        for species, freq_counts in species_freq_counts_nonzero_out.items():
            freq_str = ','.join(map(str, freq_counts))
            f.write(f"{species}_nonzero_out\t{freq_str}\n")

    print(f"结果已保存到 {output_filename}")


def main():
    # 手动指定输入文件名
    bed_file = "allpops.refTra.all.syn.bed"  # 等位基因型计数文件
    sample_file = "sample_name.txt"  # 个体物种记录文件
    regions_file = "lines.txt"  # 区域坐标文件
    output_file = "allele_frequency_results_by_region_syn.txt"  # 输出文件

    print("开始处理文件...")

    try:
        # 读取文件
        print(f"读取基因型数据文件: {bed_file}")
        position_info, genotype_data = read_bed_file(bed_file)
        print(f"  - 读取到 {len(genotype_data)} 个等位基因位点")
        print(f"  - 包含 {len(genotype_data.columns)} 个个体")

        print(f"\n读取样本信息文件: {sample_file}")
        sample_df = read_sample_file(sample_file)
        print(f"  - 读取到 {len(sample_df)} 个个体信息")

        print(f"\n读取区域坐标文件: {regions_file}")
        regions_df = read_regions_file(regions_file)
        print(f"  - 读取到 {len(regions_df)} 个目标区域")

        # 统计有多少位点在区域内
        in_region_count = 0
        for row_idx in range(len(position_info)):
            chrom = position_info.iloc[row_idx]['chrom']
            pos = position_info.iloc[row_idx]['start']
            if is_in_regions(chrom, pos, regions_df):
                in_region_count += 1

        print(f"  - {in_region_count} 个位点在区域内")
        print(f"  - {len(position_info) - in_region_count} 个位点在区域外")

        # 检查个体数量是否匹配
        if len(genotype_data.columns) != len(sample_df):
            print(
                f"警告：基因型数据中的个体数({len(genotype_data.columns)})与样本文件中的个体数({len(sample_df)})不匹配")

        # 计算等位基因频率（包括原始统计、物种特异性统计和自身非零统计）
        print("\n计算等位基因频率（包括原始统计、物种特异性统计和自身非零统计）...")
        (species_freq_counts_in, species_freq_counts_out,
         species_freq_counts_specific_in, species_freq_counts_specific_out,
         species_freq_counts_nonzero_in, species_freq_counts_nonzero_out) = \
            calculate_allele_frequencies_by_region(
                position_info, genotype_data, sample_df, regions_df)

        # 输出统计信息
        print("\n统计结果:")
        print("-" * 50)
        print("原始统计：")
        for species in species_freq_counts_in.keys():
            n_individuals = len([s for s in sample_df['species'] if s == species])
            print(f"  物种 {species}: {n_individuals} 个个体, 最大等位基因数 {n_individuals * 2}")
            print(f"    - 区域内: {species_freq_counts_in[species][-1]} 个等位基因位点")
            print(f"    - 区域外: {species_freq_counts_out[species][-1]} 个等位基因位点")

        print("\n物种特异性统计（其他物种均为0或NA）：")
        for species in species_freq_counts_specific_in.keys():
            print(f"  物种 {species}:")
            print(f"    - 区域内: {species_freq_counts_specific_in[species][-1]} 个物种特异性位点")
            print(f"    - 区域外: {species_freq_counts_specific_out[species][-1]} 个物种特异性位点")

        print("\n物种自身非零统计（当前物种基因计数>0）：")
        for species in species_freq_counts_nonzero_in.keys():
            print(f"  物种 {species}:")
            print(f"    - 区域内: {species_freq_counts_nonzero_in[species][-1]} 个非零位点")
            print(f"    - 区域外: {species_freq_counts_nonzero_out[species][-1]} 个非零位点")

        # 写入结果
        print(f"\n保存结果到文件: {output_file}")
        write_results_by_region(species_freq_counts_in, species_freq_counts_out,
                                species_freq_counts_specific_in, species_freq_counts_specific_out,
                                species_freq_counts_nonzero_in, species_freq_counts_nonzero_out,
                                output_file)

        print("\n处理完成！")

    except FileNotFoundError as e:
        print(f"错误：找不到文件 - {e}")
    except Exception as e:
        print(f"错误：{e}")
        import traceback
        traceback.print_exc()


# 如果需要作为独立脚本运行
if __name__ == "__main__":
    main()