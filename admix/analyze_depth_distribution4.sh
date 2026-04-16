#!/bin/bash
# 文件名：depth_analysis_from_stats.sh

# 步骤1: 生成深度分布统计
bcftools stats -d 0,3000,1 "$1" > depth_stats.txt

# 检查是否提供了手动峰值参数
MANUAL_PEAKS=""
if [ -n "$2" ]; then
    MANUAL_PEAKS="$2"
    echo "Using manual peaks: $MANUAL_PEAKS"
fi

# 步骤2: 使用Python处理统计数据和绘图
python3 << 'EOF'
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import re
import os

# 读取手动峰值参数（如果提供）
manual_peaks = []
if 'MANUAL_PEAKS' in os.environ and os.environ['MANUAL_PEAKS']:
    try:
        manual_peaks = [float(p) for p in os.environ['MANUAL_PEAKS'].split(',')]
        print(f"Manual peaks provided: {manual_peaks}")
    except Exception as e:
        print(f"Error parsing manual peaks: {str(e)}")
        manual_peaks = []

# 读取并解析bcftools stats输出
print("Parsing depth distribution data...")
depth_data = []
depth_data_full = []  # 用于完整统计的数据
total_sites = 0
total_sites_full = 0  # 完整数据的总位点数
max_depth = 800  # 绘图的最大深度

with open("depth_stats.txt", "r") as f:
    for line in f:
        if line.startswith("DP\t0\t"):
            parts = line.strip().split('\t')
            depth_str = parts[2]  # 深度值字符串
            
            # 处理">3000"这种特殊情况
            if depth_str.startswith('>'):
                # 对于深度大于3000的汇总行，我们将其深度设为3001
                depth = 3001
            else:
                depth = int(depth_str)
            
            sites = int(parts[5])  # 该深度的位点数
            
            # 完整统计数据集（包含所有深度）
            depth_data_full.append((depth, sites))
            total_sites_full += sites
            
            # 绘图数据集（只包含0-800）
            if depth <= max_depth:
                depth_data.append((depth, sites))
                total_sites += sites

# 转换为DataFrame
df_full = pd.DataFrame(depth_data_full, columns=['Depth', 'Sites'])
df_plot = pd.DataFrame(depth_data, columns=['Depth', 'Sites'])

print(f"Total sites in full dataset: {total_sites_full:,}")
print(f"Total sites in 0-{max_depth} DP range: {total_sites:,}")

# 如果没有数据，则报错退出
if total_sites_full == 0:
    raise ValueError("No depth data found in the input")

# 使用完整数据集进行统计计算
depths_full = df_full['Depth'].values
frequencies_full = df_full['Sites'].values

# 处理手动峰值输入
if manual_peaks:
    print("Using manually specified peaks")
    # 按测序深度（X轴位置）排序
    sorted_manual_peaks = sorted(manual_peaks)
    
    if len(sorted_manual_peaks) >= 2:
        print(f"Peak 1 (leftmost): {sorted_manual_peaks[0]}")
        print(f"Peak 2: {sorted_manual_peaks[1]}")
        print("Using second peak as boundary")
        boundary = sorted_manual_peaks[1]
        peak_type = "Manual Secondary Peak"
    else:
        print(f"Only one peak provided: {sorted_manual_peaks[0]}")
        boundary = sorted_manual_peaks[0]
        peak_type = "Manual Primary Peak"
else:
    # 自动检测峰值 - 使用完整数据集
    print("Detecting peaks in depth distribution...")
    # 设置更宽松的参数以检测平缓峰
    min_height = 0.001 * np.max(frequencies_full)  # 非常低的阈值
    peaks, properties = find_peaks(
        frequencies_full, 
        height=min_height,  
        distance=5,  # 更小的峰间距
        prominence=0.01 * np.max(frequencies_full)  # 更低的显著度
    )
    
    print(f"Found {len(peaks)} potential peaks")
    
    # 如果没有检测到峰值，使用全局最大值
    if len(peaks) == 0:
        print("Warning: No peaks detected, using global maximum as boundary")
        main_peak_idx = np.argmax(frequencies_full)
        boundary = depths_full[main_peak_idx]
        peak_type = "Global Max"
        auto_peaks = [boundary]
    else:
        # 获取峰值对应的深度值
        auto_peaks_unsorted = depths_full[peaks]
        
        # 按测序深度（X轴位置）排序
        auto_peaks = sorted(auto_peaks_unsorted)
        print(f"Auto-detected peaks (by depth order): {auto_peaks}")
        
        # 处理不同峰数量的情况
        if len(auto_peaks) >= 2:
            print(f"Peak 1 (leftmost): {auto_peaks[0]}")
            print(f"Peak 2: {auto_peaks[1]}")
            print("Using second peak as boundary")
            boundary = auto_peaks[1]
            peak_type = "Auto Secondary Peak"
        else:
            print(f"Only one peak detected: {auto_peaks[0]}")
            boundary = auto_peaks[0]
            peak_type = "Auto Primary Peak"

print(f"Boundary position ({peak_type}): {boundary}")

# 计算边界右侧的95%分位数
print("Calculating 95th percentile for right-side data...")
# 找到边界在深度数组中的位置
boundary_idx = np.searchsorted(depths_full, boundary)

# 提取边界右侧的数据
right_depths = depths_full[boundary_idx:]
right_frequencies = frequencies_full[boundary_idx:]
right_total = np.sum(right_frequencies)

# 计算累积分布
cum_right = np.cumsum(right_frequencies)
percentile_95_pos = 0.95 * right_total

# 找到95%分位点的位置
percentile_idx = np.searchsorted(cum_right, percentile_95_pos)
percentile_95 = right_depths[percentile_idx]

right_percent = right_total / total_sites_full * 100

# 输出统计结果
print("\n===== Depth Distribution Statistics =====")
print(f"Boundary position ({peak_type}): {boundary}")
print(f"95th percentile of right-side data: {percentile_95}")
print(f"Right-side data points: {right_total:,} ({right_percent:.1f}% of total)")

with open("depth_summary.txt", "w") as f:
    f.write(f"Boundary position ({peak_type}): {boundary}\n")
    f.write(f"95th percentile of right-side data: {percentile_95}\n")
    f.write(f"Right-side data proportion: {right_percent:.1f}%\n")

# 可视化 - 只使用0-800范围的数据
print("Generating visualization...")
plt.figure(figsize=(12, 8))

# 提取绘图数据
depths_plot = df_plot['Depth'].values
frequencies_plot = df_plot['Sites'].values

# 创建主图和子图
ax1 = plt.subplot(2, 1, 1)
ax2 = plt.subplot(2, 1, 2)

# 在第一个子图上绘制原始深度分布
ax1.step(depths_plot, frequencies_plot / 1000, where='mid', 
         color="steelblue", label="Depth Distribution (thousands)")

# 标注关键点
if boundary <= max_depth:
    ax1.axvline(boundary, color="red", linestyle="--", 
                label=f"Boundary ({boundary}, {peak_type})")
else:
    ax1.axvline(max_depth, color="red", linestyle="--", 
                label=f"Boundary (>{max_depth}, actual: {boundary})")

if percentile_95 <= max_depth:
    ax1.axvline(percentile_95, color="purple", linestyle="-",
                label=f"95th Percentile ({percentile_95})")
else:
    ax1.axvline(max_depth, color="purple", linestyle="-",
                label=f"95th Percentile (>{max_depth}, actual: {percentile_95})")

# 标记检测到的所有峰值（按深度顺序）
if manual_peaks:
    # 使用排序后的手动峰值
    for i, peak in enumerate(sorted_manual_peaks):
        if peak <= max_depth:
            # 在原始数据中找到最接近的索引
            idx = np.abs(depths_plot - peak).argmin()
            peak_freq = frequencies_plot[idx]
            ax1.plot(peak, peak_freq / 1000, "ro", markersize=8)
            ax1.text(peak, peak_freq / 1000 + 5, f"Peak {i+1}: {peak}", 
                     fontsize=9, ha='center')
        else:
            ax1.plot(max_depth, 5, "ro", markersize=8)
            ax1.text(max_depth, 10, f"Peak {i+1}: {peak} (out of range)", 
                     fontsize=9, ha='right')
else:
    # 标记自动检测的峰值（已按深度排序）
    for i, peak in enumerate(auto_peaks):
        if peak <= max_depth:
            # 在原始数据中找到最接近的索引
            idx = np.abs(depths_plot - peak).argmin()
            peak_freq = frequencies_plot[idx]
            ax1.plot(peak, peak_freq / 1000, "ro", markersize=8)
            ax1.text(peak, peak_freq / 1000 + 5, f"Peak {i+1}: {peak}", 
                     fontsize=9, ha='center')
        else:
            ax1.plot(max_depth, 5, "ro", markersize=8)
            ax1.text(max_depth, 10, f"Peak {i+1}: {peak} (out of range)", 
                     fontsize=9, ha='right')

ax1.set_xlim(0, max_depth)
ax1.set_ylim(bottom=0)
ax1.set_xlabel("Sequencing Depth (DP)", fontsize=12)
ax1.set_ylabel("Number of Sites (thousands)", fontsize=12)
ax1.set_title(f"Depth Distribution (0-{max_depth} DP)", fontsize=14)
ax1.legend(fontsize=10, loc='upper right')
ax1.grid(alpha=0.3, linestyle="--")

# 在第二个子图上绘制累积分布图
ax2.plot(depths_plot, np.cumsum(frequencies_plot) / total_sites * 100, 
         color="green", linewidth=2, label="Cumulative Distribution")

# 标注边界和95%分位数
if boundary <= max_depth:
    ax2.axvline(boundary, color="red", linestyle="--")
    ax2.text(boundary, 50, f"Boundary: {boundary}", 
             fontsize=9, ha='right', backgroundcolor='white')

if percentile_95 <= max_depth:
    ax2.axvline(percentile_95, color="purple", linestyle="-")
    ax2.text(percentile_95, 30, f"95%: {percentile_95}", 
             fontsize=9, ha='left', backgroundcolor='white')

ax2.set_xlim(0, max_depth)
ax2.set_ylim(0, 100)
ax2.set_xlabel("Sequencing Depth (DP)", fontsize=12)
ax2.set_ylabel("Cumulative Percentage (%)", fontsize=12)
ax2.set_title("Cumulative Distribution", fontsize=14)
ax2.grid(alpha=0.3, linestyle="--")
ax2.legend(fontsize=10)

# 添加主标题
plt.suptitle(f"Sequencing Depth Distribution Analysis\n{total_sites_full:,} variants (plot shows 0-{max_depth} DP)", 
             fontsize=16, y=0.98)

# 调整布局
plt.tight_layout(rect=[0, 0, 1, 0.96])  # 为suptitle留出空间

# 保存并显示
plt.savefig("depth_distribution.png", dpi=300)
plt.savefig("depth_distribution.pdf")
print("Analysis complete! Results saved to depth_summary.txt and depth_distribution.png")

# 保存峰值检测报告
with open("peak_detection_report.txt", "w") as f:
    f.write("Peak Detection Report\n")
    f.write("=====================\n")
    f.write(f"Total sites: {total_sites_full:,}\n")
    f.write(f"Detection method: {'Manual' if manual_peaks else 'Auto'}\n")
    
    if manual_peaks:
        f.write(f"Manual peaks (by depth order): {', '.join(map(str, sorted_manual_peaks))}\n")
    else:
        f.write(f"Auto-detected peaks (by depth order): {', '.join(map(str, auto_peaks))}\n")
    
    f.write(f"Boundary position: {boundary} ({peak_type})\n")
    f.write(f"95th percentile: {percentile_95}\n")
    f.write(f"Right-side proportion: {right_percent:.1f}%\n")
EOF
