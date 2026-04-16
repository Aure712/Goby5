#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ===== 输入文件（直接在这里指定）=====
normal_file <- "all_tra_subinvert.bed"
invert_file <- "syn_invert.bed"
dxy_file    <- "Tba.Tra.dxy.bed"  # 实际文件名

# ===== 读取数据 =====
normal <- fread(normal_file, header = FALSE, col.names = c("chr", "start", "end"))
invert <- fread(invert_file, header = FALSE, col.names = c("chr", "start", "end"))

dxy <- fread(dxy_file, header = FALSE)
dxy <- dxy[, .(chr = V1, start = V2, end = V3, dxy = V4)]

# ===== 计算中点并判定归属区域 =====
dxy[, mid := (start + end) / 2]
points <- dxy[, .(chr, start = mid, end = mid, dxy)]

normal[, type := "Normal"]
invert[, type := "Inversion"]
regions <- rbindlist(list(normal, invert), use.names = TRUE)

setkey(regions, chr, start, end)
setkey(points, chr, start, end)

# 中点落在区域内的记录
over <- foverlaps(points, regions, nomatch = 0L)
over <- over[, .(dxy, type)]

if (length(unique(over$type)) < 2) {
  stop("有效数据不足，无法比较两类区域。")
}

# ===== 通用绘图+统计函数 =====
analyze_and_plot <- function(dat, tag) {
  # 统计检验
  test_res <- wilcox.test(dxy ~ type, data = dat, exact = FALSE)
  pval <- test_res$p.value
  
  # 统计汇总表
  stats <- dat[, .(
    n = .N,
    mean = mean(dxy, na.rm = TRUE),
    median = median(dxy, na.rm = TRUE),
    sd = sd(dxy, na.rm = TRUE),
    iqr = IQR(dxy, na.rm = TRUE),
    min = min(dxy, na.rm = TRUE),
    max = max(dxy, na.rm = TRUE)
  ), by = type]
  stats[, p_value_wilcox := pval]
  
  stats_file <- sprintf("dxy_normal_vs_invert_stats_%s.tsv", tag)
  fwrite(stats, stats_file, sep = "\t")
  
  # 画图（小提琴图 + 盒须图 + 抖动散点）
  ymax <- max(dat$dxy, na.rm = TRUE)
  ymin <- min(dat$dxy, na.rm = TRUE)
  ypos <- ymax + 0.08 * (ymax - ymin)
  
  p_label <- if (pval < 0.001) "p < 0.001" else sprintf("p = %.3g", pval)
  
  g <- ggplot(dat, aes(x = type, y = dxy, fill = type)) +
    geom_violin(trim = FALSE, color = "#2b2b2b", linewidth = 0.6, alpha = 0.75) +
    geom_boxplot(width = 0.12, outlier.shape = NA, color = "#2b2b2b",
                 linewidth = 0.5, fill = "white") +
    geom_jitter(width = 0.12, alpha = 0.4, size = 0.6, color = "#1f1f1f") +
    annotate("text", x = 1.5, y = ypos, label = p_label, size = 4.2) +
    scale_fill_manual(values = c("Normal" = "#4DBBD5", "Inversion" = "#E64B35")) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12),
      axis.text = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
      plot.margin = margin(10, 12, 10, 12)
    ) +
    ylab("dxy") +
    ggtitle(sprintf("Dxy distribution in Normal vs Inversion (%s)", tag))
  
  plot_file <- sprintf("dxy_normal_vs_invert_violin_%s.png", tag)
  ggsave(plot_file, g, width = 6.5, height = 4.8, dpi = 300)
  
  list(pval = pval, stats_file = stats_file, plot_file = plot_file)
}

# ===== 原始数据分析 =====
res_raw <- analyze_and_plot(over, "raw")

# ===== 去除每个区域内前后2.5%极端值后分析 =====
over_trim <- over[, {
  q <- quantile(dxy, probs = c(0.025, 0.975), na.rm = TRUE)
  .SD[dxy >= q[1] & dxy <= q[2]]
}, by = type]

res_trim <- analyze_and_plot(over_trim, "trim_2.5pct")

# ===== 输出结果 =====
cat("Raw Wilcoxon p-value:", res_raw$pval, "\n")
cat("Raw plot:", res_raw$plot_file, "\n")
cat("Raw stats:", res_raw$stats_file, "\n\n")

cat("Trimmed Wilcoxon p-value:", res_trim$pval, "\n")
cat("Trimmed plot:", res_trim$plot_file, "\n")
cat("Trimmed stats:", res_trim$stats_file, "\n")
