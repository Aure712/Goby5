#!/usr/bin/env Rscript

# -------- 用户可直接修改的参数 --------
lines_file <- "lines.txt"
fai_file <- "Tra_sort.fa.fai"
bed_file <- "tidk_default_tra_onlychr.bed"

max_chroms <- 11        # 只读取前多少条染色体，设为 NA 表示全部
threshold_bp <- 10e6     # 特定阈值，单位 bp
prefer_breakpoint <- TRUE  # 若同时命中断点区和端区，是否优先断点区
# ----------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# 1) 读取 fai：染色体名与长度
fai <- read.table(fai_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(fai) <- c("chr", "len", "offset", "linebases", "linewidth")

if (!is.na(max_chroms)) {
  max_chroms <- max(1, min(max_chroms, nrow(fai)))
  fai <- fai[seq_len(max_chroms), , drop = FALSE]
}
chr_len <- setNames(fai$len, fai$chr)

# 2) 读取断点（lines.txt 第二列）
lines <- read.table(lines_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
colnames(lines)[1:2] <- c("chr", "breakpoint")
lines <- lines[lines$chr %in% names(chr_len), , drop = FALSE]
bp_list <- split(lines$breakpoint, lines$chr)

# 3) 读取 BED 窗口
bed <- read.table(bed_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(bed) <- c("chr", "start", "end", "count")
bed <- bed[bed$chr %in% names(chr_len), , drop = FALSE]
if (nrow(bed) == 0) {
  stop("过滤后没有可用的 BED 窗口，请检查输入文件和 max_chroms 设置。")
}
bed$mid <- floor((bed$start + bed$end) / 2)

# 区域判定函数
get_region <- function(chr, mid) {
  len <- chr_len[[chr]]
  if (is.na(len)) return(NA_character_)
  
  in_end <- (mid <= threshold_bp) || (mid >= (len - threshold_bp))
  
  bps <- bp_list[[chr]]
  if (is.null(bps)) bps <- numeric(0)
  in_bp <- if (length(bps) > 0) any(abs(mid - bps) <= threshold_bp) else FALSE
  
  if (prefer_breakpoint && in_bp) return("breakpoint")
  if (in_end) return("chrom_end")
  if (in_bp) return("breakpoint")
  return("other")
}

bed$region <- mapply(get_region, bed$chr, bed$mid, USE.NAMES = FALSE)
bed$region <- factor(bed$region, levels = c("chrom_end", "breakpoint", "other"))

# -------- 区域坐标构建（0-based, end 为开区间）--------
make_intervals <- function(starts, ends, len) {
  if (length(starts) == 0) {
    return(data.frame(start = numeric(0), end = numeric(0)))
  }
  s <- pmax(0, starts)
  e <- pmin(len, ends)
  keep <- s < e
  data.frame(start = s[keep], end = e[keep])
}

merge_intervals <- function(df) {
  if (nrow(df) == 0) return(df)
  df <- df[order(df$start, df$end), , drop = FALSE]
  res <- list()
  cur_start <- df$start[1]
  cur_end <- df$end[1]
  if (nrow(df) > 1) {
    for (i in 2:nrow(df)) {
      s <- df$start[i]
      e <- df$end[i]
      if (s <= cur_end) {
        cur_end <- max(cur_end, e)
      } else {
        res[[length(res) + 1]] <- c(cur_start, cur_end)
        cur_start <- s
        cur_end <- e
      }
    }
  }
  res[[length(res) + 1]] <- c(cur_start, cur_end)
  res_df <- as.data.frame(do.call(rbind, res))
  colnames(res_df) <- c("start", "end")
  res_df
}

subtract_intervals <- function(base, sub) {
  if (nrow(base) == 0) return(base)
  if (nrow(sub) == 0) return(base)
  base <- merge_intervals(base)
  sub <- merge_intervals(sub)
  
  res <- list()
  for (i in seq_len(nrow(base))) {
    b_start <- base$start[i]
    b_end <- base$end[i]
    cur <- b_start
    
    for (j in seq_len(nrow(sub))) {
      s <- sub$start[j]
      e <- sub$end[j]
      if (e <= cur) next
      if (s >= b_end) break
      if (s > cur) {
        res[[length(res) + 1]] <- c(cur, min(s, b_end))
      }
      cur <- max(cur, e)
      if (cur >= b_end) break
    }
    
    if (cur < b_end) {
      res[[length(res) + 1]] <- c(cur, b_end)
    }
  }
  
  if (length(res) == 0) {
    return(data.frame(start = numeric(0), end = numeric(0)))
  }
  res_df <- as.data.frame(do.call(rbind, res))
  colnames(res_df) <- c("start", "end")
  res_df
}

region_coords_list <- list()

for (chr in names(chr_len)) {
  len <- chr_len[[chr]]
  full <- data.frame(start = 0, end = len)
  
  end_intervals <- make_intervals(
    c(0, len - threshold_bp),
    c(min(threshold_bp, len), len),
    len
  )
  end_intervals <- merge_intervals(end_intervals)
  
  bps <- bp_list[[chr]]
  if (is.null(bps)) bps <- numeric(0)
  bp_intervals <- make_intervals(bps - threshold_bp, bps + threshold_bp, len)
  bp_intervals <- merge_intervals(bp_intervals)
  
  if (prefer_breakpoint) {
    bp_regions <- bp_intervals
    end_regions <- subtract_intervals(end_intervals, bp_regions)
  } else {
    end_regions <- end_intervals
    bp_regions <- subtract_intervals(bp_intervals, end_regions)
  }
  
  occupied <- rbind(end_regions, bp_regions)
  if (nrow(occupied) > 0) {
    occupied <- merge_intervals(occupied)
  }
  other_regions <- subtract_intervals(full, occupied)
  
  if (nrow(end_regions) > 0) {
    region_coords_list[[length(region_coords_list) + 1]] <- data.frame(
      chr = chr, region = "chrom_end", start = end_regions$start, end = end_regions$end
    )
  }
  if (nrow(bp_regions) > 0) {
    region_coords_list[[length(region_coords_list) + 1]] <- data.frame(
      chr = chr, region = "breakpoint", start = bp_regions$start, end = bp_regions$end
    )
  }
  if (nrow(other_regions) > 0) {
    region_coords_list[[length(region_coords_list) + 1]] <- data.frame(
      chr = chr, region = "other", start = other_regions$start, end = other_regions$end
    )
  }
}

region_coords <- bind_rows(region_coords_list)
write.table(region_coords, "region_coordinates.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)

# 4) 统计不同区域内端粒计数分布
region_summary <- bed %>%
  group_by(region) %>%
  summarise(
    n_windows = n(),
    mean_count = mean(count),
    median_count = median(count),
    sd_count = sd(count),
    zero_windows = sum(count == 0),
    zero_pct = 100 * zero_windows / n_windows,
    .groups = "drop"
  )

write.table(region_summary, "region_summary.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)

# 5) 显著性检验 + p 值表格
stats_file <- "region_stats.txt"
pvals_file <- "region_pvalues.tsv"

valid_regions <- unique(na.omit(bed$region))
if (length(valid_regions) >= 2) {
  kw <- kruskal.test(count ~ region, data = bed)
  pw <- pairwise.wilcox.test(bed$count, bed$region, p.adjust.method = "BH")
  
  sink(stats_file)
  cat("Kruskal-Wallis test\n")
  print(kw)
  cat("\nPairwise Wilcoxon test (BH adjusted)\n")
  print(pw)
  sink()
  
  pvals <- data.frame(
    test = "Kruskal-Wallis",
    group1 = "all",
    group2 = "all",
    p_value = kw$p.value,
    p_adj_method = NA_character_,
    stringsAsFactors = FALSE
  )
  
  pw_mat <- pw$p.value
  if (!is.null(pw_mat)) {
    idx <- which(!is.na(pw_mat), arr.ind = TRUE)
    if (nrow(idx) > 0) {
      pw_df <- data.frame(
        test = "Wilcoxon",
        group1 = rownames(pw_mat)[idx[, 1]],
        group2 = colnames(pw_mat)[idx[, 2]],
        p_value = pw_mat[idx],
        p_adj_method = pw$p.adjust.method,
        stringsAsFactors = FALSE
      )
      pvals <- bind_rows(pvals, pw_df)
    }
  }
  
  write.table(pvals, pvals_file, sep = "\t",
              row.names = FALSE, quote = FALSE)
} else {
  writeLines("Not enough regions for statistical test.", stats_file)
  pvals <- data.frame(
    test = "Not enough regions",
    group1 = NA_character_,
    group2 = NA_character_,
    p_value = NA_real_,
    p_adj_method = NA_character_,
    stringsAsFactors = FALSE
  )
  write.table(pvals, pvals_file, sep = "\t",
              row.names = FALSE, quote = FALSE)
}

# 6) 作图
p <- ggplot(bed, aes(x = region, y = count)) +
  geom_violin(trim = FALSE, fill = "gray85", color = "gray40") +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black") +
  theme_classic() +
  labs(
    x = "Region",
    y = "Telomere count per window",
    title = "Telomere counts by region",
    subtitle = paste("threshold_bp =", threshold_bp)
  )

ggsave("telomere_counts_by_region.png", p, width = 7, height = 4, dpi = 300)

cat("Done. Outputs:\n",
    "  region_coordinates.tsv\n",
    "  region_summary.tsv\n",
    "  region_stats.txt\n",
    "  region_pvalues.tsv\n",
    "  telomere_counts_by_region.png\n")
