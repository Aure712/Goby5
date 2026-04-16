# === R代码文件: plot_abs_stats_scan_window_with_tables.R ===

library(dplyr)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(data.table)

# ================== 核心参数设置 ==================
plot_window <- 5e5  # 绘图窗口大小 (0.5MB) - 用于高精度图像绘制
scan_window <- 2e6  # 扫描窗口大小 (2MB) - 用于p值和概率的统计检验
scan_step <- 2e6    # 扫描步长 (2MB)
x1 <- 0.007300749
x2 <- 0.005889914

summary_extension <- 5e7  # 断点两侧各延伸50Mb

# 差异分析阈值参数 (20Mb)
threshold_dist <- 20e6    

# 平滑度控制参数
smoothing_params <- list(
  enabled = FALSE,        
  method = "none",        
  span = 0.1,             
  degree = 1,             
  spar = 0.5,             
  ma_window = 3           
)

# 发表级别的颜色主题
publication_colors <- list(
  tba = "#2E86AB",        
  tra = "#A23B72",        
  marker = "#C73E1D",     
  scan = "#6A994E",       
  threshold = "#F5A623",  # 20Mb 阈值分割线颜色 (橙色)
  background = "#FAFAFA", 
  text = "#2D3748"        
)

# 自定义主题函数
create_publication_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.background = element_rect(fill = publication_colors$background, color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "white", size = 0.8),
      panel.grid.minor = element_line(color = "white", size = 0.4),
      axis.line = element_line(color = publication_colors$text, size = 0.6),
      axis.ticks = element_line(color = publication_colors$text, size = 0.5),
      axis.text = element_text(color = publication_colors$text, size = 10),
      axis.title = element_text(color = publication_colors$text, size = 12, face = "bold"),
      plot.title = element_text(color = publication_colors$text, size = 16, face = "bold", hjust = 0.5, margin = margin(b = 20)),
      plot.subtitle = element_text(color = publication_colors$text, size = 12, hjust = 0.5, margin = margin(b = 15)),
      legend.position = "bottom",
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16, face = "bold"),
      legend.margin = margin(t = 15),
      plot.margin = margin(20, 20, 20, 20)
    )
}

# 优化的读取函数
read_bed_efficient <- function(filepath, source_name) {
  first_line <- readLines(filepath, n = 1)
  n_cols <- length(strsplit(first_line, "\t")[[1]])
  
  if (n_cols == 4) {
    col_classes <- c("character", "numeric", "numeric", "numeric")
    dt <- fread(filepath, sep = "\t", header = FALSE, colClasses = col_classes, select = 1:4)
    setnames(dt, c("V1", "V2", "V3", "V4"))
    dt[, `:=`(V5 = V1, V6 = V2, V7 = (V2 + V3) / 2, V8 = V4)]
    dt[, c("V1", "V2", "V3", "V4") := NULL]
  } else {
    dt <- fread(filepath, sep = "\t", header = FALSE, select = c(5, 7, 8))
    setnames(dt, c("V5", "V7", "V8"))
  }
  dt[, source := source_name]
  return(dt)
}

get_unique_chromosomes <- function() {
  chr_tba <- fread("../../tba.bed", sep = "\t", header = FALSE, select = 5)[[1]]
  chr_tra <- fread("../../tra.bed", sep = "\t", header = FALSE, select = 5)[[1]]
  all_chr <- unique(c(chr_tba, chr_tra))
  rm(chr_tba, chr_tra); gc()
  return(all_chr)
}

create_plot_data <- function(data, window_size) {
  if (nrow(data) == 0) return(data.frame())
  dt <- as.data.table(data)
  result <- dt[, .(
    pos = mean(floor((V7 - min(V7)) / window_size) * window_size + min(V7)) + window_size/2,
    mean_y = mean(adjusted_y, na.rm = TRUE),
    se_y = sd(adjusted_y, na.rm = TRUE) / sqrt(sum(!is.na(adjusted_y))),
    n_points = .N
  ), by = .(chr = V5, window = floor((V7 - min(V7)) / window_size) * window_size + min(V7))]
  result <- result[!is.na(mean_y) & !is.na(pos)]
  setorder(result, pos)
  return(as.data.frame(result))
}

apply_smoothing <- function(data, params = smoothing_params) {
  if (nrow(data) < 4 || !params$enabled || params$method == "none") {
    data$smooth_y <- data$mean_y
    return(data)
  }
  tryCatch({
    if (params$method == "loess") {
      smooth_fit <- loess(mean_y ~ pos, data = data, span = params$span, degree = params$degree)
      data$smooth_y <- predict(smooth_fit)
    } else if (params$method == "spline") {
      smooth_fit <- smooth.spline(data$pos, data$mean_y, spar = params$spar)
      data$smooth_y <- predict(smooth_fit, data$pos)$y
    } else if (params$method == "moving_average") {
      n <- nrow(data); window <- params$ma_window
      data$smooth_y <- data$mean_y
      for (i in 1:n) {
        start_idx <- max(1, i - floor(window/2))
        end_idx <- min(n, i + floor(window/2))
        data$smooth_y[i] <- mean(data$mean_y[start_idx:end_idx], na.rm = TRUE)
      }
    } else {
      data$smooth_y <- data$mean_y
    }
    data$smooth_y[is.na(data$smooth_y)] <- data$mean_y[is.na(data$smooth_y)]
    return(data)
  }, error = function(e) {
    data$smooth_y <- data$mean_y
    return(data)
  })
}

# 竖线坐标定义
vertical_lines <- list(
  "1" = c(40492542, 80104410), "2" = 45654190, "3" = 53945679,
  "4" = 44435307, "5" = 42365773, "6" = 38964430, "7" = 39901015,
  "8" = 36517876, "9" = 41362729, "10" = 33661800
)

# 助手函数：按绝对距离聚合
aggregate_by_abs_bin <- function(df, window_size) {
  if (is.null(df) || nrow(df) == 0) return(data.frame())
  dt <- as.data.table(df)[!is.na(relative_pos) & !is.na(smooth_y)]
  dt[, abs_pos := abs(relative_pos)]
  dt[, abs_bin := round(abs_pos / window_size) * window_size]
  res <- dt[, .(
    pos = abs_bin[1],  
    mean_y = mean(smooth_y, na.rm = TRUE),
    se_y = { n <- .N; if (n <= 1) 0 else sd(smooth_y, na.rm = TRUE) / sqrt(n) },
    n_breakpoints = .N
  ), by = abs_bin]
  setorder(res, pos)
  as.data.frame(res)
}

# 助手函数：计算绝对距离短板范围
compute_common_ranges <- function(bp_list) {
  if (length(bp_list) == 0) return(NULL)
  ranges <- lapply(bp_list, function(bp) {
    rel <- c()
    if (!is.null(bp$tba) && nrow(bp$tba) > 0) rel <- c(rel, bp$tba$relative_pos)
    if (!is.null(bp$tra) && nrow(bp$tra) > 0) rel <- c(rel, bp$tra$relative_pos)
    rel <- rel[is.finite(rel)]
    if (length(rel) == 0) return(NULL)
    list(max_abs = max(abs(rel)))
  })
  ranges <- ranges[!sapply(ranges, is.null)]
  if (length(ranges) == 0) return(NULL)
  list(max_abs = min(sapply(ranges, `[[`, "max_abs")))
}

# ================== 新增：统计检验模块及表格输出 ==================
calculate_stats <- function(df_tba, df_tra, threshold, label, output_filename) {
  # 根据绝对坐标 pos 匹配对应窗口的值
  df_merged <- merge(df_tba[, c("pos", "mean_y")], 
                     df_tra[, c("pos", "mean_y")], 
                     by="pos", suffixes=c("_TBA", "_TRA"))
  
  df_in <- df_merged[df_merged$pos <= threshold, ]
  df_out <- df_merged[df_merged$pos > threshold, ]
  
  cat(sprintf("\n================= %s 统计分析 =================\n", label))
  cat(sprintf("分析设定：阈值距离为 %s Mb\n", threshold/1e6))
  
  p_val_in <- NA
  p_val_out <- NA
  
  # 阈值内的分析
  if (nrow(df_in) > 0) {
    prob_tba_higher <- mean(df_in$mean_y_TBA > df_in$mean_y_TRA, na.rm=TRUE)
    # 使用成对 Wilcoxon 符号秩检验计算 p-value (备择假设：TBA > TRA)
    p_val_in <- tryCatch({
      wilcox.test(df_in$mean_y_TBA, df_in$mean_y_TRA, paired=TRUE, alternative="greater")$p.value
    }, error = function(e) NA)
    
    cat(sprintf("\n[ 阈值内: 0 到 %sMb ] (共基于 %d 个扫描窗口进行配对检验)\n", threshold/1e6, nrow(df_in)))
    cat(sprintf("  ▶ TBA 高于 TRA 的概率: %.2f%%\n", prob_tba_higher * 100))
    cat(sprintf("  ▶ TBA 显著高于 TRA (Wilcoxon检验 p-value): %s\n", format(p_val_in, scientific=TRUE, digits=4)))
  } else {
    cat(sprintf("\n[ 阈值内: 0 到 %sMb ] 无数据！\n", threshold/1e6))
  }
  
  # 阈值外的分析
  if (nrow(df_out) > 0) {
    prob_tra_higher <- mean(df_out$mean_y_TRA > df_out$mean_y_TBA, na.rm=TRUE)
    # 使用成对 Wilcoxon 符号秩检验计算 p-value (备择假设：TRA > TBA)
    p_val_out <- tryCatch({
      wilcox.test(df_out$mean_y_TRA, df_out$mean_y_TBA, paired=TRUE, alternative="greater")$p.value
    }, error = function(e) NA)
    
    cat(sprintf("\n[ 阈值外: 大于 %sMb ] (共基于 %d 个扫描窗口进行配对检验)\n", threshold/1e6, nrow(df_out)))
    cat(sprintf("  ▶ TRA 高于 TBA 的概率: %.2f%%\n", prob_tra_higher * 100))
    cat(sprintf("  ▶ TRA 显著高于 TBA (Wilcoxon检验 p-value): %s\n", format(p_val_out, scientific=TRUE, digits=4)))
  } else {
    cat(sprintf("\n[ 阈值外: 大于 %sMb ] 无数据！\n", threshold/1e6))
  }
  cat("=================================================================\n")
  
  # ---------------- 导出统计数据表格 ----------------
  # 补充详细的统计信息列
  df_merged$Region <- ifelse(df_merged$pos <= threshold, 
                             sprintf("<= %sMb", threshold/1e6), 
                             sprintf("> %sMb", threshold/1e6))
  df_merged$Tested_Hypothesis <- ifelse(df_merged$pos <= threshold, "TBA > TRA", "TRA > TBA")
  df_merged$P_value <- ifelse(df_merged$pos <= threshold, p_val_in, p_val_out)
  
  # 重新排列一下列顺序使其更易读
  df_merged <- df_merged[, c("pos", "mean_y_TBA", "mean_y_TRA", "Region", "Tested_Hypothesis", "P_value")]
  
  # 写入 CSV 文件
  write.csv(df_merged, output_filename, row.names = FALSE, quote = FALSE)
  cat(sprintf("=> 已成功导出配套统计表格: %s\n", output_filename))
}


# ================== 数据提取与处理主循环 ==================
all_chromosomes <- get_unique_chromosomes()
all_breakpoint_data <- list()

cat("开始读取数据并提取断点区域 (跳过单染色体绘图)...\n")
for (chr in all_chromosomes) {
  chr_fullname <- chr
  chr_num <- gsub("TraScf_", "", chr_fullname)
  
  # 仅读取必需数据以加速
  tba_data <- tryCatch({
    dt <- read_bed_efficient("../../tba.bed", "TBA")
    dt <- dt[V5 == chr_fullname]
    dt[, adjusted_y := ifelse(V8 >= 0, V8 / x1, NA_real_)]
    dt[!is.na(adjusted_y) & !is.na(V7) & is.finite(adjusted_y)]
  }, error = function(e) data.table())
  
  tra_data <- tryCatch({
    dt <- read_bed_efficient("../../tra.bed", "TRA")
    dt <- dt[V5 == chr_fullname]
    dt[, adjusted_y := ifelse(V8 >= 0, V8 / x2, NA_real_)]
    dt[!is.na(adjusted_y) & !is.na(V7) & is.finite(adjusted_y)]
  }, error = function(e) data.table())
  
  chr_data <- rbindlist(list(tba_data, tra_data), use.names = TRUE, fill = TRUE)
  if (nrow(chr_data) == 0) next
  
  # 使用 plot_window(0.5MB) 生成基础序列
  plot_data_tba <- if(nrow(tba_data) > 0) apply_smoothing(create_plot_data(as.data.frame(tba_data), plot_window)) else data.frame()
  plot_data_tra <- if(nrow(tra_data) > 0) apply_smoothing(create_plot_data(as.data.frame(tra_data), plot_window)) else data.frame()
  
  if (chr_num %in% names(vertical_lines)) {
    marker_positions <- vertical_lines[[chr_num]]
    min_pos <- min(chr_data$V7, na.rm = TRUE)
    max_pos <- max(chr_data$V7, na.rm = TRUE)
    
    for (bp_idx in seq_along(marker_positions)) {
      bp_pos <- marker_positions[bp_idx]
      left_bound <- max(min_pos, bp_pos - summary_extension)
      right_bound <- min(max_pos, bp_pos + summary_extension)
      
      if (length(marker_positions) > 1) {
        if (bp_idx > 1) left_bound <- max(left_bound, (marker_positions[bp_idx-1] + bp_pos) / 2)
        if (bp_idx < length(marker_positions)) right_bound <- min(right_bound, (bp_pos + marker_positions[bp_idx+1]) / 2)
      }
      
      region_data_tba <- plot_data_tba %>% filter(pos >= left_bound & pos <= right_bound)
      region_data_tra <- plot_data_tra %>% filter(pos >= left_bound & pos <= right_bound)
      
      if (nrow(region_data_tba) > 0 || nrow(region_data_tra) > 0) {
        if (nrow(region_data_tba) > 0) { region_data_tba$relative_pos <- region_data_tba$pos - bp_pos }
        if (nrow(region_data_tra) > 0) { region_data_tra$relative_pos <- region_data_tra$pos - bp_pos }
        
        all_breakpoint_data[[length(all_breakpoint_data) + 1]] <- list(
          tba = region_data_tba, tra = region_data_tra
        )
      }
    }
  }
  rm(chr_data, plot_data_tba, plot_data_tra); gc()
}

# ==================== 开始绝对距离聚合分析与绘图 ====================
cat("\n提取完成！开始进行绝对距离相关的绘图与统计检验...\n")

if (length(all_breakpoint_data) > 0) {
  all_tba_summary <- do.call(rbind, lapply(all_breakpoint_data, function(x) x$tba))
  all_tra_summary <- do.call(rbind, lapply(all_breakpoint_data, function(x) x$tra))
  
  # -------------------- 1. 全局绝对距离聚合 (Aggregated_Abs) --------------------
  
  # 用于高精度绘图的聚合 (0.5MB 窗口)
  agg_tba_abs_plot <- aggregate_by_abs_bin(all_tba_summary, plot_window)
  agg_tra_abs_plot <- aggregate_by_abs_bin(all_tra_summary, plot_window)
  
  # 用于统计检验的聚合 (2MB 扫描窗口)
  stats_tba_abs <- aggregate_by_abs_bin(all_tba_summary, scan_window)
  stats_tra_abs <- aggregate_by_abs_bin(all_tra_summary, scan_window)
  
  if ((nrow(stats_tba_abs) + nrow(stats_tra_abs)) > 0) {
    # >>> 调用统计模块，传入 2MB 聚合的数据，并输出 CSV <<<
    calculate_stats(stats_tba_abs, stats_tra_abs, threshold_dist, 
                    "全局绝对距离汇总 (Aggregated_Abs)", 
                    "Stats_All_Breakpoints_Summary_Aggregated_Abs.csv")
    
    # 绘图逻辑 (修复图注问题：使用 aes(color=) 和 aes(fill=) 正确映射)
    y_min_abs <- min(agg_tba_abs_plot$mean_y - agg_tba_abs_plot$se_y, agg_tra_abs_plot$mean_y - agg_tra_abs_plot$se_y, na.rm=TRUE)
    y_max_abs <- max(agg_tba_abs_plot$mean_y + agg_tba_abs_plot$se_y, agg_tra_abs_plot$mean_y + agg_tra_abs_plot$se_y, na.rm=TRUE)
    y_lim_abs <- c(max(0, y_min_abs - (y_max_abs-y_min_abs)*0.1), y_max_abs + (y_max_abs-y_min_abs)*0.1)
    
    p_abs <- ggplot() +
      create_publication_theme() +
      labs(title = "Summary Aggregated by Absolute Distance", subtitle = "With Statistical Threshold Highlight", x = "Distance to Breakpoint (Mb)", y = "RHO Value") +
      coord_cartesian(ylim = y_lim_abs)
    
    if(nrow(agg_tba_abs_plot)>0) {
      p_abs <- p_abs + 
        geom_ribbon(data = agg_tba_abs_plot, aes(x=pos/1e6, ymin=mean_y-se_y, ymax=mean_y+se_y, fill="TBA"), alpha=0.12) +
        geom_line(data = agg_tba_abs_plot, aes(x=pos/1e6, y=mean_y, color="TBA"), size=0.8)
    }
    if(nrow(agg_tra_abs_plot)>0) {
      p_abs <- p_abs + 
        geom_ribbon(data = agg_tra_abs_plot, aes(x=pos/1e6, ymin=mean_y-se_y, ymax=mean_y+se_y, fill="TRA"), alpha=0.12) +
        geom_line(data = agg_tra_abs_plot, aes(x=pos/1e6, y=mean_y, color="TRA"), size=0.8)
    }
    
    # 增加断点 0 处标记、20Mb 处阈值标记，并定义正确的颜色图注
    p_abs <- p_abs +
      geom_vline(xintercept = 0, color = publication_colors$marker, linetype = "dashed", size = 1.0) +
      geom_vline(xintercept = threshold_dist/1e6, color = publication_colors$threshold, linetype = "longdash", size = 0.8) +
      annotate("text", x = threshold_dist/1e6, y = Inf, label = paste0(threshold_dist/1e6, "Mb Threshold"), 
               color = publication_colors$threshold, vjust = 1.5, hjust = -0.1, size = 4, fontface = "bold") +
      scale_y_continuous(labels = scales::scientific_format(digits = 2)) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      scale_color_manual(name="Sample", values=c("TBA"=publication_colors$tba, "TRA"=publication_colors$tra)) +
      scale_fill_manual(name="Sample", values=c("TBA"=publication_colors$tba, "TRA"=publication_colors$tra))
    
    ggsave("All_Breakpoints_Summary_Aggregated_Abs.png", plot=p_abs, width=14, height=8, dpi=600, bg="white")
    ggsave("All_Breakpoints_Summary_Aggregated_Abs.pdf", plot=p_abs, width=14, height=8, device="pdf", bg="white")
  }
  
  # -------------------- 2. 公共边界绝对距离聚合 (Aggregated_Abs_Common) --------------------
  common_ranges <- compute_common_ranges(all_breakpoint_data)
  if (!is.null(common_ranges) && is.finite(common_ranges$max_abs) && common_ranges$max_abs > 0) {
    all_tba_abs_common <- all_tba_summary[abs(all_tba_summary$relative_pos) <= common_ranges$max_abs, ]
    all_tra_abs_common <- all_tra_summary[abs(all_tra_summary$relative_pos) <= common_ranges$max_abs, ]
    
    # 用于高精度绘图的聚合 (0.5MB 窗口)
    agg_tba_abs_common_plot <- aggregate_by_abs_bin(all_tba_abs_common, plot_window)
    agg_tra_abs_common_plot <- aggregate_by_abs_bin(all_tra_abs_common, plot_window)
    
    # 用于统计检验的聚合 (2MB 扫描窗口)
    stats_tba_abs_common <- aggregate_by_abs_bin(all_tba_abs_common, scan_window)
    stats_tra_abs_common <- aggregate_by_abs_bin(all_tra_abs_common, scan_window)
    
    if ((nrow(stats_tba_abs_common) + nrow(stats_tra_abs_common)) > 0) {
      # >>> 调用统计模块，传入 2MB 聚合的数据，并输出 CSV <<<
      calculate_stats(stats_tba_abs_common, stats_tra_abs_common, threshold_dist, 
                      "公共边界绝对距离汇总 (Aggregated_Abs_Common)",
                      "Stats_All_Breakpoints_Summary_Aggregated_Abs_Common.csv")
      
      # 绘图逻辑
      y_min_comm <- min(agg_tba_abs_common_plot$mean_y - agg_tba_abs_common_plot$se_y, agg_tra_abs_common_plot$mean_y - agg_tra_abs_common_plot$se_y, na.rm=TRUE)
      y_max_comm <- max(agg_tba_abs_common_plot$mean_y + agg_tba_abs_common_plot$se_y, agg_tra_abs_common_plot$mean_y + agg_tra_abs_common_plot$se_y, na.rm=TRUE)
      y_lim_comm <- c(max(0, y_min_comm - (y_max_comm-y_min_comm)*0.1), y_max_comm + (y_max_comm-y_min_comm)*0.1)
      
      p_abs_common <- ggplot() +
        create_publication_theme() +
        labs(title = "Summary Aggregated by Absolute Distance (Common Range)", subtitle = "With Statistical Threshold Highlight", x = "Distance to Breakpoint (Mb)", y = "RHO Value") +
        coord_cartesian(ylim = y_lim_comm)
      
      if(nrow(agg_tba_abs_common_plot)>0) {
        p_abs_common <- p_abs_common + 
          geom_ribbon(data = agg_tba_abs_common_plot, aes(x=pos/1e6, ymin=mean_y-se_y, ymax=mean_y+se_y, fill="TBA"), alpha=0.12) + 
          geom_line(data = agg_tba_abs_common_plot, aes(x=pos/1e6, y=mean_y, color="TBA"), size=0.8)
      }
      if(nrow(agg_tra_abs_common_plot)>0) {
        p_abs_common <- p_abs_common + 
          geom_ribbon(data = agg_tra_abs_common_plot, aes(x=pos/1e6, ymin=mean_y-se_y, ymax=mean_y+se_y, fill="TRA"), alpha=0.12) + 
          geom_line(data = agg_tra_abs_common_plot, aes(x=pos/1e6, y=mean_y, color="TRA"), size=0.8)
      }
      
      p_abs_common <- p_abs_common +
        geom_vline(xintercept = 0, color = publication_colors$marker, linetype = "dashed", size = 1.0) +
        geom_vline(xintercept = threshold_dist/1e6, color = publication_colors$threshold, linetype = "longdash", size = 0.8) +
        annotate("text", x = threshold_dist/1e6, y = Inf, label = paste0(threshold_dist/1e6, "Mb Threshold"), 
                 color = publication_colors$threshold, vjust = 1.5, hjust = -0.1, size = 4, fontface = "bold") +
        scale_y_continuous(labels = scales::scientific_format(digits = 2)) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
        scale_color_manual(name="Sample", values=c("TBA"=publication_colors$tba, "TRA"=publication_colors$tra)) +
        scale_fill_manual(name="Sample", values=c("TBA"=publication_colors$tba, "TRA"=publication_colors$tra))
      
      ggsave("All_Breakpoints_Summary_Aggregated_Abs_Common.png", plot=p_abs_common, width=14, height=8, dpi=600, bg="white")
      ggsave("All_Breakpoints_Summary_Aggregated_Abs_Common.pdf", plot=p_abs_common, width=14, height=8, device="pdf", bg="white")
    }
  }
} else {
  cat("未检测到有效断点数据。\n")
}

cat("\n所有任务执行完毕。\n")
