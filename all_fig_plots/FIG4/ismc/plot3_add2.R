library(dplyr)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(data.table)  # 使用data.table处理大数据更高效

# 参数设置
plot_window <- 5e5  # 绘图窗口大小 (0.5MB)
scan_window <- 2e6  # 扫描窗口大小 (2MB)
scan_step <- 2e6    # 扫描步长 (2MB)
x1 <- 0.007300749
x2 <- 0.005889914

# 新增:汇总图扩展范围参数
summary_extension <- 5e7  # 断点两侧各延伸50Mb

# 新增:平滑度控制参数
smoothing_params <- list(
  enabled = FALSE,        # 是否启用平滑
  method = "none",        # 平滑方法: "loess", "spline", "moving_average", "none"
  span = 0.1,             # loess平滑度 (0.1-1.0,值越小越不平滑)
  degree = 1,             # loess多项式次数 (1或2)
  spar = 0.5,             # spline平滑参数 (0-1)
  ma_window = 3           # 移动平均窗口大小
)

# 定义发表级别的颜色主题
publication_colors <- list(
  tba = "#2E86AB",        # 专业蓝色
  tra = "#A23B72",        # 深玫瑰色
  marker = "#C73E1D",     # 深红色标记
  scan = "#6A994E",       # 深绿色扫描线
  background = "#FAFAFA", # 浅灰背景
  text = "#2D3748"        # 深灰文本
)

# 自定义主题函数 - 修改图例设置
create_publication_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      # 背景设置
      panel.background = element_rect(fill = publication_colors$background, color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      
      # 网格线设置
      panel.grid.major = element_line(color = "white", size = 0.8),
      panel.grid.minor = element_line(color = "white", size = 0.4),
      
      # 轴设置
      axis.line = element_line(color = publication_colors$text, size = 0.6),
      axis.ticks = element_line(color = publication_colors$text, size = 0.5),
      axis.text = element_text(color = publication_colors$text, size = 10),
      axis.title = element_text(color = publication_colors$text, size = 12, face = "bold"),
      
      # 标题设置
      plot.title = element_text(
        color = publication_colors$text, 
        size = 16, 
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 20)
      ),
      plot.subtitle = element_text(
        color = publication_colors$text, 
        size = 12,
        hjust = 0.5,
        margin = margin(b = 15)
      ),
      
      # 图例设置
      legend.position = "bottom",
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18, face = "bold"),
      legend.margin = margin(t = 15),
      legend.key.size = unit(1.5, "lines"),
      legend.spacing.x = unit(1.0, "cm"),
      
      # 边距设置
      plot.margin = margin(20, 20, 20, 20),
      
      # 刻面设置
      strip.background = element_rect(fill = publication_colors$background, color = publication_colors$text),
      strip.text = element_text(color = publication_colors$text, size = 11, face = "bold")
    )
}

# 优化的读取函数 - 只读取需要的列
read_bed_efficient <- function(filepath, source_name) {
  # 首先确定文件的列数
  first_line <- readLines(filepath, n = 1)
  n_cols <- length(strsplit(first_line, "\t")[[1]])
  
  if (n_cols == 4) {
    # 对于4列的文件,只读取需要的列
    col_classes <- c("character", "numeric", "numeric", "numeric")
    dt <- fread(filepath, sep = "\t", header = FALSE, 
                colClasses = col_classes,
                select = 1:4)
    setnames(dt, c("V1", "V2", "V3", "V4"))
    
    # 标准化为统一格式
    dt[, `:=`(V5 = V1, V6 = V2, V7 = (V2 + V3) / 2, V8 = V4)]
    dt[, c("V1", "V2", "V3", "V4") := NULL]  # 删除不需要的列
  } else {
    # 对于其他文件,只读取需要的列 (5, 7, 8)
    dt <- fread(filepath, sep = "\t", header = FALSE,
                select = c(5, 7, 8))
    setnames(dt, c("V5", "V7", "V8"))
  }
  
  dt[, source := source_name]
  return(dt)
}

# 获取唯一染色体列表
get_unique_chromosomes <- function() {
  # 分别读取每个文件的染色体信息
  chr_tba <- fread("../../tba.bed", sep = "\t", header = FALSE, select = 5)[[1]]
  chr_tra <- fread("../../tra.bed", sep = "\t", header = FALSE, select = 5)[[1]]
  
  # 合并并去重
  all_chr <- unique(c(chr_tba, chr_tra))
  
  # 清理内存
  rm(chr_tba, chr_tra)
  gc()
  
  return(all_chr)
}

# 优化的创建绘图数据函数
create_plot_data <- function(data, window_size) {
  if (nrow(data) == 0) return(data.frame())
  
  # 使用data.table进行高效分组计算
  dt <- as.data.table(data)
  
  result <- dt[, .(
    pos = mean(floor((V7 - min(V7)) / window_size) * window_size + min(V7)) + window_size/2,
    mean_y = mean(adjusted_y, na.rm = TRUE),
    se_y = sd(adjusted_y, na.rm = TRUE) / sqrt(sum(!is.na(adjusted_y))),
    n_points = .N
  ), by = .(chr = V5, 
            window = floor((V7 - min(V7)) / window_size) * window_size + min(V7))]
  
  result <- result[!is.na(mean_y) & !is.na(pos)]
  setorder(result, pos)
  
  return(as.data.frame(result))
}

# 改进的平滑函数
apply_smoothing <- function(data, params = smoothing_params) {
  if (nrow(data) < 4) {
    data$smooth_y <- data$mean_y
    return(data)
  }
  
  if (!params$enabled || params$method == "none") {
    data$smooth_y <- data$mean_y
    return(data)
  }
  
  tryCatch({
    if (params$method == "loess") {
      smooth_fit <- loess(mean_y ~ pos, data = data, 
                         span = params$span, 
                         degree = params$degree)
      data$smooth_y <- predict(smooth_fit)
      
    } else if (params$method == "spline") {
      smooth_fit <- smooth.spline(data$pos, data$mean_y, spar = params$spar)
      data$smooth_y <- predict(smooth_fit, data$pos)$y
      
    } else if (params$method == "moving_average") {
      n <- nrow(data)
      window <- params$ma_window
      data$smooth_y <- data$mean_y
      
      for (i in 1:n) {
        start_idx <- max(1, i - floor(window/2))
        end_idx <- min(n, i + floor(window/2))
        data$smooth_y[i] <- mean(data$mean_y[start_idx:end_idx], na.rm = TRUE)
      }
      
    } else {
      data$smooth_y <- data$mean_y
    }
    
    na_indices <- is.na(data$smooth_y)
    if (any(na_indices)) {
      data$smooth_y[na_indices] <- data$mean_y[na_indices]
    }
    
    return(data)
  }, error = function(e) {
    cat("平滑处理失败,使用原始数据:", e$message, "\n")
    data$smooth_y <- data$mean_y
    return(data)
  })
}

# 竖线坐标定义
vertical_lines <- list(
  "1" = c(40492542, 80104410),
  "2" = 45654190,
  "3" = 53945679,
  "4" = 44435307,
  "5" = 42365773,
  "6" = 38964430,
  "7" = 39901015,
  "8" = 36517876,
  "9" = 41362729,
  "10" = 33661800
)

# 显示当前平滑参数
cat("\n当前平滑参数设置:\n")
cat("启用平滑:", smoothing_params$enabled, "\n")
cat("平滑方法:", smoothing_params$method, "\n")
if (smoothing_params$method == "loess") {
  cat("Loess span:", smoothing_params$span, "(值越小越不平滑)\n")
  cat("Loess degree:", smoothing_params$degree, "\n")
} else if (smoothing_params$method == "spline") {
  cat("Spline spar:", smoothing_params$spar, "\n")
} else if (smoothing_params$method == "moving_average") {
  cat("移动平均窗口:", smoothing_params$ma_window, "\n")
}

# 获取所有染色体
all_chromosomes <- get_unique_chromosomes()

# 用于存储所有断点信息的列表
all_breakpoint_data <- list()

# 主处理循环 - 逐个染色体处理,避免一次性加载所有数据
for (chr in all_chromosomes) {
  chr_fullname <- chr
  chr_num <- gsub("TraScf_", "", chr_fullname)
  
  cat("\n正在处理染色体:", chr_fullname, "\n")
  
  # 逐个文件读取并过滤当前染色体的数据
  # TBA数据
  tba_data <- tryCatch({
    dt <- read_bed_efficient("../../tba.bed", "TBA")
    dt <- dt[V5 == chr_fullname]
    dt[, adjusted_y := ifelse(V8 >= 0, V8 / x1, NA_real_)]
    dt[!is.na(adjusted_y) & !is.na(V7) & is.finite(adjusted_y)]
  }, error = function(e) {
    cat("读取TBA数据失败:", e$message, "\n")
    data.table()
  })
  
  # TRA数据
  tra_data <- tryCatch({
    dt <- read_bed_efficient("../../tra.bed", "TRA")
    dt <- dt[V5 == chr_fullname]
    dt[, adjusted_y := ifelse(V8 >= 0, V8 / x2, NA_real_)]
    dt[!is.na(adjusted_y) & !is.na(V7) & is.finite(adjusted_y)]
  }, error = function(e) {
    cat("读取TRA数据失败:", e$message, "\n")
    data.table()
  })
  
  # 合并数据
  chr_data <- rbindlist(list(tba_data, tra_data), use.names = TRUE, fill = TRUE)
  
  if (nrow(chr_data) == 0) {
    cat("染色体", chr_fullname, "无有效数据,跳过\n")
    next
  }
  
  cat("数据点数量 - TBA:", nrow(tba_data), ", TRA:", nrow(tra_data), "\n")
  
  # 准备绘图数据
  plot_data_tba <- data.frame()
  plot_data_tra <- data.frame()
  
  # TBA数据
  if (nrow(tba_data) > 0) {
    tba_df <- as.data.frame(tba_data)
    plot_data_tba <- create_plot_data(tba_df, plot_window)
    plot_data_tba <- apply_smoothing(plot_data_tba)
    cat("TBA窗口数据点:", nrow(plot_data_tba), "\n")
  }
  
  # TRA数据
  if (nrow(tra_data) > 0) {
    tra_df <- as.data.frame(tra_data)
    plot_data_tra <- create_plot_data(tra_df, plot_window)
    plot_data_tra <- apply_smoothing(plot_data_tra)
    cat("TRA窗口数据点:", nrow(plot_data_tra), "\n")
  }
  
  # 合并所有绘图数据
  plot_data <- rbind(
    if (nrow(plot_data_tba) > 0) {
      data.frame(plot_data_tba, source = "TBA")
    } else NULL,
    if (nrow(plot_data_tra) > 0) {
      data.frame(plot_data_tra, source = "TRA")
    } else NULL
  )
  
  if (nrow(plot_data) == 0) {
    cat("无足够数据绘图,跳过染色体", chr_fullname, "\n")
    next
  }
  
  # 确定Y轴范围
  y_min <- min(
    if (nrow(plot_data_tba) > 0) min(plot_data_tba$smooth_y - plot_data_tba$se_y, na.rm = TRUE) else Inf,
    if (nrow(plot_data_tra) > 0) min(plot_data_tra$smooth_y - plot_data_tra$se_y, na.rm = TRUE) else Inf,
    na.rm = TRUE
  )
  
  y_max <- max(
    if (nrow(plot_data_tba) > 0) max(plot_data_tba$smooth_y + plot_data_tba$se_y, na.rm = TRUE) else -Inf,
    if (nrow(plot_data_tra) > 0) max(plot_data_tra$smooth_y + plot_data_tra$se_y, na.rm = TRUE) else -Inf,
    na.rm = TRUE
  )
  
  y_margin <- (y_max - y_min) * 0.1
  y_lim <- c(max(0, y_min - y_margin), y_max + y_margin)
  
  # ================== 为汇总图收集断点数据 ==================
  if (chr_num %in% names(vertical_lines)) {
    marker_positions <- vertical_lines[[chr_num]]
    
    # 获取染色体数据范围
    min_pos <- min(chr_data$V7, na.rm = TRUE)
    max_pos <- max(chr_data$V7, na.rm = TRUE)
    
    # 对每个断点收集数据
    for (bp_idx in seq_along(marker_positions)) {
      bp_pos <- marker_positions[bp_idx]
      
      # 确定左右边界
      left_bound <- max(min_pos, bp_pos - summary_extension)
      right_bound <- min(max_pos, bp_pos + summary_extension)
      
      # 如果有多个断点,避免重叠
      if (length(marker_positions) > 1) {
        if (bp_idx > 1) {
          left_bound <- max(left_bound, (marker_positions[bp_idx-1] + bp_pos) / 2)
        }
        if (bp_idx < length(marker_positions)) {
          right_bound <- min(right_bound, (bp_pos + marker_positions[bp_idx+1]) / 2)
        }
      }
      
      # 提取该区域的数据
      region_data_tba <- plot_data_tba %>% filter(pos >= left_bound & pos <= right_bound)
      region_data_tra <- plot_data_tra %>% filter(pos >= left_bound & pos <= right_bound)
      
      if (nrow(region_data_tba) > 0 || nrow(region_data_tra) > 0) {
        # 将位置转换为相对于断点的距离
        if (nrow(region_data_tba) > 0) {
          region_data_tba$relative_pos <- region_data_tba$pos - bp_pos
          region_data_tba$chr <- chr_fullname
          region_data_tba$breakpoint <- bp_pos
          region_data_tba$bp_label <- paste0("Chr", chr_num, "_BP", 
                                             if(length(marker_positions) > 1) bp_idx else "")
        }
        
        if (nrow(region_data_tra) > 0) {
          region_data_tra$relative_pos <- region_data_tra$pos - bp_pos
          region_data_tra$chr <- chr_fullname
          region_data_tra$breakpoint <- bp_pos
          region_data_tra$bp_label <- paste0("Chr", chr_num, "_BP", 
                                             if(length(marker_positions) > 1) bp_idx else "")
        }
        
        # 计算扫描线位置
        dt_scan <- as.data.table(chr_data)
        setkey(dt_scan, V7)
        
        scan_positions_bp <- c()
        
        # 左向扫描
        current_window_start <- bp_pos - scan_window
        while(current_window_start >= left_bound) {
          window_end <- current_window_start + scan_window
          window_data <- dt_scan[V7 >= current_window_start & V7 < window_end & source %in% c("TBA", "TRA")]
          
          tba_mean <- window_data[source == "TBA", mean(adjusted_y, na.rm = TRUE)]
          tra_mean <- window_data[source == "TRA", mean(adjusted_y, na.rm = TRUE)]
          
          if (!is.na(tra_mean) && !is.na(tba_mean) && tra_mean > tba_mean) {
            median_pos <- current_window_start + scan_window/2
            scan_positions_bp <- c(scan_positions_bp, median_pos)
            break
          }
          current_window_start <- current_window_start - scan_step
        }
        
        # 右向扫描
        current_window_start <- bp_pos
        while(current_window_start + scan_window <= right_bound) {
          window_end <- current_window_start + scan_window
          window_data <- dt_scan[V7 >= current_window_start & V7 < window_end & source %in% c("TBA", "TRA")]
          
          tba_mean <- window_data[source == "TBA", mean(adjusted_y, na.rm = TRUE)]
          tra_mean <- window_data[source == "TRA", mean(adjusted_y, na.rm = TRUE)]
          
          if (!is.na(tra_mean) && !is.na(tba_mean) && tra_mean > tba_mean) {
            median_pos <- current_window_start + scan_window/2
            scan_positions_bp <- c(scan_positions_bp, median_pos)
            break
          }
          current_window_start <- current_window_start + scan_step
        }
        
        # 保存断点数据
        bp_data <- list(
          tba = region_data_tba,
          tra = region_data_tra,
          scan_lines = scan_positions_bp - bp_pos,  # 转换为相对位置
          bp_label = paste0("Chr", chr_num, "_BP", if(length(marker_positions) > 1) bp_idx else "")
        )
        
        all_breakpoint_data[[length(all_breakpoint_data) + 1]] <- bp_data
      }
    }
  }
  # ================== 断点数据收集结束 ==================
  
  # 创建单个染色体的图
  p <- ggplot() +
    create_publication_theme() +
    labs(
      title = paste("Chromosome", chr_num),
      subtitle = "ISMC Analysis Comparison",
      x = "Position (Mb)",
      y = "RHO Value"
    ) +
    coord_cartesian(ylim = y_lim)
  
  # 添加标准误差带
  if (nrow(plot_data_tba) > 0) {
    p <- p + 
      geom_ribbon(
        data = plot_data_tba,
        aes(x = pos/1e6, ymin = smooth_y - se_y, ymax = smooth_y + se_y),
        fill = publication_colors$tba,
        alpha = 0.2
      )
  }
  
  if (nrow(plot_data_tra) > 0) {
    p <- p + 
      geom_ribbon(
        data = plot_data_tra,
        aes(x = pos/1e6, ymin = smooth_y - se_y, ymax = smooth_y + se_y),
        fill = publication_colors$tra,
        alpha = 0.2
      )
  }
  
  # 绘制主要数据线
  if (nrow(plot_data_tba) > 0) {
    p <- p + 
      geom_line(
        data = plot_data_tba,
        aes(x = pos/1e6, y = smooth_y),
        color = publication_colors$tba,
        size = 1.2,
        alpha = 0.9
      )
  }
  
  if (nrow(plot_data_tra) > 0) {
    p <- p + 
      geom_line(
        data = plot_data_tra,
        aes(x = pos/1e6, y = smooth_y),
        color = publication_colors$tra,
        size = 1.2,
        alpha = 0.9
      )
  }
  
  # 美化Y轴
  p <- p + 
    scale_y_continuous(
      name = "RHO Value",
      labels = scales::scientific_format(digits = 2)
    )
  
  # 美化X轴
  p <- p + 
    scale_x_continuous(
      labels = scales::number_format(accuracy = 0.1),
      breaks = scales::pretty_breaks(n = 8)
    )
  
  # 创建图例
  legend_colors <- c()
  legend_labels <- c()
  
  if (nrow(plot_data_tba) > 0) {
    legend_colors <- c(legend_colors, publication_colors$tba)
    legend_labels <- c(legend_labels, "TBA")
  }
  
  if (nrow(plot_data_tra) > 0) {
    legend_colors <- c(legend_colors, publication_colors$tra)
    legend_labels <- c(legend_labels, "TRA")
  }
  
  if (length(legend_colors) > 0) {
    # 创建虚拟数据用于图例
    dummy_data <- data.frame(
      x = rep(-Inf, length(legend_labels)),
      y = rep(-Inf, length(legend_labels)),
      Method = legend_labels
    )
    
    p <- p + 
      geom_point(
        data = dummy_data,
        aes(x = x, y = y, color = Method),
        size = 3, alpha = 0
      ) +
      scale_color_manual(
        name = "Sample",
        values = setNames(legend_colors, legend_labels),
        labels = setNames(c("ISMC Analysis of TBA", "ISMC Analysis of TRA")[1:length(legend_labels)], legend_labels)
      ) +
      guides(color = guide_legend(
        override.aes = list(
          alpha = 1, 
          size = 6,
          linewidth = 2
        ),
        title.position = "top",
        byrow = TRUE
      ))
  }
  
  # 添加标记线
  if (chr_num %in% names(vertical_lines)) {
    marker_positions <- vertical_lines[[chr_num]]
    
    p <- p + 
      geom_vline(
        xintercept = marker_positions/1e6,
        color = publication_colors$marker,
        linetype = "dashed",
        size = 0.8,
        alpha = 0.7
      )
    
    # 添加标记线标签
    for (i in seq_along(marker_positions)) {
      p <- p + 
        annotate(
          "text",
          x = marker_positions[i]/1e6,
          y = Inf,
          label = paste("Break Point"),
          color = publication_colors$marker,
          vjust = 1.2,
          size = 3.5,
          fontface = "bold"
        )
    }
    
    # 扫描逻辑 - 使用data.table进行高效扫描
    min_pos <- min(chr_data$V7, na.rm = TRUE)
    max_pos <- max(chr_data$V7, na.rm = TRUE)
    
    scan_positions <- c()
    
    # 将chr_data转换为data.table以提高扫描效率
    dt_scan <- as.data.table(chr_data)
    setkey(dt_scan, V7)  # 设置键以加速范围查询
    
    for (mark_pos in marker_positions) {
      # 左向扫描
      current_window_start <- mark_pos - scan_window
      while(current_window_start >= min_pos) {
        window_end <- current_window_start + scan_window
        
        # 使用data.table的范围查询
        window_data <- dt_scan[V7 >= current_window_start & V7 < window_end & source %in% c("TBA", "TRA")]
        
        tba_mean <- window_data[source == "TBA", mean(adjusted_y, na.rm = TRUE)]
        tra_mean <- window_data[source == "TRA", mean(adjusted_y, na.rm = TRUE)]
        
        if (!is.na(tra_mean) && !is.na(tba_mean) && tra_mean > tba_mean) {
          median_pos <- current_window_start + scan_window/2
          scan_positions <- c(scan_positions, median_pos)
          break
        }
        current_window_start <- current_window_start - scan_step
      }
      
      # 右向扫描
      current_window_start <- mark_pos
      while(current_window_start + scan_window <= max_pos) {
        window_end <- current_window_start + scan_window
        
        window_data <- dt_scan[V7 >= current_window_start & V7 < window_end & source %in% c("TBA", "TRA")]
        
        tba_mean <- window_data[source == "TBA", mean(adjusted_y, na.rm = TRUE)]
        tra_mean <- window_data[source == "TRA", mean(adjusted_y, na.rm = TRUE)]
        
        if (!is.na(tra_mean) && !is.na(tba_mean) && tra_mean > tba_mean) {
          median_pos <- current_window_start + scan_window/2
          scan_positions <- c(scan_positions, median_pos)
          break
        }
        current_window_start <- current_window_start + scan_step
      }
    }
    
    # 添加扫描结果线
    if (length(scan_positions) > 0) {
      p <- p + 
        geom_vline(
          xintercept = scan_positions/1e6,
          color = publication_colors$scan,
          linetype = "dotted",
          size = 0.6,
          alpha = 0.8
        )
    }
  }
  
  # 添加统计信息文本框
  stats_text <- paste0(
    "Window: ", format(plot_window/1e6, digits = 1), " Mb\n",
    "Scan: ", format(scan_window/1e6, digits = 1), " Mb\n",
    "Data points: ", nrow(chr_data)
  )
  
  p <- p + 
    annotate(
      "text",
      x = Inf, y = Inf,
      label = stats_text,
      hjust = 1.02, vjust = 1.02,
      size = 3,
      color = publication_colors$text,
      fontface = "italic",
      alpha = 0.8
    )
  
  # 保存高质量图形
  # PNG格式 - 高分辨率
  ggsave(
    filename = paste0(chr_fullname, "_publication.png"), 
    plot = p, 
    width = 12, height = 8, 
    dpi = 600,
    bg = "white"
  )
  
  # PDF格式 - 矢量图
  ggsave(
    filename = paste0(chr_fullname, "_publication.pdf"), 
    plot = p, 
    width = 12, height = 8, 
    device = "pdf",
    bg = "white"
  )
  
  # TIFF格式 - 印刷质量
  ggsave(
    filename = paste0(chr_fullname, "_publication.tiff"), 
    plot = p, 
    width = 12, height = 8, 
    dpi = 300,
    compression = "lzw",
    bg = "white"
  )
  
  cat("高质量图形已保存:\n")
  cat("- PNG:", paste0(chr_fullname, "_publication.png"), "\n")
  cat("- PDF:", paste0(chr_fullname, "_publication.pdf"), "\n")
  cat("- TIFF:", paste0(chr_fullname, "_publication.tiff"), "\n\n")
  
  # 清理内存 - 处理完每个染色体后清理
  rm(chr_data, plot_data_tba, plot_data_tra, plot_data, p)
  gc()
}

# ==================== 创建断点汇总图 ====================
cat("\n==================== 开始创建断点汇总图 ====================\n")

if (length(all_breakpoint_data) > 0) {
  # 合并所有断点的TBA和TRA数据
  all_tba_summary <- do.call(rbind, lapply(all_breakpoint_data, function(x) x$tba))
  all_tra_summary <- do.call(rbind, lapply(all_breakpoint_data, function(x) x$tra))
  
  # 确定Y轴范围
  y_min_summary <- min(
    if (!is.null(all_tba_summary) && nrow(all_tba_summary) > 0) 
      min(all_tba_summary$smooth_y - all_tba_summary$se_y, na.rm = TRUE) else Inf,
    if (!is.null(all_tra_summary) && nrow(all_tra_summary) > 0)
      min(all_tra_summary$smooth_y - all_tra_summary$se_y, na.rm = TRUE) else Inf,
    na.rm = TRUE
  )
  
  y_max_summary <- max(
    if (!is.null(all_tba_summary) && nrow(all_tba_summary) > 0)
      max(all_tba_summary$smooth_y + all_tba_summary$se_y, na.rm = TRUE) else -Inf,
    if (!is.null(all_tra_summary) && nrow(all_tra_summary) > 0)
      max(all_tra_summary$smooth_y + all_tra_summary$se_y, na.rm = TRUE) else -Inf,
    na.rm = TRUE
  )
  
  y_margin_summary <- (y_max_summary - y_min_summary) * 0.1
  y_lim_summary <- c(max(0, y_min_summary - y_margin_summary), y_max_summary + y_margin_summary)
  
  # 创建汇总图
  p_summary <- ggplot() +
    create_publication_theme() +
    labs(
      title = "All Breakpoints Summary",
      subtitle = "ISMC Analysis Around Chromosomal Breakpoints",
      x = "Relative Position to Breakpoint (Mb)",
      y = "RHO Value"
    ) +
    coord_cartesian(ylim = y_lim_summary)
  
  # 为每个断点添加一个facet标签
  if (!is.null(all_tba_summary) && nrow(all_tba_summary) > 0) {
    all_tba_summary$facet_label <- all_tba_summary$bp_label
  }
  if (!is.null(all_tra_summary) && nrow(all_tra_summary) > 0) {
    all_tra_summary$facet_label <- all_tra_summary$bp_label
  }
  
  # 添加标准误差带
  if (!is.null(all_tba_summary) && nrow(all_tba_summary) > 0) {
    p_summary <- p_summary + 
      geom_ribbon(
        data = all_tba_summary,
        aes(x = relative_pos/1e6, ymin = smooth_y - se_y, ymax = smooth_y + se_y, group = bp_label),
        fill = publication_colors$tba,
        alpha = 0.15
      )
  }
  
  if (!is.null(all_tra_summary) && nrow(all_tra_summary) > 0) {
    p_summary <- p_summary + 
      geom_ribbon(
        data = all_tra_summary,
        aes(x = relative_pos/1e6, ymin = smooth_y - se_y, ymax = smooth_y + se_y, group = bp_label),
        fill = publication_colors$tra,
        alpha = 0.15
      )
  }
  
  # 绘制主要数据线
  if (!is.null(all_tba_summary) && nrow(all_tba_summary) > 0) {
    p_summary <- p_summary + 
      geom_line(
        data = all_tba_summary,
        aes(x = relative_pos/1e6, y = smooth_y, group = bp_label),
        color = publication_colors$tba,
        size = 0.8,
        alpha = 0.7
      )
  }
  
  if (!is.null(all_tra_summary) && nrow(all_tra_summary) > 0) {
    p_summary <- p_summary + 
      geom_line(
        data = all_tra_summary,
        aes(x = relative_pos/1e6, y = smooth_y, group = bp_label),
        color = publication_colors$tra,
        size = 0.8,
        alpha = 0.7
      )
  }
  
  # 添加断点中心线(x=0)
  p_summary <- p_summary + 
    geom_vline(
      xintercept = 0,
      color = publication_colors$marker,
      linetype = "dashed",
      size = 1.0,
      alpha = 0.8
    ) +
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "Breakpoint",
      color = publication_colors$marker,
      vjust = 1.2,
      size = 4,
      fontface = "bold"
    )
  
  # 添加所有扫描线(绿色虚线) - 不变
  for (bp_data in all_breakpoint_data) {
    if (length(bp_data$scan_lines) > 0) {
      p_summary <- p_summary + 
        geom_vline(
          xintercept = bp_data$scan_lines/1e6,
          color = publication_colors$scan,
          linetype = "dotted",
          size = 0.5,
          alpha = 0.6
        )
    }
  }
  
  # 美化Y轴
  p_summary <- p_summary + 
    scale_y_continuous(
      name = "RHO Value",
      labels = scales::scientific_format(digits = 2)
    )
  
  # 美化X轴
  p_summary <- p_summary + 
    scale_x_continuous(
      labels = scales::number_format(accuracy = 0.1),
      breaks = scales::pretty_breaks(n = 10)
    )
  
  # 创建图例
  dummy_data_summary <- data.frame(
    x = rep(-Inf, 2),
    y = rep(-Inf, 2),
    Method = c("TBA", "TRA")
  )
  
  p_summary <- p_summary + 
    geom_point(
      data = dummy_data_summary,
      aes(x = x, y = y, color = Method),
      size = 3, alpha = 0
    ) +
    scale_color_manual(
      name = "Sample",
      values = c("TBA" = publication_colors$tba, "TRA" = publication_colors$tra),
      labels = c("TBA" = "ISMC Analysis of TBA", "TRA" = "ISMC Analysis of TRA")
    ) +
    guides(color = guide_legend(
      override.aes = list(
        alpha = 1, 
        size = 6,
        linewidth = 2
      ),
      title.position = "top",
      byrow = TRUE
    ))
  
  # 添加统计信息
  stats_text_summary <- paste0(
    "Total Breakpoints: ", length(all_breakpoint_data), "\n",
    "Window: ", format(plot_window/1e6, digits = 1), " Mb\n",
    "Extension: ", format(summary_extension/1e6, digits = 1), " Mb"
  )
  
  p_summary <- p_summary + 
    annotate(
      "text",
      x = Inf, y = Inf,
      label = stats_text_summary,
      hjust = 1.02, vjust = 1.02,
      size = 3,
      color = publication_colors$text,
      fontface = "italic",
      alpha = 0.8
    )
  
  # 保存汇总图
  ggsave(
    filename = "All_Breakpoints_Summary.png", 
    plot = p_summary, 
    width = 14, height = 8, 
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    filename = "All_Breakpoints_Summary.pdf", 
    plot = p_summary, 
    width = 14, height = 8, 
    device = "pdf",
    bg = "white"
  )
  
  ggsave(
    filename = "All_Breakpoints_Summary.tiff", 
    plot = p_summary, 
    width = 14, height = 8, 
    dpi = 300,
    compression = "lzw",
    bg = "white"
  )
  
  cat("\n断点汇总图已保存:\n")
  cat("- PNG: All_Breakpoints_Summary.png\n")
  cat("- PDF: All_Breakpoints_Summary.pdf\n")
  cat("- TIFF: All_Breakpoints_Summary.tiff\n")
  
  # ==================== 新增: 样本单线汇总图（每个Sample仅一条线） ====================
  cat("\n==================== 开始创建样本单线汇总图 ====================\n")
  
  # 助手函数: 将相对位置按与plot_window一致的网格归并到中心点上
  aggregate_by_relative_bin <- function(df, window_size) {
    if (is.null(df) || nrow(df) == 0) return(data.frame())
    dt <- as.data.table(df)[!is.na(relative_pos) & !is.na(smooth_y)]
    # 使用round保证0附近对称的网格中心
    dt[, rel_bin := round(relative_pos / window_size) * window_size]
    res <- dt[, .(
      pos = rel_bin[1],
      mean_y = mean(smooth_y, na.rm = TRUE),
      se_y = {
        n <- .N
        if (n <= 1) 0 else sd(smooth_y, na.rm = TRUE) / sqrt(n)
      },
      n_breakpoints = .N
    ), by = rel_bin]
    setorder(res, pos)
    as.data.frame(res)
  }
  
  # 新增助手函数: 将相对位置按“绝对距离”折叠后再按网格归并
  aggregate_by_abs_bin <- function(df, window_size) {
    if (is.null(df) || nrow(df) == 0) return(data.frame())
    dt <- as.data.table(df)[!is.na(relative_pos) & !is.na(smooth_y)]
    dt[, abs_pos := abs(relative_pos)]
    dt[, abs_bin := round(abs_pos / window_size) * window_size]
    res <- dt[, .(
      pos = abs_bin[1],  # 使用折叠后的距离作为x
      mean_y = mean(smooth_y, na.rm = TRUE),
      se_y = {
        n <- .N
        if (n <= 1) 0 else sd(smooth_y, na.rm = TRUE) / sqrt(n)
      },
      n_breakpoints = .N
    ), by = abs_bin]
    setorder(res, pos)
    as.data.frame(res)
  }
  
  # 新增助手函数: 计算所有断点共同覆盖的相对/绝对范围
  compute_common_ranges <- function(bp_list) {
    if (length(bp_list) == 0) return(NULL)
    ranges <- lapply(bp_list, function(bp) {
      rel <- c()
      if (!is.null(bp$tba) && nrow(bp$tba) > 0) rel <- c(rel, bp$tba$relative_pos)
      if (!is.null(bp$tra) && nrow(bp$tra) > 0) rel <- c(rel, bp$tra$relative_pos)
      rel <- rel[is.finite(rel)]
      if (length(rel) == 0) return(NULL)
      list(
        min = min(rel),
        max = max(rel),
        max_abs = max(abs(rel))
      )
    })
    ranges <- ranges[!sapply(ranges, is.null)]
    if (length(ranges) == 0) return(NULL)
    list(
      left = max(sapply(ranges, `[[`, "min")),
      right = min(sapply(ranges, `[[`, "max")),
      max_abs = min(sapply(ranges, `[[`, "max_abs"))
    )
  }
  
  agg_tba <- aggregate_by_relative_bin(all_tba_summary, plot_window)
  agg_tra <- aggregate_by_relative_bin(all_tra_summary, plot_window)
  
  if ((nrow(agg_tba) + nrow(agg_tra)) == 0) {
    cat("警告: 无可用于样本单线汇总的数据\n")
  } else {
    # 计算Y轴范围（与原汇总图一致的方式）
    y_min_single <- min(
      if (nrow(agg_tba) > 0) min(agg_tba$mean_y - agg_tba$se_y, na.rm = TRUE) else Inf,
      if (nrow(agg_tra) > 0) min(agg_tra$mean_y - agg_tra$se_y, na.rm = TRUE) else Inf,
      na.rm = TRUE
    )
    y_max_single <- max(
      if (nrow(agg_tba) > 0) max(agg_tba$mean_y + agg_tba$se_y, na.rm = TRUE) else -Inf,
      if (nrow(agg_tra) > 0) max(agg_tra$mean_y + agg_tra$se_y, na.rm = TRUE) else -Inf,
      na.rm = TRUE
    )
    y_margin_single <- (y_max_single - y_min_single) * 0.1
    y_lim_single <- c(max(0, y_min_single - y_margin_single), y_max_single + y_margin_single)
    
    # 创建样本单线汇总图
    p_summary_single <- ggplot() +
      create_publication_theme() +
      labs(
        title = "All Breakpoints Summary (Aggregated by Sample)",
        subtitle = "Mean ISMC Across Breakpoints With SE Ribbon",
        x = "Relative Position to Breakpoint (Mb)",
        y = "RHO Value"
      ) +
      coord_cartesian(ylim = y_lim_single)
    
    # 更高透明度的带状区域展示偏差（±SE）
    agg_ribbon_alpha <- 0.12  # 比原0.15更透明
    
    if (nrow(agg_tba) > 0) {
      p_summary_single <- p_summary_single +
        geom_ribbon(
          data = agg_tba,
          aes(x = pos/1e6, ymin = mean_y - se_y, ymax = mean_y + se_y),
          fill = publication_colors$tba,
          alpha = agg_ribbon_alpha
        ) +
        geom_line(
          data = agg_tba,
          aes(x = pos/1e6, y = mean_y),
          color = publication_colors$tba,
          size = 0.8,
          alpha = 0.9
        )
    }
    if (nrow(agg_tra) > 0) {
      p_summary_single <- p_summary_single +
        geom_ribbon(
          data = agg_tra,
          aes(x = pos/1e6, ymin = mean_y - se_y, ymax = mean_y + se_y),
          fill = publication_colors$tra,
          alpha = agg_ribbon_alpha
        ) +
        geom_line(
          data = agg_tra,
          aes(x = pos/1e6, y = mean_y),
          color = publication_colors$tra,
          size = 0.8,
          alpha = 0.9
        )
    }
    
    # 断点中心线(x=0) - 保持一致
    p_summary_single <- p_summary_single +
      geom_vline(
        xintercept = 0,
        color = publication_colors$marker,
        linetype = "dashed",
        size = 1.0,
        alpha = 0.8
      ) +
      annotate(
        "text",
        x = 0,
        y = Inf,
        label = "Breakpoint",
        color = publication_colors$marker,
        vjust = 1.2,
        size = 4,
        fontface = "bold"
      )
    
    # 绿色扫描虚线 - 完全保持与原图一致
    for (bp_data in all_breakpoint_data) {
      if (length(bp_data$scan_lines) > 0) {
        p_summary_single <- p_summary_single +
          geom_vline(
            xintercept = bp_data$scan_lines/1e6,
            color = publication_colors$scan,
            linetype = "dotted",
            size = 0.5,
            alpha = 0.6
          )
      }
    }
    
    # 轴格式
    p_summary_single <- p_summary_single +
      scale_y_continuous(
        name = "RHO Value",
        labels = scales::scientific_format(digits = 2)
      ) +
      scale_x_continuous(
        labels = scales::number_format(accuracy = 0.1),
        breaks = scales::pretty_breaks(n = 10)
      )
    
    # 图例（与原汇总图一致的风格）
    dummy_data_summary_single <- data.frame(
      x = rep(-Inf, 2),
      y = rep(-Inf, 2),
      Method = c("TBA", "TRA")
    )
    p_summary_single <- p_summary_single +
      geom_point(
        data = dummy_data_summary_single,
        aes(x = x, y = y, color = Method),
        size = 3, alpha = 0
      ) +
      scale_color_manual(
        name = "Sample",
        values = c("TBA" = publication_colors$tba, "TRA" = publication_colors$tra),
        labels = c("TBA" = "ISMC Analysis of TBA", "TRA" = "ISMC Analysis of TRA")
      ) +
      guides(color = guide_legend(
        override.aes = list(
          alpha = 1,
          size = 6,
          linewidth = 2
        ),
        title.position = "top",
        byrow = TRUE
      ))
    
    # 统计信息文本（说明聚合逻辑）
    stats_text_single <- paste0(
      "Total Breakpoints: ", length(all_breakpoint_data), "\n",
      "Window: ", format(plot_window/1e6, digits = 1), " Mb\n",
      "Aggregation: mean ± SE across breakpoints"
    )
    p_summary_single <- p_summary_single +
      annotate(
        "text",
        x = Inf, y = Inf,
        label = stats_text_single,
        hjust = 1.02, vjust = 1.02,
        size = 3,
        color = publication_colors$text,
        fontface = "italic",
        alpha = 0.8
      )
    
    # 保存新增图（与原汇总图一致的高质量设置）
    ggsave(
      filename = "All_Breakpoints_Summary_Aggregated.png",
      plot = p_summary_single,
      width = 14, height = 8,
      dpi = 600,
      bg = "white"
    )
    ggsave(
      filename = "All_Breakpoints_Summary_Aggregated.pdf",
      plot = p_summary_single,
      width = 14, height = 8,
      device = "pdf",
      bg = "white"
    )
    ggsave(
      filename = "All_Breakpoints_Summary_Aggregated.tiff",
      plot = p_summary_single,
      width = 14, height = 8,
      dpi = 300,
      compression = "lzw",
      bg = "white"
    )
    
    cat("\n样本单线汇总图已保存:\n")
    cat("- PNG: All_Breakpoints_Summary_Aggregated.png\n")
    cat("- PDF: All_Breakpoints_Summary_Aggregated.pdf\n")
    cat("- TIFF: All_Breakpoints_Summary_Aggregated.tiff\n")
  }
  
  # ==================== 新增: 样本单线汇总图（最短范围截止） ====================
  cat("\n==================== 开始创建样本单线汇总图(最短范围截止) ====================\n")
  
  common_ranges <- compute_common_ranges(all_breakpoint_data)
  if (!is.null(common_ranges) && is.finite(common_ranges$left) && is.finite(common_ranges$right) &&
      common_ranges$left < common_ranges$right) {
    all_tba_common <- all_tba_summary[all_tba_summary$relative_pos >= common_ranges$left &
                                       all_tba_summary$relative_pos <= common_ranges$right, ]
    all_tra_common <- all_tra_summary[all_tra_summary$relative_pos >= common_ranges$left &
                                       all_tra_summary$relative_pos <= common_ranges$right, ]
    
    agg_tba_common <- aggregate_by_relative_bin(all_tba_common, plot_window)
    agg_tra_common <- aggregate_by_relative_bin(all_tra_common, plot_window)
    
    if ((nrow(agg_tba_common) + nrow(agg_tra_common)) == 0) {
      cat("警告: 最短范围内无可用于样本单线汇总的数据\n")
    } else {
      y_min_common <- min(
        if (nrow(agg_tba_common) > 0) min(agg_tba_common$mean_y - agg_tba_common$se_y, na.rm = TRUE) else Inf,
        if (nrow(agg_tra_common) > 0) min(agg_tra_common$mean_y - agg_tra_common$se_y, na.rm = TRUE) else Inf,
        na.rm = TRUE
      )
      y_max_common <- max(
        if (nrow(agg_tba_common) > 0) max(agg_tba_common$mean_y + agg_tba_common$se_y, na.rm = TRUE) else -Inf,
        if (nrow(agg_tra_common) > 0) max(agg_tra_common$mean_y + agg_tra_common$se_y, na.rm = TRUE) else -Inf,
        na.rm = TRUE
      )
      y_margin_common <- (y_max_common - y_min_common) * 0.1
      y_lim_common <- c(max(0, y_min_common - y_margin_common), y_max_common + y_margin_common)
      
      p_summary_common <- ggplot() +
        create_publication_theme() +
        labs(
          title = "All Breakpoints Summary (Aggregated by Sample, Common Range)",
          subtitle = "Stop When Any Chromosome Ends",
          x = "Relative Position to Breakpoint (Mb)",
          y = "RHO Value"
        ) +
        coord_cartesian(ylim = y_lim_common)
      
      agg_ribbon_alpha_common <- 0.12
      
      if (nrow(agg_tba_common) > 0) {
        p_summary_common <- p_summary_common +
          geom_ribbon(
            data = agg_tba_common,
            aes(x = pos/1e6, ymin = mean_y - se_y, ymax = mean_y + se_y),
            fill = publication_colors$tba,
            alpha = agg_ribbon_alpha_common
          ) +
          geom_line(
            data = agg_tba_common,
            aes(x = pos/1e6, y = mean_y),
            color = publication_colors$tba,
            size = 0.8,
            alpha = 0.9
          )
      }
      if (nrow(agg_tra_common) > 0) {
        p_summary_common <- p_summary_common +
          geom_ribbon(
            data = agg_tra_common,
            aes(x = pos/1e6, ymin = mean_y - se_y, ymax = mean_y + se_y),
            fill = publication_colors$tra,
            alpha = agg_ribbon_alpha_common
          ) +
          geom_line(
            data = agg_tra_common,
            aes(x = pos/1e6, y = mean_y),
            color = publication_colors$tra,
            size = 0.8,
            alpha = 0.9
          )
      }
      
      p_summary_common <- p_summary_common +
        geom_vline(
          xintercept = 0,
          color = publication_colors$marker,
          linetype = "dashed",
          size = 1.0,
          alpha = 0.8
        ) +
        annotate(
          "text",
          x = 0,
          y = Inf,
          label = "Breakpoint",
          color = publication_colors$marker,
          vjust = 1.2,
          size = 4,
          fontface = "bold"
        )
      
      for (bp_data in all_breakpoint_data) {
        if (length(bp_data$scan_lines) > 0) {
          p_summary_common <- p_summary_common +
            geom_vline(
              xintercept = bp_data$scan_lines/1e6,
              color = publication_colors$scan,
              linetype = "dotted",
              size = 0.5,
              alpha = 0.6
            )
        }
      }
      
      p_summary_common <- p_summary_common +
        scale_y_continuous(
          name = "RHO Value",
          labels = scales::scientific_format(digits = 2)
        ) +
        scale_x_continuous(
          labels = scales::number_format(accuracy = 0.1),
          breaks = scales::pretty_breaks(n = 10)
        )
      
      dummy_data_summary_common <- data.frame(
        x = rep(-Inf, 2),
        y = rep(-Inf, 2),
        Method = c("TBA", "TRA")
      )
      p_summary_common <- p_summary_common +
        geom_point(
          data = dummy_data_summary_common,
          aes(x = x, y = y, color = Method),
          size = 3, alpha = 0
        ) +
        scale_color_manual(
          name = "Sample",
          values = c("TBA" = publication_colors$tba, "TRA" = publication_colors$tra),
          labels = c("TBA" = "ISMC Analysis of TBA", "TRA" = "ISMC Analysis of TRA")
        ) +
        guides(color = guide_legend(
          override.aes = list(
            alpha = 1,
            size = 6,
            linewidth = 2
          ),
          title.position = "top",
          byrow = TRUE
        ))
      
      stats_text_common <- paste0(
        "Total Breakpoints: ", length(all_breakpoint_data), "\n",
        "Window: ", format(plot_window/1e6, digits = 1), " Mb\n",
        "Common range: [", 
        format(common_ranges$left/1e6, digits = 2), ", ",
        format(common_ranges$right/1e6, digits = 2), "] Mb"
      )
      p_summary_common <- p_summary_common +
        annotate(
          "text",
          x = Inf, y = Inf,
          label = stats_text_common,
          hjust = 1.02, vjust = 1.02,
          size = 3,
          color = publication_colors$text,
          fontface = "italic",
          alpha = 0.8
        )
      
      ggsave(
        filename = "All_Breakpoints_Summary_Aggregated_Common.png",
        plot = p_summary_common,
        width = 14, height = 8,
        dpi = 600,
        bg = "white"
      )
      ggsave(
        filename = "All_Breakpoints_Summary_Aggregated_Common.pdf",
        plot = p_summary_common,
        width = 14, height = 8,
        device = "pdf",
        bg = "white"
      )
      ggsave(
        filename = "All_Breakpoints_Summary_Aggregated_Common.tiff",
        plot = p_summary_common,
        width = 14, height = 8,
        dpi = 300,
        compression = "lzw",
        bg = "white"
      )
      
      cat("\n最短范围截止样本单线汇总图已保存:\n")
      cat("- PNG: All_Breakpoints_Summary_Aggregated_Common.png\n")
      cat("- PDF: All_Breakpoints_Summary_Aggregated_Common.pdf\n")
      cat("- TIFF: All_Breakpoints_Summary_Aggregated_Common.tiff\n")
    }
  } else {
    cat("警告: 无法计算公共范围,跳过最短范围截止汇总图\n")
  }
  
  # ==================== 新增: 按“绝对距离”聚合的样本单线汇总图 ====================
  cat("\n==================== 开始创建按绝对距离聚合汇总图 ====================\n")
  
  agg_tba_abs <- aggregate_by_abs_bin(all_tba_summary, plot_window)
  agg_tra_abs <- aggregate_by_abs_bin(all_tra_summary, plot_window)
  
  if ((nrow(agg_tba_abs) + nrow(agg_tra_abs)) == 0) {
    cat("警告: 无可用于按绝对距离聚合汇总的数据\n")
  } else {
    # Y轴范围
    y_min_abs <- min(
      if (nrow(agg_tba_abs) > 0) min(agg_tba_abs$mean_y - agg_tba_abs$se_y, na.rm = TRUE) else Inf,
      if (nrow(agg_tra_abs) > 0) min(agg_tra_abs$mean_y - agg_tra_abs$se_y, na.rm = TRUE) else Inf,
      na.rm = TRUE
    )
    y_max_abs <- max(
      if (nrow(agg_tba_abs) > 0) max(agg_tba_abs$mean_y + agg_tba_abs$se_y, na.rm = TRUE) else -Inf,
      if (nrow(agg_tra_abs) > 0) max(agg_tra_abs$mean_y + agg_tra_abs$se_y, na.rm = TRUE) else -Inf,
      na.rm = TRUE
    )
    y_margin_abs <- (y_max_abs - y_min_abs) * 0.1
    y_lim_abs <- c(max(0, y_min_abs - y_margin_abs), y_max_abs + y_margin_abs)
    
    # 汇总图(绝对距离)
    p_summary_abs <- ggplot() +
      create_publication_theme() +
      labs(
        title = "All Breakpoints Summary (Aggregated by Absolute Distance)",
        subtitle = "Mean ISMC vs |Distance to Breakpoint| With SE Ribbon",
        x = "Distance to Breakpoint (Mb)",
        y = "RHO Value"
      ) +
      coord_cartesian(ylim = y_lim_abs)
    
    agg_ribbon_alpha_abs <- 0.12
    
    if (nrow(agg_tba_abs) > 0) {
      p_summary_abs <- p_summary_abs +
        geom_ribbon(
          data = agg_tba_abs,
          aes(x = pos/1e6, ymin = mean_y - se_y, ymax = mean_y + se_y),
          fill = publication_colors$tba,
          alpha = agg_ribbon_alpha_abs
        ) +
        geom_line(
          data = agg_tba_abs,
          aes(x = pos/1e6, y = mean_y),
          color = publication_colors$tba,
          size = 0.8,
          alpha = 0.9
        )
    }
    if (nrow(agg_tra_abs) > 0) {
      p_summary_abs <- p_summary_abs +
        geom_ribbon(
          data = agg_tra_abs,
          aes(x = pos/1e6, ymin = mean_y - se_y, ymax = mean_y + se_y),
          fill = publication_colors$tra,
          alpha = agg_ribbon_alpha_abs
        ) +
        geom_line(
          data = agg_tra_abs,
          aes(x = pos/1e6, y = mean_y),
          color = publication_colors$tra,
          size = 0.8,
          alpha = 0.9
        )
    }
    
    # 中心线(x=0) - 表示断点位置（距离为0）
    p_summary_abs <- p_summary_abs +
      geom_vline(
        xintercept = 0,
        color = publication_colors$marker,
        linetype = "dashed",
        size = 1.0,
        alpha = 0.8
      ) +
      annotate(
        "text",
        x = 0,
        y = Inf,
        label = "Breakpoint",
        color = publication_colors$marker,
        vjust = 1.2,
        size = 4,
        fontface = "bold"
      )
    
    # 扫描线: 使用绝对距离折叠(避免左右重复)
    scan_lines_abs <- unique(unlist(lapply(
      all_breakpoint_data,
      function(x) if (length(x$scan_lines) > 0) abs(x$scan_lines) else numeric(0)
    )))
    if (length(scan_lines_abs) > 0) {
      p_summary_abs <- p_summary_abs +
        geom_vline(
          xintercept = scan_lines_abs/1e6,
          color = publication_colors$scan,
          linetype = "dotted",
          size = 0.5,
          alpha = 0.6
        )
    }
    
    # 轴
    p_summary_abs <- p_summary_abs +
      scale_y_continuous(
        name = "RHO Value",
        labels = scales::scientific_format(digits = 2)
      ) +
      scale_x_continuous(
        labels = scales::number_format(accuracy = 0.1),
        breaks = scales::pretty_breaks(n = 10)
      )
    
    # 图例
    dummy_data_summary_abs <- data.frame(
      x = rep(-Inf, 2),
      y = rep(-Inf, 2),
      Method = c("TBA", "TRA")
    )
    p_summary_abs <- p_summary_abs +
      geom_point(
        data = dummy_data_summary_abs,
        aes(x = x, y = y, color = Method),
        size = 3, alpha = 0
      ) +
      scale_color_manual(
        name = "Sample",
        values = c("TBA" = publication_colors$tba, "TRA" = publication_colors$tra),
        labels = c("TBA" = "ISMC Analysis of TBA", "TRA" = "ISMC Analysis of TRA")
      ) +
      guides(color = guide_legend(
        override.aes = list(
          alpha = 1,
          size = 6,
          linewidth = 2
        ),
        title.position = "top",
        byrow = TRUE
      ))
    
    # 统计信息文本
    stats_text_abs <- paste0(
      "Total Breakpoints: ", length(all_breakpoint_data), "\n",
      "Window: ", format(plot_window/1e6, digits = 1), " Mb\n",
      "Aggregation: mean ± SE by |distance| (left+right combined)"
    )
    p_summary_abs <- p_summary_abs +
      annotate(
        "text",
        x = Inf, y = Inf,
        label = stats_text_abs,
        hjust = 1.02, vjust = 1.02,
        size = 3,
        color = publication_colors$text,
        fontface = "italic",
        alpha = 0.8
      )
    
    # 保存图
    ggsave(
      filename = "All_Breakpoints_Summary_Aggregated_Abs.png",
      plot = p_summary_abs,
      width = 14, height = 8,
      dpi = 600,
      bg = "white"
    )
    ggsave(
      filename = "All_Breakpoints_Summary_Aggregated_Abs.pdf",
      plot = p_summary_abs,
      width = 14, height = 8,
      device = "pdf",
      bg = "white"
    )
    ggsave(
      filename = "All_Breakpoints_Summary_Aggregated_Abs.tiff",
      plot = p_summary_abs,
      width = 14, height = 8,
      dpi = 300,
      compression = "lzw",
      bg = "white"
    )
    
    cat("\n按绝对距离聚合汇总图已保存:\n")
    cat("- PNG: All_Breakpoints_Summary_Aggregated_Abs.png\n")
    cat("- PDF: All_Breakpoints_Summary_Aggregated_Abs.pdf\n")
    cat("- TIFF: All_Breakpoints_Summary_Aggregated_Abs.tiff\n")
  }
  
  # ==================== 新增: 按“绝对距离”最短范围截止的样本单线汇总图 ====================
  cat("\n==================== 开始创建按绝对距离最短范围截止汇总图 ====================\n")
  
  common_ranges <- compute_common_ranges(all_breakpoint_data)
  if (!is.null(common_ranges) && is.finite(common_ranges$max_abs) && common_ranges$max_abs > 0) {
    all_tba_abs_common <- all_tba_summary[abs(all_tba_summary$relative_pos) <= common_ranges$max_abs, ]
    all_tra_abs_common <- all_tra_summary[abs(all_tra_summary$relative_pos) <= common_ranges$max_abs, ]
    
    agg_tba_abs_common <- aggregate_by_abs_bin(all_tba_abs_common, plot_window)
    agg_tra_abs_common <- aggregate_by_abs_bin(all_tra_abs_common, plot_window)
    
    if ((nrow(agg_tba_abs_common) + nrow(agg_tra_abs_common)) == 0) {
      cat("警告: 最短范围内无可用于按绝对距离聚合汇总的数据\n")
    } else {
      y_min_abs_common <- min(
        if (nrow(agg_tba_abs_common) > 0) min(agg_tba_abs_common$mean_y - agg_tba_abs_common$se_y, na.rm = TRUE) else Inf,
        if (nrow(agg_tra_abs_common) > 0) min(agg_tra_abs_common$mean_y - agg_tra_abs_common$se_y, na.rm = TRUE) else Inf,
        na.rm = TRUE
      )
      y_max_abs_common <- max(
        if (nrow(agg_tba_abs_common) > 0) max(agg_tba_abs_common$mean_y + agg_tba_abs_common$se_y, na.rm = TRUE) else -Inf,
        if (nrow(agg_tra_abs_common) > 0) max(agg_tra_abs_common$mean_y + agg_tra_abs_common$se_y, na.rm = TRUE) else -Inf,
        na.rm = TRUE
      )
      y_margin_abs_common <- (y_max_abs_common - y_min_abs_common) * 0.1
      y_lim_abs_common <- c(max(0, y_min_abs_common - y_margin_abs_common), y_max_abs_common + y_margin_abs_common)
      
      p_summary_abs_common <- ggplot() +
        create_publication_theme() +
        labs(
          title = "All Breakpoints Summary (Aggregated by Absolute Distance, Common Range)",
          subtitle = "Stop When Any Chromosome Ends",
          x = "Distance to Breakpoint (Mb)",
          y = "RHO Value"
        ) +
        coord_cartesian(ylim = y_lim_abs_common)
      
      agg_ribbon_alpha_abs_common <- 0.12
      
      if (nrow(agg_tba_abs_common) > 0) {
        p_summary_abs_common <- p_summary_abs_common +
          geom_ribbon(
            data = agg_tba_abs_common,
            aes(x = pos/1e6, ymin = mean_y - se_y, ymax = mean_y + se_y),
            fill = publication_colors$tba,
            alpha = agg_ribbon_alpha_abs_common
          ) +
          geom_line(
            data = agg_tba_abs_common,
            aes(x = pos/1e6, y = mean_y),
            color = publication_colors$tba,
            size = 0.8,
            alpha = 0.9
          )
      }
      if (nrow(agg_tra_abs_common) > 0) {
        p_summary_abs_common <- p_summary_abs_common +
          geom_ribbon(
            data = agg_tra_abs_common,
            aes(x = pos/1e6, ymin = mean_y - se_y, ymax = mean_y + se_y),
            fill = publication_colors$tra,
            alpha = agg_ribbon_alpha_abs_common
          ) +
          geom_line(
            data = agg_tra_abs_common,
            aes(x = pos/1e6, y = mean_y),
            color = publication_colors$tra,
            size = 0.8,
            alpha = 0.9
          )
      }
      
      p_summary_abs_common <- p_summary_abs_common +
        geom_vline(
          xintercept = 0,
          color = publication_colors$marker,
          linetype = "dashed",
          size = 1.0,
          alpha = 0.8
        ) +
        annotate(
          "text",
          x = 0,
          y = Inf,
          label = "Breakpoint",
          color = publication_colors$marker,
          vjust = 1.2,
          size = 4,
          fontface = "bold"
        )
      
      scan_lines_abs_common <- unique(unlist(lapply(
        all_breakpoint_data,
        function(x) if (length(x$scan_lines) > 0) abs(x$scan_lines) else numeric(0)
      )))
      if (length(scan_lines_abs_common) > 0) {
        p_summary_abs_common <- p_summary_abs_common +
          geom_vline(
            xintercept = scan_lines_abs_common/1e6,
            color = publication_colors$scan,
            linetype = "dotted",
            size = 0.5,
            alpha = 0.6
          )
      }
      
      p_summary_abs_common <- p_summary_abs_common +
        scale_y_continuous(
          name = "RHO Value",
          labels = scales::scientific_format(digits = 2)
        ) +
        scale_x_continuous(
          labels = scales::number_format(accuracy = 0.1),
          breaks = scales::pretty_breaks(n = 10)
        )
      
      dummy_data_summary_abs_common <- data.frame(
        x = rep(-Inf, 2),
        y = rep(-Inf, 2),
        Method = c("TBA", "TRA")
      )
      p_summary_abs_common <- p_summary_abs_common +
        geom_point(
          data = dummy_data_summary_abs_common,
          aes(x = x, y = y, color = Method),
          size = 3, alpha = 0
        ) +
        scale_color_manual(
          name = "Sample",
          values = c("TBA" = publication_colors$tba, "TRA" = publication_colors$tra),
          labels = c("TBA" = "ISMC Analysis of TBA", "TRA" = "ISMC Analysis of TRA")
        ) +
        guides(color = guide_legend(
          override.aes = list(
            alpha = 1,
            size = 6,
            linewidth = 2
          ),
          title.position = "top",
          byrow = TRUE
        ))
      
      stats_text_abs_common <- paste0(
        "Total Breakpoints: ", length(all_breakpoint_data), "\n",
        "Window: ", format(plot_window/1e6, digits = 1), " Mb\n",
        "Common max |distance|: ", format(common_ranges$max_abs/1e6, digits = 2), " Mb"
      )
      p_summary_abs_common <- p_summary_abs_common +
        annotate(
          "text",
          x = Inf, y = Inf,
          label = stats_text_abs_common,
          hjust = 1.02, vjust = 1.02,
          size = 3,
          color = publication_colors$text,
          fontface = "italic",
          alpha = 0.8
        )
      
      ggsave(
        filename = "All_Breakpoints_Summary_Aggregated_Abs_Common.png",
        plot = p_summary_abs_common,
        width = 14, height = 8,
        dpi = 600,
        bg = "white"
      )
      ggsave(
        filename = "All_Breakpoints_Summary_Aggregated_Abs_Common.pdf",
        plot = p_summary_abs_common,
        width = 14, height = 8,
        device = "pdf",
        bg = "white"
      )
      ggsave(
        filename = "All_Breakpoints_Summary_Aggregated_Abs_Common.tiff",
        plot = p_summary_abs_common,
        width = 14, height = 8,
        dpi = 300,
        compression = "lzw",
        bg = "white"
      )
      
      cat("\n按绝对距离最短范围截止汇总图已保存:\n")
      cat("- PNG: All_Breakpoints_Summary_Aggregated_Abs_Common.png\n")
      cat("- PDF: All_Breakpoints_Summary_Aggregated_Abs_Common.pdf\n")
      cat("- TIFF: All_Breakpoints_Summary_Aggregated_Abs_Common.tiff\n")
    }
  } else {
    cat("警告: 无法计算绝对距离公共范围,跳过最短范围截止绝对距离汇总图\n")
  }
  
} else {
  cat("\n警告:未找到任何断点数据,无法创建汇总图\n")
}

# 创建综合统计报告
cat("\n=== 分析完成 ===\n")
cat("所有染色体的高质量图形已生成\n")
cat("断点汇总图已生成\n")
cat("新增样本单线汇总图已生成\n")
cat("新增按绝对距离聚合的汇总图已生成\n")
cat("新增最短范围截止汇总图已生成\n")
cat("建议使用PDF格式用于出版物投稿\n")
cat("建议使用TIFF格式用于印刷质量要求\n")
cat("PNG格式适合网络展示和演示\n")

# 使用说明
cat("\n=== 平滑度调整说明 ===\n")
cat("1. 调整 smoothing_params$span 值来控制平滑程度:\n")
cat("   - span = 0.05: 非常不平滑,贴近原始数据\n")
cat("   - span = 0.1:  较不平滑 (当前设置)\n")
cat("   - span = 0.3:  中等平滑\n")
cat("   - span = 0.5:  较平滑\n")
cat("   - span = 1.0:  非常平滑\n")
cat("2. 或者设置 smoothing_params$method = 'none' 完全不平滑\n")
cat("3. 其他平滑方法: 'spline', 'moving_average'\n")
cat("\n=== 汇总图说明 ===\n")
cat("- 汇总图显示所有断点周围的ISMC模式\n")
cat("- X轴=0表示断点位置\n")
cat("- 每个断点的数据向两侧延伸至", format(summary_extension/1e6, digits=1), "Mb\n")
cat("- 绿色虚线为扫描结果,两种汇总图中保持一致\n")
cat("- 样本单线汇总图将所有断点在相对位置上聚合为每个Sample的一条曲线,带状区域为跨断点的±SE\n")
cat("- 新增的“按绝对距离聚合”图将左右两侧折叠后按|距离|聚合,展示与断点距离的整体趋势\n")
cat("- 新增的“最短范围截止”图在任一染色体数据终止时不再延伸\n")
