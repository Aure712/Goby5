# === R代码文件: abs_bp_stats_no_midpoint.R ===
# 修改说明：
#   本版本已去掉“同一 TraScf 上多个断点时用相邻断点中点截断分析范围”的逻辑。
#   每个断点独立使用 bp_pos ± summary_extension 作为分析范围，仍受染色体实际边界限制。
#   其余 P 值计算、breakpoint/direction 筛选、聚合统计和画图逻辑保持不变。


library(dplyr)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(data.table)

# ================== 核心参数设置 ==================
# ---------- 输入文件与样本命名设置 ----------
# primary_sample 可在 "tbi" 与 "tba" 之间切换：
#   primary_sample <- "tbi"  # 图例、统计表、控制台输出中显示为 TBI
#   primary_sample <- "tba"  # 图例、统计表、控制台输出中显示为 TBA
primary_sample <- "tba"

input_files <- list(
  tbi = "/public4/group_crf/home/g21shaoy23/ismc/tbi2tra/output2.bed",
  tba = "/fast3/group_crf/home/g21shaoy23/goby5/ismc/tba/20inds/trans/output2.bed",
  tra = "/fast3/group_crf/home/g21shaoy23/goby5/ismc/tba/20inds/trans/tra/new_output2.bed" #for tba
)

# 各样本独立归一化系数
# 说明：
#   - tbi：TBI 样本专用归一化系数
#   - tba：TBA 样本专用归一化系数
#   - tra：TRA 样本专用归一化系数
#
# 注意：原脚本中的 x1 = 0.007300749 原本用于 TBA。
# 如果 TBI 与 TBA 的真实归一化系数不同，请在这里分别修改 tbi 和 tba 的数值。
normalization_factors <- list(
  tbi = 0.011549063,  # TODO: 请替换为 TBI 的真实归一化系数
  tba = 0.007300749,  # 原脚本 TBA/x1 的归一化系数；如有新值请修改
  tra = 0.005889914   # 原脚本 TRA/x2 的归一化系数
)

# 输出目录；默认当前工作目录。如需统一输出到指定文件夹，可改为例如 "./results"
output_dir <- "."

# 根据开关自动生成后续所有命名
primary_sample <- tolower(primary_sample)
if (!primary_sample %in% c("tbi", "tba")) {
  stop("primary_sample 只能设置为 'tbi' 或 'tba'。")
}
if (is.null(input_files[[primary_sample]]) || !nzchar(input_files[[primary_sample]])) {
  stop(sprintf("input_files 中未设置 %s 的输入文件路径。", primary_sample))
}
if (is.null(input_files$tra) || !nzchar(input_files$tra)) {
  stop("input_files 中未设置 tra 的输入文件路径。")
}

primary_bed_file <- input_files[[primary_sample]]
tra_bed_file <- input_files$tra
primary_label <- toupper(primary_sample)
tra_label <- "TRA"
primary_suffix <- primary_label
tra_suffix <- tra_label

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 检查并读取归一化系数
if (is.null(normalization_factors[[primary_sample]]) || !is.finite(normalization_factors[[primary_sample]])) {
  stop(sprintf("normalization_factors 中未设置 %s 的有效归一化系数。", primary_sample))
}
if (is.null(normalization_factors$tra) || !is.finite(normalization_factors$tra)) {
  stop("normalization_factors 中未设置 tra 的有效归一化系数。")
}

plot_window <- 5e5  # 绘图窗口大小 (0.5MB) - 用于高精度图像绘制
scan_window <- 2e6  # 扫描窗口大小 (2MB) - 用于p值和概率的统计检验
scan_step <- 2e6    # 扫描步长 (2MB)
x1 <- normalization_factors[[primary_sample]]  # primary_sample = "tbi" 时取 tbi 系数；= "tba" 时取 tba 系数
x2 <- normalization_factors$tra

summary_extension <- 5e7  # 断点两侧各延伸50Mb

# 差异分析阈值参数 (20Mb)
threshold_dist <- 20e6

# P 值计算方法选择
#   - "wilcox"：配对 Wilcoxon 符号秩检验，保持原脚本默认行为
#   - "t_test"：配对 t-test
#
# 注意：
#   1) 这里只切换 P 值计算方法。
#   2) 用于检验的数据、breakpoint 配对逻辑、<=threshold/>threshold 的方向设置、
#      全局均值/概率统计以及画图逻辑都保持不变。
p_value_method <- "t_test"   # 可改为 "t_test" "wilcox"

p_value_method <- tolower(gsub("-", "_", p_value_method))
if (p_value_method %in% c("ttest", "t.test")) {
  p_value_method <- "t_test"
}
if (!p_value_method %in% c("wilcox", "t_test")) {
  stop("p_value_method 只能设置为 'wilcox' 或 't_test'。")
}

# ================== 新增：TraScf / 断点编号 / 方向筛选设置 ==================
# 说明：
#   1) 该模块只决定哪些断点、哪些方向进入后续 all_breakpoint_data。
#   2) 后续绝对距离聚合、breakpoint 配对 P 值、region P 值和画图逻辑均保持不变。
#   3) 默认不筛选；即 selected_trascf <- NULL 且 selected_breakpoint_targets <- NULL 时，
#      行为与上一版脚本一致。

# ---- 指定要进入分析的 TraScf 染色体 ----
# 可写成 c("TraScf_1", "TraScf_5") 或 c("1", "5")。
# 不筛选时设为 NULL 或 character(0)。
selected_trascf <- NULL
# 示例：
#selected_trascf <- c("TraScf_1", "TraScf_5", "TraScf_8", "TraScf_10")

# ---- 指定要进入分析的断点编号和方向 ----
# chr：TraScf 染色体名称，可写 "TraScf_1" 或 "1"。
# bp_index：该 TraScf 上按断点坐标从小到大排序后的断点编号。
# direction：可选 "left"、"right" 或 "both"。
#   - "left"：只保留 relative_pos <= 0 的一侧。
#   - "right"：只保留 relative_pos >= 0 的一侧。
#   - "both"：左右两侧都保留，相当于不按方向删该断点的数据。
# 不筛选断点和方向时设为 NULL 或 data.frame()。
selected_breakpoint_targets <- NULL
# 示例：
# selected_breakpoint_targets <- data.frame(
#   chr = c("TraScf_1", "TraScf_5", "TraScf_8", "TraScf_10"),
#   bp_index = c(1L, 1L, 1L, 1L),
#   direction = c("right", "right", "right", "right"),
#   stringsAsFactors = FALSE
# )


normalize_trascf_name <- function(x) {
  x <- as.character(x)
  x <- gsub("^TraScf_", "", x)
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  x
}

is_selected_trascf <- function(chr_name, selected_trascf) {
  if (is.null(selected_trascf) || length(selected_trascf) == 0) {
    return(TRUE)
  }

  chr_clean <- normalize_trascf_name(chr_name)
  selected_clean <- normalize_trascf_name(selected_trascf)

  chr_clean %in% selected_clean
}

is_breakpoint_direction_filter_enabled <- function(selected_breakpoint_targets) {
  !is.null(selected_breakpoint_targets) && nrow(as.data.frame(selected_breakpoint_targets)) > 0
}

normalize_breakpoint_direction <- function(direction) {
  tolower(trimws(as.character(direction)))
}

format_selected_breakpoint_targets <- function(selected_breakpoint_targets) {
  if (!is_breakpoint_direction_filter_enabled(selected_breakpoint_targets)) {
    return("ALL breakpoints and directions")
  }

  targets <- as.data.frame(selected_breakpoint_targets, stringsAsFactors = FALSE)
  paste(
    paste0(
      targets$chr,
      "_BP", targets$bp_index,
      "_", normalize_breakpoint_direction(targets$direction)
    ),
    collapse = ", "
  )
}

get_selected_breakpoint_directions <- function(chr_name, bp_index, selected_breakpoint_targets) {
  all_dirs <- c("left", "right")

  if (!is_breakpoint_direction_filter_enabled(selected_breakpoint_targets)) {
    return(all_dirs)
  }

  targets <- as.data.frame(selected_breakpoint_targets, stringsAsFactors = FALSE)
  chr_clean <- normalize_trascf_name(chr_name)
  target_chr_clean <- normalize_trascf_name(targets$chr)
  target_bp_index <- suppressWarnings(as.integer(targets$bp_index))

  keep <- target_chr_clean == chr_clean & target_bp_index == as.integer(bp_index)
  if (!any(keep, na.rm = TRUE)) {
    return(character(0))
  }

  dirs <- normalize_breakpoint_direction(targets$direction[keep])
  if (any(dirs == "both")) {
    dirs <- c(dirs[dirs != "both"], all_dirs)
  }

  unique(dirs[dirs %in% all_dirs])
}

is_selected_breakpoint_target <- function(chr_name, bp_index, selected_breakpoint_targets) {
  length(get_selected_breakpoint_directions(chr_name, bp_index, selected_breakpoint_targets)) > 0
}

filter_region_data_by_selected_directions <- function(region_data, selected_directions) {
  if (is.null(region_data) || nrow(region_data) == 0) {
    return(region_data)
  }

  selected_directions <- unique(normalize_breakpoint_direction(selected_directions))
  selected_directions <- selected_directions[selected_directions %in% c("left", "right")]

  if (length(selected_directions) == 0) {
    return(region_data[0, ])
  }

  if (all(c("left", "right") %in% selected_directions)) {
    return(region_data)
  }

  dt <- as.data.table(region_data)
  keep <- rep(FALSE, nrow(dt))
  if ("left" %in% selected_directions) {
    keep <- keep | dt$relative_pos <= 0
  }
  if ("right" %in% selected_directions) {
    keep <- keep | dt$relative_pos >= 0
  }

  dt[keep]
}

validate_breakpoint_direction_filter <- function(selected_breakpoint_targets, vertical_lines) {
  if (!is_breakpoint_direction_filter_enabled(selected_breakpoint_targets)) {
    return(invisible(TRUE))
  }

  targets <- tryCatch(
    as.data.frame(selected_breakpoint_targets, stringsAsFactors = FALSE),
    error = function(e) stop("selected_breakpoint_targets 必须是 data.frame 或可转换为 data.frame 的对象。")
  )

  required_cols <- c("chr", "bp_index", "direction")
  missing_cols <- setdiff(required_cols, names(targets))
  if (length(missing_cols) > 0) {
    stop("selected_breakpoint_targets 缺少列: ", paste(missing_cols, collapse = ", "))
  }

  target_bp_index <- suppressWarnings(as.integer(targets$bp_index))
  if (any(is.na(target_bp_index) | target_bp_index < 1)) {
    stop("selected_breakpoint_targets$bp_index 必须是大于等于 1 的整数。")
  }

  target_direction <- normalize_breakpoint_direction(targets$direction)
  invalid_direction <- setdiff(target_direction, c("left", "right", "both"))
  if (length(invalid_direction) > 0) {
    stop("selected_breakpoint_targets$direction 只能包含 'left'、'right' 或 'both'。")
  }

  target_chr_clean <- normalize_trascf_name(targets$chr)
  missing_chr <- setdiff(unique(target_chr_clean), names(vertical_lines))
  if (length(missing_chr) > 0) {
    stop("selected_breakpoint_targets 中以下染色体没有预设断点: ", paste(missing_chr, collapse = ", "))
  }

  for (i in seq_len(nrow(targets))) {
    chr_i <- target_chr_clean[i]
    bp_i <- target_bp_index[i]
    n_bp <- length(sort(as.numeric(vertical_lines[[chr_i]])))
    if (bp_i > n_bp) {
      stop(
        "selected_breakpoint_targets 中 ", targets$chr[i], " 的 bp_index = ", bp_i,
        " 超出该染色体断点数量。该染色体共有 ", n_bp, " 个断点。"
      )
    }
  }

  invisible(TRUE)
}

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
  tbi = "#2E86AB",
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

# 输出文件路径辅助函数
output_path <- function(filename) {
  file.path(output_dir, filename)
}

# 样本颜色辅助函数
sample_color_values <- function() {
  values <- c(publication_colors[[primary_sample]], publication_colors$tra)
  names(values) <- c(primary_label, tra_label)
  values
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

read_chromosomes_efficient <- function(filepath) {
  first_line <- readLines(filepath, n = 1)
  n_cols <- length(strsplit(first_line, "\t")[[1]])
  chr_col <- if (n_cols == 4) 1 else 5
  fread(filepath, sep = "\t", header = FALSE, select = chr_col)[[1]]
}

get_unique_chromosomes <- function(primary_file, tra_file) {
  chr_primary <- read_chromosomes_efficient(primary_file)
  chr_tra <- read_chromosomes_efficient(tra_file)
  all_chr <- unique(c(chr_primary, chr_tra))
  rm(chr_primary, chr_tra); gc()
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

# 检查 selected_breakpoint_targets 的 chr、bp_index 和 direction 是否有效。
validate_breakpoint_direction_filter(selected_breakpoint_targets, vertical_lines)

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

# 助手函数：按 breakpoint + 绝对距离窗口聚合
# 说明：
#   - 该函数专用于 P 值计算。
#   - 每个 breakpoint 在每个 abs_bin 中保留一个 mean_y。
#   - 后续在同一个 abs_bin 内，将多个 breakpoint 作为配对重复进行 Wilcoxon 检验。
aggregate_by_abs_bin_per_breakpoint <- function(df, window_size) {
  empty_bp_agg <- function() {
    data.frame(
      breakpoint_id = character(),
      pos = numeric(),
      mean_y = numeric(),
      n_points_in_breakpoint_bin = integer()
    )
  }

  if (is.null(df) || nrow(df) == 0) return(empty_bp_agg())
  if (!"breakpoint_id" %in% colnames(df)) {
    stop("用于 per-bin P 值计算的数据缺少 breakpoint_id 列。")
  }

  dt <- as.data.table(df)[!is.na(relative_pos) & !is.na(smooth_y) & !is.na(breakpoint_id)]
  if (nrow(dt) == 0) return(empty_bp_agg())

  dt[, abs_pos := abs(relative_pos)]
  dt[, abs_bin := round(abs_pos / window_size) * window_size]

  res <- dt[, .(
    pos = abs_bin[1],
    mean_y = mean(smooth_y, na.rm = TRUE),
    n_points_in_breakpoint_bin = .N
  ), by = .(breakpoint_id, abs_bin)]

  setorder(res, pos, breakpoint_id)
  as.data.frame(res[, .(breakpoint_id, pos, mean_y, n_points_in_breakpoint_bin)])
}

# 助手函数：计算绝对距离短板范围
compute_common_ranges <- function(bp_list) {
  if (length(bp_list) == 0) return(NULL)
  ranges <- lapply(bp_list, function(bp) {
    rel <- c()
    if (!is.null(bp$primary) && nrow(bp$primary) > 0) rel <- c(rel, bp$primary$relative_pos)
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
pvalue_method_label <- function(method) {
  method <- tolower(gsub("-", "_", method))
  if (method %in% c("ttest", "t.test")) method <- "t_test"
  if (method == "wilcox") return("paired Wilcoxon")
  if (method == "t_test") return("paired t-test")
  return(method)
}

calculate_paired_p_value <- function(x, y, alternative, method, min_pairs_for_test) {
  method <- tolower(gsub("-", "_", method))
  if (method %in% c("ttest", "t.test")) method <- "t_test"

  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  d <- x - y

  n_pairs <- length(d)
  n_nonzero <- sum(d != 0, na.rm = TRUE)

  if (n_pairs < min_pairs_for_test) {
    return(list(
      p_value = NA_real_,
      n_pairs = n_pairs,
      n_nonzero = n_nonzero,
      note = sprintf("not_tested: fewer than %d paired breakpoints", min_pairs_for_test)
    ))
  }

  if (method == "wilcox") {
    if (n_nonzero == 0) {
      return(list(
        p_value = NA_real_,
        n_pairs = n_pairs,
        n_nonzero = n_nonzero,
        note = "not_tested: all paired differences are zero"
      ))
    }

    p_val <- tryCatch({
      suppressWarnings(wilcox.test(x, y, paired = TRUE, alternative = alternative)$p.value)
    }, error = function(e) NA_real_)

    return(list(
      p_value = p_val,
      n_pairs = n_pairs,
      n_nonzero = n_nonzero,
      note = ifelse(is.na(p_val), "test_failed", "ok")
    ))
  }

  if (method == "t_test") {
    diff_sd <- suppressWarnings(sd(d, na.rm = TRUE))
    if (!is.finite(diff_sd) || diff_sd == 0) {
      return(list(
        p_value = NA_real_,
        n_pairs = n_pairs,
        n_nonzero = n_nonzero,
        note = "not_tested: paired differences have zero or undefined variance"
      ))
    }

    p_val <- tryCatch({
      suppressWarnings(t.test(x, y, paired = TRUE, alternative = alternative)$p.value)
    }, error = function(e) NA_real_)

    return(list(
      p_value = p_val,
      n_pairs = n_pairs,
      n_nonzero = n_nonzero,
      note = ifelse(is.na(p_val), "test_failed", "ok")
    ))
  }

  stop("未知的 p_value_method。")
}

calculate_stats <- function(df_primary, df_tra,
                            df_primary_bp, df_tra_bp,
                            threshold, label, output_filename,
                            min_pairs_for_test = 3,
                            p_value_method_arg = p_value_method) {
  # 根据绝对坐标 pos 匹配对应窗口的全局均值。
  # 注意：这部分仍保持原逻辑，用于输出 mean_y 与概率统计；绘图逻辑也不受影响。
  df_merged <- merge(df_primary[, c("pos", "mean_y")],
                     df_tra[, c("pos", "mean_y")],
                     by = "pos", suffixes = c(paste0("_", primary_suffix), paste0("_", tra_suffix)))

  primary_col <- paste0("mean_y_", primary_suffix)
  tra_col <- paste0("mean_y_", tra_suffix)

  p_value_method_arg <- tolower(gsub("-", "_", p_value_method_arg))
  if (p_value_method_arg %in% c("ttest", "t.test")) {
    p_value_method_arg <- "t_test"
  }
  if (!p_value_method_arg %in% c("wilcox", "t_test")) {
    stop("p_value_method_arg 只能设置为 'wilcox' 或 't_test'。")
  }
  selected_test_label <- pvalue_method_label(p_value_method_arg)

  df_in <- df_merged[df_merged$pos <= threshold, ]
  df_out <- df_merged[df_merged$pos > threshold, ]

  region_in_label <- sprintf("<= %sMb", threshold/1e6)
  region_out_label <- sprintf("> %sMb", threshold/1e6)

  cat(sprintf("\n================= %s 统计分析 =================\n", label))
  cat(sprintf("分析设定：阈值距离为 %s Mb；P 值按 breakpoint 配对重复计算；当前检验方法：%s\n",
              threshold/1e6, selected_test_label))
  cat("说明：\n")
  cat("  1) per-bin P 值：每个 abs bin 内，以 breakpoint_id 作为配对重复。\n")
  cat("  2) region P 值：每个 breakpoint 先在 <=threshold 或 >threshold 区域内求均值，再以 breakpoint_id 作为配对重复。\n")
  cat("  3) p_value_method 只改变 P 值检验方法，不改变进入检验的数值。\n")

  # ---------------- 每个 abs bin 的配对检验 ----------------
  # 新逻辑：在同一个 pos/abs_bin 内，用 breakpoint_id 进行 primary 与 TRA 配对，
  # 每个 pos 单独得到一个 P_value。
  df_bp_merged <- merge(df_primary_bp[, c("breakpoint_id", "pos", "mean_y")],
                        df_tra_bp[, c("breakpoint_id", "pos", "mean_y")],
                        by = c("breakpoint_id", "pos"),
                        suffixes = c(paste0("_", primary_suffix), paste0("_", tra_suffix)))

  primary_bp_col <- paste0("mean_y_", primary_suffix)
  tra_bp_col <- paste0("mean_y_", tra_suffix)

  if (nrow(df_bp_merged) > 0) {
    bp_dt <- as.data.table(df_bp_merged)
    pval_by_pos <- bp_dt[, {
      x <- get(primary_bp_col)
      y <- get(tra_bp_col)
      ok <- is.finite(x) & is.finite(y)
      x <- x[ok]
      y <- y[ok]
      d <- x - y
      n_pairs <- length(d)
      n_nonzero <- sum(d != 0, na.rm = TRUE)

      test_alt <- if (pos[1] <= threshold) "greater" else "less"
      # pos <= threshold: 检验 primary > TRA。
      # pos > threshold: 沿用原脚本方向，检验 TRA > primary；等价于 primary < TRA，alternative = "less"。
      test_res <- calculate_paired_p_value(
        x, y,
        alternative = test_alt,
        method = p_value_method_arg,
        min_pairs_for_test = min_pairs_for_test
      )
      p_val <- test_res$p_value

      .(
        N_pairs_for_test = n_pairs,
        N_nonzero_pairs = n_nonzero,
        Mean_primary_by_breakpoint = if (n_pairs > 0) mean(x, na.rm = TRUE) else NA_real_,
        Mean_TRA_by_breakpoint = if (n_pairs > 0) mean(y, na.rm = TRUE) else NA_real_,
        Mean_diff_primary_minus_TRA = if (n_pairs > 0) mean(d, na.rm = TRUE) else NA_real_,
        Prob_primary_higher_by_breakpoint = if (n_pairs > 0) mean(d > 0, na.rm = TRUE) else NA_real_,
        P_value_method = selected_test_label,
        P_value = p_val,
        P_value_note = test_res$note
      )
    }, by = pos]
    setorder(pval_by_pos, pos)
  } else {
    pval_by_pos <- data.table(
      pos = numeric(),
      N_pairs_for_test = integer(),
      N_nonzero_pairs = integer(),
      Mean_primary_by_breakpoint = numeric(),
      Mean_TRA_by_breakpoint = numeric(),
      Mean_diff_primary_minus_TRA = numeric(),
      Prob_primary_higher_by_breakpoint = numeric(),
      P_value_method = character(),
      P_value = numeric(),
      P_value_note = character()
    )
  }

  # ---------------- <=threshold 和 >threshold 整段的配对检验 ----------------
  # 这里仍然坚持“breakpoint 作为重复”的原则：
  #   Step 1: 对每个 breakpoint，在 region 内把多个 abs bin 的 mean_y 求平均。
  #   Step 2: 用每个 breakpoint 的 primary region mean 与 TRA region mean 做配对检验。
  if (nrow(df_bp_merged) > 0) {
    region_dt <- as.data.table(df_bp_merged)
    region_dt[, Region := ifelse(pos <= threshold, region_in_label, region_out_label)]

    bp_region_mean <- region_dt[, .(
      Primary_region_mean = mean(get(primary_bp_col), na.rm = TRUE),
      TRA_region_mean = mean(get(tra_bp_col), na.rm = TRUE),
      N_abs_bins_in_region_for_breakpoint = .N
    ), by = .(breakpoint_id, Region)]

    region_pvals <- bp_region_mean[, {
      x <- Primary_region_mean
      y <- TRA_region_mean
      ok <- is.finite(x) & is.finite(y)
      x <- x[ok]
      y <- y[ok]
      d <- x - y
      n_pairs <- length(d)
      n_nonzero <- sum(d != 0, na.rm = TRUE)

      test_alt <- if (Region[1] == region_in_label) "greater" else "less"
      # <=threshold: 检验 primary > TRA。
      # >threshold: 检验 TRA > primary，等价于 primary < TRA。
      test_res <- calculate_paired_p_value(
        x, y,
        alternative = test_alt,
        method = p_value_method_arg,
        min_pairs_for_test = min_pairs_for_test
      )
      p_val <- test_res$p_value

      .(
        Tested_Hypothesis = if (Region[1] == region_in_label) sprintf("%s > %s", primary_label, tra_label) else sprintf("%s > %s", tra_label, primary_label),
        N_breakpoint_pairs_region = n_pairs,
        N_nonzero_pairs_region = n_nonzero,
        Mean_primary_region_by_breakpoint = if (n_pairs > 0) mean(x, na.rm = TRUE) else NA_real_,
        Mean_TRA_region_by_breakpoint = if (n_pairs > 0) mean(y, na.rm = TRUE) else NA_real_,
        Mean_diff_primary_minus_TRA_region = if (n_pairs > 0) mean(d, na.rm = TRUE) else NA_real_,
        Prob_primary_higher_by_breakpoint_region = if (n_pairs > 0) mean(d > 0, na.rm = TRUE) else NA_real_,
        Median_abs_bins_per_breakpoint_region = if (n_pairs > 0) as.numeric(median(N_abs_bins_in_region_for_breakpoint, na.rm = TRUE)) else NA_real_,
        Region_P_value_method = selected_test_label,
        P_value_region = p_val,
        Region_P_value_note = test_res$note
      )
    }, by = Region]

    # 固定输出顺序：先 <=threshold，再 >threshold。
    region_pvals[, Region_order := ifelse(Region == region_in_label, 1, 2)]
    setorder(region_pvals, Region_order)
    region_pvals[, Region_order := NULL]
  } else {
    bp_region_mean <- data.table(
      breakpoint_id = character(),
      Region = character(),
      Primary_region_mean = numeric(),
      TRA_region_mean = numeric(),
      N_abs_bins_in_region_for_breakpoint = integer()
    )
    region_pvals <- data.table(
      Region = character(),
      Tested_Hypothesis = character(),
      N_breakpoint_pairs_region = integer(),
      N_nonzero_pairs_region = integer(),
      Mean_primary_region_by_breakpoint = numeric(),
      Mean_TRA_region_by_breakpoint = numeric(),
      Mean_diff_primary_minus_TRA_region = numeric(),
      Prob_primary_higher_by_breakpoint_region = numeric(),
      Median_abs_bins_per_breakpoint_region = numeric(),
      Region_P_value_method = character(),
      P_value_region = numeric(),
      Region_P_value_note = character()
    )
  }

  # ---------------- 控制台输出：per-bin P 值 + 整段 region P 值 ----------------
  if (nrow(df_in) > 0) {
    prob_primary_higher <- mean(df_in[[primary_col]] > df_in[[tra_col]], na.rm = TRUE)
    p_in <- pval_by_pos[pos <= threshold & !is.na(P_value)]
    region_in <- region_pvals[Region == region_in_label]

    cat(sprintf("\n[ 阈值内: 0 到 %sMb ]\n", threshold/1e6))
    cat(sprintf("  ▶ 基于全局窗口均值，%s 高于 %s 的窗口比例: %.2f%%\n", primary_label, tra_label, prob_primary_higher * 100))
    cat(sprintf("  ▶ per-bin %s：共 %d 个 abs bin 有可计算 P 值；每个 bin 使用 breakpoint 作为配对重复\n",
                selected_test_label, nrow(p_in)))
    if (nrow(p_in) > 0) {
      cat(sprintf("  ▶ per-bin P 值范围: min=%s, median=%s, max=%s\n",
                  format(min(p_in$P_value, na.rm = TRUE), scientific = TRUE, digits = 4),
                  format(median(p_in$P_value, na.rm = TRUE), scientific = TRUE, digits = 4),
                  format(max(p_in$P_value, na.rm = TRUE), scientific = TRUE, digits = 4)))
    }
    if (nrow(region_in) > 0) {
      cat(sprintf("  ▶ 整段 region %s：%s，breakpoint 配对数=%d，p-value=%s\n",
                  selected_test_label,
                  region_in$Tested_Hypothesis[1],
                  region_in$N_breakpoint_pairs_region[1],
                  format(region_in$P_value_region[1], scientific = TRUE, digits = 4)))
    }
  } else {
    cat(sprintf("\n[ 阈值内: 0 到 %sMb ] 无数据！\n", threshold/1e6))
  }

  if (nrow(df_out) > 0) {
    prob_tra_higher <- mean(df_out[[tra_col]] > df_out[[primary_col]], na.rm = TRUE)
    p_out <- pval_by_pos[pos > threshold & !is.na(P_value)]
    region_out <- region_pvals[Region == region_out_label]

    cat(sprintf("\n[ 阈值外: 大于 %sMb ]\n", threshold/1e6))
    cat(sprintf("  ▶ 基于全局窗口均值，%s 高于 %s 的窗口比例: %.2f%%\n", tra_label, primary_label, prob_tra_higher * 100))
    cat(sprintf("  ▶ per-bin %s：共 %d 个 abs bin 有可计算 P 值；每个 bin 使用 breakpoint 作为配对重复\n",
                selected_test_label, nrow(p_out)))
    if (nrow(p_out) > 0) {
      cat(sprintf("  ▶ per-bin P 值范围: min=%s, median=%s, max=%s\n",
                  format(min(p_out$P_value, na.rm = TRUE), scientific = TRUE, digits = 4),
                  format(median(p_out$P_value, na.rm = TRUE), scientific = TRUE, digits = 4),
                  format(max(p_out$P_value, na.rm = TRUE), scientific = TRUE, digits = 4)))
    }
    if (nrow(region_out) > 0) {
      cat(sprintf("  ▶ 整段 region %s：%s，breakpoint 配对数=%d，p-value=%s\n",
                  selected_test_label,
                  region_out$Tested_Hypothesis[1],
                  region_out$N_breakpoint_pairs_region[1],
                  format(region_out$P_value_region[1], scientific = TRUE, digits = 4)))
    }
  } else {
    cat(sprintf("\n[ 阈值外: 大于 %sMb ] 无数据！\n", threshold/1e6))
  }
  cat("=================================================================\n")

  # ---------------- 导出 per-bin 统计数据表格 ----------------
  df_merged$Region <- ifelse(df_merged$pos <= threshold, region_in_label, region_out_label)
  df_merged$Tested_Hypothesis <- ifelse(df_merged$pos <= threshold,
                                        sprintf("%s > %s", primary_label, tra_label),
                                        sprintf("%s > %s", tra_label, primary_label))

  df_merged <- merge(df_merged, as.data.frame(pval_by_pos), by = "pos", all.x = TRUE)

  # 把整段 region P 值也合并到每个 pos 行，方便在同一个 CSV 里查看。
  region_pvals_for_merge <- as.data.frame(region_pvals[, .(
    Region,
    Region_P_value_method = Region_P_value_method,
    Region_P_value = P_value_region,
    Region_P_value_note = Region_P_value_note,
    Region_N_breakpoint_pairs = N_breakpoint_pairs_region,
    Region_N_nonzero_pairs = N_nonzero_pairs_region
  )])
  df_merged <- merge(df_merged, region_pvals_for_merge, by = "Region", all.x = TRUE)

  # 重新排列列顺序。
  # P_value 是每个 pos/abs_bin 独立的 P 值；Region_P_value 是 <=threshold 或 >threshold 整段 P 值。
  df_merged <- df_merged[, c("pos", primary_col, tra_col,
                             "Region", "Tested_Hypothesis",
                             "N_pairs_for_test", "N_nonzero_pairs",
                             "Mean_primary_by_breakpoint", "Mean_TRA_by_breakpoint",
                             "Mean_diff_primary_minus_TRA",
                             "Prob_primary_higher_by_breakpoint",
                             "P_value_method", "P_value", "P_value_note",
                             "Region_N_breakpoint_pairs", "Region_N_nonzero_pairs",
                             "Region_P_value_method", "Region_P_value", "Region_P_value_note")]

  write.csv(df_merged, output_path(output_filename), row.names = FALSE, quote = FALSE)
  cat(sprintf("=> 已成功导出 per-bin 配套统计表格: %s\n", output_path(output_filename)))

  # ---------------- 单独导出整段 region P 值表格 ----------------
  region_output_filename <- sub("\\.csv$", "_Region_P_values.csv", output_filename)
  region_pvals_out <- as.data.frame(region_pvals)
  write.csv(region_pvals_out, output_path(region_output_filename), row.names = FALSE, quote = FALSE)
  cat(sprintf("=> 已成功导出整段 region P 值表格: %s\n", output_path(region_output_filename)))

  # 额外导出每个 breakpoint 的 region 均值，便于检查整段 P 值是由哪些配对数据得到的。
  bp_region_output_filename <- sub("\\.csv$", "_Breakpoint_Region_Means.csv", output_filename)
  bp_region_mean_out <- as.data.frame(bp_region_mean)
  write.csv(bp_region_mean_out, output_path(bp_region_output_filename), row.names = FALSE, quote = FALSE)
  cat(sprintf("=> 已成功导出 breakpoint-level region 均值表格: %s\n", output_path(bp_region_output_filename)))
}

# ================== 数据提取与处理主循环 ==================
cat(sprintf("当前 primary 样本: %s\n", primary_label))
cat(sprintf("primary 输入文件: %s\n", primary_bed_file))
cat(sprintf("TRA 输入文件: %s\n", tra_bed_file))
cat(sprintf("%s 归一化系数: %g\n", primary_label, x1))
cat(sprintf("%s 归一化系数: %g\n", tra_label, x2))
cat(sprintf("P 值计算方法: %s\n", pvalue_method_label(p_value_method)))
cat("selected_trascf: ",
    ifelse(is.null(selected_trascf) || length(selected_trascf) == 0,
           "ALL",
           paste(selected_trascf, collapse = ", ")),
    "\n", sep = "")
cat("selected_breakpoint_targets: ",
    format_selected_breakpoint_targets(selected_breakpoint_targets),
    "\n", sep = "")

all_chromosomes <- get_unique_chromosomes(primary_bed_file, tra_bed_file)

# ================== 按指定 TraScf 染色体过滤 ==================
# 只改变进入后续断点提取的染色体；后续统计和画图逻辑不改变。
all_chromosomes_original <- all_chromosomes
all_chromosomes <- all_chromosomes[
  vapply(
    all_chromosomes,
    is_selected_trascf,
    logical(1),
    selected_trascf = selected_trascf
  )
]

cat("\nTraScf 染色体筛选设置:\n")
cat("  筛选前染色体数量:", length(all_chromosomes_original), "\n")
cat("  筛选后染色体数量:", length(all_chromosomes), "\n")
cat("  实际进入统计和绘图的染色体:",
    ifelse(length(all_chromosomes) == 0, "NONE", paste(all_chromosomes, collapse = ", ")),
    "\n")

if (length(all_chromosomes) == 0) {
  stop("没有任何染色体匹配 selected_trascf。请检查输入文件中的染色体名称是否类似 TraScf_7 或 7。")
}

all_breakpoint_data <- list()

cat("开始读取数据并提取断点区域 (跳过单染色体绘图)...\n")
for (chr in all_chromosomes) {
  chr_fullname <- chr
  chr_num <- normalize_trascf_name(chr_fullname)

  # 仅读取必需数据以加速
  primary_data <- tryCatch({
    dt <- read_bed_efficient(primary_bed_file, primary_label)
    dt <- dt[V5 == chr_fullname]
    dt[, adjusted_y := ifelse(V8 >= 0, V8 / x1, NA_real_)]
    dt[!is.na(adjusted_y) & !is.na(V7) & is.finite(adjusted_y)]
  }, error = function(e) data.table())

  tra_data <- tryCatch({
    dt <- read_bed_efficient(tra_bed_file, tra_label)
    dt <- dt[V5 == chr_fullname]
    dt[, adjusted_y := ifelse(V8 >= 0, V8 / x2, NA_real_)]
    dt[!is.na(adjusted_y) & !is.na(V7) & is.finite(adjusted_y)]
  }, error = function(e) data.table())

  chr_data <- rbindlist(list(primary_data, tra_data), use.names = TRUE, fill = TRUE)
  if (nrow(chr_data) == 0) next

  # 使用 plot_window(0.5MB) 生成基础序列。
  # 注意：后续方向筛选发生在 relative_pos 生成之后，只改变进入 all_breakpoint_data 的行。
  plot_data_primary <- if(nrow(primary_data) > 0) apply_smoothing(create_plot_data(as.data.frame(primary_data), plot_window)) else data.frame()
  plot_data_tra <- if(nrow(tra_data) > 0) apply_smoothing(create_plot_data(as.data.frame(tra_data), plot_window)) else data.frame()

  if (chr_num %in% names(vertical_lines)) {
    # 断点编号按该 TraScf 上断点坐标从小到大排序后确定，与 plot3_selected_TraScf_BP_direction.R 保持一致。
    marker_positions <- sort(as.numeric(vertical_lines[[chr_num]]))
    min_pos <- min(chr_data$V7, na.rm = TRUE)
    max_pos <- max(chr_data$V7, na.rm = TRUE)

    for (bp_idx in seq_along(marker_positions)) {
      selected_directions <- get_selected_breakpoint_directions(
        chr_fullname,
        bp_idx,
        selected_breakpoint_targets
      )

      if (length(selected_directions) == 0) {
        cat(
          "跳过断点:", chr_fullname,
          "BP", bp_idx,
          "坐标", marker_positions[bp_idx],
          "（未在 selected_breakpoint_targets 中指定）\n"
        )
        next
      }

      bp_pos <- marker_positions[bp_idx]
      left_bound <- max(min_pos, bp_pos - summary_extension)
      right_bound <- min(max_pos, bp_pos + summary_extension)

      # 每个断点独立分析：
      # 不再使用相邻两个断点的中点截断 left_bound / right_bound。
      # 因此，即使同一条 TraScf 上有多个断点，每个断点仍按 bp_pos ± summary_extension 取范围；
      # 只受当前染色体实际 min_pos / max_pos 边界限制。
      # left_bound <- max(min_pos, bp_pos - summary_extension)
      # right_bound <- min(max_pos, bp_pos + summary_extension)

      cat(
        "进入分析的断点:", chr_fullname,
        "BP", bp_idx,
        "坐标", bp_pos,
        "方向", paste(selected_directions, collapse = ","),
        "\n"
      )

      region_data_primary <- plot_data_primary %>% filter(pos >= left_bound & pos <= right_bound)
      region_data_tra <- plot_data_tra %>% filter(pos >= left_bound & pos <= right_bound)

      if (nrow(region_data_primary) > 0 || nrow(region_data_tra) > 0) {
        # 为每个断点创建唯一 ID，用于后续在同一个 abs_bin 内做 primary 与 TRA 的配对检验。
        # 注意：这里仍以 breakpoint 为重复单位；方向只用于筛选该 breakpoint 的保留侧，不额外把方向当作重复。
        bp_id <- paste0(chr_fullname, "_bp", bp_idx, "_", bp_pos)

        if (nrow(region_data_primary) > 0) {
          region_data_primary$relative_pos <- region_data_primary$pos - bp_pos
          region_data_primary$breakpoint_id <- bp_id
          region_data_primary <- filter_region_data_by_selected_directions(
            region_data_primary,
            selected_directions
          )
        }
        if (nrow(region_data_tra) > 0) {
          region_data_tra$relative_pos <- region_data_tra$pos - bp_pos
          region_data_tra$breakpoint_id <- bp_id
          region_data_tra <- filter_region_data_by_selected_directions(
            region_data_tra,
            selected_directions
          )
        }

        if (nrow(region_data_primary) > 0 || nrow(region_data_tra) > 0) {
          all_breakpoint_data[[length(all_breakpoint_data) + 1]] <- list(
            primary = region_data_primary,
            tra = region_data_tra
          )
        } else {
          cat(
            "  该断点在指定方向筛选后没有可用窗口数据，跳过:",
            chr_fullname, "BP", bp_idx, "\n"
          )
        }
      }
    }
  }
  rm(chr_data, plot_data_primary, plot_data_tra); gc()
}

# ==================== 开始绝对距离聚合分析与绘图 ====================
cat("\n提取完成！开始进行绝对距离相关的绘图与统计检验...\n")

if (length(all_breakpoint_data) > 0) {
  all_primary_summary <- do.call(rbind, lapply(all_breakpoint_data, function(x) x$primary))
  all_tra_summary <- do.call(rbind, lapply(all_breakpoint_data, function(x) x$tra))

  # -------------------- 1. 全局绝对距离聚合 (Aggregated_Abs) --------------------

  # 用于高精度绘图的聚合 (0.5MB 窗口)
  agg_primary_abs_plot <- aggregate_by_abs_bin(all_primary_summary, plot_window)
  agg_tra_abs_plot <- aggregate_by_abs_bin(all_tra_summary, plot_window)

  # 用于统计检验的全局聚合均值 (2MB 扫描窗口)：保留原有均值/概率输出逻辑
  stats_primary_abs <- aggregate_by_abs_bin(all_primary_summary, scan_window)
  stats_tra_abs <- aggregate_by_abs_bin(all_tra_summary, scan_window)

  # 用于 P 值计算的 per-breakpoint 聚合 (2MB 扫描窗口)：每个 breakpoint 在每个 abs_bin 保留一个值
  stats_primary_abs_bp <- aggregate_by_abs_bin_per_breakpoint(all_primary_summary, scan_window)
  stats_tra_abs_bp <- aggregate_by_abs_bin_per_breakpoint(all_tra_summary, scan_window)

  if ((nrow(stats_primary_abs) + nrow(stats_tra_abs)) > 0) {
    # >>> 调用统计模块，P 值按每个 abs_bin 内的 breakpoint 配对重复计算，检验方法由 p_value_method 指定 <<<
    calculate_stats(stats_primary_abs, stats_tra_abs,
                    stats_primary_abs_bp, stats_tra_abs_bp,
                    threshold_dist,
                    "全局绝对距离汇总 (Aggregated_Abs)",
                    "Stats_All_Breakpoints_Summary_Aggregated_Abs.csv")

    # 绘图逻辑 (修复图注问题：使用 aes(color=) 和 aes(fill=) 正确映射)
    y_min_abs <- min(agg_primary_abs_plot$mean_y - agg_primary_abs_plot$se_y, agg_tra_abs_plot$mean_y - agg_tra_abs_plot$se_y, na.rm = TRUE)
    y_max_abs <- max(agg_primary_abs_plot$mean_y + agg_primary_abs_plot$se_y, agg_tra_abs_plot$mean_y + agg_tra_abs_plot$se_y, na.rm = TRUE)
    y_lim_abs <- c(max(0, y_min_abs - (y_max_abs-y_min_abs)*0.1), y_max_abs + (y_max_abs-y_min_abs)*0.1)

    p_abs <- ggplot() +
      create_publication_theme() +
      labs(title = "Summary Aggregated by Absolute Distance", subtitle = "With Statistical Threshold Highlight", x = "Distance to Breakpoint (Mb)", y = "RHO Value") +
      coord_cartesian(ylim = y_lim_abs)

    if(nrow(agg_primary_abs_plot)>0) {
      p_abs <- p_abs +
        geom_ribbon(data = agg_primary_abs_plot, aes(x = pos/1e6, ymin = mean_y-se_y, ymax = mean_y+se_y, fill = primary_label), alpha = 0.12) +
        geom_line(data = agg_primary_abs_plot, aes(x = pos/1e6, y = mean_y, color = primary_label), size = 0.8)
    }
    if(nrow(agg_tra_abs_plot)>0) {
      p_abs <- p_abs +
        geom_ribbon(data = agg_tra_abs_plot, aes(x = pos/1e6, ymin = mean_y-se_y, ymax = mean_y+se_y, fill = tra_label), alpha = 0.12) +
        geom_line(data = agg_tra_abs_plot, aes(x = pos/1e6, y = mean_y, color = tra_label), size = 0.8)
    }

    # 增加断点 0 处标记、20Mb 处阈值标记，并定义正确的颜色图注
    p_abs <- p_abs +
      geom_vline(xintercept = 0, color = publication_colors$marker, linetype = "dashed", size = 1.0) +
      geom_vline(xintercept = threshold_dist/1e6, color = publication_colors$threshold, linetype = "longdash", size = 0.8) +
      annotate("text", x = threshold_dist/1e6, y = Inf, label = paste0(threshold_dist/1e6, "Mb Threshold"),
               color = publication_colors$threshold, vjust = 1.5, hjust = -0.1, size = 4, fontface = "bold") +
      scale_y_continuous(labels = scales::scientific_format(digits = 2)) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      scale_color_manual(name = "Sample", values = sample_color_values()) +
      scale_fill_manual(name = "Sample", values = sample_color_values())

    ggsave(output_path("All_Breakpoints_Summary_Aggregated_Abs.png"), plot = p_abs, width = 14, height = 8, dpi = 600, bg = "white")
    ggsave(output_path("All_Breakpoints_Summary_Aggregated_Abs.pdf"), plot = p_abs, width = 14, height = 8, device = "pdf", bg = "white")
  }

  # -------------------- 2. 公共边界绝对距离聚合 (Aggregated_Abs_Common) --------------------
  common_ranges <- compute_common_ranges(all_breakpoint_data)
  if (!is.null(common_ranges) && is.finite(common_ranges$max_abs) && common_ranges$max_abs > 0) {
    all_primary_abs_common <- all_primary_summary[abs(all_primary_summary$relative_pos) <= common_ranges$max_abs, ]
    all_tra_abs_common <- all_tra_summary[abs(all_tra_summary$relative_pos) <= common_ranges$max_abs, ]

    # 用于高精度绘图的聚合 (0.5MB 窗口)
    agg_primary_abs_common_plot <- aggregate_by_abs_bin(all_primary_abs_common, plot_window)
    agg_tra_abs_common_plot <- aggregate_by_abs_bin(all_tra_abs_common, plot_window)

    # 用于统计检验的全局聚合均值 (2MB 扫描窗口)：保留原有均值/概率输出逻辑
    stats_primary_abs_common <- aggregate_by_abs_bin(all_primary_abs_common, scan_window)
    stats_tra_abs_common <- aggregate_by_abs_bin(all_tra_abs_common, scan_window)

    # 用于 P 值计算的 per-breakpoint 聚合 (2MB 扫描窗口)：每个 breakpoint 在每个 abs_bin 保留一个值
    stats_primary_abs_common_bp <- aggregate_by_abs_bin_per_breakpoint(all_primary_abs_common, scan_window)
    stats_tra_abs_common_bp <- aggregate_by_abs_bin_per_breakpoint(all_tra_abs_common, scan_window)

    if ((nrow(stats_primary_abs_common) + nrow(stats_tra_abs_common)) > 0) {
      # >>> 调用统计模块，P 值按每个 abs_bin 内的 breakpoint 配对重复计算，检验方法由 p_value_method 指定 <<<
      calculate_stats(stats_primary_abs_common, stats_tra_abs_common,
                      stats_primary_abs_common_bp, stats_tra_abs_common_bp,
                      threshold_dist,
                      "公共边界绝对距离汇总 (Aggregated_Abs_Common)",
                      "Stats_All_Breakpoints_Summary_Aggregated_Abs_Common.csv")

      # 绘图逻辑
      y_min_comm <- min(agg_primary_abs_common_plot$mean_y - agg_primary_abs_common_plot$se_y, agg_tra_abs_common_plot$mean_y - agg_tra_abs_common_plot$se_y, na.rm = TRUE)
      y_max_comm <- max(agg_primary_abs_common_plot$mean_y + agg_primary_abs_common_plot$se_y, agg_tra_abs_common_plot$mean_y + agg_tra_abs_common_plot$se_y, na.rm = TRUE)
      y_lim_comm <- c(max(0, y_min_comm - (y_max_comm-y_min_comm)*0.1), y_max_comm + (y_max_comm-y_min_comm)*0.1)

      p_abs_common <- ggplot() +
        create_publication_theme() +
        labs(title = "Summary Aggregated by Absolute Distance (Common Range)", subtitle = "With Statistical Threshold Highlight", x = "Distance to Breakpoint (Mb)", y = "RHO Value") +
        coord_cartesian(ylim = y_lim_comm)

      if(nrow(agg_primary_abs_common_plot)>0) {
        p_abs_common <- p_abs_common +
          geom_ribbon(data = agg_primary_abs_common_plot, aes(x = pos/1e6, ymin = mean_y-se_y, ymax = mean_y+se_y, fill = primary_label), alpha = 0.12) +
          geom_line(data = agg_primary_abs_common_plot, aes(x = pos/1e6, y = mean_y, color = primary_label), size = 0.8)
      }
      if(nrow(agg_tra_abs_common_plot)>0) {
        p_abs_common <- p_abs_common +
          geom_ribbon(data = agg_tra_abs_common_plot, aes(x = pos/1e6, ymin = mean_y-se_y, ymax = mean_y+se_y, fill = tra_label), alpha = 0.12) +
          geom_line(data = agg_tra_abs_common_plot, aes(x = pos/1e6, y = mean_y, color = tra_label), size = 0.8)
      }

      p_abs_common <- p_abs_common +
        geom_vline(xintercept = 0, color = publication_colors$marker, linetype = "dashed", size = 1.0) +
        geom_vline(xintercept = threshold_dist/1e6, color = publication_colors$threshold, linetype = "longdash", size = 0.8) +
        annotate("text", x = threshold_dist/1e6, y = Inf, label = paste0(threshold_dist/1e6, "Mb Threshold"),
                 color = publication_colors$threshold, vjust = 1.5, hjust = -0.1, size = 4, fontface = "bold") +
        scale_y_continuous(labels = scales::scientific_format(digits = 2)) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
        scale_color_manual(name = "Sample", values = sample_color_values()) +
        scale_fill_manual(name = "Sample", values = sample_color_values())

      ggsave(output_path("All_Breakpoints_Summary_Aggregated_Abs_Common.png"), plot = p_abs_common, width = 14, height = 8, dpi = 600, bg = "white")
      ggsave(output_path("All_Breakpoints_Summary_Aggregated_Abs_Common.pdf"), plot = p_abs_common, width = 14, height = 8, device = "pdf", bg = "white")
    }
  }
} else {
  cat("未检测到有效断点数据。\n")
}

cat("\n所有任务执行完毕。\n")
