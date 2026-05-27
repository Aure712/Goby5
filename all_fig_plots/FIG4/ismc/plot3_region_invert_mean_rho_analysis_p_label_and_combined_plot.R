library(dplyr)
library(ggplot2)
library(scales)
library(data.table)

# ============================================================
# 0. 参数设置区：运行前主要修改这里
# ============================================================

# ---- 原始4个BED/ISMC输入文件，保持不变 ----
input_files <- list(
  tba_bed          = "/fast3/group_crf/home/g21shaoy23/goby5/ismc/tba/20inds/trans/output2.bed",
  tra_for_tba_bed = "/fast3/group_crf/home/g21shaoy23/goby5/ismc/tba/20inds/trans/tra/new_output2.bed",
  tbi_bed          = "/public4/group_crf/home/g21shaoy23/ismc/tbi2tra/output2.bed",
  tra_for_tbi_bed = "/public4/group_crf/home/g21shaoy23/ismc/tbi2tra/tra/new_output2.bed"
)

# ---- 新增：invert BED文件；这两个BED中的坐标都是TRA坐标 ----
invert_files <- list(
  invert_tba = "invert_tba.bed",
  invert_tbi = "invert_tbi.bed"
)

# ---- 每个样本/群体对应的x值：RHO = BED数值 / x ----
x_values <- list(
  tba = 0.007300749,
  tra = 0.005889914,
  tbi = 0.011549063
)

# ---- 指定要统计的TRA染色体区间 ----
# 每一行是一个目标区间；同一染色体可以写多行。
# chr名称需要和tra_for_tba_bed、tra_for_tbi_bed、invert_tba、invert_tbi中的染色体名称一致。
target_tra_regions <- data.table(
  chr = c("TraScf_1","TraScf_1","TraScf_2","TraScf_2","TraScf_3","TraScf_3","TraScf_4","TraScf_4","TraScf_5","TraScf_5","TraScf_6","TraScf_6","TraScf_7","TraScf_7","TraScf_8","TraScf_8","TraScf_9","TraScf_9","TraScf_10","TraScf_10"),
  start = c(0,106968322,0,81585602,0,82758597,0,74874537,0,70273994,0,68458668,0,67720776,0,65834950,0,67186786,0,62435997),
  end = c(13497514,120400279,15218063,99551308,17981893,97165056,14811769,90094152,14121924,84228105,12988143,83205787,13300338,81630657,12172625,80493487,13787576,80098815,11220600,76823096),
  region_id = c("TraScf_1_L1","TraScf_1_R1","TraScf_2_L1","TraScf_2_R1","TraScf_3_L1","TraScf_3_R1","TraScf_4_L1","TraScf_4_R1","TraScf_5_L1","TraScf_5_R1","TraScf_6_L1","TraScf_6_R1","TraScf_7_L1","TraScf_7_R1","TraScf_8_L1","TraScf_8_R1","TraScf_9_L1","TraScf_9_R1","TraScf_10_L1","TraScf_10_R1")
)

# ---- 输出设置 ----
output_dir <- "region_invert_rho_results"
output_prefix <- "region_invert_rho"
save_png <- TRUE
save_pdf <- TRUE
save_tiff <- TRUE

# ---- 低内存设置 ----
# TRUE：使用awk按染色体从BED/ISMC中筛选行，只把当前染色体读入R。
use_awk_chr_filter <- TRUE
data_table_threads <- 1
setDTthreads(data_table_threads)
run_gc_each_chr <- TRUE

# ---- 显著性检验设置 ----
# within_source：同一Source内部 perfect vs background，默认非配对Wilcoxon。
# diff_of_diff：比较 TRA(perfect-background) 与 REF(perfect-background) 的差异，按target_region_id配对。
significance_test_method <- "wilcox"   # 可选 "wilcox" 或 "t.test"
significance_p_adjust_method <- "BH"
# 控制图上显著性括号标注使用哪一种P值：
#   "p_value"：标注原始P值；
#   "p_adj"  ：标注经过 significance_p_adjust_method 校正后的P值。
# 注意：输出表格中会同时保留 p_value 和 p_adj 两列。
plot_p_value_type <- "p_value"
min_values_for_test <- 2

# ---- 图形设置 ----
bar_point_alpha <- 0.65
bar_point_size <- 2.0
bar_width <- 0.65

publication_colors <- list(
  tba = "#2E86AB",
  tra = "#A23B72",
  tbi = "#F18F01",
  perfect = "#A23B72",
  background = "#9CA3AF",
  text = "#2D3748",
  bracket = "#2D3748",
  background_panel = "#FAFAFA"
)

batch_configs <- list(
  TBA_batch = list(
    batch = "TBA_batch",
    ref_source = "TBA",
    ref_file_key = "tba_bed",
    tra_file_key = "tra_for_tba_bed",
    invert_key = "invert_tba",
    ref_x_key = "tba",
    tra_x_key = "tra"
  ),
  TBI_batch = list(
    batch = "TBI_batch",
    ref_source = "TBI",
    ref_file_key = "tbi_bed",
    tra_file_key = "tra_for_tbi_bed",
    invert_key = "invert_tbi",
    ref_x_key = "tbi",
    tra_x_key = "tra"
  )
)

# ============================================================
# 1. 通用辅助函数
# ============================================================

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

create_publication_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.background = element_rect(fill = publication_colors$background_panel, color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "white", linewidth = 0.8),
      panel.grid.minor = element_line(color = "white", linewidth = 0.4),
      axis.line = element_line(color = publication_colors$text, linewidth = 0.6),
      axis.ticks = element_line(color = publication_colors$text, linewidth = 0.5),
      axis.text = element_text(color = publication_colors$text, size = 10),
      axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
      axis.title = element_text(color = publication_colors$text, size = 12, face = "bold"),
      plot.title = element_text(color = publication_colors$text, size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(color = publication_colors$text, size = 11, hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      plot.margin = margin(20, 25, 20, 25)
    )
}

check_file_exists <- function(filepath, label) {
  if (!file.exists(filepath)) {
    stop(paste0("找不到输入文件: ", label, " = ", filepath))
  }
}

check_config <- function() {
  for (nm in names(input_files)) check_file_exists(input_files[[nm]], nm)
  for (nm in names(invert_files)) check_file_exists(invert_files[[nm]], nm)

  for (nm in names(x_values)) {
    if (!is.numeric(x_values[[nm]]) || length(x_values[[nm]]) != 1 ||
        is.na(x_values[[nm]]) || x_values[[nm]] == 0) {
      stop(paste0("x_values$", nm, " 必须是非零数值"))
    }
  }

  required_cols <- c("chr", "start", "end", "region_id")
  missing_cols <- setdiff(required_cols, names(target_tra_regions))
  if (length(missing_cols) > 0) {
    stop("target_tra_regions 缺少列: ", paste(missing_cols, collapse = ", "))
  }

  target_tra_regions[, `:=`(
    chr = as.character(chr),
    start = as.numeric(start),
    end = as.numeric(end),
    region_id = as.character(region_id)
  )]

  if (any(is.na(target_tra_regions$chr)) || any(target_tra_regions$chr == "")) {
    stop("target_tra_regions$chr 不能有NA或空值")
  }
  if (any(!is.finite(target_tra_regions$start)) || any(!is.finite(target_tra_regions$end))) {
    stop("target_tra_regions$start/end 必须是有效数字")
  }
  if (any(target_tra_regions$end <= target_tra_regions$start)) {
    stop("target_tra_regions 中存在 end <= start 的区间")
  }
  if (any(is.na(target_tra_regions$region_id)) || any(target_tra_regions$region_id == "")) {
    stop("target_tra_regions$region_id 不能有NA或空值")
  }
  if (anyDuplicated(target_tra_regions$region_id) > 0) {
    stop("target_tra_regions$region_id 不能重复")
  }

  if (!significance_test_method %in% c("wilcox", "t.test")) {
    stop("significance_test_method 只能是 'wilcox' 或 't.test'")
  }
  if (!plot_p_value_type %in% c("p_value", "p_adj")) {
    stop("plot_p_value_type 只能是 'p_value' 或 'p_adj'")
  }
}

read_first_line <- function(filepath) {
  con <- NULL
  on.exit({
    if (!is.null(con)) close(con)
  }, add = TRUE)

  if (grepl("\\.gz$", filepath, ignore.case = TRUE)) {
    con <- gzfile(filepath, open = "rt")
  } else {
    con <- file(filepath, open = "rt")
  }
  line <- readLines(con, n = 1, warn = FALSE)
  if (length(line) == 0) stop(paste0("文件为空: ", filepath))
  line
}

detect_bed_format <- function(filepath, label) {
  first_line <- read_first_line(filepath)
  n_cols <- length(strsplit(first_line, "\t", fixed = TRUE)[[1]])
  if (n_cols == 4) {
    list(label = label, path = filepath, n_cols = n_cols, format = "bed4", chr_col = 1)
  } else if (n_cols >= 8) {
    list(label = label, path = filepath, n_cols = n_cols, format = "ismc8", chr_col = 5)
  } else {
    stop(paste0(
      "文件 ", label, " = ", filepath, " 的列数为 ", n_cols,
      "。当前脚本只支持4列BED，或至少8列且第5/7/8列分别为chr/pos/value的文件。"
    ))
  }
}

make_reader_cmd <- function(filepath) {
  if (grepl("\\.gz$", filepath, ignore.case = TRUE)) {
    paste("gzip -dc", shQuote(filepath))
  } else {
    paste("cat", shQuote(filepath))
  }
}

read_chr_bed <- function(file_info, chr_name, source_name, x_value, file_id) {
  if (use_awk_chr_filter) {
    if (file_info$format == "bed4") {
      awk_prog <- "BEGIN{FS=OFS=\"\\t\"} $1==C {print $1,$2,$3,$4}"
    } else {
      awk_prog <- "BEGIN{FS=OFS=\"\\t\"} $5==C {print $5,$7,$8}"
    }
    cmd <- paste0(
      make_reader_cmd(file_info$path),
      " | awk -v C=", shQuote(chr_name), " '", awk_prog, "'"
    )
    raw_dt <- tryCatch(
      fread(cmd = cmd, sep = "\t", header = FALSE, showProgress = FALSE),
      error = function(e) data.table()
    )
  } else {
    if (file_info$format == "bed4") {
      raw_dt <- fread(file_info$path, header = FALSE, select = 1:4, showProgress = FALSE)
      raw_dt <- raw_dt[as.character(V1) == chr_name]
    } else {
      raw_dt <- fread(file_info$path, header = FALSE, select = c(5, 7, 8), showProgress = FALSE)
      setnames(raw_dt, c("V1", "V2", "V3"))
      raw_dt <- raw_dt[as.character(V1) == chr_name]
    }
  }

  if (nrow(raw_dt) == 0) {
    return(data.table(chr = character(), pos = numeric(), adjusted_y = numeric(), source = character(), file_id = character()))
  }

  if (file_info$format == "bed4") {
    if (ncol(raw_dt) < 4) stop(paste0("读取 ", file_info$label, " 时列数不足。"))
    setnames(raw_dt, c("V1", "V2", "V3", "V4"))
    dt <- raw_dt[, .(
      chr = as.character(V1),
      pos = (as.numeric(V2) + as.numeric(V3)) / 2,
      raw_y = as.numeric(V4)
    )]
  } else {
    if (ncol(raw_dt) < 3) stop(paste0("读取 ", file_info$label, " 时列数不足。"))
    setnames(raw_dt, c("V1", "V2", "V3"))
    dt <- raw_dt[, .(
      chr = as.character(V1),
      pos = as.numeric(V2),
      raw_y = as.numeric(V3)
    )]
  }

  rm(raw_dt)

  dt <- dt[!is.na(chr) & chr == chr_name & !is.na(pos) & is.finite(pos)]
  dt[, adjusted_y := ifelse(!is.na(raw_y) & raw_y >= 0, raw_y / x_value, NA_real_)]
  dt <- dt[!is.na(adjusted_y) & is.finite(adjusted_y)]
  dt[, `:=`(source = source_name, file_id = file_id)]
  dt[, raw_y := NULL]
  dt[]
}

read_invert_bed <- function(filepath, label) {
  dt <- fread(filepath, header = FALSE, showProgress = FALSE)
  if (ncol(dt) < 3) {
    stop(paste0(label, " 至少需要3列: chr, start, end"))
  }
  dt <- dt[, .(
    chr = as.character(V1),
    start = as.numeric(V2),
    end = as.numeric(V3)
  )]
  dt <- dt[!is.na(chr) & chr != "" & is.finite(start) & is.finite(end) & end > start]
  dt[, invert_id := paste0(label, "_INV", sprintf("%06d", seq_len(.N)))]
  setorder(dt, chr, start, end)
  dt[]
}

merge_segments <- function(seg_dt) {
  if (is.null(seg_dt) || nrow(seg_dt) == 0) {
    return(data.table(start = numeric(), end = numeric()))
  }
  dt <- copy(as.data.table(seg_dt))[, .(start = as.numeric(start), end = as.numeric(end))]
  dt <- dt[is.finite(start) & is.finite(end) & end > start]
  if (nrow(dt) == 0) return(data.table(start = numeric(), end = numeric()))
  setorder(dt, start, end)

  starts <- numeric()
  ends <- numeric()
  cur_start <- dt$start[1]
  cur_end <- dt$end[1]

  if (nrow(dt) >= 2) {
    for (i in 2:nrow(dt)) {
      if (dt$start[i] <= cur_end) {
        cur_end <- max(cur_end, dt$end[i])
      } else {
        starts <- c(starts, cur_start)
        ends <- c(ends, cur_end)
        cur_start <- dt$start[i]
        cur_end <- dt$end[i]
      }
    }
  }
  starts <- c(starts, cur_start)
  ends <- c(ends, cur_end)

  data.table(start = starts, end = ends)
}

subtract_segments_from_target <- function(target_start, target_end, remove_dt) {
  remove_dt <- merge_segments(remove_dt)
  if (nrow(remove_dt) == 0) {
    return(data.table(start = target_start, end = target_end))
  }

  out <- data.table(start = numeric(), end = numeric())
  cur <- target_start

  for (i in seq_len(nrow(remove_dt))) {
    rs <- max(remove_dt$start[i], target_start)
    re <- min(remove_dt$end[i], target_end)
    if (!is.finite(rs) || !is.finite(re) || re <= rs) next

    if (rs > cur) {
      out <- rbindlist(list(out, data.table(start = cur, end = rs)), use.names = TRUE)
    }
    cur <- max(cur, re)
  }

  if (cur < target_end) {
    out <- rbindlist(list(out, data.table(start = cur, end = target_end)), use.names = TRUE)
  }

  out <- out[end > start]
  out[]
}

make_analysis_intervals <- function(batch_name, invert_dt, targets_dt) {
  intervals <- data.table()
  excluded_partial <- data.table()

  targets <- copy(targets_dt)
  targets[, `:=`(
    chr = as.character(chr),
    start = as.numeric(start),
    end = as.numeric(end),
    region_id = as.character(region_id)
  )]

  for (i in seq_len(nrow(targets))) {
    tr <- targets[i]
    inv_chr <- invert_dt[chr == tr$chr]

    overlap <- inv_chr[start < tr$end & end > tr$start]
    perfect <- overlap[start >= tr$start & end <= tr$end]

    partial <- overlap[!(start >= tr$start & end <= tr$end)]
    if (nrow(partial) > 0) {
      partial_out <- copy(partial)
      partial_out[, `:=`(
        batch = batch_name,
        target_region_id = tr$region_id,
        target_chr = tr$chr,
        target_start = tr$start,
        target_end = tr$end,
        overlap_start = pmax(start, tr$start),
        overlap_end = pmin(end, tr$end)
      )]
      partial_out <- partial_out[overlap_end > overlap_start]
      excluded_partial <- rbindlist(list(excluded_partial, partial_out), use.names = TRUE, fill = TRUE)
    }

    if (nrow(perfect) > 0) {
      perfect_out <- copy(perfect)
      perfect_out[, `:=`(
        batch = batch_name,
        target_region_id = tr$region_id,
        target_chr = tr$chr,
        target_start = tr$start,
        target_end = tr$end,
        interval_type = "perfect_invert_inside_target",
        interval_source = "invert_fully_inside_target"
      )]
      perfect_out[, interval_local_index := seq_len(.N)]
      perfect_out[, interval_id := paste(batch, target_region_id, "perfect", interval_local_index, sep = "__")]
      intervals <- rbindlist(list(intervals, perfect_out), use.names = TRUE, fill = TRUE)
    }

    remove_segments <- overlap[, .(
      start = pmax(start, tr$start),
      end = pmin(end, tr$end)
    )]
    remove_segments <- remove_segments[end > start]

    background <- subtract_segments_from_target(tr$start, tr$end, remove_segments)
    if (nrow(background) > 0) {
      background[, `:=`(
        chr = tr$chr,
        batch = batch_name,
        target_region_id = tr$region_id,
        target_chr = tr$chr,
        target_start = tr$start,
        target_end = tr$end,
        interval_type = "target_minus_all_invert_overlap",
        interval_source = "target_after_removing_all_invert_overlap",
        invert_id = NA_character_
      )]
      background[, interval_local_index := seq_len(.N)]
      background[, interval_id := paste(batch, target_region_id, "background", interval_local_index, sep = "__")]
      intervals <- rbindlist(list(intervals, background), use.names = TRUE, fill = TRUE)
    }
  }

  if (nrow(intervals) > 0) {
    intervals[, `:=`(
      chr = as.character(chr),
      start = as.numeric(start),
      end = as.numeric(end),
      length = as.numeric(end - start),
      interval_type = factor(
        interval_type,
        levels = c("perfect_invert_inside_target", "target_minus_all_invert_overlap")
      )
    )]
    setcolorder(intervals, c(
      "batch", "target_region_id", "chr", "start", "end", "length",
      "interval_type", "interval_id", "interval_source", "invert_id",
      "target_chr", "target_start", "target_end"
    ))
    setorder(intervals, batch, chr, start, end, interval_type)
  }

  list(intervals = intervals, excluded_partial = excluded_partial)
}

compute_interval_means_for_source <- function(intervals_chr, data_chr, source_name) {
  intervals_chr <- copy(as.data.table(intervals_chr))
  if (nrow(intervals_chr) == 0) return(data.table())

  # 关键修复点：
  # 不再通过merge后猜测 n_points.x / n_points.y 这样的后缀列名；
  # 直接逐区间显式生成 mean_rho 与 n_points，避免 object 'n_points.y' not found。
  res <- vector("list", nrow(intervals_chr))

  for (i in seq_len(nrow(intervals_chr))) {
    int <- intervals_chr[i]
    vals <- data_chr[pos >= int$start & pos <= int$end, adjusted_y]
    vals <- vals[!is.na(vals) & is.finite(vals)]

    res[[i]] <- data.table(
      batch = int$batch,
      target_region_id = int$target_region_id,
      chr = int$chr,
      start = int$start,
      end = int$end,
      length = int$length,
      interval_type = as.character(int$interval_type),
      interval_id = int$interval_id,
      source = source_name,
      mean_rho = if (length(vals) == 0) NA_real_ else mean(vals),
      n_points = as.integer(length(vals))
    )
  }

  rbindlist(res, use.names = TRUE, fill = TRUE)
}

compute_all_interval_means_for_batch <- function(batch_cfg, file_infos, intervals) {
  batch_name <- batch_cfg$batch
  ref_source <- batch_cfg$ref_source

  if (is.null(intervals) || nrow(intervals) == 0) {
    return(data.table())
  }

  chr_vec <- sort(unique(intervals$chr))
  all_means <- list()

  for (chr_name in chr_vec) {
    cat("  - ", batch_name, ": reading chromosome ", chr_name, "\n", sep = "")
    intervals_chr <- intervals[chr == chr_name]
    if (nrow(intervals_chr) == 0) next

    tra_dt <- read_chr_bed(
      file_infos[[batch_cfg$tra_file_key]],
      chr_name,
      "TRA",
      x_values[[batch_cfg$tra_x_key]],
      batch_cfg$tra_file_key
    )

    ref_dt <- read_chr_bed(
      file_infos[[batch_cfg$ref_file_key]],
      chr_name,
      ref_source,
      x_values[[batch_cfg$ref_x_key]],
      batch_cfg$ref_file_key
    )

    tra_means <- compute_interval_means_for_source(intervals_chr, tra_dt, "TRA")
    ref_means <- compute_interval_means_for_source(intervals_chr, ref_dt, ref_source)

    all_means[[length(all_means) + 1]] <- rbindlist(
      list(tra_means, ref_means),
      use.names = TRUE,
      fill = TRUE
    )

    rm(tra_dt, ref_dt, tra_means, ref_means)
    if (run_gc_each_chr) gc()
  }

  out <- rbindlist(all_means, use.names = TRUE, fill = TRUE)
  if (nrow(out) > 0) {
    out[, interval_type := factor(
      interval_type,
      levels = c("perfect_invert_inside_target", "target_minus_all_invert_overlap")
    )]
    setorder(out, batch, source, interval_type, chr, start, end)
  }
  out[]
}

summarize_categories <- function(mean_dt) {
  if (is.null(mean_dt) || nrow(mean_dt) == 0) return(data.table())
  dt <- as.data.table(mean_dt)[!is.na(mean_rho) & is.finite(mean_rho)]
  if (nrow(dt) == 0) return(data.table())

  dt[, .(
    n_intervals = .N,
    n_intervals_with_points = sum(n_points > 0),
    total_points = sum(n_points, na.rm = TRUE),
    mean_rho = mean(mean_rho, na.rm = TRUE),
    median_rho = median(mean_rho, na.rm = TRUE),
    sd_rho = if (.N <= 1) NA_real_ else sd(mean_rho, na.rm = TRUE),
    se_rho = if (.N <= 1) NA_real_ else sd(mean_rho, na.rm = TRUE) / sqrt(.N),
    min_rho = min(mean_rho, na.rm = TRUE),
    max_rho = max(mean_rho, na.rm = TRUE)
  ), by = .(batch, source, interval_type)]
}

safe_p_label <- function(p, prefix = "p") {
  if (is.na(p) || !is.finite(p)) return(paste0(prefix, "=NA"))
  if (p < 2.2e-16) return(paste0(prefix, "<2.2e-16"))
  paste0(prefix, "=", formatC(p, format = "e", digits = 2))
}

p_to_stars <- function(p) {
  if (is.na(p) || !is.finite(p)) return("NA")
  if (p <= 0.001) return("***")
  if (p <= 0.01) return("**")
  if (p <= 0.05) return("*")
  "ns"
}

plot_p_label_prefix <- function() {
  if (plot_p_value_type == "p_adj") {
    return("p adj")
  }
  "p value"
}

run_two_group_test <- function(x, y, paired = FALSE, method = significance_test_method) {
  x <- x[!is.na(x) & is.finite(x)]
  y <- y[!is.na(y) & is.finite(y)]

  if (paired && length(x) != length(y)) {
    stop("paired test requires equal length vectors after alignment")
  }

  if (length(x) < min_values_for_test || length(y) < min_values_for_test) {
    return(list(p_value = NA_real_, statistic = NA_real_, method = method, note = "not_enough_values"))
  }

  out <- tryCatch({
    if (method == "wilcox") {
      test <- suppressWarnings(wilcox.test(x, y, paired = paired, exact = FALSE))
    } else {
      test <- suppressWarnings(t.test(x, y, paired = paired))
    }
    list(
      p_value = as.numeric(test$p.value),
      statistic = unname(as.numeric(test$statistic[1])),
      method = method,
      note = "ok"
    )
  }, error = function(e) {
    list(p_value = NA_real_, statistic = NA_real_, method = method, note = paste0("test_error: ", e$message))
  })
  out
}

compute_significance_tests_for_batch <- function(mean_dt, batch_name, ref_source) {
  dt <- as.data.table(mean_dt)[batch == batch_name & !is.na(mean_rho) & is.finite(mean_rho)]
  if (nrow(dt) == 0) return(data.table())

  perfect_type <- "perfect_invert_inside_target"
  bg_type <- "target_minus_all_invert_overlap"

  tests <- list()

  # 1) 同一Source内部：perfect vs background
  for (src in c("TRA", ref_source)) {
    x <- dt[source == src & interval_type == perfect_type, mean_rho]
    y <- dt[source == src & interval_type == bg_type, mean_rho]
    tt <- run_two_group_test(x, y, paired = FALSE)

    tests[[length(tests) + 1]] <- data.table(
      batch = batch_name,
      comparison_type = "within_source",
      comparison = paste0(src, ": perfect vs background"),
      source_1 = src,
      interval_type_1 = perfect_type,
      source_2 = src,
      interval_type_2 = bg_type,
      paired = FALSE,
      n_1 = length(x[!is.na(x) & is.finite(x)]),
      n_2 = length(y[!is.na(y) & is.finite(y)]),
      estimate_1 = if (length(x[!is.na(x) & is.finite(x)]) == 0) NA_real_ else mean(x, na.rm = TRUE),
      estimate_2 = if (length(y[!is.na(y) & is.finite(y)]) == 0) NA_real_ else mean(y, na.rm = TRUE),
      estimate_difference = if (
        length(x[!is.na(x) & is.finite(x)]) == 0 ||
        length(y[!is.na(y) & is.finite(y)]) == 0
      ) NA_real_ else mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
      statistic = tt$statistic,
      p_value = tt$p_value,
      method = tt$method,
      note = tt$note
    )
  }

  # 2) 差异的差异：
  #    delta_TRA = mean(TRA perfect within target_region_id) - mean(TRA background within target_region_id)
  #    delta_REF = mean(REF perfect within target_region_id) - mean(REF background within target_region_id)
  #    检验 delta_TRA 和 delta_REF 是否有显著差异。
  delta_by_region <- dt[, .(
    region_mean_rho = mean(mean_rho, na.rm = TRUE)
  ), by = .(target_region_id, source, interval_type)]

  wide <- dcast(
    delta_by_region,
    target_region_id + source ~ interval_type,
    value.var = "region_mean_rho"
  )

  if (!perfect_type %in% names(wide)) wide[, (perfect_type) := NA_real_]
  if (!bg_type %in% names(wide)) wide[, (bg_type) := NA_real_]

  wide[, delta := get(perfect_type) - get(bg_type)]
  delta_wide <- dcast(
    wide[, .(target_region_id, source, delta)],
    target_region_id ~ source,
    value.var = "delta"
  )

  if (!"TRA" %in% names(delta_wide)) delta_wide[, TRA := NA_real_]
  if (!ref_source %in% names(delta_wide)) delta_wide[, (ref_source) := NA_real_]

  paired_delta <- delta_wide[
    is.finite(TRA) & is.finite(get(ref_source)),
    .(
      target_region_id,
      delta_TRA = TRA,
      delta_REF = get(ref_source),
      diff_of_diff = TRA - get(ref_source)
    )
  ]

  if (nrow(paired_delta) >= min_values_for_test) {
    tt <- run_two_group_test(paired_delta$delta_TRA, paired_delta$delta_REF, paired = TRUE)
  } else {
    tt <- list(p_value = NA_real_, statistic = NA_real_, method = significance_test_method, note = "not_enough_paired_target_regions")
  }

  tests[[length(tests) + 1]] <- data.table(
    batch = batch_name,
    comparison_type = "difference_of_differences",
    comparison = paste0("(TRA perfect-background) vs (", ref_source, " perfect-background)"),
    source_1 = "TRA",
    interval_type_1 = "perfect_minus_background",
    source_2 = ref_source,
    interval_type_2 = "perfect_minus_background",
    paired = TRUE,
    n_1 = nrow(paired_delta),
    n_2 = nrow(paired_delta),
    estimate_1 = if (nrow(paired_delta) == 0) NA_real_ else mean(paired_delta$delta_TRA, na.rm = TRUE),
    estimate_2 = if (nrow(paired_delta) == 0) NA_real_ else mean(paired_delta$delta_REF, na.rm = TRUE),
    estimate_difference = if (nrow(paired_delta) == 0) NA_real_ else mean(paired_delta$diff_of_diff, na.rm = TRUE),
    statistic = tt$statistic,
    p_value = tt$p_value,
    method = tt$method,
    note = tt$note
  )

  res <- rbindlist(tests, use.names = TRUE, fill = TRUE)
  res[, p_adj := p.adjust(p_value, method = significance_p_adjust_method)]
  res[, significance := vapply(p_adj, p_to_stars, character(1))]
  res[]

  attr(res, "paired_delta") <- paired_delta
  res
}

save_plot_all_formats <- function(plot_obj, filename_base, width = 11, height = 7) {
  if (save_png) {
    ggsave(
      filename = file.path(output_dir, paste0(filename_base, ".png")),
      plot = plot_obj,
      width = width,
      height = height,
      dpi = 600,
      bg = "white"
    )
  }
  if (save_pdf) {
    ggsave(
      filename = file.path(output_dir, paste0(filename_base, ".pdf")),
      plot = plot_obj,
      width = width,
      height = height,
      device = "pdf",
      bg = "white"
    )
  }
  if (save_tiff) {
    ggsave(
      filename = file.path(output_dir, paste0(filename_base, ".tiff")),
      plot = plot_obj,
      width = width,
      height = height,
      dpi = 300,
      compression = "lzw",
      bg = "white"
    )
  }
}

add_sig_bracket <- function(p, x1, x2, y, label, tick_height) {
  p +
    annotate("segment", x = x1, xend = x2, y = y, yend = y, color = publication_colors$bracket, linewidth = 0.45) +
    annotate("segment", x = x1, xend = x1, y = y, yend = y - tick_height, color = publication_colors$bracket, linewidth = 0.45) +
    annotate("segment", x = x2, xend = x2, y = y, yend = y - tick_height, color = publication_colors$bracket, linewidth = 0.45) +
    annotate("text", x = (x1 + x2) / 2, y = y + tick_height * 0.35, label = label, size = 3.1, color = publication_colors$bracket)
}

make_bar_distribution_plot <- function(mean_dt, summary_dt, sig_dt, batch_name, ref_source) {
  dt <- as.data.table(mean_dt)[batch == batch_name & !is.na(mean_rho) & is.finite(mean_rho)]
  if (nrow(dt) == 0) {
    cat("警告:", batch_name, "无可用于绘图的区间均值，跳过柱状图。\n")
    return(invisible(NULL))
  }

  perfect_type <- "perfect_invert_inside_target"
  bg_type <- "target_minus_all_invert_overlap"

  dt[, interval_type_short := fifelse(interval_type == perfect_type, "Perfect", "Background")]
  dt[, plot_group := paste(source, interval_type_short, sep = "_")]
  group_order <- c("TRA_Perfect", "TRA_Background", paste0(ref_source, "_Perfect"), paste0(ref_source, "_Background"))
  group_labels <- c(
    "TRA_Perfect" = "TRA perfect",
    "TRA_Background" = "TRA background"
  )
  group_labels[paste0(ref_source, "_Perfect")] <- paste(ref_source, "perfect")
  group_labels[paste0(ref_source, "_Background")] <- paste(ref_source, "background")
  dt[, plot_group := factor(plot_group, levels = group_order)]

  plot_summary <- dt[, .(
    bar_mean = mean(mean_rho, na.rm = TRUE),
    bar_se = if (.N <= 1) 0 else sd(mean_rho, na.rm = TRUE) / sqrt(.N),
    n = .N
  ), by = plot_group]
  plot_summary <- plot_summary[!is.na(plot_group)]
  plot_summary[, x_num := as.numeric(plot_group)]

  dt <- dt[!is.na(plot_group)]
  dt[, x_num := as.numeric(plot_group)]

  fill_values <- c(
    "TRA_Perfect" = publication_colors$tra,
    "TRA_Background" = publication_colors$tra
  )
  fill_values[paste0(ref_source, "_Perfect")] <- ifelse(ref_source == "TBA", publication_colors$tba, publication_colors$tbi)
  fill_values[paste0(ref_source, "_Background")] <- ifelse(ref_source == "TBA", publication_colors$tba, publication_colors$tbi)

  alpha_values <- c(
    "TRA_Perfect" = 0.95,
    "TRA_Background" = 0.45
  )
  alpha_values[paste0(ref_source, "_Perfect")] <- 0.95
  alpha_values[paste0(ref_source, "_Background")] <- 0.45

  y_max_data <- max(c(dt$mean_rho, plot_summary$bar_mean + plot_summary$bar_se), na.rm = TRUE)
  if (!is.finite(y_max_data) || y_max_data <= 0) y_max_data <- 1
  y_step <- y_max_data * 0.12
  tick_height <- y_max_data * 0.035
  y_limit <- y_max_data + y_step * 4.0

  p <- ggplot() +
    create_publication_theme() +
    geom_col(
      data = plot_summary,
      aes(x = plot_group, y = bar_mean, fill = plot_group, alpha = plot_group),
      width = bar_width,
      color = "grey25",
      linewidth = 0.25
    ) +
    geom_errorbar(
      data = plot_summary,
      aes(x = plot_group, ymin = pmax(0, bar_mean - bar_se), ymax = bar_mean + bar_se),
      width = 0.20,
      linewidth = 0.45
    ) +
    geom_jitter(
      data = dt,
      aes(x = plot_group, y = mean_rho),
      width = 0.12,
      height = 0,
      alpha = bar_point_alpha,
      size = bar_point_size
    ) +
    scale_fill_manual(values = fill_values, guide = "none") +
    scale_alpha_manual(values = alpha_values, guide = "none") +
    scale_x_discrete(labels = group_labels, drop = FALSE) +
    scale_y_continuous(labels = scales::scientific_format(digits = 2), limits = c(0, y_limit), expand = expansion(mult = c(0, 0.02))) +
    labs(
      title = paste0(batch_name, ": Mean RHO by interval category"),
      subtitle = "Bars show mean ± SE across analysis intervals; points show per-interval mean RHO",
      x = "Category",
      y = "Mean normalized RHO"
    )

  # 显著性标注：同一Source内部 perfect vs background
  sig <- as.data.table(sig_dt)[batch == batch_name]
  get_p_for <- function(comp_type, src) {
    tmp <- sig[comparison_type == comp_type & source_1 == src]
    if (nrow(tmp) == 0 || !plot_p_value_type %in% names(tmp)) return(NA_real_)
    tmp[[plot_p_value_type]][1]
  }

  p_tra <- get_p_for("within_source", "TRA")
  p_ref <- get_p_for("within_source", ref_source)
  p_did <- sig[comparison_type == "difference_of_differences", ][[plot_p_value_type]][1]
  if (length(p_did) == 0) p_did <- NA_real_

  p <- add_sig_bracket(
    p,
    x1 = 1,
    x2 = 2,
    y = y_max_data + y_step * 0.7,
    label = paste0("TRA: ", safe_p_label(p_tra, plot_p_label_prefix()), " ", p_to_stars(p_tra)),
    tick_height = tick_height
  )

  p <- add_sig_bracket(
    p,
    x1 = 3,
    x2 = 4,
    y = y_max_data + y_step * 1.5,
    label = paste0(ref_source, ": ", safe_p_label(p_ref, plot_p_label_prefix()), " ", p_to_stars(p_ref)),
    tick_height = tick_height
  )

  p <- add_sig_bracket(
    p,
    x1 = 1.5,
    x2 = 3.5,
    y = y_max_data + y_step * 2.6,
    label = paste0("Diff-in-diff: ", safe_p_label(p_did, plot_p_label_prefix()), " ", p_to_stars(p_did)),
    tick_height = tick_height
  )

  filename_base <- paste0(output_prefix, "_", batch_name, "_bar_distribution")
  save_plot_all_formats(p, filename_base, width = 11, height = 7)
  invisible(p)
}

make_combined_bar_distribution_plot <- function(mean_dt, sig_dt) {
  dt <- as.data.table(mean_dt)[!is.na(mean_rho) & is.finite(mean_rho)]
  if (nrow(dt) == 0) {
    cat("警告: 无可用于绘图的区间均值，跳过合并柱状图。\n")
    return(invisible(NULL))
  }

  batch_ref_map <- rbindlist(lapply(batch_configs, function(cfg) {
    data.table(batch = cfg$batch, ref_source = cfg$ref_source)
  }), use.names = TRUE, fill = TRUE)

  dt <- merge(dt, batch_ref_map, by = "batch", all.x = TRUE)
  dt <- dt[!is.na(ref_source)]
  if (nrow(dt) == 0) {
    cat("警告: 合并柱状图没有匹配到batch配置，跳过。\n")
    return(invisible(NULL))
  }

  perfect_type <- "perfect_invert_inside_target"
  bg_type <- "target_minus_all_invert_overlap"

  dt[, interval_type_short := fifelse(interval_type == perfect_type, "Perfect", "Background")]
  dt[, x_num := fifelse(
    source == "TRA" & interval_type_short == "Perfect", 1,
    fifelse(
      source == "TRA" & interval_type_short == "Background", 2,
      fifelse(
        source == ref_source & interval_type_short == "Perfect", 3,
        fifelse(source == ref_source & interval_type_short == "Background", 4, NA_real_)
      )
    )
  )]
  dt <- dt[!is.na(x_num)]

  dt[, batch_label := paste0(batch, "\nReference: ", ref_source)]
  dt[, fill_color := fifelse(
    source == "TRA",
    publication_colors$tra,
    fifelse(ref_source == "TBA", publication_colors$tba, publication_colors$tbi)
  )]
  dt[, alpha_value := fifelse(interval_type_short == "Perfect", 0.95, 0.45)]

  plot_summary <- dt[, .(
    bar_mean = mean(mean_rho, na.rm = TRUE),
    bar_se = if (.N <= 1) 0 else sd(mean_rho, na.rm = TRUE) / sqrt(.N),
    n = .N,
    fill_color = fill_color[1],
    alpha_value = alpha_value[1]
  ), by = .(batch, batch_label, ref_source, x_num)]

  sig_all <- as.data.table(sig_dt)

  bracket_rows <- list()
  for (bn in unique(dt$batch)) {
    dt_bn <- dt[batch == bn]
    sum_bn <- plot_summary[batch == bn]
    ref_src <- unique(dt_bn$ref_source)[1]
    batch_lab <- unique(dt_bn$batch_label)[1]

    y_max_data <- max(c(dt_bn$mean_rho, sum_bn$bar_mean + sum_bn$bar_se), na.rm = TRUE)
    if (!is.finite(y_max_data) || y_max_data <= 0) y_max_data <- 1
    y_step <- y_max_data * 0.12
    tick_height <- y_max_data * 0.035

    sig_bn <- sig_all[batch == bn]

    get_p_for <- function(comp_type, src) {
      tmp <- sig_bn[comparison_type == comp_type & source_1 == src]
      if (nrow(tmp) == 0 || !plot_p_value_type %in% names(tmp)) return(NA_real_)
      tmp[[plot_p_value_type]][1]
    }

    p_tra <- get_p_for("within_source", "TRA")
    p_ref <- get_p_for("within_source", ref_src)
    p_did <- sig_bn[comparison_type == "difference_of_differences", ][[plot_p_value_type]][1]
    if (length(p_did) == 0) p_did <- NA_real_

    bracket_rows[[length(bracket_rows) + 1]] <- data.table(
      batch_label = batch_lab,
      x1 = c(1, 3, 1.5),
      x2 = c(2, 4, 3.5),
      y = c(
        y_max_data + y_step * 0.7,
        y_max_data + y_step * 1.5,
        y_max_data + y_step * 2.6
      ),
      tick_height = tick_height,
      label = c(
        paste0("TRA: ", safe_p_label(p_tra, plot_p_label_prefix()), " ", p_to_stars(p_tra)),
        paste0(ref_src, ": ", safe_p_label(p_ref, plot_p_label_prefix()), " ", p_to_stars(p_ref)),
        paste0("Diff-in-diff: ", safe_p_label(p_did, plot_p_label_prefix()), " ", p_to_stars(p_did))
      )
    )
  }

  bracket_dt <- rbindlist(bracket_rows, use.names = TRUE, fill = TRUE)

  p <- ggplot() +
    create_publication_theme() +
    geom_col(
      data = plot_summary,
      aes(x = x_num, y = bar_mean, fill = fill_color, alpha = alpha_value),
      width = bar_width,
      color = "grey25",
      linewidth = 0.25
    ) +
    geom_errorbar(
      data = plot_summary,
      aes(x = x_num, ymin = pmax(0, bar_mean - bar_se), ymax = bar_mean + bar_se),
      width = 0.20,
      linewidth = 0.45
    ) +
    geom_jitter(
      data = dt,
      aes(x = x_num, y = mean_rho),
      width = 0.12,
      height = 0,
      alpha = bar_point_alpha,
      size = bar_point_size
    ) +
    geom_segment(
      data = bracket_dt,
      aes(x = x1, xend = x2, y = y, yend = y),
      inherit.aes = FALSE,
      color = publication_colors$bracket,
      linewidth = 0.45
    ) +
    geom_segment(
      data = bracket_dt,
      aes(x = x1, xend = x1, y = y, yend = y - tick_height),
      inherit.aes = FALSE,
      color = publication_colors$bracket,
      linewidth = 0.45
    ) +
    geom_segment(
      data = bracket_dt,
      aes(x = x2, xend = x2, y = y, yend = y - tick_height),
      inherit.aes = FALSE,
      color = publication_colors$bracket,
      linewidth = 0.45
    ) +
    geom_text(
      data = bracket_dt,
      aes(x = (x1 + x2) / 2, y = y + tick_height * 0.35, label = label),
      inherit.aes = FALSE,
      size = 3.0,
      color = publication_colors$bracket
    ) +
    facet_wrap(~ batch_label, nrow = 1, scales = "free_y") +
    scale_fill_identity() +
    scale_alpha_identity() +
    scale_x_continuous(
      breaks = c(1, 2, 3, 4),
      labels = c("TRA\nperfect", "TRA\nbackground", "REF\nperfect", "REF\nbackground")
    ) +
    scale_y_continuous(
      labels = scales::scientific_format(digits = 2),
      expand = expansion(mult = c(0, 0.18))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = "Mean RHO by interval category across batches",
      subtitle = paste0(
        "Bars show mean ± SE across analysis intervals; points show per-interval mean RHO; labels show ",
        plot_p_label_prefix()
      ),
      x = "Category",
      y = "Mean normalized RHO"
    )

  filename_base <- paste0(output_prefix, "_ALL_batches_bar_distribution_combined")
  save_plot_all_formats(p, filename_base, width = 14, height = 7)
  invisible(p)
}


analyze_batch <- function(batch_cfg, file_infos) {
  batch_name <- batch_cfg$batch
  ref_source <- batch_cfg$ref_source

  cat("\n============================================================\n")
  cat("开始分析 batch:", batch_name, "\n")
  cat("参考Source:", ref_source, "\n")

  invert_dt <- read_invert_bed(invert_files[[batch_cfg$invert_key]], batch_cfg$invert_key)
  cat("invert区间数:", nrow(invert_dt), "\n")

  interval_info <- make_analysis_intervals(batch_name, invert_dt, target_tra_regions)
  intervals <- interval_info$intervals
  excluded_partial <- interval_info$excluded_partial

  cat("最终分析区间数:", nrow(intervals), "\n")
  if (nrow(intervals) > 0) {
    print(intervals[, .(n_intervals = .N, total_length = sum(length)), by = interval_type])
  }
  if (nrow(excluded_partial) > 0) {
    cat("partial-overlap invert区间被排除数:", nrow(excluded_partial), "\n")
  }

  interval_means <- compute_all_interval_means_for_batch(batch_cfg, file_infos, intervals)
  category_summary <- summarize_categories(interval_means)

  significance_tests <- compute_significance_tests_for_batch(interval_means, batch_name, ref_source)
  paired_delta <- attr(significance_tests, "paired_delta")
  if (is.null(paired_delta)) paired_delta <- data.table()

  make_bar_distribution_plot(interval_means, category_summary, significance_tests, batch_name, ref_source)

  list(
    intervals = intervals,
    excluded_partial = excluded_partial,
    interval_means = interval_means,
    category_summary = category_summary,
    significance_tests = significance_tests,
    paired_delta = paired_delta
  )
}

# ============================================================
# 2. 主流程
# ============================================================

cat("\n当前关键参数设置:\n")
cat("output_dir:", output_dir, "\n")
cat("use_awk_chr_filter:", use_awk_chr_filter, "\n")
cat("data.table threads:", getDTthreads(), "\n")
cat("significance_test_method:", significance_test_method, "\n")
cat("significance_p_adjust_method:", significance_p_adjust_method, "\n")
cat("plot_p_value_type:", plot_p_value_type, "\n")

check_config()

file_infos <- list(
  tba_bed = detect_bed_format(input_files$tba_bed, "tba_bed"),
  tra_for_tba_bed = detect_bed_format(input_files$tra_for_tba_bed, "tra_for_tba_bed"),
  tbi_bed = detect_bed_format(input_files$tbi_bed, "tbi_bed"),
  tra_for_tbi_bed = detect_bed_format(input_files$tra_for_tbi_bed, "tra_for_tbi_bed")
)

cat("\n输入文件格式识别结果:\n")
for (nm in names(file_infos)) {
  cat("-", nm, ":", file_infos[[nm]]$format, ", n_cols =", file_infos[[nm]]$n_cols, "\n")
}

all_results <- list()
for (bn in names(batch_configs)) {
  all_results[[bn]] <- analyze_batch(batch_configs[[bn]], file_infos)
}

all_intervals <- rbindlist(lapply(all_results, `[[`, "intervals"), use.names = TRUE, fill = TRUE)
all_excluded_partial <- rbindlist(lapply(all_results, `[[`, "excluded_partial"), use.names = TRUE, fill = TRUE)
all_interval_means <- rbindlist(lapply(all_results, `[[`, "interval_means"), use.names = TRUE, fill = TRUE)
all_category_summary <- rbindlist(lapply(all_results, `[[`, "category_summary"), use.names = TRUE, fill = TRUE)
all_significance_tests <- rbindlist(lapply(all_results, `[[`, "significance_tests"), use.names = TRUE, fill = TRUE)

# ---- 新增：将两个batch原本分别输出的柱状分布图合并到同一张图中 ----
make_combined_bar_distribution_plot(all_interval_means, all_significance_tests)

all_paired_delta <- rbindlist(lapply(names(all_results), function(nm) {
  x <- all_results[[nm]]$paired_delta
  if (is.null(x) || nrow(x) == 0) return(data.table())
  x[, batch := nm]
  x
}), use.names = TRUE, fill = TRUE)

# ---- 输出明细表 ----
analysis_intervals_path <- file.path(output_dir, paste0(output_prefix, "_ALL_analysis_intervals.tsv"))
excluded_partial_path <- file.path(output_dir, paste0(output_prefix, "_ALL_excluded_partial_invert_overlaps.tsv"))
interval_mean_path <- file.path(output_dir, paste0(output_prefix, "_ALL_interval_mean_rho.tsv"))
category_summary_path <- file.path(output_dir, paste0(output_prefix, "_ALL_category_summary.tsv"))
significance_all_path <- file.path(output_dir, paste0(output_prefix, "_ALL_significance_tests.tsv"))
paired_delta_path <- file.path(output_dir, paste0(output_prefix, "_ALL_target_region_deltas.tsv"))

fwrite(all_intervals, analysis_intervals_path, sep = "\t", quote = FALSE, na = "NA")
fwrite(all_excluded_partial, excluded_partial_path, sep = "\t", quote = FALSE, na = "NA")
fwrite(all_interval_means, interval_mean_path, sep = "\t", quote = FALSE, na = "NA")
fwrite(all_category_summary, category_summary_path, sep = "\t", quote = FALSE, na = "NA")
fwrite(all_significance_tests, significance_all_path, sep = "\t", quote = FALSE, na = "NA")
fwrite(all_paired_delta, paired_delta_path, sep = "\t", quote = FALSE, na = "NA")

for (bn in names(all_results)) {
  sig_path <- file.path(output_dir, paste0(output_prefix, "_", bn, "_significance_tests.tsv"))
  fwrite(all_results[[bn]]$significance_tests, sig_path, sep = "\t", quote = FALSE, na = "NA")
}

cat("\n=== 分析完成 ===\n")
cat("输出目录:", output_dir, "\n")
cat("分析区间表:", analysis_intervals_path, "\n")
cat("partial-overlap排除表:", excluded_partial_path, "\n")
cat("区间平均RHO表:", interval_mean_path, "\n")
cat("类别汇总表:", category_summary_path, "\n")
cat("显著性检验总表:", significance_all_path, "\n")
cat("target-region delta表:", paired_delta_path, "\n")
cat("合并柱状图:", file.path(output_dir, paste0(output_prefix, "_ALL_batches_bar_distribution_combined.[png/pdf/tiff]")), "\n")
cat("\n说明:\n")
cat("1. perfect_invert_inside_target 表示 invert 区间完整落在目标TRA区间内。\n")
cat("2. target_minus_all_invert_overlap 表示目标TRA区间去除所有invert重叠后的剩余部分。\n")
cat("3. partial-overlap invert 不进入 perfect 类别，其与目标区间的重叠部分也会从 background 中剔除。\n")
cat("4. within_source 检验同一Source内部 perfect vs background。\n")
cat("5. difference_of_differences 检验 TRA(perfect-background) 与 REF(perfect-background) 的差异是否显著不同。\n")
