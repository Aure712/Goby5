library(dplyr)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(data.table)

# ============================================================
# 0. 统一参数设置区：运行前主要修改这里
# ============================================================

# ---- 输入文件：4个BED文件 ----
# 说明：
# 1) tba_bed 与 tra_for_tba_bed 的坐标一致
# 2) tbi_bed 与 tra_for_tbi_bed 的坐标一致
# 3) 上述两组坐标之间可以不一致
input_files <- list(
  tba_bed         = "/fast3/group_crf/home/g21shaoy23/goby5/ismc/tba/20inds/trans/output2.bed",
  tra_for_tba_bed = "/fast3/group_crf/home/g21shaoy23/goby5/ismc/tba/20inds/trans/tra/new_output2.bed",
  tbi_bed         = "/public4/group_crf/home/g21shaoy23/ismc/tbi2tra/output2.bed",
  tra_for_tbi_bed = "/public4/group_crf/home/g21shaoy23/ismc/tbi2tra/tra/new_output2.bed"
)

# ---- 每个样本/群体对应的x值：RHO = BED数值 / x ----
x_values <- list(
  tba = 0.007300749,
  tra = 0.005889914,
  tbi = 0.011549063   
)

# ---- 坐标处理模式 ----
# "union_tra"：TBA和TBI各自BED不变；两个TRA BED按坐标取并集；每条线使用自己的坐标窗口绘图。
#              若两个TRA文件有相同坐标，默认对TRA数值取平均，避免重复坐标被加权两次。
# "intersection_4"：4个BED文件取坐标交集；只在4个文件都存在的坐标上绘图。
coordinate_mode <- "union_tra"   # 可改为 "intersection_4"

# ---- 指定要统计和绘图的TraScf染色体 ----
# 只处理这里列出的TraScf；其余TraScf染色体会在主循环前跳过。
# 示例：
#   selected_trascf <- c("TraScf_7", "TraScf_9", "TraScf_10")
#   selected_trascf <- c("7", "9", "10")
#
# 如果希望恢复处理所有染色体，可设置为：
#   selected_trascf <- NULL
# 或：
#   selected_trascf <- character(0)
selected_trascf <- c("TraScf_1", "TraScf_5", "TraScf_8", "TraScf_10")

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

# ---- 指定要分析的断点编号和方向 ----
# 说明：
# 1) chr：TraScf染色体名称，可写成 "TraScf_7" 或 "7"。
# 2) bp_index：该TraScf染色体上按断点坐标从小到大排序后的断点编号。
#    例如 TraScf_1 有 40492542 和 80104410 两个断点，
#    则 40492542 是 bp_index = 1，80104410 是 bp_index = 2。
# 3) direction：只分析该断点的哪个方向，可选 "left"、"right" 或 "both"。
#    如果同一个断点需要同时分析左右两侧，也可以写一行 direction = "both"，
#    或写两行分别为 "left" 和 "right"。
#
# 如果希望恢复分析 selected_trascf 中所有断点和所有方向，可设置为：
#   selected_breakpoint_targets <- NULL
# 或：
#   selected_breakpoint_targets <- data.frame()
#
# 当前示例：只分析 TraScf_7 的1号断点左侧、TraScf_9 的1号断点右侧、
#          TraScf_10 的1号断点左侧。运行前请按实际需求修改这里。
selected_breakpoint_targets <- data.frame(
  chr = c("TraScf_1", "TraScf_5", "TraScf_8", "TraScf_10"),
  bp_index = c(1L, 1L, 1L, 1L),
  direction = c("right", "right", "right", "right"),
  stringsAsFactors = FALSE
)

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
    return(region_data[0])
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


# ---- 输出设置 ----
output_dir <- paste0("plots_", coordinate_mode)
output_prefix <- coordinate_mode
save_png <- TRUE
save_pdf <- TRUE
save_tiff <- TRUE

# 如果希望尽量保持原始文件名，例如 TraScf_1_publication.png，而不是 union_tra_TraScf_1_publication.png，
# 可以改为 FALSE。默认 TRUE 是为了避免两种坐标模式的结果互相覆盖。
add_mode_prefix_to_filenames <- TRUE

# ---- 低内存设置 ----
# TRUE：使用 awk 按染色体从BED中筛选行，只把当前染色体读入R。
# 这是降低内存占用的关键设置。HPC/Linux/Slurm 环境通常都支持 awk。
use_awk_chr_filter <- TRUE

# data.table 多线程有时会增加瞬时内存峰值。内存紧张时建议保持1。
data_table_threads <- 1
setDTthreads(data_table_threads)

# 每处理完一个染色体是否主动gc。
run_gc_each_chr <- TRUE

# ---- 绘图窗口和扫描窗口 ----
plot_window <- 5e5       # 绘图窗口大小，例如 5e5 = 0.5 Mb
scan_window <- 2e6       # 扫描窗口大小，例如 2e6 = 2 Mb
scan_step <- 2e6         # 扫描步长，例如 2e6 = 2 Mb
summary_extension <- 5e7 # 断点汇总图中，断点两侧各延伸范围，例如 5e7 = 50 Mb

# ---- 扫描线逻辑 ----
# 原代码逻辑是 TRA > TBA 时记录扫描线。
# 如果希望TRA同时大于TBA和TBI才记录扫描线，可改成 c("TBA", "TBI")。
scan_reference_sources <- c("TBA")

# 是否在图上绘制绿色扫描虚线。
# TRUE：绘制；FALSE：不绘制，但仍会计算并可输出绿色虚线对应的详细数值。
draw_scan_lines <- TRUE

# 是否输出绿色扫描虚线的详细数值表。
# 输出字段包括：染色体、断点、左右方向、扫描线绝对位置、相对断点距离、
# 扫描窗口起止坐标、TRA/TBA/TBI窗口均值、窗口内点数等。
write_scan_line_details <- TRUE
scan_line_details_file <- paste0(output_prefix, "_Green_Scan_Lines_Details.tsv")

# ---- 平滑度控制参数 ----
smoothing_params <- list(
  enabled = FALSE,        # 是否启用平滑
  method = "none",        # "loess", "spline", "moving_average", "none"
  span = 0.1,             # loess平滑度，越小越贴近原始数据
  degree = 1,             # loess多项式次数，1或2
  spar = 0.5,             # spline平滑参数
  ma_window = 3           # 移动平均窗口大小
)

# ---- 竖线坐标定义 ----
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

# ---- 颜色主题 ----
publication_colors <- list(
  tba = "#2E86AB",        # TBA：蓝色
  tra = "#A23B72",        # TRA：玫瑰色
  tbi = "#F18F01",        # TBI：橙色
  marker = "#C73E1D",     # 断点标记线
  scan = "#6A994E",       # 扫描线
  background = "#FAFAFA",
  text = "#2D3748"
)

source_levels <- c("TBA", "TRA", "TBI")
source_color_values <- c(
  "TBA" = publication_colors$tba,
  "TRA" = publication_colors$tra,
  "TBI" = publication_colors$tbi
)
source_labels <- c(
  "TBA" = "ISMC Analysis of TBA",
  "TRA" = "ISMC Analysis of TRA",
  "TBI" = "ISMC Analysis of TBI"
)

# ============================================================
# 1. 通用辅助函数
# ============================================================

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

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
      legend.position = "bottom",
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18, face = "bold"),
      legend.margin = margin(t = 15),
      legend.key.size = unit(1.5, "lines"),
      legend.spacing.x = unit(1.0, "cm"),
      plot.margin = margin(20, 20, 20, 20),
      strip.background = element_rect(fill = publication_colors$background, color = publication_colors$text),
      strip.text = element_text(color = publication_colors$text, size = 11, face = "bold")
    )
}

clean_chr_name <- function(chr_name) {
  gsub("^TraScf_", "", chr_name)
}

sort_chromosomes <- function(chr_vec) {
  chr_vec <- unique(chr_vec)
  chr_vec <- chr_vec[!is.na(chr_vec) & chr_vec != ""]
  chr_num <- suppressWarnings(as.numeric(clean_chr_name(chr_vec)))
  chr_vec[order(is.na(chr_num), chr_num, chr_vec)]
}

check_file_exists <- function(filepath, label) {
  if (!file.exists(filepath)) {
    stop(paste0("找不到输入文件: ", label, " = ", filepath))
  }
}

check_config <- function() {
  if (!coordinate_mode %in% c("union_tra", "intersection_4")) {
    stop("coordinate_mode 只能是 'union_tra' 或 'intersection_4'")
  }
  for (nm in names(input_files)) check_file_exists(input_files[[nm]], nm)
  for (nm in names(x_values)) {
    if (!is.numeric(x_values[[nm]]) || length(x_values[[nm]]) != 1 || is.na(x_values[[nm]]) || x_values[[nm]] == 0) {
      stop(paste0("x_values$", nm, " 必须是非零数值"))
    }
  }

  if (is_breakpoint_direction_filter_enabled(selected_breakpoint_targets)) {
    targets <- tryCatch(
      as.data.frame(selected_breakpoint_targets, stringsAsFactors = FALSE),
      error = function(e) stop("selected_breakpoint_targets 必须是 data.frame 或可转换为 data.frame 的对象")
    )

    required_cols <- c("chr", "bp_index", "direction")
    missing_cols <- setdiff(required_cols, names(targets))
    if (length(missing_cols) > 0) {
      stop("selected_breakpoint_targets 缺少列: ", paste(missing_cols, collapse = ", "))
    }

    target_bp_index <- suppressWarnings(as.integer(targets$bp_index))
    if (any(is.na(target_bp_index) | target_bp_index < 1)) {
      stop("selected_breakpoint_targets$bp_index 必须是大于等于1的整数")
    }

    target_direction <- normalize_breakpoint_direction(targets$direction)
    invalid_direction <- setdiff(target_direction, c("left", "right", "both"))
    if (length(invalid_direction) > 0) {
      stop("selected_breakpoint_targets$direction 只能包含 'left'、'right' 或 'both'")
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
  }

  invalid_scan_ref <- setdiff(scan_reference_sources, c("TBA", "TBI"))
  if (length(invalid_scan_ref) > 0) {
    stop("scan_reference_sources 只能包含 'TBA' 和/或 'TBI'")
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

# 兼容两类输入：
# 1) 4列BED: chr, start, end, value；位置使用(start+end)/2；坐标key使用start:end
# 2) 原脚本中的多列文件: 使用第5列为chr、第7列为pos、第8列为value；坐标key使用pos
# 注意：如果某些文件是4列、某些文件是8列，intersection_4 可能因为坐标key定义不同而无交集。
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

get_chromosomes_from_file <- function(file_info) {
  if (use_awk_chr_filter) {
    cmd <- paste0(
      make_reader_cmd(file_info$path),
      " | awk -F '\\t' '{print $", file_info$chr_col, "}' | sort -u"
    )
    chrs <- tryCatch(system(cmd, intern = TRUE), error = function(e) character())
    chrs <- chrs[!is.na(chrs) & chrs != ""]
    return(chrs)
  }

  # 备用方案：会读取该文件的染色体列，内存占用高于awk方案。
  if (file_info$format == "bed4") {
    tmp <- fread(file_info$path, header = FALSE, select = 1, showProgress = FALSE)
  } else {
    tmp <- fread(file_info$path, header = FALSE, select = 5, showProgress = FALSE)
  }
  unique(as.character(tmp[[1]]))
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
    # 备用方案：读取必要列后在R中筛选，内存占用较高。
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
    return(data.table(chr = character(), pos = numeric(), pos_key = character(), adjusted_y = numeric(), source = character()))
  }

  if (file_info$format == "bed4") {
    if (ncol(raw_dt) < 4) stop(paste0("读取 ", file_info$label, " 时列数不足。"))
    setnames(raw_dt, c("V1", "V2", "V3", "V4"))
    dt <- raw_dt[, .(
      chr = as.character(V1),
      pos = (as.numeric(V2) + as.numeric(V3)) / 2,
      pos_key = paste0(as.character(V2), ":", as.character(V3)),
      raw_y = as.numeric(V4)
    )]
  } else {
    if (ncol(raw_dt) < 3) stop(paste0("读取 ", file_info$label, " 时列数不足。"))
    setnames(raw_dt, c("V1", "V2", "V3"))
    dt <- raw_dt[, .(
      chr = as.character(V1),
      pos = as.numeric(V2),
      pos_key = as.character(V2),
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

empty_plot_dt <- function() {
  data.table(chr = character(), pos = numeric(), pos_key = character(), adjusted_y = numeric(), source = character())
}

empty_scan_detail_dt <- function() {
  data.table(
    chr = character(),
    bp_label = character(),
    breakpoint = numeric(),
    direction = character(),
    scan_position = numeric(),
    relative_pos = numeric(),
    scan_position_Mb = numeric(),
    relative_pos_Mb = numeric(),
    window_start = numeric(),
    window_end = numeric(),
    window_start_Mb = numeric(),
    window_end_Mb = numeric(),
    tra_mean = numeric(),
    tba_mean = numeric(),
    tbi_mean = numeric(),
    n_total = integer(),
    n_TRA = integer(),
    n_TBA = integer(),
    n_TBI = integer(),
    scan_window = numeric(),
    scan_step = numeric(),
    scan_reference_sources = character(),
    coordinate_mode = character()
  )
}

# 同一source、同一坐标重复时取平均，避免重复坐标被重复加权。
deduplicate_by_pos <- function(dt, source_name, chr_name) {
  if (is.null(dt) || nrow(dt) == 0) return(empty_plot_dt())
  dt <- as.data.table(dt)
  dt[, source := source_name]
  res <- dt[, .(
    chr = chr_name,
    pos = mean(pos, na.rm = TRUE),
    adjusted_y = mean(adjusted_y, na.rm = TRUE)
  ), by = .(source, pos_key)]
  res <- res[!is.na(pos) & is.finite(pos) & !is.na(adjusted_y) & is.finite(adjusted_y)]
  setcolorder(res, c("chr", "pos", "pos_key", "adjusted_y", "source"))
  res[]
}

get_common_pos_keys <- function(...) {
  dts <- list(...)
  key_tables <- lapply(dts, function(x) {
    if (is.null(x) || nrow(x) == 0) return(data.table(pos_key = character()))
    unique(as.data.table(x)[, .(pos_key)])
  })
  common_dt <- Reduce(function(a, b) fintersect(a, b), key_tables)
  common_dt$pos_key
}

prepare_chr_plot_table <- function(chr_name, file_infos) {
  cat("\n读取当前染色体数据:", chr_name, "\n")

  tba_dt <- read_chr_bed(file_infos$tba_bed, chr_name, "TBA", x_values$tba, "tba_bed")
  tra_tba_dt <- read_chr_bed(file_infos$tra_for_tba_bed, chr_name, "TRA", x_values$tra, "tra_for_tba_bed")
  tbi_dt <- read_chr_bed(file_infos$tbi_bed, chr_name, "TBI", x_values$tbi, "tbi_bed")
  tra_tbi_dt <- read_chr_bed(file_infos$tra_for_tbi_bed, chr_name, "TRA", x_values$tra, "tra_for_tbi_bed")

  cat("当前染色体原始有效记录数 - TBA:", nrow(tba_dt),
      ", TRA_for_TBA:", nrow(tra_tba_dt),
      ", TBI:", nrow(tbi_dt),
      ", TRA_for_TBI:", nrow(tra_tbi_dt), "\n")

  if (coordinate_mode == "union_tra") {
    tba_keep <- deduplicate_by_pos(tba_dt, "TBA", chr_name)
    tbi_keep <- deduplicate_by_pos(tbi_dt, "TBI", chr_name)
    tra_keep <- deduplicate_by_pos(
      rbindlist(list(tra_tba_dt, tra_tbi_dt), use.names = TRUE, fill = TRUE),
      "TRA",
      chr_name
    )
    plot_dt <- rbindlist(list(tba_keep, tra_keep, tbi_keep), use.names = TRUE, fill = TRUE)

  } else if (coordinate_mode == "intersection_4") {
    common_keys <- get_common_pos_keys(tba_dt, tra_tba_dt, tbi_dt, tra_tbi_dt)
    cat("当前染色体4个BED共有坐标数:", length(common_keys), "\n")

    if (length(common_keys) == 0) {
      plot_dt <- empty_plot_dt()
    } else {
      tba_keep <- deduplicate_by_pos(tba_dt[pos_key %chin% common_keys], "TBA", chr_name)
      tbi_keep <- deduplicate_by_pos(tbi_dt[pos_key %chin% common_keys], "TBI", chr_name)
      tra_keep <- deduplicate_by_pos(
        rbindlist(list(
          tra_tba_dt[pos_key %chin% common_keys],
          tra_tbi_dt[pos_key %chin% common_keys]
        ), use.names = TRUE, fill = TRUE),
        "TRA",
        chr_name
      )
      plot_dt <- rbindlist(list(tba_keep, tra_keep, tbi_keep), use.names = TRUE, fill = TRUE)
    }
  }

  rm(tba_dt, tra_tba_dt, tbi_dt, tra_tbi_dt)

  if (nrow(plot_dt) > 0) {
    plot_dt <- plot_dt[source %in% source_levels]
    plot_dt[, source := factor(as.character(source), levels = source_levels)]
    setorder(plot_dt, pos, source)
  }
  plot_dt[]
}

create_plot_data <- function(data, window_size, bin_origin = NULL) {
  if (is.null(data) || nrow(data) == 0) return(data.frame())

  dt <- as.data.table(data)
  if (is.null(bin_origin)) bin_origin <- min(dt$pos, na.rm = TRUE)

  dt <- dt[!is.na(pos) & !is.na(adjusted_y) & is.finite(pos) & is.finite(adjusted_y)]
  if (nrow(dt) == 0) return(data.frame())

  dt[, window_id := floor((pos - bin_origin) / window_size)]
  result <- dt[, .(
    pos = bin_origin + window_id[1] * window_size + window_size / 2,
    mean_y = mean(adjusted_y, na.rm = TRUE),
    se_y = {
      n_non_na <- sum(!is.na(adjusted_y))
      if (n_non_na <= 1) 0 else sd(adjusted_y, na.rm = TRUE) / sqrt(n_non_na)
    },
    n_points = .N
  ), by = .(chr, source, window_id)]

  result <- result[!is.na(mean_y) & !is.na(pos)]
  setorder(result, source, pos)
  as.data.frame(result)
}

apply_smoothing <- function(data, params = smoothing_params) {
  if (is.null(data) || nrow(data) == 0) return(data)
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
      smooth_fit <- loess(mean_y ~ pos, data = data, span = params$span, degree = params$degree)
      data$smooth_y <- predict(smooth_fit)
    } else if (params$method == "spline") {
      smooth_fit <- smooth.spline(data$pos, data$mean_y, spar = params$spar)
      data$smooth_y <- predict(smooth_fit, data$pos)$y
    } else if (params$method == "moving_average") {
      n <- nrow(data)
      window <- params$ma_window
      data$smooth_y <- data$mean_y
      for (i in seq_len(n)) {
        start_idx <- max(1, i - floor(window / 2))
        end_idx <- min(n, i + floor(window / 2))
        data$smooth_y[i] <- mean(data$mean_y[start_idx:end_idx], na.rm = TRUE)
      }
    } else {
      data$smooth_y <- data$mean_y
    }

    na_indices <- is.na(data$smooth_y)
    if (any(na_indices)) data$smooth_y[na_indices] <- data$mean_y[na_indices]
    data
  }, error = function(e) {
    cat("平滑处理失败，使用原始数据:", e$message, "\n")
    data$smooth_y <- data$mean_y
    data
  })
}

create_smoothed_plot_data <- function(chr_dt, window_size) {
  if (is.null(chr_dt) || nrow(chr_dt) == 0) return(data.table())
  chr_dt <- as.data.table(chr_dt)
  bin_origin <- min(chr_dt$pos, na.rm = TRUE)

  res_list <- lapply(source_levels, function(src) {
    tmp <- chr_dt[as.character(source) == src]
    if (nrow(tmp) == 0) return(data.frame())
    pd <- create_plot_data(tmp, window_size, bin_origin = bin_origin)
    if (nrow(pd) == 0) return(data.frame())
    pd <- apply_smoothing(pd)
    pd$source <- factor(src, levels = source_levels)
    pd
  })

  res <- rbindlist(res_list, use.names = TRUE, fill = TRUE)
  if (nrow(res) > 0) setorder(res, source, pos)
  res[]
}

get_y_limits <- function(plot_dt) {
  if (is.null(plot_dt) || nrow(plot_dt) == 0) return(c(0, 1))
  y_low <- plot_dt$smooth_y - plot_dt$se_y
  y_high <- plot_dt$smooth_y + plot_dt$se_y
  y_min <- suppressWarnings(min(y_low, na.rm = TRUE))
  y_max <- suppressWarnings(max(y_high, na.rm = TRUE))
  if (!is.finite(y_min) || !is.finite(y_max)) return(c(0, 1))
  y_margin <- (y_max - y_min) * 0.1
  if (!is.finite(y_margin) || y_margin == 0) y_margin <- max(abs(y_max), 1) * 0.05
  c(max(0, y_min - y_margin), y_max + y_margin)
}

add_main_series <- function(p, plot_dt, ribbon_alpha = 0.18, line_alpha = 0.9, line_size = 1.2) {
  if (is.null(plot_dt) || nrow(plot_dt) == 0) return(p)
  plot_dt$source <- factor(as.character(plot_dt$source), levels = source_levels)

  p +
    geom_ribbon(
      data = plot_dt,
      aes(x = pos / 1e6, ymin = smooth_y - se_y, ymax = smooth_y + se_y, fill = source, group = source),
      alpha = ribbon_alpha,
      color = NA
    ) +
    geom_line(
      data = plot_dt,
      aes(x = pos / 1e6, y = smooth_y, color = source, group = source),
      size = line_size,
      alpha = line_alpha
    ) +
    scale_color_manual(name = "Sample", values = source_color_values, labels = source_labels, drop = FALSE) +
    scale_fill_manual(name = "Sample", values = source_color_values, labels = source_labels, drop = FALSE) +
    guides(
      fill = "none",
      color = guide_legend(
        override.aes = list(alpha = 1, size = 2.5),
        title.position = "top",
        byrow = TRUE
      )
    )
}

scan_breakpoint <- function(chr_dt, chr_name, bp_label, bp_pos, left_bound, right_bound, scan_directions = c("left", "right")) {
  if (is.null(chr_dt) || nrow(chr_dt) == 0) return(empty_scan_detail_dt())

  scan_directions <- unique(normalize_breakpoint_direction(scan_directions))
  scan_directions <- scan_directions[scan_directions %in% c("left", "right")]
  if (length(scan_directions) == 0) return(empty_scan_detail_dt())

  dt_scan <- as.data.table(chr_dt)
  setkey(dt_scan, pos)

  scan_details <- empty_scan_detail_dt()

  safe_mean <- function(x) {
    x <- x[!is.na(x) & is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    mean(x)
  }

  evaluate_window <- function(window_start, window_end, direction) {
    window_data <- dt_scan[pos >= window_start & pos < window_end]
    if (nrow(window_data) == 0) return(NULL)

    src <- as.character(window_data$source)

    tra_mean <- safe_mean(window_data$adjusted_y[src == "TRA"])
    tba_mean <- safe_mean(window_data$adjusted_y[src == "TBA"])
    tbi_mean <- safe_mean(window_data$adjusted_y[src == "TBI"])

    ref_mean_map <- c("TBA" = tba_mean, "TBI" = tbi_mean)

    condition_met <- is.finite(tra_mean)
    for (ref_src in scan_reference_sources) {
      ref_mean <- ref_mean_map[[ref_src]]
      if (!is.finite(ref_mean) || !(tra_mean > ref_mean)) {
        condition_met <- FALSE
        break
      }
    }

    if (!condition_met) return(NULL)

    scan_position <- window_start + scan_window / 2

    data.table(
      chr = chr_name,
      bp_label = bp_label,
      breakpoint = bp_pos,
      direction = direction,
      scan_position = scan_position,
      relative_pos = scan_position - bp_pos,
      scan_position_Mb = scan_position / 1e6,
      relative_pos_Mb = (scan_position - bp_pos) / 1e6,
      window_start = window_start,
      window_end = window_end,
      window_start_Mb = window_start / 1e6,
      window_end_Mb = window_end / 1e6,
      tra_mean = tra_mean,
      tba_mean = tba_mean,
      tbi_mean = tbi_mean,
      n_total = nrow(window_data),
      n_TRA = sum(src == "TRA"),
      n_TBA = sum(src == "TBA"),
      n_TBI = sum(src == "TBI"),
      scan_window = scan_window,
      scan_step = scan_step,
      scan_reference_sources = paste(scan_reference_sources, collapse = ","),
      coordinate_mode = coordinate_mode
    )
  }

  # 左向扫描：从断点左侧开始，逐步向更左侧扫描，记录首次满足条件的窗口中点。
  if ("left" %in% scan_directions) {
    current_window_start <- bp_pos - scan_window
    while (current_window_start >= left_bound) {
      window_end <- current_window_start + scan_window
      hit <- evaluate_window(current_window_start, window_end, direction = "left")
      if (!is.null(hit) && nrow(hit) > 0) {
        scan_details <- rbindlist(list(scan_details, hit), use.names = TRUE, fill = TRUE)
        break
      }
      current_window_start <- current_window_start - scan_step
    }
  }

  # 右向扫描：从断点右侧开始，逐步向更右侧扫描，记录首次满足条件的窗口中点。
  if ("right" %in% scan_directions) {
    current_window_start <- bp_pos
    while (current_window_start + scan_window <= right_bound) {
      window_end <- current_window_start + scan_window
      hit <- evaluate_window(current_window_start, window_end, direction = "right")
      if (!is.null(hit) && nrow(hit) > 0) {
        scan_details <- rbindlist(list(scan_details, hit), use.names = TRUE, fill = TRUE)
        break
      }
      current_window_start <- current_window_start + scan_step
    }
  }

  scan_details[]
}

save_plot_all_formats <- function(plot_obj, filename_base, width = 12, height = 8) {
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

make_output_basename <- function(base_name) {
  if (add_mode_prefix_to_filenames) paste0(output_prefix, "_", base_name) else base_name
}

# ============================================================
# 2. 断点汇总图辅助函数
# ============================================================

aggregate_by_relative_bin <- function(df, window_size, use_abs = FALSE) {
  if (is.null(df) || nrow(df) == 0) return(data.frame())
  dt <- as.data.table(df)[!is.na(relative_pos) & !is.na(smooth_y) & is.finite(relative_pos) & is.finite(smooth_y)]
  if (nrow(dt) == 0) return(data.frame())

  if (use_abs) {
    dt[, bin_pos := round(abs(relative_pos) / window_size) * window_size]
  } else {
    dt[, bin_pos := round(relative_pos / window_size) * window_size]
  }

  res <- dt[, .(
    pos = bin_pos[1],
    mean_y = mean(smooth_y, na.rm = TRUE),
    se_y = {
      n <- sum(!is.na(smooth_y))
      if (n <= 1) 0 else sd(smooth_y, na.rm = TRUE) / sqrt(n)
    },
    n_points = .N,
    n_breakpoints = uniqueN(bp_label)
  ), by = .(source, bin_pos)]

  res[, source := factor(as.character(source), levels = source_levels)]
  setorder(res, source, pos)
  as.data.frame(res)
}

compute_common_ranges <- function(summary_dt) {
  if (is.null(summary_dt) || nrow(summary_dt) == 0) return(NULL)
  dt <- as.data.table(summary_dt)[is.finite(relative_pos)]
  if (nrow(dt) == 0) return(NULL)

  ranges <- dt[, .(
    left = min(relative_pos, na.rm = TRUE),
    right = max(relative_pos, na.rm = TRUE),
    max_abs = max(abs(relative_pos), na.rm = TRUE)
  ), by = bp_label]

  if (nrow(ranges) == 0) return(NULL)
  list(
    left = max(ranges$left, na.rm = TRUE),
    right = min(ranges$right, na.rm = TRUE),
    max_abs = min(ranges$max_abs, na.rm = TRUE)
  )
}

add_breakpoint_and_scan_lines <- function(p, scan_dt = NULL, use_abs_scan = FALSE, x_min = NULL, x_max = NULL) {
  p <- p +
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

  if (draw_scan_lines && !is.null(scan_dt) && nrow(scan_dt) > 0) {
    scan_plot <- as.data.table(scan_dt)
    if (use_abs_scan) {
      scan_plot[, plot_pos := abs(relative_pos)]
    } else {
      scan_plot[, plot_pos := relative_pos]
    }
    scan_plot <- unique(scan_plot[is.finite(plot_pos), .(plot_pos)])
    if (!is.null(x_min)) scan_plot <- scan_plot[plot_pos >= x_min]
    if (!is.null(x_max)) scan_plot <- scan_plot[plot_pos <= x_max]

    if (nrow(scan_plot) > 0) {
      p <- p +
        geom_vline(
          data = scan_plot,
          aes(xintercept = plot_pos / 1e6),
          color = publication_colors$scan,
          linetype = "dotted",
          size = 0.5,
          alpha = 0.6,
          inherit.aes = FALSE
        )
    }
  }

  p
}

make_individual_breakpoint_summary_plot <- function(summary_dt, scan_dt) {
  y_lim <- get_y_limits(summary_dt)
  summary_dt$source <- factor(as.character(summary_dt$source), levels = source_levels)

  p <- ggplot() +
    create_publication_theme() +
    labs(
      title = "All Breakpoints Summary",
      subtitle = paste0("ISMC Analysis Around Chromosomal Breakpoints (mode: ", coordinate_mode, ")"),
      x = "Relative Position to Breakpoint (Mb)",
      y = "RHO Value"
    ) +
    coord_cartesian(ylim = y_lim)

  p <- p +
    geom_ribbon(
      data = summary_dt,
      aes(
        x = relative_pos / 1e6,
        ymin = smooth_y - se_y,
        ymax = smooth_y + se_y,
        fill = source,
        group = interaction(source, bp_label)
      ),
      alpha = 0.12,
      color = NA
    ) +
    geom_line(
      data = summary_dt,
      aes(
        x = relative_pos / 1e6,
        y = smooth_y,
        color = source,
        group = interaction(source, bp_label)
      ),
      size = 0.8,
      alpha = 0.7
    ) +
    scale_color_manual(name = "Sample", values = source_color_values, labels = source_labels, drop = FALSE) +
    scale_fill_manual(name = "Sample", values = source_color_values, labels = source_labels, drop = FALSE) +
    guides(
      fill = "none",
      color = guide_legend(
        override.aes = list(alpha = 1, size = 2.5),
        title.position = "top",
        byrow = TRUE
      )
    )

  p <- add_breakpoint_and_scan_lines(p, scan_dt, use_abs_scan = FALSE)

  p <- p +
    scale_y_continuous(name = "RHO Value", labels = scales::scientific_format(digits = 2)) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.1), breaks = scales::pretty_breaks(n = 10)) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = paste0(
        "Total Breakpoints: ", uniqueN(as.data.table(summary_dt)$bp_label), "\n",
        "Window: ", format(plot_window / 1e6, digits = 1), " Mb\n",
        "Extension: ", format(summary_extension / 1e6, digits = 1), " Mb"
      ),
      hjust = 1.02,
      vjust = 1.02,
      size = 3,
      color = publication_colors$text,
      fontface = "italic",
      alpha = 0.8
    )

  p
}

make_aggregated_summary_plot <- function(agg_dt, scan_dt, filename_base, title, subtitle, x_label, use_abs_scan = FALSE, stats_extra = NULL) {
  if (is.null(agg_dt) || nrow(agg_dt) == 0) {
    cat("警告:", title, "无可用数据，跳过。\n")
    return(invisible(NULL))
  }

  agg_dt <- as.data.table(agg_dt)
  agg_dt[, source := factor(as.character(source), levels = source_levels)]
  y_lim <- get_y_limits(data.table(
    smooth_y = agg_dt$mean_y,
    se_y = agg_dt$se_y
  ))

  x_min <- min(agg_dt$pos, na.rm = TRUE)
  x_max <- max(agg_dt$pos, na.rm = TRUE)

  p <- ggplot() +
    create_publication_theme() +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = "RHO Value"
    ) +
    coord_cartesian(ylim = y_lim)

  p <- p +
    geom_ribbon(
      data = agg_dt,
      aes(x = pos / 1e6, ymin = mean_y - se_y, ymax = mean_y + se_y, fill = source, group = source),
      alpha = 0.12,
      color = NA
    ) +
    geom_line(
      data = agg_dt,
      aes(x = pos / 1e6, y = mean_y, color = source, group = source),
      size = 0.9,
      alpha = 0.9
    ) +
    scale_color_manual(name = "Sample", values = source_color_values, labels = source_labels, drop = FALSE) +
    scale_fill_manual(name = "Sample", values = source_color_values, labels = source_labels, drop = FALSE) +
    guides(
      fill = "none",
      color = guide_legend(
        override.aes = list(alpha = 1, size = 2.5),
        title.position = "top",
        byrow = TRUE
      )
    )

  p <- add_breakpoint_and_scan_lines(p, scan_dt, use_abs_scan = use_abs_scan, x_min = x_min, x_max = x_max)

  if (is.null(stats_extra)) stats_extra <- "Aggregation: mean ± SE across breakpoints"
  p <- p +
    scale_y_continuous(name = "RHO Value", labels = scales::scientific_format(digits = 2)) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.1), breaks = scales::pretty_breaks(n = 10)) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = paste0(
        "Window: ", format(plot_window / 1e6, digits = 1), " Mb\n",
        stats_extra
      ),
      hjust = 1.02,
      vjust = 1.02,
      size = 3,
      color = publication_colors$text,
      fontface = "italic",
      alpha = 0.8
    )

  save_plot_all_formats(p, filename_base, width = 14, height = 8)
  invisible(p)
}

# ============================================================
# 3. 主流程：低内存逐染色体处理
# ============================================================

cat("\n当前关键参数设置:\n")
cat("coordinate_mode:", coordinate_mode, "\n")
cat("selected_trascf:", ifelse(is.null(selected_trascf) || length(selected_trascf) == 0, "ALL", paste(selected_trascf, collapse = ", ")), "\n")
cat("selected_breakpoint_targets:", format_selected_breakpoint_targets(selected_breakpoint_targets), "\n")
cat("plot_window:", plot_window, "\n")
cat("scan_window:", scan_window, "\n")
cat("scan_step:", scan_step, "\n")
cat("summary_extension:", summary_extension, "\n")
cat("scan_reference_sources:", paste(scan_reference_sources, collapse = ", "), "\n")
cat("draw_scan_lines:", draw_scan_lines, "\n")
cat("write_scan_line_details:", write_scan_line_details, "\n")
cat("smoothing enabled:", smoothing_params$enabled, "\n")
cat("smoothing method:", smoothing_params$method, "\n")
cat("use_awk_chr_filter:", use_awk_chr_filter, "\n")
cat("data.table threads:", getDTthreads(), "\n")

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

cat("\n正在获取染色体列表。低内存模式下只提取chr列，不读入完整BED。\n")
chr_lists <- lapply(file_infos, get_chromosomes_from_file)
all_chromosomes <- sort_chromosomes(unique(unlist(chr_lists, use.names = FALSE)))

cat("染色体数量:", length(all_chromosomes), "\n")
cat("染色体列表:", paste(all_chromosomes, collapse = ", "), "\n")

# ================== 按指定TraScf染色体过滤 ==================
# 注意：这里只过滤进入统计、画图和断点扫描的染色体；
# 其余统计函数、绘图函数和扫描函数内部逻辑不做任何修改。
all_chromosomes_original <- all_chromosomes

all_chromosomes <- all_chromosomes[
  vapply(
    all_chromosomes,
    is_selected_trascf,
    logical(1),
    selected_trascf = selected_trascf
  )
]

cat("\nTraScf染色体筛选设置:\n")
cat("selected_trascf:", ifelse(is.null(selected_trascf) || length(selected_trascf) == 0, "ALL", paste(selected_trascf, collapse = ", ")), "\n")
cat("筛选前染色体数量:", length(all_chromosomes_original), "\n")
cat("筛选后染色体数量:", length(all_chromosomes), "\n")
cat("实际进入统计和绘图的染色体:", paste(all_chromosomes, collapse = ", "), "\n")

if (length(all_chromosomes) == 0) {
  stop(
    "没有任何染色体匹配 selected_trascf。请检查输入文件中的染色体名称是否类似 TraScf_7 或 7。"
  )
}

all_breakpoint_summary <- list()
all_scan_summary <- list()

for (chr_fullname in all_chromosomes) {
  chr_num <- clean_chr_name(chr_fullname)
  cat("\n============================================================\n")
  cat("正在处理染色体:", chr_fullname, "\n")

  chr_dt <- prepare_chr_plot_table(chr_fullname, file_infos)
  if (nrow(chr_dt) == 0) {
    cat("染色体", chr_fullname, "在当前坐标模式下无有效数据，跳过。\n")
    if (run_gc_each_chr) gc()
    next
  }

  cat("当前染色体处理后记录数:\n")
  print(chr_dt[, .(n_records = .N, n_coords = uniqueN(pos_key)), by = source])

  n_data_points_chr <- nrow(chr_dt)

  # 后续绘图和扫描不再需要pos_key，删除该列降低当前染色体内存占用。
  chr_dt[, pos_key := NULL]

  plot_data <- create_smoothed_plot_data(chr_dt, plot_window)
  if (nrow(plot_data) == 0) {
    cat("染色体", chr_fullname, "无足够窗口数据绘图，跳过。\n")
    rm(chr_dt, plot_data)
    if (run_gc_each_chr) gc()
    next
  }

  cat("当前染色体窗口数据点:\n")
  print(as.data.table(plot_data)[, .(n_windows = .N), by = source])

  y_lim <- get_y_limits(plot_data)

  # ================== 为汇总图收集断点数据 ==================
  scan_positions_chr <- empty_scan_detail_dt()
  if (chr_num %in% names(vertical_lines)) {
    # 断点编号按该TraScf上断点坐标从小到大排序后确定。
    marker_positions <- sort(as.numeric(vertical_lines[[chr_num]]))
    min_pos <- min(chr_dt$pos, na.rm = TRUE)
    max_pos <- max(chr_dt$pos, na.rm = TRUE)

    for (bp_idx in seq_along(marker_positions)) {
      selected_directions <- get_selected_breakpoint_directions(chr_fullname, bp_idx, selected_breakpoint_targets)
      bp_pos <- marker_positions[bp_idx]

      if (length(selected_directions) == 0) {
        cat(
          "跳过断点:", chr_fullname,
          "BP", bp_idx,
          "坐标", bp_pos,
          "（未在 selected_breakpoint_targets 中指定）\n"
        )
        next
      }

      left_bound <- max(min_pos, bp_pos - summary_extension)
      right_bound <- min(max_pos, bp_pos + summary_extension)


      bp_label <- paste0("Chr", chr_num, "_BP", if (length(marker_positions) > 1) bp_idx else "")

      cat(
        "进入分析的断点:", bp_label,
        "坐标", bp_pos,
        "方向", paste(selected_directions, collapse = ","),
        "\n"
      )

      # 汇总图只保存窗口级数据，不保存原始BED点，降低内存占用。
      region_data <- as.data.table(plot_data)[pos >= left_bound & pos <= right_bound]
      if (nrow(region_data) > 0) {
        region_data[, `:=`(
          relative_pos = pos - bp_pos,
          breakpoint = bp_pos,
          bp_label = bp_label,
          chr = chr_fullname
        )]
        region_data <- filter_region_data_by_selected_directions(region_data, selected_directions)
        if (nrow(region_data) > 0) {
          all_breakpoint_summary[[length(all_breakpoint_summary) + 1]] <- region_data
        }
      }

      scan_details_bp <- scan_breakpoint(
        chr_dt,
        chr_fullname,
        bp_label,
        bp_pos,
        left_bound,
        right_bound,
        scan_directions = selected_directions
      )
      if (!is.null(scan_details_bp) && nrow(scan_details_bp) > 0) {
        scan_positions_chr <- rbindlist(list(scan_positions_chr, scan_details_bp), use.names = TRUE, fill = TRUE)
        all_scan_summary[[length(all_scan_summary) + 1]] <- scan_details_bp

        cat("绿色扫描虚线命中:", bp_label, "\n")
        print(scan_details_bp[, .(
          direction,
          scan_position,
          relative_pos,
          scan_position_Mb,
          relative_pos_Mb,
          window_start,
          window_end,
          tra_mean,
          tba_mean,
          tbi_mean,
          n_total,
          n_TRA,
          n_TBA,
          n_TBI
        )])
      }
    }
  }
  # ================== 断点数据收集结束 ==================

  # 创建单个染色体图
  p <- ggplot() +
    create_publication_theme() +
    labs(
      title = paste("Chromosome", chr_num),
      subtitle = paste0("ISMC Analysis Comparison (mode: ", coordinate_mode, ")"),
      x = "Position (Mb)",
      y = "RHO Value"
    ) +
    coord_cartesian(ylim = y_lim)

  p <- add_main_series(p, plot_data)

  p <- p +
    scale_y_continuous(name = "RHO Value", labels = scales::scientific_format(digits = 2)) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.1), breaks = scales::pretty_breaks(n = 8))

  # 添加断点线和扫描线
  if (chr_num %in% names(vertical_lines)) {
    # 只绘制进入分析的断点；断点编号按坐标从小到大排序。
    marker_positions_all <- sort(as.numeric(vertical_lines[[chr_num]]))
    marker_keep <- vapply(
      seq_along(marker_positions_all),
      function(bp_idx) is_selected_breakpoint_target(chr_fullname, bp_idx, selected_breakpoint_targets),
      logical(1)
    )
    marker_positions <- marker_positions_all[marker_keep]

    if (length(marker_positions) > 0) {
      p <- p +
        geom_vline(
          xintercept = marker_positions / 1e6,
          color = publication_colors$marker,
          linetype = "dashed",
          size = 0.8,
          alpha = 0.7
        )

      for (i in seq_along(marker_positions)) {
        p <- p +
          annotate(
            "text",
            x = marker_positions[i] / 1e6,
            y = Inf,
            label = "Break Point",
            color = publication_colors$marker,
            vjust = 1.2,
            size = 3.5,
            fontface = "bold"
          )
      }
    }

    if (draw_scan_lines && nrow(scan_positions_chr) > 0) {
      p <- p +
        geom_vline(
          xintercept = unique(scan_positions_chr$scan_position) / 1e6,
          color = publication_colors$scan,
          linetype = "dotted",
          size = 0.6,
          alpha = 0.8
        )
    }
  }

  stats_text <- paste0(
    "Mode: ", coordinate_mode, "\n",
    "Window: ", format(plot_window / 1e6, digits = 1), " Mb\n",
    "Scan: ", format(scan_window / 1e6, digits = 1), " Mb\n",
    "Data points: ", n_data_points_chr
  )

  p <- p +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = stats_text,
      hjust = 1.02,
      vjust = 1.02,
      size = 3,
      color = publication_colors$text,
      fontface = "italic",
      alpha = 0.8
    )

  filename_base <- make_output_basename(paste0(chr_fullname, "_publication"))
  save_plot_all_formats(p, filename_base, width = 12, height = 8)

  cat("图形已保存到目录:", output_dir, "\n")
  cat("文件前缀:", filename_base, "\n")

  rm(chr_dt, plot_data, p)
  if (run_gc_each_chr) gc()
}

# ============================================================
# 4. 创建断点汇总图
# ============================================================

cat("\n==================== 开始创建断点汇总图 ====================\n")

if (length(all_breakpoint_summary) > 0) {
  breakpoint_summary <- rbindlist(all_breakpoint_summary, use.names = TRUE, fill = TRUE)
  breakpoint_summary[, source := factor(as.character(source), levels = source_levels)]

  scan_summary <- if (length(all_scan_summary) > 0) {
    rbindlist(all_scan_summary, use.names = TRUE, fill = TRUE)
  } else {
    empty_scan_detail_dt()
  }

  if (write_scan_line_details) {
    scan_line_details_path <- file.path(output_dir, scan_line_details_file)
    fwrite(scan_summary, scan_line_details_path, sep = "\t", quote = FALSE, na = "NA")
    cat("绿色扫描虚线详细数值已输出:", scan_line_details_path, "\n")
  }

  # 4.1 每个断点分别一条线的汇总图
  p_summary <- make_individual_breakpoint_summary_plot(breakpoint_summary, scan_summary)
  save_plot_all_formats(
    p_summary,
    make_output_basename("All_Breakpoints_Summary"),
    width = 14,
    height = 8
  )
  cat("断点汇总图已保存。\n")

  # 4.2 样本单线汇总图：按相对位置聚合
  cat("\n==================== 开始创建样本单线汇总图 ====================\n")
  agg_relative <- aggregate_by_relative_bin(breakpoint_summary, plot_window, use_abs = FALSE)
  make_aggregated_summary_plot(
    agg_relative,
    scan_summary,
    make_output_basename("All_Breakpoints_Summary_Aggregated"),
    title = "All Breakpoints Summary (Aggregated by Sample)",
    subtitle = paste0("Mean ISMC Across Breakpoints With SE Ribbon (mode: ", coordinate_mode, ")"),
    x_label = "Relative Position to Breakpoint (Mb)",
    use_abs_scan = FALSE,
    stats_extra = paste0("Aggregation: mean ± SE across breakpoints\nMode: ", coordinate_mode)
  )

  # 4.3 样本单线汇总图：相对位置公共范围
  cat("\n==================== 开始创建样本单线汇总图(最短范围截止) ====================\n")
  common_ranges <- compute_common_ranges(breakpoint_summary)
  if (!is.null(common_ranges) && is.finite(common_ranges$left) && is.finite(common_ranges$right) && common_ranges$left < common_ranges$right) {
    breakpoint_common <- breakpoint_summary[relative_pos >= common_ranges$left & relative_pos <= common_ranges$right]
    agg_common <- aggregate_by_relative_bin(breakpoint_common, plot_window, use_abs = FALSE)
    make_aggregated_summary_plot(
      agg_common,
      scan_summary,
      make_output_basename("All_Breakpoints_Summary_Aggregated_Common"),
      title = "All Breakpoints Summary (Aggregated by Sample, Common Range)",
      subtitle = paste0("Stop When Any Chromosome Ends (mode: ", coordinate_mode, ")"),
      x_label = "Relative Position to Breakpoint (Mb)",
      use_abs_scan = FALSE,
      stats_extra = paste0(
        "Common range: [",
        format(common_ranges$left / 1e6, digits = 2), ", ",
        format(common_ranges$right / 1e6, digits = 2), "] Mb\nMode: ", coordinate_mode
      )
    )
    rm(breakpoint_common, agg_common)
  } else {
    cat("警告: 无法计算公共范围，跳过最短范围截止汇总图。\n")
  }

  # 4.4 按绝对距离聚合
  cat("\n==================== 开始创建按绝对距离聚合汇总图 ====================\n")
  agg_abs <- aggregate_by_relative_bin(breakpoint_summary, plot_window, use_abs = TRUE)
  make_aggregated_summary_plot(
    agg_abs,
    scan_summary,
    make_output_basename("All_Breakpoints_Summary_Aggregated_Abs"),
    title = "All Breakpoints Summary (Aggregated by Absolute Distance)",
    subtitle = paste0("Mean ISMC vs |Distance to Breakpoint| With SE Ribbon (mode: ", coordinate_mode, ")"),
    x_label = "Distance to Breakpoint (Mb)",
    use_abs_scan = TRUE,
    stats_extra = paste0("Aggregation: mean ± SE by |distance|\nMode: ", coordinate_mode)
  )

  # 4.5 按绝对距离聚合：公共最大距离
  cat("\n==================== 开始创建按绝对距离最短范围截止汇总图 ====================\n")
  if (!is.null(common_ranges) && is.finite(common_ranges$max_abs) && common_ranges$max_abs > 0) {
    breakpoint_abs_common <- breakpoint_summary[abs(relative_pos) <= common_ranges$max_abs]
    agg_abs_common <- aggregate_by_relative_bin(breakpoint_abs_common, plot_window, use_abs = TRUE)
    make_aggregated_summary_plot(
      agg_abs_common,
      scan_summary,
      make_output_basename("All_Breakpoints_Summary_Aggregated_Abs_Common"),
      title = "All Breakpoints Summary (Aggregated by Absolute Distance, Common Range)",
      subtitle = paste0("Stop When Any Chromosome Ends (mode: ", coordinate_mode, ")"),
      x_label = "Distance to Breakpoint (Mb)",
      use_abs_scan = TRUE,
      stats_extra = paste0(
        "Common max |distance|: ", format(common_ranges$max_abs / 1e6, digits = 2), " Mb\nMode: ", coordinate_mode
      )
    )
    rm(breakpoint_abs_common, agg_abs_common)
  } else {
    cat("警告: 无法计算绝对距离公共范围，跳过最短范围截止绝对距离汇总图。\n")
  }

  rm(breakpoint_summary, scan_summary, agg_relative, agg_abs)
  if (run_gc_each_chr) gc()

} else {
  cat("\n警告: 未找到任何断点数据，无法创建汇总图。\n")
}

# ============================================================
# 5. 结束提示
# ============================================================

cat("\n=== 分析完成 ===\n")
cat("低内存版本已完成。\n")
cat("主要内存优化：\n")
cat("1. 不再一次性读入4个完整BED文件。\n")
cat("2. 每次只用awk筛选并读入一个染色体的数据。\n")
cat("3. 每个染色体处理完立即释放原始点数据。\n")
cat("4. 断点汇总图只保存窗口级数据，不保存全基因组原始BED点。\n")
cat("输出目录:", output_dir, "\n")
if (write_scan_line_details) {
  cat("绿色扫描虚线详细数值文件:", file.path(output_dir, scan_line_details_file), "\n")
}
cat("建议使用PDF格式用于出版物投稿；TIFF用于印刷质量要求；PNG适合网络展示和演示。\n")

cat("\n=== 坐标模式说明 ===\n")
cat("- union_tra: TBA和TBI的坐标各自保留，两个TRA文件坐标取并集，TRA重复坐标取平均。\n")
cat("- intersection_4: 只保留4个BED文件共有坐标，两个TRA文件在共有坐标上的值取平均。\n")

cat("\n=== 绿色扫描虚线说明 ===\n")
cat("- draw_scan_lines = TRUE/FALSE 控制是否在图上绘制绿色扫描虚线。\n")
cat("- write_scan_line_details = TRUE/FALSE 控制是否输出绿色扫描虚线详细数值表。\n")
cat("- 详细数值表中的scan_position为绿色虚线的绝对坐标，relative_pos为相对断点坐标。\n")
cat("- direction表示从断点左侧或右侧扫描时首次满足条件的窗口。\n")
cat("- selected_breakpoint_targets 可控制只分析指定TraScf、指定断点编号和指定方向。\n")
cat("- bp_index 按每个TraScf中断点坐标从小到大排序后编号。\n")

cat("\n=== 平滑度调整说明 ===\n")
cat("1. smoothing_params$enabled = TRUE 后可启用平滑。\n")
cat("2. method可选: 'loess', 'spline', 'moving_average', 'none'。\n")
cat("3. loess的span越小越贴近原始窗口均值，越大越平滑。\n")
