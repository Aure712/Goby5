#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
  library(grid) # unit()
})

# =========================
# 0) USER CONTROLS
# =========================
CATEGORY_ORDER <- c(
  "non_invert_non_ismc",
  "invert_non_ismc",
  "non_invert_ismc",
  "invert_ismc"
)

BASELINE_CATEGORY <- "non_invert_non_ismc"

FIG_DPI <- 600

# 散点控制：
# = 0  -> 完全不画散点
# < 0  -> 画全部散点
# > 0  -> 每组最多画这么多个散点
MAX_POINTS_PER_GROUP <- 0

POINT_ALPHA <- 0.10
POINT_SIZE <- 0.40
JITTER_W <- 0.10

# 更紧凑、更饱满
BOXPLOT_WIDTH <- 0.24
VIOLIN_WIDTH <- 0.95

PALETTE <- c(
  "non_invert_non_ismc" = "#0072B2",
  "invert_non_ismc" = "#009E73",
  "non_invert_ismc" = "#E69F00",
  "invert_ismc" = "#D55E00"
)

LABEL_MAP <- c(
  "non_invert_non_ismc" = "non-inv & non-FSR",
  "invert_non_ismc" = "inv & non-FSR",
  "non_invert_ismc" = "non-inv & FSR",
  "invert_ismc" = "inv & FSR"
)

# =========================
# 视觉强调设置
# =========================
STAT_LABEL_DIGITS <- 4
SHOW_MEDIAN_LINE <- TRUE
SHOW_MEDIAN_TEXT <- TRUE
SHOW_MEAN_LINE <- FALSE
SHOW_MEAN_TEXT <- FALSE
DRAW_STAT_POINTS <- FALSE

STAT_TEXT_SIZE <- 3.0
STAT_TEXT_COLOR_MEDIAN <- "#6D4C41"
STAT_TEXT_COLOR_MEAN <- "#37474F"
STAT_LABEL_NUDGE_X <- 0.18
MEDIAN_TEXT_OFFSET_Y <- 0.00
MEAN_TEXT_OFFSET_Y <- 0.00
MIN_SEP_BETWEEN_MED_MEAN_TEXT <- 0.06

SHOW_STAT_LINES <- TRUE
STAT_LINE_HALF_HEIGHT <- BOXPLOT_WIDTH / 2
STAT_LINE_MED_COLOR <- "#9edae5"
STAT_LINE_MEAN_COLOR <- "#00838F"
STAT_LINE_MED_LWD <- 0.95
STAT_LINE_MEAN_LWD <- 0.95
STAT_LINE_MED_LTY <- "solid"
STAT_LINE_MEAN_LTY <- "solid"
STAT_LINE_END <- "butt"

# =========================
# Peak（密度峰值）标记
# =========================
SHOW_PEAK_LINE <- FALSE
SHOW_PEAK_TEXT <- FALSE
PEAK_USE_DISPLAY_RANGE <- FALSE

PEAK_LABEL_DIGITS <- 4
PEAK_TEXT_SIZE <- 3.0
PEAK_TEXT_COLOR <- "#7B1FA2"
PEAK_LABEL_NUDGE_X <- 0.18
PEAK_TEXT_OFFSET_Y <- 0.00

PEAK_LINE_HALF_HEIGHT <- BOXPLOT_WIDTH / 2
PEAK_LINE_COLOR <- "#7B1FA2"
PEAK_LINE_LWD <- 1.00
PEAK_LINE_LTY <- "solid"

# =========================
# 当 SHOW_MEAN_LINE = TRUE 且 SHOW_MEDIAN_LINE = FALSE 时，
# 自动隐藏 boxplot 自带的中位数线（用箱体同色覆盖）
# =========================
HIDE_BOXPLOT_MEDIAN_WHEN_MEAN_ONLY <- TRUE
BOXPLOT_MEDIAN_COVER_LWD <- 1.20
BOXPLOT_FILL_ALPHA <- 0.75

# =========================
# 对照组参考虚线（贯穿所有组别）
# =========================
CONTROL_REF_MODE <- "median"   # "none", "median", "mean", "both"
CONTROL_REF_USE_DISPLAY_RANGE <- FALSE

CONTROL_REF_MEDIAN_COLOR <- "#8E24AA"
CONTROL_REF_MEAN_COLOR <- "#3949AB"
CONTROL_REF_MEDIAN_LTY <- "dashed"
CONTROL_REF_MEAN_LTY <- "dashed"
CONTROL_REF_MEDIAN_LWD <- 0.50
CONTROL_REF_MEAN_LWD <- 0.50

# =========================
# Display-range controls
# =========================
# 1.00 = full range
# 0.95 = middle 95%
# 0.90 = middle 90%
VIOLIN_MIDDLE_PROP <- 0.90
BOXPLOT_MIDDLE_PROP <- 0.90

POINTS_FOLLOW_VIOLIN_RANGE <- TRUE
STATS_USE_DISPLAY_RANGE <- FALSE

# 两两比较的统计检验是否基于显示后的范围
# FALSE = 基于完整数据（推荐）
# TRUE  = 基于显示后的数据
PAIRWISE_USE_DISPLAY_RANGE <- FALSE

# =========================
# Pairwise significance controls
# =========================
SHOW_PAIRWISE_SIGNIF <- TRUE
PAIRWISE_TEST <- "wilcox"
PAIRWISE_P_ADJUST <- "BH"
PAIRWISE_PLOT_MODE <- "baseline_only"
PAIRWISE_SHOW_NS <- TRUE
PAIRWISE_LABEL_MODE <- "p"

PAIRWISE_P_DIGITS <- 3
PAIRWISE_P_FORMAT <- "g"
PAIRWISE_P_SCI_DIGITS <- 2
PAIRWISE_P_MAX_CHARS <- 6
PAIRWISE_LABEL_MULTILINE <- TRUE

PAIRWISE_LINE_COLOR <- "#333333"
PAIRWISE_LINE_LWD <- 0.35
PAIRWISE_TEXT_COLOR <- "#333333"
PAIRWISE_TEXT_SIZE <- 2.2

# 关键：统计可基于完整数据，但括号布局基于当前图真正显示的数据范围
PAIRWISE_POSITION_USE_DISPLAY_RANGE <- TRUE

# p 值括号区域：放在显示数据范围之外，但只占一个很小的专用右侧带
PAIRWISE_TOP_MARGIN_FRAC <- 0.025
PAIRWISE_STEP_FRAC <- 0.045
PAIRWISE_BRACKET_DEPTH_FRAC <- 0.010
PAIRWISE_LABEL_GAP_FRAC <- 0.006

# 同一“中点”附近标签轻微错开
PAIRWISE_TEXT_STAGGER <- 0.08

# 画图区留白：非常小
PLOT_EXPAND_LEFT <- 0.005
PLOT_EXPAND_RIGHT <- 0.005

# =========================
# 1) INPUT FILES
# =========================
bed_non_invert <- "all-invert.bed"
bed_non_ismc <- "all-fusion.bed"
bed_ismc <- "ismc.bed"
bed_invert <- "syn_invert3.bed"

pi_tba_file <- "Tba.refTra.pi.bed"
pi_tra_file <- "Tra.refTra.pi.bed"
dxy_file <- "Tba.Tra.dxy.bed"
fst_file <- "Tba.Tra.fst.bed"

out_windows <- "window_stats_by_region.tsv"
out_summary <- "summary_by_region.tsv"
out_tests <- "tests_by_region.txt"
out_dropped <- "dropped_windows.tsv"
out_pairwise_table <- "pairwise_all_comparisons.tsv"
out_plot_pairwise_table <- "pairwise_plot_labels.tsv"

out_plot_pi_tra_png <- "fig_mbe_pi_tra.png"
out_plot_pi_tra_pdf <- "fig_mbe_pi_tra.pdf"
out_plot_pi_tba_png <- "fig_mbe_pi_tba.png"
out_plot_pi_tba_pdf <- "fig_mbe_pi_tba.pdf"
out_plot_pi_both_png <- "fig_mbe_pi_both.png"
out_plot_pi_both_pdf <- "fig_mbe_pi_both.pdf"
out_plot_dxy_png <- "fig_mbe_dxy.png"
out_plot_dxy_pdf <- "fig_mbe_dxy.pdf"
out_plot_fst_png <- "fig_mbe_fst.png"
out_plot_fst_pdf <- "fig_mbe_fst.pdf"

# =========================
# 2) Helpers
# =========================
read_region_bed <- function(path) {
  dt <- fread(path, header = FALSE)
  if (ncol(dt) < 3) stop("Region BED needs >=3 columns: chr start end: ", path)
  setnames(dt, 1:3, c("chr", "start", "end"))
  dt <- dt[, .(
    chr = as.character(chr),
    start = as.integer(start),
    end = as.integer(end)
  )]
  dt <- dt[end > start]
  GRanges(seqnames = dt$chr, ranges = IRanges(start = dt$start + 1L, end = dt$end))
}

read_window_bed <- function(path, value_name) {
  dt <- fread(path, header = FALSE)
  if (ncol(dt) < 4) stop("Window BED needs >=4 columns: chr start end value: ", path)
  setnames(dt, 1:4, c("chr", "start", "end", value_name))
  dt <- dt[, .(
    chr = as.character(chr),
    start = as.integer(start),
    end = as.integer(end),
    v = as.numeric(get(value_name))
  )]
  setnames(dt, "v", value_name)
  dt <- dt[end > start]
  dt
}

label_by_midpoint <- function(mid_gr, region_gr) {
  hits <- findOverlaps(mid_gr, region_gr, ignore.strand = TRUE)
  flag <- rep(FALSE, length(mid_gr))
  flag[queryHits(hits)] <- TRUE
  flag
}

cliffs_delta <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  nx <- length(x)
  ny <- length(y)
  if (nx == 0 || ny == 0) return(NA_real_)
  r <- rank(c(x, y), ties.method = "average")
  rx <- sum(r[seq_len(nx)])
  ux <- rx - nx * (nx + 1) / 2
  (2 * ux) / (nx * ny) - 1
}

block_jackknife_median_diff <- function(dt, groupA, groupB, value_col, block_col = "chr") {
  blocks <- unique(dt[[block_col]])
  diffs <- numeric(0)
  
  for (b in blocks) {
    sub <- dt[dt[[block_col]] != b, ]
    a <- sub[sub$category == groupA, ][[value_col]]
    bb <- sub[sub$category == groupB, ][[value_col]]
    if (length(a) > 0 && length(bb) > 0) {
      diffs <- c(diffs, median(a, na.rm = TRUE) - median(bb, na.rm = TRUE))
    }
  }
  
  if (length(diffs) < 5) {
    return(list(est = NA_real_, se = NA_real_, n_blocks = length(diffs)))
  }
  
  est <- median(dt[category == groupA][[value_col]], na.rm = TRUE) -
    median(dt[category == groupB][[value_col]], na.rm = TRUE)
  se <- sqrt((length(diffs) - 1) * mean((diffs - mean(diffs))^2))
  list(est = est, se = se, n_blocks = length(diffs))
}

theme_mbe <- function() {
  theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11, color = "black"),
      axis.line = element_line(linewidth = 0.6, color = "black"),
      axis.ticks = element_line(linewidth = 0.6, color = "black"),
      axis.ticks.length = unit(2.2, "pt"),
      panel.border = element_rect(fill = NA, linewidth = 0.6, color = "black"),
      plot.margin = margin(4, 4, 4, 4),
      legend.position = "none"
    )
}

subsample_points <- function(d, group_col) {
  set.seed(1)
  
  if (MAX_POINTS_PER_GROUP == 0) {
    return(d[0])
  }
  
  d[, {
    n <- .N
    if (MAX_POINTS_PER_GROUP < 0 || n <= MAX_POINTS_PER_GROUP) {
      .SD
    } else {
      .SD[sample.int(n, MAX_POINTS_PER_GROUP)]
    }
  }, by = group_col]
}

middle_prop_to_probs <- function(middle_prop, arg_name = "middle_prop") {
  if (length(middle_prop) != 1 || !is.finite(middle_prop) || middle_prop <= 0 || middle_prop > 1) {
    stop(arg_name, " must be a single number in (0, 1].")
  }
  lower <- (1 - middle_prop) / 2
  c(lower, 1 - lower)
}

filter_middle_range_by_group <- function(d, value_col, group_cols, middle_prop = 1.00) {
  x <- copy(d)
  if (nrow(x) == 0) return(x)
  if (middle_prop >= 1) return(x)
  
  probs <- middle_prop_to_probs(middle_prop, "middle_prop")
  
  qdt <- x[is.finite(get(value_col)), .(
    q_low = as.numeric(quantile(get(value_col), probs[1], na.rm = TRUE, names = FALSE, type = 7)),
    q_high = as.numeric(quantile(get(value_col), probs[2], na.rm = TRUE, names = FALSE, type = 7))
  ), by = group_cols]
  
  x <- merge(x, qdt, by = group_cols, all.x = TRUE, sort = FALSE)
  x <- x[is.finite(get(value_col)) & get(value_col) >= q_low & get(value_col) <= q_high]
  x[, c("q_low", "q_high") := NULL]
  x[]
}

build_stat_text_positions <- function(stat_dt, log10_scale = FALSE) {
  out <- copy(stat_dt)
  
  if (log10_scale) {
    out[, `:=`(
      med_log = log10(med),
      mean_log = log10(mean)
    )]
    out[, med_label_log := med_log + MEDIAN_TEXT_OFFSET_Y]
    out[, mean_label_log := mean_log + MEAN_TEXT_OFFSET_Y]
    
    if (MIN_SEP_BETWEEN_MED_MEAN_TEXT > 0) {
      idx <- abs(out$mean_label_log - out$med_label_log) < MIN_SEP_BETWEEN_MED_MEAN_TEXT
      if (any(idx)) {
        mid <- (out$mean_label_log[idx] + out$med_label_log[idx]) / 2
        out$med_label_log[idx] <- mid - MIN_SEP_BETWEEN_MED_MEAN_TEXT / 2
        out$mean_label_log[idx] <- mid + MIN_SEP_BETWEEN_MED_MEAN_TEXT / 2
      }
    }
    
    out[, `:=`(
      med_label_y = 10^med_label_log,
      mean_label_y = 10^mean_label_log
    )]
    out[, c("med_log", "mean_log", "med_label_log", "mean_label_log") := NULL]
  } else {
    out[, med_label_y := med + MEDIAN_TEXT_OFFSET_Y]
    out[, mean_label_y := mean + MEAN_TEXT_OFFSET_Y]
    
    if (MIN_SEP_BETWEEN_MED_MEAN_TEXT > 0) {
      idx <- abs(out$mean_label_y - out$med_label_y) < MIN_SEP_BETWEEN_MED_MEAN_TEXT
      if (any(idx)) {
        mid <- (out$mean_label_y[idx] + out$med_label_y[idx]) / 2
        out$med_label_y[idx] <- mid - MIN_SEP_BETWEEN_MED_MEAN_TEXT / 2
        out$mean_label_y[idx] <- mid + MIN_SEP_BETWEEN_MED_MEAN_TEXT / 2
      }
    }
  }
  
  out
}

build_peak_text_positions <- function(peak_dt, log10_scale = FALSE) {
  out <- copy(peak_dt)
  
  if (log10_scale) {
    out[, peak_log := log10(peak)]
    out[, peak_label_log := peak_log + PEAK_TEXT_OFFSET_Y]
    out[, peak_label_y := 10^peak_label_log]
    out[, c("peak_log", "peak_label_log") := NULL]
  } else {
    out[, peak_label_y := peak + PEAK_TEXT_OFFSET_Y]
  }
  
  out
}

p_to_stars <- function(p) {
  if (!is.finite(p)) return("ns")
  if (p <= 0.001) "***" else if (p <= 0.01) "**" else if (p <= 0.05) "*" else "ns"
}

format_p_value_short <- function(p) {
  if (!is.finite(p)) return("p=NA")
  if (p == 0) return("p=0")
  s <- formatC(p, digits = PAIRWISE_P_DIGITS, format = PAIRWISE_P_FORMAT)
  if (nchar(s) > PAIRWISE_P_MAX_CHARS) {
    s <- formatC(p, digits = PAIRWISE_P_SCI_DIGITS, format = "e")
  }
  s <- sub("^0\\.", ".", s)
  s <- sub("^(-)0\\.", "\\1.", s)
  paste0("p=", s)
}

build_pairwise_label <- function(p, stars) {
  mode <- match.arg(PAIRWISE_LABEL_MODE, c("p", "stars", "both"))
  ptxt <- format_p_value_short(p)
  
  if (mode == "p") {
    return(ptxt)
  } else if (mode == "stars") {
    return(stars)
  } else {
    if (PAIRWISE_LABEL_MULTILINE) {
      return(paste0(stars, "\n", ptxt))
    } else {
      return(paste0(stars, " (", ptxt, ")"))
    }
  }
}

add_pairwise_text_stagger <- function(pair_dt) {
  if (nrow(pair_dt) == 0) return(pair_dt)
  
  out <- copy(pair_dt)
  out[, xmid_base := (x + xend) / 2]
  out[, `:=`(dup_id = seq_len(.N), dup_n = .N), by = xmid_base]
  out[, x_label := xmid_base]
  out[dup_n > 1, x_label := xmid_base + (dup_id - (dup_n + 1) / 2) * PAIRWISE_TEXT_STAGGER]
  out[, c("dup_id", "dup_n") := NULL]
  out
}

should_hide_boxplot_median <- function() {
  isTRUE(HIDE_BOXPLOT_MEDIAN_WHEN_MEAN_ONLY) &&
    isTRUE(SHOW_MEAN_LINE) &&
    !isTRUE(SHOW_MEDIAN_LINE)
}

build_box_median_cover <- function(d_box, value_col, facet_col = NULL) {
  if (nrow(d_box) == 0) return(data.table())
  
  by_cols <- c("category_lab", "category")
  if (!is.null(facet_col)) by_cols <- c(by_cols, facet_col)
  
  out <- d_box[, .(
    box_med = median(get(value_col), na.rm = TRUE)
  ), by = by_cols]
  
  out[, x_num := as.numeric(category_lab)]
  out[, `:=`(
    x_bottom = x_num - BOXPLOT_WIDTH / 2,
    x_top = x_num + BOXPLOT_WIDTH / 2
  )]
  
  out[, cover_color := grDevices::adjustcolor(
    PALETTE[as.character(category)],
    alpha.f = BOXPLOT_FILL_ALPHA
  )]
  
  out[]
}

compute_density_peak <- function(x, log10_scale = FALSE, adjust = 1) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  
  if (log10_scale) {
    x <- x[x > 0]
    if (length(x) == 0) return(NA_real_)
    tx <- log10(x)
  } else {
    tx <- x
  }
  
  tx <- tx[is.finite(tx)]
  if (length(tx) == 0) return(NA_real_)
  
  u <- unique(tx)
  if (length(u) == 1) {
    return(if (log10_scale) 10^u[1] else u[1])
  }
  
  den <- density(tx, adjust = adjust, na.rm = TRUE)
  peak_t <- den$x[which.max(den$y)]
  
  if (log10_scale) 10^peak_t else peak_t
}

build_peak_table <- function(d, value_col, log10_scale = FALSE, facet_col = NULL) {
  if (nrow(d) == 0) return(data.table())
  
  by_cols <- c("category_lab", "category")
  if (!is.null(facet_col)) by_cols <- c(by_cols, facet_col)
  
  out <- d[, .(
    peak = compute_density_peak(get(value_col), log10_scale = log10_scale)
  ), by = by_cols]
  
  out <- out[is.finite(peak)]
  if (nrow(out) == 0) return(out)
  
  out[, x_num := as.numeric(category_lab)]
  out[, `:=`(
    x_bottom = x_num - PEAK_LINE_HALF_HEIGHT,
    x_top = x_num + PEAK_LINE_HALF_HEIGHT
  )]
  out[, peak_label := sprintf(paste0("peak=%.", PEAK_LABEL_DIGITS, "g"), peak)]
  out[]
}

compute_pairwise_table_core <- function(d, value_col, group_col, group_levels, log10_scale = FALSE) {
  sub <- copy(d[is.finite(get(value_col))])
  if (log10_scale) {
    sub <- sub[get(value_col) > 0]
  }
  if (nrow(sub) == 0) return(data.table())
  
  group_levels <- as.character(group_levels)
  present <- group_levels[group_levels %in% unique(as.character(sub[[group_col]]))]
  if (length(present) < 2) return(data.table())
  
  combs <- t(combn(present, 2))
  out_list <- vector("list", nrow(combs))
  
  for (i in seq_len(nrow(combs))) {
    g1 <- combs[i, 1]
    g2 <- combs[i, 2]
    
    x <- sub[as.character(get(group_col)) == g1, get(value_col)]
    y <- sub[as.character(get(group_col)) == g2, get(value_col)]
    
    p_raw <- NA_real_
    if (length(x) > 0 && length(y) > 0) {
      if (PAIRWISE_TEST == "wilcox") {
        p_raw <- suppressWarnings(wilcox.test(x, y)$p.value)
      } else if (PAIRWISE_TEST == "t") {
        p_raw <- suppressWarnings(t.test(x, y)$p.value)
      } else {
        stop("Unknown PAIRWISE_TEST: ", PAIRWISE_TEST)
      }
    }
    
    out_list[[i]] <- data.table(
      group1 = g1,
      group2 = g2,
      n1 = sum(is.finite(x)),
      n2 = sum(is.finite(y)),
      median1 = median(x, na.rm = TRUE),
      median2 = median(y, na.rm = TRUE),
      mean1 = mean(x, na.rm = TRUE),
      mean2 = mean(y, na.rm = TRUE),
      p_raw = p_raw
    )
  }
  
  out <- rbindlist(out_list, fill = TRUE)
  out[, p_adj := p.adjust(p_raw, method = PAIRWISE_P_ADJUST)]
  out[, stars := vapply(p_adj, p_to_stars, character(1))]
  out[, `:=`(
    median_diff = median1 - median2,
    mean_diff = mean1 - mean2
  )]
  out[, label := mapply(build_pairwise_label, p_adj, stars, USE.NAMES = FALSE)]
  out[]
}

compute_pairwise_table <- function(d, value_col, group_col, group_levels,
                                   log10_scale = FALSE, facet_col = NULL, metric_name = NULL) {
  if (is.null(facet_col)) {
    out <- compute_pairwise_table_core(d, value_col, group_col, group_levels, log10_scale = log10_scale)
  } else {
    out <- d[, compute_pairwise_table_core(.SD, value_col, group_col, group_levels, log10_scale = log10_scale), by = facet_col]
  }
  
  if (!is.null(metric_name) && nrow(out) > 0) {
    out[, metric := metric_name]
  }
  out[]
}

filter_pairwise_for_plot <- function(pair_dt,
                                     plot_mode = PAIRWISE_PLOT_MODE,
                                     baseline_group = NULL) {
  mode <- match.arg(plot_mode, c("baseline_only", "all"))
  out <- copy(pair_dt)
  
  if (mode == "baseline_only") {
    if (is.null(baseline_group)) {
      stop("baseline_group must not be NULL when plot_mode = 'baseline_only'.")
    }
    out <- out[group1 == baseline_group | group2 == baseline_group]
  }
  
  if (!PAIRWISE_SHOW_NS) {
    out <- out[stars != "ns"]
  }
  
  out[]
}

# stats_d：用于统计检验
# position_d：用于图上括号布局
build_pairwise_annotations <- function(stats_d, position_d, value_col, group_col, group_levels,
                                       log10_scale,
                                       plot_mode = PAIRWISE_PLOT_MODE,
                                       baseline_group = NULL) {
  sub_stats <- copy(stats_d[is.finite(get(value_col))])
  sub_pos <- copy(position_d[is.finite(get(value_col))])
  
  if (log10_scale) {
    sub_stats <- sub_stats[get(value_col) > 0]
    sub_pos <- sub_pos[get(value_col) > 0]
  }
  
  if (nrow(sub_stats) == 0 || nrow(sub_pos) == 0) return(data.table())
  
  if (log10_scale) {
    sub_pos[, v_trans := log10(get(value_col))]
  } else {
    sub_pos[, v_trans := get(value_col)]
  }
  
  group_levels <- as.character(group_levels)
  
  out <- compute_pairwise_table_core(
    sub_stats,
    value_col = value_col,
    group_col = group_col,
    group_levels = group_levels,
    log10_scale = FALSE
  )
  
  out <- filter_pairwise_for_plot(
    out,
    plot_mode = plot_mode,
    baseline_group = baseline_group
  )
  
  if (nrow(out) == 0) return(data.table())
  
  out[, `:=`(
    x = match(group1, group_levels),
    xend = match(group2, group_levels)
  )]
  
  out[, span := abs(xend - x)]
  setorder(out, span, x, xend)
  
  rng <- range(sub_pos$v_trans, na.rm = TRUE)
  span_val <- diff(rng)
  if (!is.finite(span_val) || span_val == 0) span_val <- 1
  
  top_margin <- PAIRWISE_TOP_MARGIN_FRAC * span_val
  step <- PAIRWISE_STEP_FRAC * span_val
  depth <- PAIRWISE_BRACKET_DEPTH_FRAC * span_val
  gap <- PAIRWISE_LABEL_GAP_FRAC * span_val
  
  out[, y_trans := max(rng) + top_margin + (seq_len(.N) - 1) * step]
  
  if (log10_scale) {
    out[, `:=`(
      y = 10^y_trans,
      y_tip = 10^(y_trans - depth),
      y_text = 10^(y_trans + gap)
    )]
  } else {
    out[, `:=`(
      y = y_trans,
      y_tip = y_trans - depth,
      y_text = y_trans + gap
    )]
  }
  
  out <- add_pairwise_text_stagger(out)
  out[]
}

pairwise_annotations <- function(stats_d, position_d, value_col, group_col, group_levels,
                                 log10_scale = FALSE, facet_col = NULL,
                                 plot_mode = PAIRWISE_PLOT_MODE,
                                 baseline_group = NULL) {
  if (is.null(facet_col)) {
    return(build_pairwise_annotations(
      stats_d = stats_d,
      position_d = position_d,
      value_col = value_col,
      group_col = group_col,
      group_levels = group_levels,
      log10_scale = log10_scale,
      plot_mode = plot_mode,
      baseline_group = baseline_group
    ))
  }
  
  stat_list <- split(stats_d, by = facet_col, keep.by = TRUE)
  pos_list <- split(position_d, by = facet_col, keep.by = TRUE)
  
  keys <- intersect(names(stat_list), names(pos_list))
  if (length(keys) == 0) return(data.table())
  
  res <- rbindlist(lapply(keys, function(k) {
    build_pairwise_annotations(
      stats_d = stat_list[[k]],
      position_d = pos_list[[k]],
      value_col = value_col,
      group_col = group_col,
      group_levels = group_levels,
      log10_scale = log10_scale,
      plot_mode = plot_mode,
      baseline_group = baseline_group
    )[, (facet_col) := k]
  }), fill = TRUE)
  
  res[]
}

add_pairwise_layers <- function(p, pair_dt) {
  if (nrow(pair_dt) == 0) return(p)
  
  p +
    geom_segment(
      data = pair_dt,
      aes(x = x, xend = xend, y = y, yend = y),
      inherit.aes = FALSE,
      colour = PAIRWISE_LINE_COLOR,
      linewidth = PAIRWISE_LINE_LWD,
      lineend = "butt"
    ) +
    geom_segment(
      data = pair_dt,
      aes(x = x, xend = x, y = y_tip, yend = y),
      inherit.aes = FALSE,
      colour = PAIRWISE_LINE_COLOR,
      linewidth = PAIRWISE_LINE_LWD,
      lineend = "butt"
    ) +
    geom_segment(
      data = pair_dt,
      aes(x = xend, xend = xend, y = y_tip, yend = y),
      inherit.aes = FALSE,
      colour = PAIRWISE_LINE_COLOR,
      linewidth = PAIRWISE_LINE_LWD,
      lineend = "butt"
    ) +
    geom_text(
      data = pair_dt,
      aes(x = x_label, y = y_text, label = label),
      inherit.aes = FALSE,
      colour = PAIRWISE_TEXT_COLOR,
      size = PAIRWISE_TEXT_SIZE,
      hjust = 0,
      vjust = 0.5,
      lineheight = 0.9
    )
}

# 直接把“图上真实使用的 pair_dt”导出成表
extract_plot_pairwise_table <- function(pair_dt, plot_name, facet_col = NULL) {
  out <- copy(pair_dt)
  
  if (nrow(out) == 0) {
    return(data.table(
      plot_name = character(),
      facet = character(),
      group1 = character(),
      group1_label = character(),
      group2 = character(),
      group2_label = character(),
      n1 = integer(),
      n2 = integer(),
      median1 = numeric(),
      median2 = numeric(),
      median_diff = numeric(),
      mean1 = numeric(),
      mean2 = numeric(),
      mean_diff = numeric(),
      p_raw = numeric(),
      p_adj = numeric(),
      stars = character(),
      label = character()
    ))
  }
  
  out[, plot_name := plot_name]
  
  if (!is.null(facet_col) && facet_col %in% names(out)) {
    out[, facet := as.character(get(facet_col))]
  } else {
    out[, facet := NA_character_]
  }
  
  out[, `:=`(
    group1_label = as.character(group1),
    group2_label = as.character(group2)
  )]
  
  keep_cols <- c(
    "plot_name", "facet",
    "group1", "group1_label",
    "group2", "group2_label",
    "n1", "n2",
    "median1", "median2", "median_diff",
    "mean1", "mean2", "mean_diff",
    "p_raw", "p_adj", "stars", "label"
  )
  
  out[, ..keep_cols]
}

build_control_reference <- function(d, value_col, group_col, baseline_group,
                                    log10_scale = FALSE, facet_col = NULL) {
  sub <- copy(d[is.finite(get(value_col))])
  if (log10_scale) {
    sub <- sub[get(value_col) > 0]
  }
  if (nrow(sub) == 0) return(data.table())
  
  sub <- sub[as.character(get(group_col)) == baseline_group]
  if (nrow(sub) == 0) return(data.table())
  
  if (is.null(facet_col)) {
    out <- sub[, .(
      ref_median = median(get(value_col), na.rm = TRUE),
      ref_mean = mean(get(value_col), na.rm = TRUE)
    )]
  } else {
    out <- sub[, .(
      ref_median = median(get(value_col), na.rm = TRUE),
      ref_mean = mean(get(value_col), na.rm = TRUE)
    ), by = facet_col]
  }
  out[]
}

add_control_reference_layers <- function(p, ref_dt) {
  if (nrow(ref_dt) == 0) return(p)
  
  mode <- match.arg(CONTROL_REF_MODE, c("none", "median", "mean", "both"))
  if (mode == "none") return(p)
  
  if (mode %in% c("median", "both")) {
    p <- p + geom_hline(
      data = ref_dt,
      aes(yintercept = ref_median),
      inherit.aes = FALSE,
      colour = CONTROL_REF_MEDIAN_COLOR,
      linetype = CONTROL_REF_MEDIAN_LTY,
      linewidth = CONTROL_REF_MEDIAN_LWD
    )
  }
  
  if (mode %in% c("mean", "both")) {
    p <- p + geom_hline(
      data = ref_dt,
      aes(yintercept = ref_mean),
      inherit.aes = FALSE,
      colour = CONTROL_REF_MEAN_COLOR,
      linetype = CONTROL_REF_MEAN_LTY,
      linewidth = CONTROL_REF_MEAN_LWD
    )
  }
  
  p
}

is_log_metric <- function(metric_name) {
  metric_name %chin% c("pi_tra", "pi_tba", "dxy")
}

plot_metric_mbe <- function(dt, value_col, title, ylab,
                            log10 = FALSE, hline0 = FALSE,
                            out_png, out_pdf,
                            width = 5.8, height = 4.2,
                            plot_name = NULL) {
  d <- copy(dt)
  d <- d[is.finite(get(value_col))]
  if (log10) d <- d[get(value_col) > 0]
  
  d <- d[category %chin% CATEGORY_ORDER]
  d[, category := factor(category, levels = CATEGORY_ORDER)]
  d[, category_lab := LABEL_MAP[as.character(category)]]
  d[, category_lab := factor(category_lab, levels = LABEL_MAP[CATEGORY_ORDER])]
  
  baseline_lab <- unname(LABEL_MAP[BASELINE_CATEGORY])
  
  d_violin <- filter_middle_range_by_group(
    d, value_col = value_col, group_cols = "category_lab",
    middle_prop = VIOLIN_MIDDLE_PROP
  )
  
  d_box <- filter_middle_range_by_group(
    d, value_col = value_col, group_cols = "category_lab",
    middle_prop = BOXPLOT_MIDDLE_PROP
  )
  
  pts_source <- if (POINTS_FOLLOW_VIOLIN_RANGE) d_violin else d
  pts <- subsample_points(
    pts_source[, .(category, category_lab, v = get(value_col))],
    "category_lab"
  )
  
  stat_source <- if (STATS_USE_DISPLAY_RANGE) d_box else d
  stat_dt <- stat_source[, .(
    med = median(get(value_col), na.rm = TRUE),
    mean = mean(get(value_col), na.rm = TRUE)
  ), by = category_lab]
  
  stat_dt[, median_label := sprintf(paste0("median=%.", STAT_LABEL_DIGITS, "g"), med)]
  stat_dt[, mean_label := sprintf(paste0("mean=%.", STAT_LABEL_DIGITS, "g"), mean)]
  
  stat_seg <- copy(stat_dt)
  stat_seg[, x_num := as.numeric(category_lab)]
  stat_seg[, `:=`(
    x_bottom = x_num - STAT_LINE_HALF_HEIGHT,
    x_top = x_num + STAT_LINE_HALF_HEIGHT
  )]
  
  stat_text <- build_stat_text_positions(stat_dt, log10_scale = log10)
  
  peak_source <- if (PEAK_USE_DISPLAY_RANGE) d_violin else d
  peak_dt <- build_peak_table(
    peak_source,
    value_col = value_col,
    log10_scale = log10
  )
  peak_text <- build_peak_text_positions(peak_dt, log10_scale = log10)
  
  box_median_cover_dt <- build_box_median_cover(
    d_box,
    value_col = value_col
  )
  
  pair_stats_source <- if (PAIRWISE_USE_DISPLAY_RANGE) d_box else d
  pair_position_source <- if (PAIRWISE_POSITION_USE_DISPLAY_RANGE) d_box else pair_stats_source
  
  control_source <- if (CONTROL_REF_USE_DISPLAY_RANGE) d_box else d
  control_ref_dt <- build_control_reference(
    control_source,
    value_col = value_col,
    group_col = "category_lab",
    baseline_group = baseline_lab,
    log10_scale = log10
  )
  
  p <- ggplot() +
    geom_violin(
      data = d_violin,
      aes(x = category_lab, y = get(value_col), fill = category, colour = category),
      width = VIOLIN_WIDTH,
      trim = TRUE,
      linewidth = 0.35,
      alpha = 0.60
    ) +
    geom_boxplot(
      data = d_box,
      aes(x = category_lab, y = get(value_col), fill = category, colour = category),
      width = BOXPLOT_WIDTH,
      outlier.shape = NA,
      linewidth = 0.35,
      alpha = BOXPLOT_FILL_ALPHA
    ) +
    scale_fill_manual(values = PALETTE) +
    scale_colour_manual(values = PALETTE) +
    labs(x = NULL, y = ylab, title = title) +
    theme_mbe() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
    coord_flip(clip = "on")
  
  if (should_hide_boxplot_median() && nrow(box_median_cover_dt) > 0) {
    p <- p + geom_segment(
      data = box_median_cover_dt,
      aes(x = x_bottom, xend = x_top, y = box_med, yend = box_med),
      inherit.aes = FALSE,
      linewidth = BOXPLOT_MEDIAN_COVER_LWD,
      lineend = "butt",
      colour = box_median_cover_dt$cover_color
    )
  }
  
  if (nrow(pts) > 0) {
    p <- p + geom_point(
      data = pts,
      aes(x = category_lab, y = v),
      inherit.aes = FALSE,
      position = position_jitter(width = JITTER_W, height = 0),
      alpha = POINT_ALPHA, size = POINT_SIZE, colour = "black"
    )
  }
  
  if (hline0) {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35)
  }
  
  if (log10) {
    p <- p + scale_y_log10(expand = expansion(mult = c(PLOT_EXPAND_LEFT, PLOT_EXPAND_RIGHT)))
  } else {
    p <- p + scale_y_continuous(expand = expansion(mult = c(PLOT_EXPAND_LEFT, PLOT_EXPAND_RIGHT)))
  }
  
  p <- add_control_reference_layers(p, control_ref_dt)
  
  if (SHOW_STAT_LINES && SHOW_MEDIAN_LINE) {
    p <- p + geom_segment(
      data = stat_seg,
      aes(x = x_bottom, xend = x_top, y = med, yend = med),
      inherit.aes = FALSE,
      colour = STAT_LINE_MED_COLOR,
      linewidth = STAT_LINE_MED_LWD,
      linetype = STAT_LINE_MED_LTY,
      lineend = STAT_LINE_END
    )
  }
  
  if (SHOW_STAT_LINES && SHOW_MEAN_LINE) {
    p <- p + geom_segment(
      data = stat_seg,
      aes(x = x_bottom, xend = x_top, y = mean, yend = mean),
      inherit.aes = FALSE,
      colour = STAT_LINE_MEAN_COLOR,
      linewidth = STAT_LINE_MEAN_LWD,
      linetype = STAT_LINE_MEAN_LTY,
      lineend = STAT_LINE_END
    )
  }
  
  if (SHOW_PEAK_LINE && nrow(peak_dt) > 0) {
    p <- p + geom_segment(
      data = peak_dt,
      aes(x = x_bottom, xend = x_top, y = peak, yend = peak),
      inherit.aes = FALSE,
      colour = PEAK_LINE_COLOR,
      linewidth = PEAK_LINE_LWD,
      linetype = PEAK_LINE_LTY,
      lineend = "butt"
    )
  }
  
  if (DRAW_STAT_POINTS && (SHOW_MEDIAN_LINE || SHOW_MEDIAN_TEXT)) {
    p <- p + geom_point(
      data = stat_dt,
      aes(x = category_lab, y = med),
      inherit.aes = FALSE,
      size = 2.0,
      colour = STAT_LINE_MED_COLOR
    )
  }
  
  if (DRAW_STAT_POINTS && (SHOW_MEAN_LINE || SHOW_MEAN_TEXT)) {
    p <- p + geom_point(
      data = stat_dt,
      aes(x = category_lab, y = mean),
      inherit.aes = FALSE,
      size = 2.0,
      colour = STAT_LINE_MEAN_COLOR,
      shape = 17
    )
  }
  
  if (SHOW_MEDIAN_TEXT) {
    p <- p + geom_text(
      data = stat_text,
      aes(x = category_lab, y = med_label_y, label = median_label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 0.5,
      nudge_x = STAT_LABEL_NUDGE_X,
      size = STAT_TEXT_SIZE,
      colour = STAT_TEXT_COLOR_MEDIAN
    )
  }
  
  if (SHOW_MEAN_TEXT) {
    p <- p + geom_text(
      data = stat_text,
      aes(x = category_lab, y = mean_label_y, label = mean_label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 0.5,
      nudge_x = STAT_LABEL_NUDGE_X,
      size = STAT_TEXT_SIZE,
      colour = STAT_TEXT_COLOR_MEAN
    )
  }
  
  if (SHOW_PEAK_TEXT && nrow(peak_text) > 0) {
    p <- p + geom_text(
      data = peak_text,
      aes(x = category_lab, y = peak_label_y, label = peak_label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 0.5,
      nudge_x = PEAK_LABEL_NUDGE_X,
      size = PEAK_TEXT_SIZE,
      colour = PEAK_TEXT_COLOR
    )
  }
  
  pair_dt <- data.table()
  
  if (SHOW_PAIRWISE_SIGNIF) {
    pair_dt <- pairwise_annotations(
      stats_d = pair_stats_source,
      position_d = pair_position_source,
      value_col = value_col,
      group_col = "category_lab",
      group_levels = levels(d$category_lab),
      log10_scale = log10,
      plot_mode = PAIRWISE_PLOT_MODE,
      baseline_group = baseline_lab
    )
    p <- add_pairwise_layers(p, pair_dt)
  }
  
  ggsave(out_png, p, width = width, height = height, dpi = FIG_DPI, bg = "white")
  ggsave(out_pdf, p, width = width, height = height, device = cairo_pdf, bg = "white")
  
  plot_pairwise_out <- extract_plot_pairwise_table(
    pair_dt = pair_dt,
    plot_name = if (is.null(plot_name)) value_col else plot_name
  )
  
  return(invisible(plot_pairwise_out))
}

plot_pi_both_mbe <- function(dt, out_png, out_pdf, width = 5.8, height = 6.0,
                             plot_name = "pi_both") {
  d <- dt[is.finite(pi_tra) & is.finite(pi_tba)]
  d <- d[pi_tra > 0 & pi_tba > 0]
  d <- d[category %chin% CATEGORY_ORDER]
  d[, category := factor(category, levels = CATEGORY_ORDER)]
  d[, category_lab := LABEL_MAP[as.character(category)]]
  d[, category_lab := factor(category_lab, levels = LABEL_MAP[CATEGORY_ORDER])]
  
  baseline_lab <- unname(LABEL_MAP[BASELINE_CATEGORY])
  
  long <- rbind(
    d[, .(category, category_lab, species = "Tra", pi = pi_tra)],
    d[, .(category, category_lab, species = "Tba", pi = pi_tba)]
  )
  
  long_violin <- filter_middle_range_by_group(
    long, value_col = "pi", group_cols = c("species", "category_lab"),
    middle_prop = VIOLIN_MIDDLE_PROP
  )
  
  long_box <- filter_middle_range_by_group(
    long, value_col = "pi", group_cols = c("species", "category_lab"),
    middle_prop = BOXPLOT_MIDDLE_PROP
  )
  
  pts_source <- if (POINTS_FOLLOW_VIOLIN_RANGE) long_violin else long
  pts <- subsample_points(
    pts_source[, .(category, category_lab, species, pi)],
    c("species", "category_lab")
  )
  
  stat_source <- if (STATS_USE_DISPLAY_RANGE) long_box else long
  stat_dt <- stat_source[, .(
    med = median(pi, na.rm = TRUE),
    mean = mean(pi, na.rm = TRUE)
  ), by = .(category_lab, species)]
  
  stat_dt[, median_label := sprintf(paste0("median=%.", STAT_LABEL_DIGITS, "g"), med)]
  stat_dt[, mean_label := sprintf(paste0("mean=%.", STAT_LABEL_DIGITS, "g"), mean)]
  
  stat_seg <- copy(stat_dt)
  stat_seg[, x_num := as.numeric(category_lab)]
  stat_seg[, `:=`(
    x_bottom = x_num - STAT_LINE_HALF_HEIGHT,
    x_top = x_num + STAT_LINE_HALF_HEIGHT
  )]
  
  stat_text <- build_stat_text_positions(stat_dt, log10_scale = TRUE)
  
  peak_source <- if (PEAK_USE_DISPLAY_RANGE) long_violin else long
  peak_dt <- build_peak_table(
    peak_source,
    value_col = "pi",
    log10_scale = TRUE,
    facet_col = "species"
  )
  peak_text <- build_peak_text_positions(peak_dt, log10_scale = TRUE)
  
  box_median_cover_dt <- build_box_median_cover(
    long_box,
    value_col = "pi",
    facet_col = "species"
  )
  
  pair_stats_source <- if (PAIRWISE_USE_DISPLAY_RANGE) long_box else long
  pair_position_source <- if (PAIRWISE_POSITION_USE_DISPLAY_RANGE) long_box else pair_stats_source
  
  control_source <- if (CONTROL_REF_USE_DISPLAY_RANGE) long_box else long
  control_ref_dt <- build_control_reference(
    control_source,
    value_col = "pi",
    group_col = "category_lab",
    baseline_group = baseline_lab,
    log10_scale = TRUE,
    facet_col = "species"
  )
  
  p <- ggplot() +
    geom_violin(
      data = long_violin,
      aes(x = category_lab, y = pi, fill = category, colour = category),
      width = VIOLIN_WIDTH,
      trim = TRUE,
      linewidth = 0.35,
      alpha = 0.60
    ) +
    geom_boxplot(
      data = long_box,
      aes(x = category_lab, y = pi, fill = category, colour = category),
      width = BOXPLOT_WIDTH,
      outlier.shape = NA,
      linewidth = 0.35,
      alpha = BOXPLOT_FILL_ALPHA
    ) +
    facet_wrap(~ species, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = PALETTE) +
    scale_colour_manual(values = PALETTE) +
    labs(x = NULL, y = "Pi (log10 scale)", title = "Nucleotide diversity (Pi) in Tra and Tba") +
    theme_mbe() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
    coord_flip(clip = "on") +
    scale_y_log10(expand = expansion(mult = c(PLOT_EXPAND_LEFT, PLOT_EXPAND_RIGHT)))
  
  if (should_hide_boxplot_median() && nrow(box_median_cover_dt) > 0) {
    p <- p + geom_segment(
      data = box_median_cover_dt,
      aes(x = x_bottom, xend = x_top, y = box_med, yend = box_med),
      inherit.aes = FALSE,
      linewidth = BOXPLOT_MEDIAN_COVER_LWD,
      lineend = "butt",
      colour = box_median_cover_dt$cover_color
    )
  }
  
  if (nrow(pts) > 0) {
    p <- p + geom_point(
      data = pts,
      aes(x = category_lab, y = pi),
      inherit.aes = FALSE,
      position = position_jitter(width = JITTER_W, height = 0),
      alpha = POINT_ALPHA, size = POINT_SIZE, colour = "black"
    )
  }
  
  p <- add_control_reference_layers(p, control_ref_dt)
  
  if (SHOW_STAT_LINES && SHOW_MEDIAN_LINE) {
    p <- p + geom_segment(
      data = stat_seg,
      aes(x = x_bottom, xend = x_top, y = med, yend = med),
      inherit.aes = FALSE,
      colour = STAT_LINE_MED_COLOR,
      linewidth = STAT_LINE_MED_LWD,
      linetype = STAT_LINE_MED_LTY,
      lineend = STAT_LINE_END
    )
  }
  
  if (SHOW_STAT_LINES && SHOW_MEAN_LINE) {
    p <- p + geom_segment(
      data = stat_seg,
      aes(x = x_bottom, xend = x_top, y = mean, yend = mean),
      inherit.aes = FALSE,
      colour = STAT_LINE_MEAN_COLOR,
      linewidth = STAT_LINE_MEAN_LWD,
      linetype = STAT_LINE_MEAN_LTY,
      lineend = STAT_LINE_END
    )
  }
  
  if (SHOW_PEAK_LINE && nrow(peak_dt) > 0) {
    p <- p + geom_segment(
      data = peak_dt,
      aes(x = x_bottom, xend = x_top, y = peak, yend = peak),
      inherit.aes = FALSE,
      colour = PEAK_LINE_COLOR,
      linewidth = PEAK_LINE_LWD,
      linetype = PEAK_LINE_LTY,
      lineend = "butt"
    )
  }
  
  if (DRAW_STAT_POINTS && (SHOW_MEDIAN_LINE || SHOW_MEDIAN_TEXT)) {
    p <- p + geom_point(
      data = stat_dt,
      aes(x = category_lab, y = med),
      inherit.aes = FALSE,
      size = 2.0,
      colour = STAT_LINE_MED_COLOR
    )
  }
  
  if (DRAW_STAT_POINTS && (SHOW_MEAN_LINE || SHOW_MEAN_TEXT)) {
    p <- p + geom_point(
      data = stat_dt,
      aes(x = category_lab, y = mean),
      inherit.aes = FALSE,
      size = 2.0,
      colour = STAT_LINE_MEAN_COLOR,
      shape = 17
    )
  }
  
  if (SHOW_MEDIAN_TEXT) {
    p <- p + geom_text(
      data = stat_text,
      aes(x = category_lab, y = med_label_y, label = median_label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 0.5,
      nudge_x = STAT_LABEL_NUDGE_X,
      size = STAT_TEXT_SIZE,
      colour = STAT_TEXT_COLOR_MEDIAN
    )
  }
  
  if (SHOW_MEAN_TEXT) {
    p <- p + geom_text(
      data = stat_text,
      aes(x = category_lab, y = mean_label_y, label = mean_label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 0.5,
      nudge_x = STAT_LABEL_NUDGE_X,
      size = STAT_TEXT_SIZE,
      colour = STAT_TEXT_COLOR_MEAN
    )
  }
  
  if (SHOW_PEAK_TEXT && nrow(peak_text) > 0) {
    p <- p + geom_text(
      data = peak_text,
      aes(x = category_lab, y = peak_label_y, label = peak_label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 0.5,
      nudge_x = PEAK_LABEL_NUDGE_X,
      size = PEAK_TEXT_SIZE,
      colour = PEAK_TEXT_COLOR
    )
  }
  
  pair_dt <- data.table()
  
  if (SHOW_PAIRWISE_SIGNIF) {
    pair_dt <- pairwise_annotations(
      stats_d = pair_stats_source,
      position_d = pair_position_source,
      value_col = "pi",
      group_col = "category_lab",
      group_levels = levels(d$category_lab),
      log10_scale = TRUE,
      facet_col = "species",
      plot_mode = PAIRWISE_PLOT_MODE,
      baseline_group = baseline_lab
    )
    p <- add_pairwise_layers(p, pair_dt)
  }
  
  ggsave(out_png, p, width = width, height = height, dpi = FIG_DPI, bg = "white")
  ggsave(out_pdf, p, width = width, height = height, device = cairo_pdf, bg = "white")
  
  plot_pairwise_out <- extract_plot_pairwise_table(
    pair_dt = pair_dt,
    plot_name = plot_name,
    facet_col = "species"
  )
  
  return(invisible(plot_pairwise_out))
}

hex_to_itemRgb <- function(hex_vec) {
  vapply(hex_vec, function(h) {
    rgb <- grDevices::col2rgb(h)
    paste(rgb[1, 1], rgb[2, 1], rgb[3, 1], sep = ",")
  }, character(1))
}

merge_touching_bed4 <- function(bed_dt) {
  x <- copy(bed_dt[, .(chr, start, end, name)])
  x <- x[end > start]
  setorder(x, chr, start, end)
  
  if (nrow(x) == 0) {
    return(x[, .(chr, start, end, name)])
  }
  
  x[, prev_end := shift(end), by = chr]
  x[, new_block := is.na(prev_end) | (start > prev_end)]
  x[, block_id := cumsum(new_block), by = chr]
  
  merged <- x[, .(
    start = min(start),
    end = max(end),
    name = name[1]
  ), by = .(chr, block_id)]
  
  merged[, block_id := NULL]
  setcolorder(merged, c("chr", "start", "end", "name"))
  merged[]
}

write_palette_bed <- function(dt, prefix = "palette_regions") {
  bed_dt <- copy(dt[, .(chr, start, end, category)])
  bed_dt <- bed_dt[category %chin% names(PALETTE)]
  bed_dt[, name := as.character(category)]
  bed_dt <- bed_dt[, .(chr, start, end, name)]
  
  for (cat_name in names(PALETTE)) {
    sub <- bed_dt[name == cat_name]
    sub <- sub[order(chr, start, end)]
    
    fwrite(
      sub[, .(chr, start, end, name)],
      file = paste0(prefix, ".", cat_name, ".bed"),
      sep = "\t",
      col.names = FALSE
    )
    
    sub_merged <- merge_touching_bed4(sub)
    
    fwrite(
      sub_merged[, .(chr, start, end, name)],
      file = paste0(prefix, ".", cat_name, ".merged.bed"),
      sep = "\t",
      col.names = FALSE
    )
  }
}

# =========================
# 4) Read and merge windows
# =========================
message("Reading window files...")
pi_tra <- read_window_bed(pi_tra_file, "pi_tra")
pi_tba <- read_window_bed(pi_tba_file, "pi_tba")
dxy <- read_window_bed(dxy_file, "dxy")
fst <- read_window_bed(fst_file, "fst")

message("Merging by chr/start/end...")
setkey(pi_tra, chr, start, end)
setkey(pi_tba, chr, start, end)
setkey(dxy, chr, start, end)
setkey(fst, chr, start, end)

dt <- merge(pi_tra, pi_tba, by = c("chr", "start", "end"), all = FALSE)
dt <- merge(dt, dxy, by = c("chr", "start", "end"), all = FALSE)
dt <- merge(dt, fst, by = c("chr", "start", "end"), all = FALSE)

dt <- dt[is.finite(pi_tra) & is.finite(pi_tba) & is.finite(dxy) & is.finite(fst)]
dt <- dt[dxy > 0]
dt[, mid0 := as.integer(floor((start + end) / 2))]
mid_gr <- GRanges(seqnames = dt$chr, ranges = IRanges(start = dt$mid0 + 1L, width = 1L))
dt[, mean_pi := (pi_tra + pi_tba) / 2]

# =========================
# 5) Strict region labeling
# =========================
message("Reading region beds...")
gr_invert <- read_region_bed(bed_invert)
gr_non_invert <- read_region_bed(bed_non_invert)
gr_ismc <- read_region_bed(bed_ismc)
gr_non_ismc <- read_region_bed(bed_non_ismc)

message("Assigning labels with overlap exclusion...")
in_invert <- label_by_midpoint(mid_gr, gr_invert)
in_non_invert <- label_by_midpoint(mid_gr, gr_non_invert)
in_ismc <- label_by_midpoint(mid_gr, gr_ismc)
in_non_ismc <- label_by_midpoint(mid_gr, gr_non_ismc)

overlap_inv <- in_invert & in_non_invert
overlap_ismc <- in_ismc & in_non_ismc

inv_unique_ok <- xor(in_invert, in_non_invert)
recomb_unique_ok <- xor(in_ismc, in_non_ismc)

unlabeled_inv <- (!in_invert) & (!in_non_invert)
unlabeled_rec <- (!in_ismc) & (!in_non_ismc)

keep <- inv_unique_ok & recomb_unique_ok

dt[, drop_reason := NA_character_]
dt[overlap_inv == TRUE, drop_reason := "overlap_invert_vs_noninvert"]
dt[overlap_ismc == TRUE, drop_reason := fifelse(
  is.na(drop_reason),
  "overlap_ismc_vs_nonismc",
  paste(drop_reason, "overlap_ismc_vs_nonismc", sep = ";")
)]
dt[is.na(drop_reason) & unlabeled_inv == TRUE & unlabeled_rec == TRUE, drop_reason := "unlabeled_both_axes"]
dt[is.na(drop_reason) & unlabeled_inv == TRUE, drop_reason := "unlabeled_invert_axis"]
dt[is.na(drop_reason) & unlabeled_rec == TRUE, drop_reason := "unlabeled_ismc_axis"]

dropped <- dt[keep == FALSE]
dropped[, `:=`(
  in_invert = in_invert[keep == FALSE],
  in_non_invert = in_non_invert[keep == FALSE],
  in_ismc = in_ismc[keep == FALSE],
  in_non_ismc = in_non_ismc[keep == FALSE]
)]
fwrite(dropped, out_dropped, sep = "\t")
message("Wrote dropped windows: ", out_dropped)

dt_strict <- dt[keep]
dt_strict[, inv_status := ifelse(in_invert[keep], "invert", "non_invert")]
dt_strict[, recomb_status := ifelse(in_ismc[keep], "ismc", "non_ismc")]
dt_strict[, category := paste(inv_status, recomb_status, sep = "_")]

message("Windows total: ", nrow(dt))
message("Windows kept (strict): ", nrow(dt_strict))

write_palette_bed(dt_strict, prefix = "palette_regions")
message("Wrote original and merged BED4 files for PALETTE regions.")

# =========================
# 6) Save tables
# =========================
setcolorder(dt_strict, c(
  "chr", "start", "end", "mid0", "category", "inv_status", "recomb_status",
  "pi_tra", "pi_tba", "mean_pi", "dxy", "fst"
))
fwrite(dt_strict, out_windows, sep = "\t")
message("Wrote per-window stats: ", out_windows)

summary_dt <- dt_strict[, .(
  n = .N,
  median_pi_tra = median(pi_tra), mean_pi_tra = mean(pi_tra),
  median_pi_tba = median(pi_tba), mean_pi_tba = mean(pi_tba),
  median_dxy = median(dxy), mean_dxy = mean(dxy),
  median_fst = median(fst), mean_fst = mean(fst)
), by = category][order(-n)]
fwrite(summary_dt, out_summary, sep = "\t")
message("Wrote summary: ", out_summary)

# =========================
# 7) Tests
# =========================
baseline <- BASELINE_CATEGORY
metrics <- list(pi_tra = "pi_tra", pi_tba = "pi_tba", dxy = "dxy", fst = "fst")
pairwise_all_results <- list()

sink(out_tests)
cat("=== Strict overlap-excluded analysis; baseline =", baseline, "===\n\n")
cat("Kept windows:", nrow(dt_strict), "\n\n")

for (mname in names(metrics)) {
  col <- metrics[[mname]]
  log_metric <- is_log_metric(mname)
  
  cat("\n\n==============================\n")
  cat("Metric:", mname, " (column:", col, ")\n")
  cat("==============================\n")
  
  sub <- dt_strict[is.finite(get(col))]
  if (log_metric) {
    sub <- sub[get(col) > 0]
  }
  
  cat("\n--- Kruskal-Wallis ---\n")
  print(kruskal.test(sub[[col]] ~ sub$category))
  
  cat("\n--- Pairwise comparisons (all pairs; adjusted across all pairs) ---\n")
  pairwise_all <- compute_pairwise_table(
    sub,
    value_col = col,
    group_col = "category",
    group_levels = CATEGORY_ORDER,
    log10_scale = FALSE,
    metric_name = mname
  )
  
  if (nrow(pairwise_all) > 0) {
    pairwise_all[, `:=`(
      group1_label = LABEL_MAP[group1],
      group2_label = LABEL_MAP[group2],
      baseline_vs_other = (group1 == baseline | group2 == baseline)
    )]
    print(pairwise_all[, .(
      metric,
      group1, group1_label,
      group2, group2_label,
      n1, n2,
      median1, median2, median_diff,
      mean1, mean2, mean_diff,
      p_raw, p_adj, stars,
      label,
      baseline_vs_other
    )])
  } else {
    cat("No pairwise comparisons available.\n")
  }
  
  pairwise_all_results[[mname]] <- copy(pairwise_all)
  
  cat("\n--- Effect sizes vs baseline ---\n")
  cats <- sort(unique(sub$category))
  cats <- cats[cats != baseline]
  
  for (cname in cats) {
    x <- sub[category == cname, get(col)]
    y <- sub[category == baseline, get(col)]
    mdiff <- median(x, na.rm = TRUE) - median(y, na.rm = TRUE)
    cd <- cliffs_delta(x, y)
    bj <- block_jackknife_median_diff(
      sub[category %in% c(baseline, cname)],
      groupA = cname,
      groupB = baseline,
      value_col = col,
      block_col = "chr"
    )
    cat("\nCategory:", cname, "vs", baseline, "\n")
    cat("  median(A)-median(B) =", mdiff, "\n")
    cat("  Cliff's delta =", cd, "\n")
    cat("  Block-jackknife median diff est =", bj$est, " SE =", bj$se, " n_blocks =", bj$n_blocks, "\n")
  }
}
sink()
message("Wrote tests: ", out_tests)

pairwise_all_dt <- rbindlist(pairwise_all_results, fill = TRUE)
if (nrow(pairwise_all_dt) > 0) {
  pairwise_all_dt[, `:=`(
    group1_label = LABEL_MAP[group1],
    group2_label = LABEL_MAP[group2],
    baseline_vs_other = (group1 == BASELINE_CATEGORY | group2 == BASELINE_CATEGORY)
  )]
  setcolorder(pairwise_all_dt, c(
    "metric",
    "group1", "group1_label",
    "group2", "group2_label",
    "n1", "n2",
    "median1", "median2", "median_diff",
    "mean1", "mean2", "mean_diff",
    "p_raw", "p_adj", "stars",
    "baseline_vs_other",
    "label"
  ))
  fwrite(pairwise_all_dt, out_pairwise_table, sep = "\t")
  message("Wrote all-pairwise table: ", out_pairwise_table)
}

# =========================
# 8) Plots
# =========================
plot_pairwise_results <- list()

plot_pairwise_results[["pi_tra"]] <- plot_metric_mbe(
  dt_strict, "pi_tra",
  title = "Nucleotide diversity in Tra",
  ylab = "Pi (Tra)",
  log10 = TRUE, hline0 = FALSE,
  out_png = out_plot_pi_tra_png, out_pdf = out_plot_pi_tra_pdf,
  plot_name = "pi_tra"
)

plot_pairwise_results[["pi_tba"]] <- plot_metric_mbe(
  dt_strict, "pi_tba",
  title = "Nucleotide diversity in Tba",
  ylab = "Pi (Tba)",
  log10 = TRUE, hline0 = FALSE,
  out_png = out_plot_pi_tba_png, out_pdf = out_plot_pi_tba_pdf,
  plot_name = "pi_tba"
)

plot_pairwise_results[["pi_both"]] <- plot_pi_both_mbe(
  dt_strict,
  out_png = out_plot_pi_both_png, out_pdf = out_plot_pi_both_pdf,
  plot_name = "pi_both"
)

plot_pairwise_results[["dxy"]] <- plot_metric_mbe(
  dt_strict, "dxy",
  title = "Interspecific divergence (Dxy)",
  ylab = "Dxy",
  log10 = TRUE, hline0 = FALSE,
  out_png = out_plot_dxy_png, out_pdf = out_plot_dxy_pdf,
  plot_name = "dxy"
)

plot_pairwise_results[["fst"]] <- plot_metric_mbe(
  dt_strict, "fst",
  title = "Genetic differentiation (FST)",
  ylab = "FST",
  log10 = FALSE, hline0 = FALSE,
  out_png = out_plot_fst_png, out_pdf = out_plot_fst_pdf,
  plot_name = "fst"
)

pairwise_plot_dt <- rbindlist(plot_pairwise_results, fill = TRUE)

if (nrow(pairwise_plot_dt) > 0) {
  setcolorder(pairwise_plot_dt, c(
    "plot_name",
    "facet",
    "group1", "group1_label",
    "group2", "group2_label",
    "n1", "n2",
    "median1", "median2", "median_diff",
    "mean1", "mean2", "mean_diff",
    "p_raw", "p_adj", "stars",
    "label"
  ))
  fwrite(pairwise_plot_dt, out_plot_pairwise_table, sep = "\t")
  message("Wrote plot-matched pairwise table: ", out_plot_pairwise_table)
}

message("Wrote figures (PNG + PDF). Done.")
