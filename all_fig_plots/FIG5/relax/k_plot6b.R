## ========== k_plot.R (ismc-group vs noismc-group；支持 P/FDR/NONE；最终比较图改为论文风格 median+IQR 浮动柱图) ==========

.libPaths("/data/projects/rcui/R/x86_64-pc-linux-gnu-library/4.1")
#setwd("/data3/projects/yshao/relax/7_genespace_GO/")

suppressPackageStartupMessages({
  library(ggplot2)
})

## ---------- 用户可调参数 ----------
nName <- "Relax_Tra"
sRefSp <- "Tridentigerradiatus"

# 原始分析文件
orig_relax_file <- "/data3/projects/yshao/relax/6_relax_Tra/sum_Tra.txt"

# 新增分析文件
new_relax_file  <- "/data3/projects/yshao/relax/6_relax_Tba/sum_Tba.txt"
new_relax_tag   <- "Tba"   # 仅用于输出文件名前缀，可自行修改

# 注释文件
pgid2ref_file <- "/data3/projects/yshao/relax/7_genespace_GO/synorthos.txt"
groupid2transid_file <- "/data3/projects/yshao/relax/4_exprotAln_genespace2/groupid2transid_tra.txt"

# 两个 mrna 列表文件
ismc_mrna_file   <- "/data3/projects/yshao/relax/7_genespace_GO/k_plot2/pure_fsr/mRNA_ids.txt"
noismc_mrna_file <- "/data3/projects/yshao/relax/7_genespace_GO/k_plot2/pure_fsr/control_mRNA_ids.txt"

# ----------------- 筛选参数 -----------------
# 可选: "P" / "FDR" / "NONE"
#   P    -> 按 V5(P值) 筛选
#   FDR  -> 按 BH-FDR 筛选
#   NONE -> 不筛选，直接用全部数据做 V9 分析
sig_filter_mode <- "NONE"

# 当 sig_filter_mode == "P" 时使用
sig_p_cutoff <- 0.05

# 当 sig_filter_mode == "FDR" 时使用
sig_fdr_cutoff <- 0.1

# FDR 的列名（由 V5 现算得到）
fdr_colname <- "FDR_BH"

# 输出前缀
orig_prefix_base <- "V9_ismcgroup_vs_noismcgroup"
new_prefix_base  <- paste0("V9_ismcgroup_vs_noismcgroup_", new_relax_tag, "_byOriginalGroup")
group_record_base <- paste0("ismc_noismc_group_record_", new_relax_tag)

# 差异分布分析输出前缀
diff_prefix_base <- paste0("V9_pairwiseDiff_noismc_minus_ismc_original_vs_", new_relax_tag)

plot_title_orig_base <- "V9: ismc-group vs noismc-group (Median + IQR, Wilcoxon)"
plot_title_new_base  <- paste0("V9: ismc-group vs noismc-group in ", new_relax_tag,
                               " by original groups (Median + IQR, Wilcoxon)")
plot_title_diff_base <- paste0("Pairwise V9 difference summaries (noismc - ismc): original vs ", new_relax_tag)

x_label      <- ""
y_label      <- "K value"
caption_lab  <- "Floating bars show IQR; center line = median; Wilcoxon rank-sum test (two-sided)."

group_order  <- c("ismc", "noismc")
group_labels <- c("ismc" = "ismc groups", "noismc" = "noismc groups")
group_colors <- c("ismc groups" = "#377eb8", "noismc groups" = "#e41a1c")

# 最终比较图标签与配色
diff_group_order <- c("original", "new")

new_analysis_label <- paste0(new_relax_tag, " analysis")

diff_group_labels <- c(
  "original" = "Original analysis",
  "new" = new_analysis_label
)

diff_group_colors <- setNames(
  c("#4daf4a", "#984ea3"),
  c("Original analysis", new_analysis_label)
)

## ---------- 通用函数 ----------
to_num <- function(x) suppressWarnings(as.numeric(x))

normalize_filter_mode <- function(x) {
  x <- toupper(trimws(x))
  if (x %in% c("NOFILTER", "ALL")) x <- "NONE"
  if (!(x %in% c("P", "FDR", "NONE"))) {
    stop("sig_filter_mode must be one of 'P', 'FDR', or 'NONE'.")
  }
  x
}

make_num_tag <- function(x) {
  x <- format(x, scientific = FALSE, trim = TRUE)
  x <- gsub("\\.", "p", x)
  x <- gsub("-", "m", x)
  x
}

star_sig <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) "****"
  else if (p < 1e-3) "***"
  else if (p < 1e-2) "**"
  else if (p < 5e-2) "*"
  else "ns"
}

calc_summary <- function(x) {
  n  <- length(x)
  m  <- mean(x)
  sdv <- if (n > 1) sd(x) else 0
  se <- if (n > 0) sdv / sqrt(n) else NA_real_
  q  <- stats::quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  data.frame(
    n = n,
    mean = m,
    sd = sdv,
    se = se,
    median = q[2],
    q1 = q[1],
    q3 = q[3],
    iqr = q[3] - q[1],
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

add_p_fdr_columns <- function(dat, p_col = "V5", fdr_col = "FDR_BH") {
  if (!(p_col %in% colnames(dat))) {
    stop(sprintf("Column %s not found.", p_col))
  }
  p_num <- to_num(dat[[p_col]])
  fdr_vec <- rep(NA_real_, length(p_num))
  idx <- which(!is.na(p_num))
  if (length(idx) > 0) {
    fdr_vec[idx] <- p.adjust(p_num[idx], method = "BH")
  }
  dat$P_numeric <- p_num
  dat[[fdr_col]] <- fdr_vec
  dat
}

flag_sig_rows <- function(dat,
                          filter_mode = "P",
                          p_cutoff = 0.05,
                          fdr_cutoff = 0.1,
                          p_col = "P_numeric",
                          fdr_col = "FDR_BH") {
  filter_mode <- normalize_filter_mode(filter_mode)

  if (!(p_col %in% colnames(dat))) stop(sprintf("Column %s not found.", p_col))
  if (!(fdr_col %in% colnames(dat))) stop(sprintf("Column %s not found.", fdr_col))

  if (filter_mode == "NONE") {
    pass_flag <- rep(TRUE, nrow(dat))
  } else if (filter_mode == "P") {
    pass_flag <- !is.na(dat[[p_col]]) & dat[[p_col]] < p_cutoff
  } else {
    pass_flag <- !is.na(dat[[fdr_col]]) & dat[[fdr_col]] < fdr_cutoff
  }

  dat$pass_sig_filter <- pass_flag
  dat
}

get_filter_desc <- function(filter_mode = "P", p_cutoff = 0.05, fdr_cutoff = 0.1) {
  filter_mode <- normalize_filter_mode(filter_mode)
  if (filter_mode == "NONE") {
    "No significance filter (use all rows)"
  } else if (filter_mode == "P") {
    paste0("P < ", format(p_cutoff, scientific = FALSE, trim = TRUE))
  } else {
    paste0("BH-FDR < ", format(fdr_cutoff, scientific = FALSE, trim = TRUE))
  }
}

get_filter_suffix <- function(filter_mode = "P", p_cutoff = 0.05, fdr_cutoff = 0.1) {
  filter_mode <- normalize_filter_mode(filter_mode)
  if (filter_mode == "NONE") {
    "NoFilter"
  } else if (filter_mode == "P") {
    paste0("Plt", make_num_tag(p_cutoff))
  } else {
    paste0("FDRlt", make_num_tag(fdr_cutoff))
  }
}

make_group_count_table <- function(dat, count_list_named) {
  if (length(count_list_named) == 0) {
    stop("count_list_named is empty.")
  }

  if (nrow(dat) == 0) {
    out <- data.frame(V1 = character(0), stringsAsFactors = FALSE)
    for (nm in names(count_list_named)) {
      out[[nm]] <- integer(0)
    }
    return(out)
  }

  tmp <- data.frame(V1 = as.character(dat$V1), stringsAsFactors = FALSE)
  for (nm in names(count_list_named)) {
    tmp[[nm]] <- count_list_named[[nm]]
  }

  aggregate(. ~ V1, data = tmp, FUN = sum)
}

fill_na_zero <- function(df, cols) {
  for (cc in cols) {
    if (cc %in% colnames(df)) {
      df[[cc]][is.na(df[[cc]])] <- 0L
    }
  }
  df
}

assign_group_class <- function(n_ismc, n_noismc) {
  out <- rep(NA_character_, length(n_ismc))
  out[n_ismc > 0 & n_noismc == 0] <- "ismc"
  out[n_ismc == 0 & n_noismc > 0] <- "noismc"
  out[n_ismc > 0 & n_noismc > 0]  <- "both"
  out[n_ismc == 0 & n_noismc == 0] <- "neither"
  out
}

assign_rule_text <- function(group_class) {
  out <- rep(NA_character_, length(group_class))
  out[group_class == "ismc"] <- "assigned_to_ismc_if_group_has_ismc_mrna_and_no_noismc_mrna"
  out[group_class == "noismc"] <- "assigned_to_noismc_if_group_has_noismc_mrna_and_no_ismc_mrna"
  out[group_class == "both"] <- "excluded_if_group_has_both_ismc_mrna_and_noismc_mrna"
  out[group_class == "neither"] <- "excluded_if_group_has_neither_ismc_mrna_nor_noismc_mrna"
  out
}

## ---------- 新增：导出与 relax_file 格式一致的分析行 ----------
extract_relax_like_rows <- function(dat, relax_colnames) {
  miss_cols <- setdiff(relax_colnames, colnames(dat))
  if (length(miss_cols) > 0) {
    stop(sprintf(
      "The following relax-file columns are missing in the data to export: %s",
      paste(miss_cols, collapse = ", ")
    ))
  }
  dat[, relax_colnames, drop = FALSE]
}

write_relax_like_rows <- function(dat, relax_colnames, outfile) {
  out_dat <- extract_relax_like_rows(dat, relax_colnames)
  write.table(out_dat,
              file = outfile,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
  cat(sprintf("Analysis rows saved in original relax-file format: %s (n = %d)\n",
              outfile, nrow(out_dat)))
}

run_v9_analysis <- function(dat_in,
                            dat_out,
                            plot_title,
                            out_prefix,
                            filter_mode,
                            p_cutoff,
                            fdr_cutoff,
                            x_label = "",
                            y_label = "K value",
                            caption_lab = "",
                            group_order = c("ismc", "noismc"),
                            group_labels = c("ismc" = "ismc groups", "noismc" = "noismc groups"),
                            group_colors = c("ismc groups" = "#377eb8", "noismc groups" = "#e41a1c")) {

  if (!("V9" %in% colnames(dat_in)) || !("V9" %in% colnames(dat_out))) {
    stop("V9 column not found in one or both input data frames.")
  }

  filter_desc <- get_filter_desc(filter_mode, p_cutoff, fdr_cutoff)
  full_caption <- paste0(caption_lab, " Significance filter: ", filter_desc, ".")

  x_in  <- to_num(dat_in$V9)
  x_out <- to_num(dat_out$V9)
  x_in  <- x_in[!is.na(x_in)]
  x_out <- x_out[!is.na(x_out)]

  if (length(x_in) == 0 || length(x_out) == 0) {
    stop(sprintf("Not enough numeric V9 data after filtering for one or both groups in analysis: %s", out_prefix))
  }

  wtest <- wilcox.test(x_in, x_out, alternative = "two.sided", exact = FALSE)
  pval  <- wtest$p.value
  p_star <- star_sig(pval)

  cat(sprintf("[%s] Wilcoxon: W = %.0f, p = %.3g, stars = %s\n",
              out_prefix, as.numeric(wtest$statistic), pval, p_star))

  sum_in  <- calc_summary(x_in);  sum_in$group  <- "ismc"
  sum_out <- calc_summary(x_out); sum_out$group <- "noismc"

  summary_df <- rbind(sum_in, sum_out)
  summary_df <- summary_df[, c("group","n","mean","sd","se","median","q1","q3","iqr","min","max")]

  plot_df <- summary_df
  plot_df$group <- factor(plot_df$group, levels = group_order, labels = group_labels[group_order])

  y_top   <- max(plot_df$q3, na.rm = TRUE)
  y_base  <- min(plot_df$q1, na.rm = TRUE)
  y_span  <- max(1e-8, y_top - y_base)
  y_brkt  <- y_top + 0.08 * y_span
  y_star  <- y_brkt + 0.04 * y_span
  x_left  <- 1
  x_right <- 2

  p <- ggplot(plot_df, aes(x = group, fill = group)) +
    geom_crossbar(aes(y = median, ymin = q1, ymax = q3),
                  width = 0.55, fatten = 2.0,
                  color = "#2b2b2b", alpha = 0.95) +
    scale_fill_manual(values = group_colors) +
    labs(title = plot_title, x = x_label, y = y_label, caption = full_caption) +
    theme_classic(base_size = 15) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0, face = "bold", size = 16),
      axis.title.x = element_text(size = 13),
      axis.title.y = element_text(size = 13),
      axis.text = element_text(color = "black", size = 12),
      plot.margin = margin(10, 16, 10, 10)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.24)))

  p <- p +
    geom_segment(aes(x = x_left, xend = x_right, y = y_brkt, yend = y_brkt),
                 inherit.aes = FALSE, linewidth = 0.6, color = "#2b2b2b") +
    geom_segment(aes(x = x_left, xend = x_left, y = y_brkt, yend = y_brkt - 0.02 * y_span),
                 inherit.aes = FALSE, linewidth = 0.6, color = "#2b2b2b") +
    geom_segment(aes(x = x_right, xend = x_right, y = y_brkt, yend = y_brkt - 0.02 * y_span),
                 inherit.aes = FALSE, linewidth = 0.6, color = "#2b2b2b") +
    annotate("text", x = 1.5, y = y_star, label = p_star, size = 5.2)

  plot_png    <- paste0(out_prefix, "_median_iqr.png")
  plot_pdf    <- paste0(out_prefix, "_median_iqr.pdf")
  summary_tsv <- paste0(out_prefix, "_summary.tsv")

  ggsave(plot_png, p, width = 5.6, height = 4.4, dpi = 600)
  ggsave(plot_pdf, p, width = 5.6, height = 4.4, dpi = 600)
  cat(sprintf("[%s] Figure saved: %s, %s\n", out_prefix, plot_png, plot_pdf))

  out_tab <- data.frame(
    group = as.character(plot_df$group),
    n = summary_df$n,
    mean = round(summary_df$mean, 6),
    sd = round(summary_df$sd, 6),
    se = round(summary_df$se, 6),
    median = round(summary_df$median, 6),
    q1 = round(summary_df$q1, 6),
    q3 = round(summary_df$q3, 6),
    iqr = round(summary_df$iqr, 6),
    min = round(summary_df$min, 6),
    max = round(summary_df$max, 6),
    wilcoxon_W = as.numeric(wtest$statistic),
    wilcoxon_p = signif(pval, 6),
    signif_stars = p_star,
    alternative = as.character(wtest$alternative),
    method = as.character(wtest$method),
    filter_mode = normalize_filter_mode(filter_mode),
    filter_desc = filter_desc,
    stringsAsFactors = FALSE
  )

  write.table(out_tab, file = summary_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("[%s] Summary table saved: %s\n", out_prefix, summary_tsv))

  invisible(list(
    summary_df = summary_df,
    out_tab = out_tab,
    wtest = wtest,
    pval = pval,
    p_star = p_star,
    plot_png = plot_png,
    plot_pdf = plot_pdf,
    summary_tsv = summary_tsv
  ))
}

## ---------- 差异分布：全量 pairwise differences ----------
## 差异定义：noismc V9 - ismc V9
make_pairwise_diff <- function(x_noismc, x_ismc) {
  x_noismc <- as.numeric(x_noismc)
  x_ismc <- as.numeric(x_ismc)

  x_noismc <- x_noismc[!is.na(x_noismc)]
  x_ismc <- x_ismc[!is.na(x_ismc)]

  if (length(x_noismc) == 0 || length(x_ismc) == 0) {
    stop("Both x_noismc and x_ismc must contain at least one numeric value.")
  }

  as.numeric(rep(x_noismc, each = length(x_ismc)) - rep(x_ismc, times = length(x_noismc)))
}

## ---------- 最终比较图：改为 summary floating bars (median + IQR) ----------
## 检验仍然使用全量 pairwise differences
run_diff_distribution_analysis <- function(orig_dat_in,
                                           orig_dat_out,
                                           new_dat_in,
                                           new_dat_out,
                                           out_prefix,
                                           plot_title,
                                           filter_mode,
                                           p_cutoff,
                                           fdr_cutoff,
                                           new_tag = "Tba",
                                           diff_group_order = c("original", "new"),
                                           diff_group_labels = c("original" = "Original analysis",
                                                                 "new" = "Tba analysis"),
                                           diff_group_colors = c("Original analysis" = "#4daf4a",
                                                                 "Tba analysis" = "#984ea3")) {

  x_orig_ismc <- to_num(orig_dat_in$V9)
  x_orig_noismc <- to_num(orig_dat_out$V9)
  x_new_ismc <- to_num(new_dat_in$V9)
  x_new_noismc <- to_num(new_dat_out$V9)

  x_orig_ismc <- x_orig_ismc[!is.na(x_orig_ismc)]
  x_orig_noismc <- x_orig_noismc[!is.na(x_orig_noismc)]
  x_new_ismc <- x_new_ismc[!is.na(x_new_ismc)]
  x_new_noismc <- x_new_noismc[!is.na(x_new_noismc)]

  if (length(x_orig_ismc) == 0 || length(x_orig_noismc) == 0) {
    stop("Original analysis does not have enough numeric V9 values to build pairwise difference distribution.")
  }
  if (length(x_new_ismc) == 0 || length(x_new_noismc) == 0) {
    stop("New analysis does not have enough numeric V9 values to build pairwise difference distribution.")
  }

  ## 全量 pairwise differences：用于统计检验与 summary
  diff_orig <- make_pairwise_diff(x_orig_noismc, x_orig_ismc)
  diff_new  <- make_pairwise_diff(x_new_noismc, x_new_ismc)

  wtest <- wilcox.test(diff_orig, diff_new, alternative = "two.sided", exact = FALSE)
  pval  <- wtest$p.value
  p_star <- star_sig(pval)

  cat(sprintf("[%s] Difference-distribution Wilcoxon: W = %.0f, p = %.3g, stars = %s\n",
              out_prefix, as.numeric(wtest$statistic), pval, p_star))

  sum_orig <- calc_summary(diff_orig)
  sum_new  <- calc_summary(diff_new)

  sum_orig$analysis <- "original"
  sum_new$analysis  <- "new"

  summary_df <- rbind(sum_orig, sum_new)
  summary_df$n_ismc_values <- c(length(x_orig_ismc), length(x_new_ismc))
  summary_df$n_noismc_values <- c(length(x_orig_noismc), length(x_new_noismc))
  summary_df$n_pairwise_diffs <- c(length(diff_orig), length(diff_new))
  summary_df$difference_definition <- "all pairwise differences: noismc V9 - ismc V9"
  summary_df$filter_mode <- normalize_filter_mode(filter_mode)
  summary_df$filter_desc <- get_filter_desc(filter_mode, p_cutoff, fdr_cutoff)

  plot_df <- summary_df
  plot_df$analysis <- factor(plot_df$analysis,
                             levels = diff_group_order,
                             labels = diff_group_labels[diff_group_order])

  y_top   <- max(plot_df$q3, na.rm = TRUE)
  y_base  <- min(plot_df$q1, na.rm = TRUE)
  y_span  <- max(1e-8, y_top - y_base)

  y_brkt  <- y_top + 0.09 * y_span
  y_star  <- y_brkt + 0.05 * y_span
  label_y <- plot_df$q3 + 0.025 * y_span

  filter_desc <- get_filter_desc(filter_mode, p_cutoff, fdr_cutoff)
  full_caption <- paste0(
    "Floating bars show IQR; center line = median. ",
    "The Wilcoxon test uses the full pairwise difference distributions (noismc V9 - ismc V9) from the two analyses. ",
    "Significance filter: ", filter_desc, "."
  )

  plot_df$median_label <- sprintf("Median = %.4f", plot_df$median)

  p <- ggplot(plot_df, aes(x = analysis, fill = analysis)) +
    geom_crossbar(aes(y = median, ymin = q1, ymax = q3),
                  width = 0.58, fatten = 2.2,
                  color = "#2b2b2b", alpha = 0.96) +
    geom_text(aes(y = label_y, label = median_label),
              vjust = 0, size = 3.6, fontface = "bold", color = "#2b2b2b") +
    scale_fill_manual(values = diff_group_colors) +
    labs(
      title = plot_title,
      x = "",
      y = "Pairwise difference in V9 (noismc - ismc)",
      caption = full_caption
    ) +
    theme_classic(base_size = 15) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0, face = "bold", size = 16),
      axis.title.x = element_text(size = 13),
      axis.title.y = element_text(size = 13),
      axis.text = element_text(color = "black", size = 12),
      plot.margin = margin(12, 18, 10, 10)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.28))) +
    geom_segment(aes(x = 1, xend = 2, y = y_brkt, yend = y_brkt),
                 inherit.aes = FALSE, linewidth = 0.65, color = "#2b2b2b") +
    geom_segment(aes(x = 1, xend = 1, y = y_brkt, yend = y_brkt - 0.022 * y_span),
                 inherit.aes = FALSE, linewidth = 0.65, color = "#2b2b2b") +
    geom_segment(aes(x = 2, xend = 2, y = y_brkt, yend = y_brkt - 0.022 * y_span),
                 inherit.aes = FALSE, linewidth = 0.65, color = "#2b2b2b") +
    annotate("text", x = 1.5, y = y_star, label = p_star, size = 5.4)

  plot_png    <- paste0(out_prefix, "_median_iqr.png")
  plot_pdf    <- paste0(out_prefix, "_median_iqr.pdf")
  summary_tsv <- paste0(out_prefix, "_summary.tsv")

  ggsave(plot_png, p, width = 6.4, height = 4.9, dpi = 600)
  ggsave(plot_pdf, p, width = 6.4, height = 4.9, dpi = 600)
  cat(sprintf("[%s] Difference-summary figure saved: %s, %s\n", out_prefix, plot_png, plot_pdf))

  out_tab <- data.frame(
    analysis = diff_group_labels[summary_df$analysis],
    n_ismc_values = summary_df$n_ismc_values,
    n_noismc_values = summary_df$n_noismc_values,
    n_pairwise_diffs = summary_df$n_pairwise_diffs,
    mean = round(summary_df$mean, 6),
    sd = round(summary_df$sd, 6),
    se = round(summary_df$se, 6),
    median = round(summary_df$median, 6),
    q1 = round(summary_df$q1, 6),
    q3 = round(summary_df$q3, 6),
    iqr = round(summary_df$iqr, 6),
    min = round(summary_df$min, 6),
    max = round(summary_df$max, 6),
    wilcoxon_W = as.numeric(wtest$statistic),
    wilcoxon_p = signif(pval, 6),
    signif_stars = p_star,
    alternative = as.character(wtest$alternative),
    method = as.character(wtest$method),
    difference_definition = summary_df$difference_definition,
    filter_mode = summary_df$filter_mode,
    filter_desc = summary_df$filter_desc,
    stringsAsFactors = FALSE
  )

  write.table(out_tab, file = summary_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("[%s] Difference-summary table saved: %s\n", out_prefix, summary_tsv))

  invisible(list(
    diff_orig = diff_orig,
    diff_new = diff_new,
    out_tab = out_tab,
    wtest = wtest,
    pval = pval,
    p_star = p_star,
    plot_png = plot_png,
    plot_pdf = plot_pdf,
    summary_tsv = summary_tsv
  ))
}

## ---------- 标准化筛选参数 ----------
sig_filter_mode <- normalize_filter_mode(sig_filter_mode)
filter_desc <- get_filter_desc(sig_filter_mode, sig_p_cutoff, sig_fdr_cutoff)
filter_suffix <- get_filter_suffix(sig_filter_mode, sig_p_cutoff, sig_fdr_cutoff)

orig_prefix <- paste0(orig_prefix_base, "_", filter_suffix)
new_prefix  <- paste0(new_prefix_base, "_", filter_suffix)
group_record_tsv <- paste0(group_record_base, "_", filter_suffix, ".tsv")
diff_prefix <- paste0(diff_prefix_base, "_", filter_suffix)

plot_title_orig <- paste0(plot_title_orig_base, " [", filter_desc, "]")
plot_title_new  <- paste0(plot_title_new_base, " [", filter_desc, "]")
plot_title_diff <- paste0(plot_title_diff_base, " [", filter_desc, "]")

## ---------- 新增：4个“参与分析行”的输出文件名 ----------
orig_ismc_rows_file    <- paste0("orig_relax_rows_ismc_", filter_suffix, ".txt")
orig_noismc_rows_file  <- paste0("orig_relax_rows_noismc_", filter_suffix, ".txt")
new_ismc_rows_file     <- paste0(new_relax_tag, "_relax_rows_ismc_", filter_suffix, ".txt")
new_noismc_rows_file   <- paste0(new_relax_tag, "_relax_rows_noismc_", filter_suffix, ".txt")

cat(sprintf("Using significance filter: %s\n", filter_desc))

## ---------- 读取原始数据 ----------
datPGID2Ref <- read.table(pgid2ref_file,
                          header = TRUE, stringsAsFactors = FALSE)

## 保留原始 relax_file 的列结构，供后续“按原格式导出”
datRelax_orig_raw <- read.table(orig_relax_file,
                                header = FALSE, sep = "\t", fill = TRUE, quote = "",
                                stringsAsFactors = FALSE)
orig_relax_colnames <- colnames(datRelax_orig_raw)

datRelax_orig <- datRelax_orig_raw

datGroupID2TransID <- read.table(groupid2transid_file,
                                 sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(datGroupID2TransID) <- c("OrthoID", "mrna")

datRelax_orig <- merge(datRelax_orig, datGroupID2TransID,
                       by.x = "V1", by.y = "OrthoID",
                       all.x = TRUE, all.y = FALSE)

datRelax_orig$V1 <- as.character(datRelax_orig$V1)

ismc_mrna <- read.table(ismc_mrna_file,
                        header = FALSE, sep = "\t", fill = TRUE, quote = "",
                        stringsAsFactors = FALSE)

noismc_mrna <- read.table(noismc_mrna_file,
                          header = FALSE, sep = "\t", fill = TRUE, quote = "",
                          stringsAsFactors = FALSE)

ismc_vec <- unique(as.character(ismc_mrna$V1))
noismc_vec <- unique(as.character(noismc_mrna$V1))

## ---------- 原始数据中标记 mrna 是否属于两个列表 ----------
datRelax_orig$is_ismc_mrna <- !is.na(datRelax_orig$mrna) & datRelax_orig$mrna %in% ismc_vec
datRelax_orig$is_noismc_mrna <- !is.na(datRelax_orig$mrna) & datRelax_orig$mrna %in% noismc_vec

## ---------- 为原始文件增加 P 和 FDR 列，并打筛选标记 ----------
datRelax_orig <- add_p_fdr_columns(datRelax_orig, p_col = "V5", fdr_col = fdr_colname)
datRelax_orig <- flag_sig_rows(datRelax_orig,
                               filter_mode = sig_filter_mode,
                               p_cutoff = sig_p_cutoff,
                               fdr_cutoff = sig_fdr_cutoff,
                               p_col = "P_numeric",
                               fdr_col = fdr_colname)

## ---------- 先在原始数据中按 V1 / OrthoID 统计 group 属于哪一类 ----------
## 规则：
## ismc    : group 中至少有一条 mrna 命中 ismc_mrna_file，且没有命中 noismc_mrna_file
## noismc  : group 中至少有一条 mrna 命中 noismc_mrna_file，且没有命中 ismc_mrna_file
## both    : 两边都命中，默认不参与比较
## neither : 两边都不命中，默认不参与比较

datRelax_orig_nonna <- datRelax_orig[!is.na(datRelax_orig$V1) & datRelax_orig$V1 != "", , drop = FALSE]
datRelax_orig_nonna$numeric_v9_flag <- as.integer(!is.na(to_num(datRelax_orig_nonna$V9)))

group_stat_orig <- make_group_count_table(
  datRelax_orig_nonna,
  list(
    n_rows_original = rep(1L, nrow(datRelax_orig_nonna)),
    n_ismc_mrna_rows_original = as.integer(datRelax_orig_nonna$is_ismc_mrna),
    n_noismc_mrna_rows_original = as.integer(datRelax_orig_nonna$is_noismc_mrna),
    n_numericV9_original = datRelax_orig_nonna$numeric_v9_flag
  )
)

group_stat_orig$group_class <- assign_group_class(
  group_stat_orig$n_ismc_mrna_rows_original,
  group_stat_orig$n_noismc_mrna_rows_original
)

group_stat_orig$assignment_rule <- assign_rule_text(group_stat_orig$group_class)
group_stat_orig$used_in_group_comparison <- ifelse(group_stat_orig$group_class %in% c("ismc", "noismc"), "yes", "no")

ismc_group_ids    <- sort(unique(group_stat_orig$V1[group_stat_orig$group_class == "ismc"]))
noismc_group_ids  <- sort(unique(group_stat_orig$V1[group_stat_orig$group_class == "noismc"]))
both_group_ids    <- sort(unique(group_stat_orig$V1[group_stat_orig$group_class == "both"]))
neither_group_ids <- sort(unique(group_stat_orig$V1[group_stat_orig$group_class == "neither"]))

cat(sprintf("[group record] unique ismc groups: %d\n", length(ismc_group_ids)))
cat(sprintf("[group record] unique noismc groups: %d\n", length(noismc_group_ids)))
cat(sprintf("[group record] unique both groups (excluded): %d\n", length(both_group_ids)))
cat(sprintf("[group record] unique neither groups (excluded): %d\n", length(neither_group_ids)))

## ---------- 原始分析：先按所选模式筛选，再比较 ismc-group vs noismc-group ----------
datRelax_orig_filtered <- datRelax_orig[datRelax_orig$pass_sig_filter, , drop = FALSE]
datRelax_orig_filtered <- merge(
  datRelax_orig_filtered,
  group_stat_orig[, c("V1", "group_class", "used_in_group_comparison")],
  by = "V1",
  all.x = TRUE,
  all.y = FALSE
)

cat(sprintf("[original] total rows before filter: %d; after filter: %d\n",
            nrow(datRelax_orig), nrow(datRelax_orig_filtered)))

datRelax_in  <- datRelax_orig_filtered[datRelax_orig_filtered$group_class == "ismc", , drop = FALSE]
datRelax_out <- datRelax_orig_filtered[datRelax_orig_filtered$group_class == "noismc", , drop = FALSE]
datRelax_excluded <- datRelax_orig_filtered[!(datRelax_orig_filtered$group_class %in% c("ismc", "noismc")), , drop = FALSE]

cat(sprintf("[original] rows after filter -> ismc-group: %d, noismc-group: %d, excluded(both/neither/unmapped): %d\n",
            nrow(datRelax_in), nrow(datRelax_out), nrow(datRelax_excluded)))

## ---------- 新增：导出原始 relax_file 中真正参与分析的两组行（保持原始格式） ----------
write_relax_like_rows(datRelax_in,  orig_relax_colnames, orig_ismc_rows_file)
write_relax_like_rows(datRelax_out, orig_relax_colnames, orig_noismc_rows_file)

orig_res <- run_v9_analysis(
  dat_in = datRelax_in,
  dat_out = datRelax_out,
  plot_title = plot_title_orig,
  out_prefix = orig_prefix,
  filter_mode = sig_filter_mode,
  p_cutoff = sig_p_cutoff,
  fdr_cutoff = sig_fdr_cutoff,
  x_label = x_label,
  y_label = y_label,
  caption_lab = caption_lab,
  group_order = group_order,
  group_labels = group_labels,
  group_colors = group_colors
)

## ---------- 记录：原始数据中，按当前筛选模式后每个 group 还剩多少行 ----------
datRelax_orig_filtered_nonna <- datRelax_orig_filtered[
  !is.na(datRelax_orig_filtered$V1) & datRelax_orig_filtered$V1 != "", , drop = FALSE
]
datRelax_orig_filtered_nonna$numeric_v9_flag <- as.integer(!is.na(to_num(datRelax_orig_filtered_nonna$V9)))

group_stat_orig_pass_sig <- make_group_count_table(
  datRelax_orig_filtered_nonna,
  list(
    n_rows_original_pass_sig = rep(1L, nrow(datRelax_orig_filtered_nonna)),
    n_numericV9_original_pass_sig = datRelax_orig_filtered_nonna$numeric_v9_flag
  )
)

## ---------- 读取新的 datRelax ----------
## 同样先保留原始列结构，供后续“按原格式导出”
datRelax_new_raw <- read.table(new_relax_file,
                               header = FALSE, sep = "\t", fill = TRUE, quote = "",
                               stringsAsFactors = FALSE)
new_relax_colnames <- colnames(datRelax_new_raw)

datRelax_new <- datRelax_new_raw
datRelax_new$V1 <- as.character(datRelax_new$V1)

## ---------- 为新文件增加 P 和 FDR 列，并打筛选标记 ----------
datRelax_new <- add_p_fdr_columns(datRelax_new, p_col = "V5", fdr_col = fdr_colname)
datRelax_new <- flag_sig_rows(datRelax_new,
                              filter_mode = sig_filter_mode,
                              p_cutoff = sig_p_cutoff,
                              fdr_cutoff = sig_fdr_cutoff,
                              p_col = "P_numeric",
                              fdr_col = fdr_colname)

## ---------- 用原始 group 分类去对应分析新的 datRelax ----------
group_map <- group_stat_orig[, c("V1", "group_class", "used_in_group_comparison")]
datRelax_new2 <- merge(datRelax_new, group_map, by = "V1", all.x = TRUE, all.y = FALSE)

## 先看新数据总体映射情况（筛选前）
datRelax_new_in_total       <- datRelax_new2[datRelax_new2$group_class == "ismc", , drop = FALSE]
datRelax_new_out_total      <- datRelax_new2[datRelax_new2$group_class == "noismc", , drop = FALSE]
datRelax_new_excluded_total <- datRelax_new2[!(datRelax_new2$group_class %in% c("ismc", "noismc")), , drop = FALSE]

cat(sprintf("[new: %s] total rows before filter -> ismc-group: %d, noismc-group: %d, excluded(both/neither/unmapped): %d, total: %d\n",
            new_relax_tag,
            nrow(datRelax_new_in_total),
            nrow(datRelax_new_out_total),
            nrow(datRelax_new_excluded_total),
            nrow(datRelax_new2)))

## 真正用于分析的，是按当前模式筛选后的行
datRelax_new2_filtered <- datRelax_new2[datRelax_new2$pass_sig_filter, , drop = FALSE]

datRelax_new_in       <- datRelax_new2_filtered[datRelax_new2_filtered$group_class == "ismc", , drop = FALSE]
datRelax_new_out      <- datRelax_new2_filtered[datRelax_new2_filtered$group_class == "noismc", , drop = FALSE]
datRelax_new_excluded <- datRelax_new2_filtered[!(datRelax_new2_filtered$group_class %in% c("ismc", "noismc")), , drop = FALSE]

cat(sprintf("[new: %s] rows after filter -> ismc-group: %d, noismc-group: %d, excluded(both/neither/unmapped): %d, total: %d\n",
            new_relax_tag,
            nrow(datRelax_new_in),
            nrow(datRelax_new_out),
            nrow(datRelax_new_excluded),
            nrow(datRelax_new2_filtered)))

## ---------- 新增：导出新 relax_file 中真正参与分析的两组行（保持原始格式） ----------
write_relax_like_rows(datRelax_new_in,  new_relax_colnames, new_ismc_rows_file)
write_relax_like_rows(datRelax_new_out, new_relax_colnames, new_noismc_rows_file)

new_res <- run_v9_analysis(
  dat_in = datRelax_new_in,
  dat_out = datRelax_new_out,
  plot_title = plot_title_new,
  out_prefix = new_prefix,
  filter_mode = sig_filter_mode,
  p_cutoff = sig_p_cutoff,
  fdr_cutoff = sig_fdr_cutoff,
  x_label = x_label,
  y_label = y_label,
  caption_lab = caption_lab,
  group_order = group_order,
  group_labels = group_labels,
  group_colors = group_colors
)

## ---------- 最终比较分析：全量 pairwise difference distribution 的 summary 对比 ----------
diff_res <- run_diff_distribution_analysis(
  orig_dat_in = datRelax_in,
  orig_dat_out = datRelax_out,
  new_dat_in = datRelax_new_in,
  new_dat_out = datRelax_new_out,
  out_prefix = diff_prefix,
  plot_title = plot_title_diff,
  filter_mode = sig_filter_mode,
  p_cutoff = sig_p_cutoff,
  fdr_cutoff = sig_fdr_cutoff,
  new_tag = new_relax_tag,
  diff_group_order = diff_group_order,
  diff_group_labels = diff_group_labels,
  diff_group_colors = diff_group_colors
)

## ---------- 输出记录 group 的表格文件 ----------
datRelax_new_nonna <- datRelax_new2[!is.na(datRelax_new2$V1) & datRelax_new2$V1 != "", , drop = FALSE]
datRelax_new_nonna$numeric_v9_flag_new <- as.integer(!is.na(to_num(datRelax_new_nonna$V9)))

new_group_stat <- make_group_count_table(
  datRelax_new_nonna,
  list(
    n_rows_new = rep(1L, nrow(datRelax_new_nonna)),
    n_numericV9_new = datRelax_new_nonna$numeric_v9_flag_new
  )
)

datRelax_new_filtered_nonna <- datRelax_new2_filtered[
  !is.na(datRelax_new2_filtered$V1) & datRelax_new2_filtered$V1 != "", , drop = FALSE
]
datRelax_new_filtered_nonna$numeric_v9_flag_new <- as.integer(!is.na(to_num(datRelax_new_filtered_nonna$V9)))

new_group_stat_pass_sig <- make_group_count_table(
  datRelax_new_filtered_nonna,
  list(
    n_rows_new_pass_sig = rep(1L, nrow(datRelax_new_filtered_nonna)),
    n_numericV9_new_pass_sig = datRelax_new_filtered_nonna$numeric_v9_flag_new
  )
)

group_record <- merge(group_stat_orig, group_stat_orig_pass_sig, by = "V1", all.x = TRUE, all.y = FALSE)
group_record <- merge(group_record, new_group_stat, by = "V1", all.x = TRUE, all.y = FALSE)
group_record <- merge(group_record, new_group_stat_pass_sig, by = "V1", all.x = TRUE, all.y = FALSE)

group_record <- fill_na_zero(group_record, c(
  "n_rows_original_pass_sig",
  "n_numericV9_original_pass_sig",
  "n_rows_new",
  "n_numericV9_new",
  "n_rows_new_pass_sig",
  "n_numericV9_new_pass_sig"
))

group_record$present_in_new_datRelax <- ifelse(group_record$n_rows_new > 0, "yes", "no")

colnames(group_record)[colnames(group_record) == "V1"] <- "OrthoID"

group_record <- group_record[, c(
  "group_class",
  "OrthoID",
  "n_rows_original",
  "n_ismc_mrna_rows_original",
  "n_noismc_mrna_rows_original",
  "n_numericV9_original",
  "n_rows_original_pass_sig",
  "n_numericV9_original_pass_sig",
  "assignment_rule",
  "used_in_group_comparison",
  "n_rows_new",
  "n_numericV9_new",
  "n_rows_new_pass_sig",
  "n_numericV9_new_pass_sig",
  "present_in_new_datRelax"
)]

group_record$filter_mode <- sig_filter_mode
group_record$filter_desc <- filter_desc

group_class_order <- c("ismc", "noismc", "both", "neither")
group_record$group_class <- factor(group_record$group_class, levels = group_class_order)
group_record <- group_record[order(group_record$group_class, group_record$OrthoID), , drop = FALSE]
group_record$group_class <- as.character(group_record$group_class)

write.table(group_record, file = group_record_tsv,
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Group record table saved: %s\n", group_record_tsv))

## ---------- 可选：分别导出原始/新文件中带 P/FDR/筛选标记/分组类别 的明细表 ----------
orig_annot_tsv <- paste0("annotated_original_relax_ismc_vs_noismc_", filter_suffix, ".tsv")
new_annot_tsv  <- paste0("annotated_", new_relax_tag, "_relax_ismc_vs_noismc_", filter_suffix, ".tsv")

write.table(datRelax_orig_filtered, file = orig_annot_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(datRelax_new2, file = new_annot_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Annotated original filtered table saved: %s\n", orig_annot_tsv))
cat(sprintf("Annotated new table saved: %s\n", new_annot_tsv))

cat("New relax-format exports finished:\n")
cat(sprintf("  %s\n", orig_ismc_rows_file))
cat(sprintf("  %s\n", orig_noismc_rows_file))
cat(sprintf("  %s\n", new_ismc_rows_file))
cat(sprintf("  %s\n", new_noismc_rows_file))

## ========== 结束 ==========
