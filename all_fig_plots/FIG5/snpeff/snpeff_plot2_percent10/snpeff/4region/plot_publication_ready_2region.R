# =========================================================
# 二项分布模拟 + FSR/noFSR 比值统计 + 论文风格柱状图 + 同频率 Tba/Tra p值
# 输出：PDF + PNG + 统计结果CSV
# =========================================================

# ---------- 0. 安装/加载包 ----------
need_pkgs <- c("dplyr", "tidyr", "purrr", "ggplot2", "readr", "scales")
to_install <- need_pkgs[!need_pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(need_pkgs, library, character.only = TRUE))

# ---------- 1. 参数 ----------
n_iter <- 10000
set.seed(20260319)

# 误差条：可选 "SD"、"SE"、"CI95"
# 建议论文图默认用 CI95（基于模拟分布的 2.5% 和 97.5% 分位数）
error_bar <- "CI95"

# p值方法：
# "wald"       = 基于 log(比值) 的 Wald 检验（推荐，结果不受模拟次数影响）
# "simulation" = 基于10000次模拟差值的双侧尾部概率
p_value_method <- "wald"

# 论文排版通常把图标题写在正文图注里；若需要标题可改成字符串
plot_title <- NULL
plot_caption <- paste0(
  "Bars show the mean simulated FSR/noFSR ratio across ", n_iter,
  " binomial simulations; error bars show ",
  ifelse(toupper(error_bar) == "CI95", "95% simulation interval",
         ifelse(toupper(error_bar) == "SE", "SE", "SD")),
  "."
)

outdir <- "binom_ratio_results"
dir.create(outdir, showWarnings = FALSE)

freq_levels <- c("<25%", "25%-50%", "50%-<100%", "100%")
type_levels <- c("Tba", "Tra")

# 输出图片参数
fig_width <- 8.6
fig_height <- 5.8
png_dpi <- 600

# 配色（色盲友好、适合论文）
#palette_fill <- c(
#  "Tba" = "#0F766E",  # deep teal
#  "Tra" = "#C46A2D"   # copper orange
#)
palette_fill <- c(
  "Tba" = "#6B6ECF",  # muted indigo
  "Tra" = "#A67C00"   # olive gold
)


#palette_fill <- c(
#  "Tba" = "#4E79A7",
#  "Tra" = "#E15759"
#)

# ---------- 2. 直接在这里填写原始计数 ----------
# 说明：
# group = test 表示 FSR
# group = control 表示 noFSR
# high  = 成功数
# syn   = 总数

raw_counts <- tibble::tribble(
  ~group,     ~measure_name, ~freq,        ~type, ~high, ~syn,
  "test",     "FSR",         "<25%",       "Tba",    307, 39497,
  "test",     "FSR",         "<25%",       "Tra",    116, 9288,
  "test",     "FSR",         "25%-50%",    "Tba",     95, 29821,
  "test",     "FSR",         "25%-50%",    "Tra",     45, 4880,
  "test",     "FSR",         "50%-<100%",  "Tba",    376, 163154,
  "test",     "FSR",         "50%-<100%",  "Tra",     39, 32763,
  "test",     "FSR",         "100%",       "Tba",   2716, 1257510,
  "test",     "FSR",         "100%",       "Tra",   1320, 727978,

  "control",  "noFSR",       "<25%",       "Tba",    229, 27891,
  "control",  "noFSR",       "<25%",       "Tra",    133, 14799,
  "control",  "noFSR",       "25%-50%",    "Tba",     83, 19939,
  "control",  "noFSR",       "25%-50%",    "Tra",     49, 9682,
  "control",  "noFSR",       "50%-<100%",  "Tba",    454, 118525,
  "control",  "noFSR",       "50%-<100%",  "Tra",     94, 34637,
  "control",  "noFSR",       "100%",       "Tba",   1586, 928192,
  "control",  "noFSR",       "100%",       "Tra",    766, 541756
)

# ---------- 3. 基本检查 ----------
stopifnot(all(raw_counts$high >= 0))
stopifnot(all(raw_counts$syn > 0))
stopifnot(all(raw_counts$high <= raw_counts$syn))
stopifnot(setequal(unique(raw_counts$group), c("test", "control")))
stopifnot(setequal(unique(raw_counts$type), type_levels))
stopifnot(setequal(unique(raw_counts$freq), freq_levels))

expected_n <- 2 * length(freq_levels) * length(type_levels)
if (nrow(raw_counts) != expected_n) {
  stop("raw_counts 行数不对，请确保 test/control × 4个频率 × 2个类型 全部填写。")
}

dup_check <- raw_counts %>%
  dplyr::count(group, freq, type) %>%
  dplyr::filter(n != 1)

if (nrow(dup_check) > 0) {
  stop("存在重复或缺失的 group/freq/type 组合，请检查 raw_counts。")
}

raw_counts <- raw_counts %>%
  dplyr::mutate(
    group = factor(group, levels = c("test", "control")),
    freq  = factor(freq, levels = freq_levels),
    type  = factor(type, levels = type_levels)
  )

# ---------- 4. 二项分布模拟函数 ----------
simulate_binom <- function(group, measure_name, freq, type, high, syn, n_iter) {
  p_obs <- high / syn
  sim_success <- rbinom(n = n_iter, size = syn, prob = p_obs)

  tibble::tibble(
    group = as.character(group),
    measure_name = as.character(measure_name),
    freq = as.character(freq),
    type = as.character(type),
    iter = seq_len(n_iter),
    high = high,
    syn = syn,
    observed_prop = p_obs,
    sim_success = sim_success,
    sim_prop = sim_success / syn
  )
}

# ---------- 5. 对每个 group/freq/type 做 n_iter 次二项分布模拟 ----------
sim_values <- purrr::pmap_dfr(raw_counts, simulate_binom, n_iter = n_iter) %>%
  dplyr::mutate(
    freq = factor(freq, levels = freq_levels),
    type = factor(type, levels = type_levels)
  )

# ---------- 6. 各组合（FSR 或 noFSR）的模拟汇总 ----------
group_summary <- sim_values %>%
  dplyr::group_by(group, measure_name, freq, type) %>%
  dplyr::summarise(
    high = dplyr::first(high),
    syn = dplyr::first(syn),
    observed_prop = dplyr::first(observed_prop),
    mean_prop = mean(sim_prop),
    sd_prop = sd(sim_prop),
    se_prop = sd_prop / sqrt(dplyr::n()),
    lower_prop_95 = stats::quantile(sim_prop, 0.025, na.rm = TRUE),
    upper_prop_95 = stats::quantile(sim_prop, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# ---------- 7. 计算 FSR / noFSR 的 n_iter 次比值 ----------
ratio_values <- sim_values %>%
  dplyr::select(group, freq, type, iter, sim_prop) %>%
  tidyr::pivot_wider(names_from = group, values_from = sim_prop) %>%
  dplyr::rename(
    test_prop = test,
    control_prop = control
  ) %>%
  dplyr::mutate(
    ratio = dplyr::if_else(control_prop == 0, NA_real_, test_prop / control_prop)
  )

# 原始观测值对应的比例和比值
observed_ratio <- raw_counts %>%
  dplyr::mutate(observed_prop = high / syn) %>%
  dplyr::select(group, freq, type, observed_prop) %>%
  tidyr::pivot_wider(names_from = group, values_from = observed_prop, names_prefix = "obs_") %>%
  dplyr::mutate(obs_ratio = obs_test / obs_control)

# ---------- 8. 比值统计汇总 ----------
ratio_summary <- ratio_values %>%
  dplyr::group_by(freq, type) %>%
  dplyr::summarise(
    n_ratio = sum(!is.na(ratio)),
    mean_test_prop = mean(test_prop, na.rm = TRUE),
    sd_test_prop = sd(test_prop, na.rm = TRUE),
    se_test_prop = sd_test_prop / sqrt(sum(!is.na(test_prop))),
    mean_control_prop = mean(control_prop, na.rm = TRUE),
    sd_control_prop = sd(control_prop, na.rm = TRUE),
    se_control_prop = sd_control_prop / sqrt(sum(!is.na(control_prop))),
    mean_ratio = mean(ratio, na.rm = TRUE),
    sd_ratio = sd(ratio, na.rm = TRUE),
    se_ratio = sd_ratio / sqrt(n_ratio),
    lower_ratio_95 = stats::quantile(ratio, 0.025, na.rm = TRUE),
    upper_ratio_95 = stats::quantile(ratio, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(observed_ratio, by = c("freq", "type"))

# ---------- 9. 计算同一频率下 Tba vs Tra 的 p 值 ----------
# 9A. Wald 检验：比较两个 log(FSR/noFSR) 是否不同
get_pvalues_wald <- function(raw_counts) {
  log_ratio_tbl <- raw_counts %>%
    dplyr::select(group, freq, type, high, syn) %>%
    tidyr::pivot_wider(names_from = group, values_from = c(high, syn), names_sep = "_") %>%
    dplyr::mutate(
      p_test = ifelse(
        high_test == 0 | high_test == syn_test,
        (high_test + 0.5) / (syn_test + 1),
        high_test / syn_test
      ),
      p_control = ifelse(
        high_control == 0 | high_control == syn_control,
        (high_control + 0.5) / (syn_control + 1),
        high_control / syn_control
      ),
      log_ratio = log(p_test / p_control),
      se_log_ratio = sqrt(
        (1 - p_test) / (syn_test * p_test) +
          (1 - p_control) / (syn_control * p_control)
      )
    ) %>%
    dplyr::select(freq, type, log_ratio, se_log_ratio)

  log_ratio_tbl %>%
    tidyr::pivot_wider(
      names_from = type,
      values_from = c(log_ratio, se_log_ratio),
      names_sep = "_"
    ) %>%
    dplyr::mutate(
      z = (log_ratio_Tba - log_ratio_Tra) /
        sqrt(se_log_ratio_Tba^2 + se_log_ratio_Tra^2),
      p_value = 2 * stats::pnorm(-abs(z))
    ) %>%
    dplyr::select(freq, p_value)
}

# 9B. 基于模拟差值的双侧尾部概率
get_pvalues_sim <- function(ratio_values) {
  ratio_values %>%
    dplyr::select(freq, type, iter, ratio) %>%
    tidyr::pivot_wider(names_from = type, values_from = ratio) %>%
    dplyr::group_by(freq) %>%
    dplyr::summarise(
      p_value = {
        d <- Tba - Tra
        d <- d[is.finite(d)]
        if (length(d) == 0) {
          NA_real_
        } else {
          2 * min(mean(d >= 0), mean(d <= 0))
        }
      },
      .groups = "drop"
    )
}

if (!tolower(p_value_method) %in% c("wald", "simulation")) {
  stop("p_value_method 只能是 'wald' 或 'simulation'")
}

p_values <- if (tolower(p_value_method) == "simulation") {
  get_pvalues_sim(ratio_values)
} else {
  get_pvalues_wald(raw_counts)
}

p_values <- p_values %>%
  dplyr::mutate(
    freq = factor(freq, levels = freq_levels),
    p_label = dplyr::case_when(
      is.na(p_value)  ~ "p = NA",
      p_value < 0.001 ~ "p < 0.001",
      TRUE            ~ sprintf("p = %.3f", p_value)
    )
  )

# ---------- 10. 为作图准备数据 ----------
plot_df <- ratio_summary %>%
  dplyr::mutate(
    freq = factor(freq, levels = freq_levels),
    type = factor(type, levels = type_levels),
    x_id = match(freq, freq_levels),
    x = x_id + dplyr::case_when(
      type == "Tba" ~ -0.18,
      type == "Tra" ~  0.18,
      TRUE ~ 0
    ),
    error = dplyr::case_when(
      toupper(error_bar) == "SE"   ~ se_ratio,
      toupper(error_bar) == "SD"   ~ sd_ratio,
      toupper(error_bar) == "CI95" ~ NA_real_,
      TRUE ~ NA_real_
    ),
    ymin = dplyr::case_when(
      toupper(error_bar) == "CI95" ~ lower_ratio_95,
      TRUE ~ pmax(mean_ratio - error, 0)
    ),
    ymax = dplyr::case_when(
      toupper(error_bar) == "CI95" ~ upper_ratio_95,
      TRUE ~ mean_ratio + error
    )
  )

max_y <- max(plot_df$ymax, na.rm = TRUE)
tip_len <- max(0.03 * max_y, 0.02)

p_annot <- plot_df %>%
  dplyr::group_by(freq, x_id) %>%
  dplyr::summarise(
    y_bracket = max(ymax, na.rm = TRUE) * 1.06,
    .groups = "drop"
  ) %>%
  dplyr::left_join(p_values, by = "freq") %>%
  dplyr::mutate(
    x_left = x_id - 0.18,
    x_right = x_id + 0.18,
    x_mid = x_id,
    y_text = y_bracket + tip_len * 0.55,
    y_tip = y_bracket - tip_len
  )

# ---------- 11. 作图（论文风格） ----------
p <- ggplot2::ggplot() +
  ggplot2::geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey50"
  ) +
  ggplot2::geom_col(
    data = plot_df,
    ggplot2::aes(x = x, y = mean_ratio, fill = type),
    width = 0.28,
    color = "white",
    linewidth = 0.5,
    alpha = 0.98
  ) +
  ggplot2::geom_errorbar(
    data = plot_df,
    ggplot2::aes(x = x, ymin = ymin, ymax = ymax),
    width = 0.06,
    linewidth = 0.55,
    color = "black"
  ) +
  ggplot2::geom_segment(
    data = p_annot,
    ggplot2::aes(x = x_left, xend = x_right, y = y_bracket, yend = y_bracket),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  ggplot2::geom_segment(
    data = p_annot,
    ggplot2::aes(x = x_left, xend = x_left, y = y_bracket, yend = y_tip),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  ggplot2::geom_segment(
    data = p_annot,
    ggplot2::aes(x = x_right, xend = x_right, y = y_bracket, yend = y_tip),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  ggplot2::geom_text(
    data = p_annot,
    ggplot2::aes(x = x_mid, y = y_text, label = p_label),
    inherit.aes = FALSE,
    vjust = 0,
    size = 3.8,
    fontface = "bold"
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq_along(freq_levels),
    labels = freq_levels,
    expand = ggplot2::expansion(mult = c(0.04, 0.04))
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(accuracy = 0.01),
    expand = ggplot2::expansion(mult = c(0, 0.15))
  ) +
  ggplot2::scale_fill_manual(values = palette_fill, name = NULL) +
  ggplot2::labs(
    title = plot_title,
    x = "Frequency",
    y = "FSR / noFSR",
    caption = plot_caption
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(color = "grey88", linewidth = 0.35),
    panel.border = ggplot2::element_rect(color = "black", linewidth = 0.7),
    axis.title = ggplot2::element_text(face = "bold", color = "black"),
    axis.text = ggplot2::element_text(color = "black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    plot.caption = ggplot2::element_text(hjust = 0, size = 9, color = "grey30"),
    plot.margin = ggplot2::margin(t = 10, r = 18, b = 10, l = 10)
  )

print(p)

# ---------- 12. 导出结果 ----------
final_summary <- ratio_summary %>%
  dplyr::left_join(p_values %>% dplyr::select(freq, p_value, p_label), by = "freq") %>%
  dplyr::arrange(freq, type)

readr::write_csv(sim_values, file.path(outdir, "binomial_simulated_values_each_group.csv"))
readr::write_csv(group_summary, file.path(outdir, "group_summary_test_control.csv"))
readr::write_csv(ratio_values, file.path(outdir, "ratio_values_FSR_div_noFSR.csv"))
readr::write_csv(p_values, file.path(outdir, "p_values_Tba_vs_Tra_by_frequency.csv"))
readr::write_csv(final_summary, file.path(outdir, "ratio_summary_with_pvalues.csv"))

# 为了在不同系统上获得更稳定的 PDF 输出，优先使用 cairo_pdf
save_pdf_plot <- function(plot_obj, filename, width, height) {
  if (capabilities("cairo")) {
    ggplot2::ggsave(
      filename = filename,
      plot = plot_obj,
      width = width,
      height = height,
      units = "in",
      device = grDevices::cairo_pdf,
      bg = "white"
    )
  } else {
    ggplot2::ggsave(
      filename = filename,
      plot = plot_obj,
      width = width,
      height = height,
      units = "in",
      device = "pdf",
      bg = "white"
    )
  }
}

save_pdf_plot(
  plot_obj = p,
  filename = file.path(outdir, "FSR_noFSR_barplot_publication.pdf"),
  width = fig_width,
  height = fig_height
)

ggplot2::ggsave(
  filename = file.path(outdir, "FSR_noFSR_barplot_publication.png"),
  plot = p,
  width = fig_width,
  height = fig_height,
  units = "in",
  dpi = png_dpi,
  bg = "white"
)

# ---------- 13. 在控制台查看核心结果 ----------
print(final_summary)
