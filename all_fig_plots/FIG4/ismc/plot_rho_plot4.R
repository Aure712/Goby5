options(scipen = 999)

############################################################
## 1. 用户需要修改的配置区
############################################################

## 每个 list 代表一个物种：
## name         : 物种名，会显示在图标题和 legend 中
## path         : 该物种 rho 结果所在的总目录
## nChr         : 该物种染色体数量
## chrom_prefix : 输出 rho 文件时染色体名前缀
## color        : 该物种画图颜色
##
## 假设每个物种目录结构如下：
## /path/to/species_A/1/out.rho.10kb.bedgraph
## /path/to/species_A/2/out.rho.10kb.bedgraph
## ...
## /path/to/species_A/nChr/out.rho.10kb.bedgraph

species_config <- list(
  list(
    name = "Tba",
    path = "/fast3/group_crf/home/g21shaoy23/goby5/ismc/tba/20inds/",
    nChr = 22,
    chrom_prefix = "Tba_chr",
    color = "blue"
  ),
  list(
    name = "Tbi",
    path = "/fast3/group_crf/home/g21shaoy23/goby5/ismc/tbi/20inds/",
    nChr = 22,
    chrom_prefix = "Tbi_chr",
    color = "red"
  ),
  list(
    name = "Tra",
    path = "/fast3/group_crf/home/g21shaoy23/goby5/ismc/tra/20inds/",
    nChr = 11,
    chrom_prefix = "Tra_chr",
    color = "darkgreen"
  )
)

## 输入文件名
## 注意：文件名虽然还是 out.rho.10kb.bedgraph，
## 但现在实际原始窗口长度由 input_window_size 控制。
input_file_name <- "out.rho.10kb.bedgraph"

############################################################
## 新增：窗口长度和平滑设置
############################################################

## 原始输入数据的窗口长度。
## 你现在的数据是 100 kbp，所以这里设置为 100000。
## 该参数用于根据 chromEnd 重新计算 chromStart。
input_window_size <- 100000

## 是否对画图曲线进行平滑。
## TRUE  : 如果 plot_window_size > input_window_size，则把多个原始窗口合并取平均后画图。
## FALSE : 不平滑，直接用原始 100 kbp 窗口画图。
smooth_curve <- TRUE

## 绘图使用的新窗口长度。
## 例如：
## 100000  = 不改变窗口，仍然使用 100 kbp；
## 500000  = 每 500 kbp 合并一次，平均其中包含的 100 kbp 窗口；
## 1000000 = 每 1 Mbp 合并一次。
plot_window_size <- 1000000

## 是否输出用于画图的窗口数据。
## 如果开启平滑，这个文件就是平滑后的数据；
## 如果不开启平滑，这个文件和原始窗口基本一致。
write_plot_window_rho_files <- TRUE

## 原脚本中样本列的选择逻辑是：
## dat19ind[, 4:(ncol(dat19ind)-3)]
## 这里保留这个逻辑，但做成参数，方便后续调整。
sample_col_start <- 4
n_tail_non_sample <- 3

## 输出目录
out_dir <- "rho_multi_species_output"

## 是否生成“同一染色体编号不同物种叠加图”
make_overlay_plot <- TRUE

## 是否额外生成：
## 所有 Tra 锚定染色体对应关系都画在同一页 PDF 中
make_tra_aligned_one_page_plot <- TRUE

############################################################
## 指定哪些染色体画图时需要翻转
############################################################

## 这里的翻转只影响画图，不改变 rho 输出文件中的真实坐标。
##
## 写法：
## flip_chromosomes <- list(
##   Tba = c(16, 15, 12),
##   Tbi = c(14, 15, 13),
##   Tra = c(1)
## )
##
## 表示：
## Tba 的 16、15、12 号染色体反向画；
## Tbi 的 14、15、13 号染色体反向画；
## Tra 的 1 号染色体反向画。
##
## 如果某个物种没有需要翻转的染色体，可以写成 c()。

flip_chromosomes <- list(
  Tba = c(16, 2, 4, 3, 11, 13, 9, 6, 18),
  Tbi = c(14, 13, 1, 7, 21, 18, 22, 5, 11),
  Tra = c()
)

## 新增图中物种行顺序
## 第一行通常建议放 Tra，因为它是锚定物种
aligned_species_row_order <- c("Tra", "Tbi", "Tba")

## 新增图的 y 轴范围设置：
## "cell"   : 每个小图独立 y 轴范围，最容易看清每条曲线
## "row"    : 同一物种一整行使用相同 y 轴范围
## "global" : 整张图所有小图使用相同 y 轴范围
aligned_y_scale <- "cell"

## 每组染色体之间拼接时的空隙比例
aligned_chr_gap_ratio <- 0.05

## Tra 锚定的染色体对应关系配置
##
## 每个 list 是图中的一列。
## group_label : 这一列的标题
## Tra         : Tra 对应染色体，通常只有 1 条
## Tbi         : 与该 Tra 染色体对应的 Tbi 染色体，可以是 1 条或多条
## Tba         : 与该 Tra 染色体对应的 Tba 染色体，可以是 1 条或多条

tra_chr_groups <- list(
  list(
    group_label = "Tra1",
    Tra = c(1),
    Tbi = c(14, 15, 13),
    Tba = c(16, 15, 12)
  ),
  list(
    group_label = "Tra2",
    Tra = c(2),
    Tbi = c(4, 2),
    Tba = c(5, 2)
  ),
  list(
    group_label = "Tra3",
    Tra = c(3),
    Tbi = c(16, 1),
    Tba = c(17, 1)
  ),
  list(
    group_label = "Tra4",
    Tra = c(4),
    Tbi = c(9, 3),
    Tba = c(4, 3)
  ),
  list(
    group_label = "Tra5",
    Tra = c(5),
    Tbi = c(7, 12),
    Tba = c(11, 13)
  ),
  list(
    group_label = "Tra6",
    Tra = c(6),
    Tbi = c(19, 8),
    Tba = c(19, 9)
  ),
  list(
    group_label = "Tra7",
    Tra = c(7),
    Tbi = c(6, 21),
    Tba = c(6, 20)
  ),
  list(
    group_label = "Tra8",
    Tra = c(8),
    Tbi = c(20, 10),
    Tba = c(21, 8)
  ),
  list(
    group_label = "Tra9",
    Tra = c(9),
    Tbi = c(18, 17),
    Tba = c(14, 18)
  ),
  list(
    group_label = "Tra10",
    Tra = c(10),
    Tbi = c(22, 5),
    Tba = c(22, 7)
  ),
  list(
    group_label = "Tra11",
    Tra = c(11),
    Tbi = c(11),
    Tba = c(10)
  )
)


############################################################
## 2. 输出文件准备
############################################################

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

rho_all_file <- file.path(out_dir, "rho.all_species.txt")
rho_plot_all_file <- file.path(out_dir, "rho.plot_window.all_species.txt")
r2_out_file <- file.path(out_dir, "rho.screen.out.txt")
pdf_file <- file.path(out_dir, "rho.plot.by_species.pdf")
overlay_pdf_file <- file.path(out_dir, "rho.plot.overlay_by_chr.pdf")
tra_aligned_pdf_file <- file.path(out_dir, "rho.plot.tra_aligned_one_page.pdf")

## 删除旧结果，避免 append 时混入旧数据
unlink(rho_all_file)
unlink(rho_plot_all_file)
unlink(r2_out_file)
unlink(pdf_file)
unlink(overlay_pdf_file)
unlink(tra_aligned_pdf_file)
unlink(file.path(out_dir, "rho.plot_window.*.txt"))

## 用于保存所有物种所有染色体的数据
all_plot_data <- list()


############################################################
## 3. 辅助函数
############################################################

safe_file_name <- function(x) {
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[^A-Za-z0-9_.-]", "_", x)
  return(x)
}

append_table_with_header <- function(dat, file) {
  write.table(
    dat,
    file = file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = !file.exists(file),
    append = file.exists(file)
  )
}

append_table_no_header <- function(dat, file) {
  write.table(
    dat,
    file = file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE
  )
}

calc_r2_for_samples <- function(dat, sample_cols) {
  arrR2 <- c()

  for (sSample in sample_cols) {
    x <- dat[[sSample]]
    y <- dat$sample_mean

    ok <- is.finite(x) & is.finite(y)

    if (sum(ok) < 2) {
      arrR2 <- c(arrR2, NA_real_)
      next
    }

    if (sd(x[ok], na.rm = TRUE) == 0 || sd(y[ok], na.rm = TRUE) == 0) {
      arrR2 <- c(arrR2, NA_real_)
      next
    }

    fit <- lm(y[ok] ~ x[ok])
    oS <- summary(fit)
    arrR2 <- c(arrR2, oS$r.squared)
  }

  return(arrR2)
}

get_species_color <- function(species_name, species_config) {
  for (sp in species_config) {
    if (sp$name == species_name) {
      return(sp$color)
    }
  }
  return("black")
}

should_flip_chr <- function(species_name, chr_id, flip_chromosomes) {
  if (is.null(flip_chromosomes[[species_name]])) {
    return(FALSE)
  }

  chr_vec <- flip_chromosomes[[species_name]]

  if (length(chr_vec) == 0) {
    return(FALSE)
  }

  return(chr_id %in% chr_vec)
}

reverse_x_coord <- function(x) {
  ok <- is.finite(x)

  if (sum(ok) < 2) {
    return(x)
  }

  xmin <- min(x[ok], na.rm = TRUE)
  xmax <- max(x[ok], na.rm = TRUE)

  x_new <- x
  x_new[ok] <- xmin + xmax - x[ok]

  return(x_new)
}

get_range_or_null <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) == 0) {
    return(NULL)
  }

  r <- range(x, na.rm = TRUE)

  if (!is.finite(r[1]) || !is.finite(r[2])) {
    return(NULL)
  }

  if (r[1] == r[2]) {
    r[1] <- r[1] - 1
    r[2] <- r[2] + 1
  }

  return(r)
}

smooth_data_by_window <- function(dat,
                                  input_window_size = 100000,
                                  plot_window_size = 100000,
                                  smooth_curve = TRUE) {
  ## dat 需要至少包含：
  ## chrom, chromStart, chromEnd, chromMid, sample_mean

  dat <- dat[order(dat$chromStart, dat$chromEnd), , drop = FALSE]

  dat$n_old_windows <- 1

  if (!isTRUE(smooth_curve)) {
    return(dat)
  }

  if (is.null(plot_window_size) || is.na(plot_window_size)) {
    return(dat)
  }

  if (plot_window_size <= input_window_size) {
    return(dat)
  }

  ok <- is.finite(dat$chromStart) & is.finite(dat$chromEnd) & is.finite(dat$sample_mean)

  if (sum(ok) == 0) {
    return(dat)
  }

  dat_ok <- dat[ok, , drop = FALSE]

  base_start <- min(dat_ok$chromStart, na.rm = TRUE)

  ## 按新窗口划分。
  ## 一个新窗口包含所有 chromStart 落入该窗口范围内的旧窗口。
  bin_id <- floor((dat_ok$chromStart - base_start) / plot_window_size)

  idx_list <- split(seq_len(nrow(dat_ok)), bin_id)

  out_list <- lapply(idx_list, function(idx) {
    chr_start <- min(dat_ok$chromStart[idx], na.rm = TRUE)
    chr_end <- max(dat_ok$chromEnd[idx], na.rm = TRUE)

    data.frame(
      chrom = dat_ok$chrom[idx[1]],
      chromStart = chr_start,
      chromEnd = chr_end,
      chromMid = (chr_start + chr_end) / 2,
      sample_mean = mean(dat_ok$sample_mean[idx], na.rm = TRUE),
      n_old_windows = length(idx),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out_list)
  rownames(out) <- NULL

  return(out)
}

get_chr_data <- function(all_plot_data, species_name, chr_id) {
  key <- paste0("Chr_", chr_id)

  if (is.null(all_plot_data[[key]])) {
    return(NULL)
  }

  if (is.null(all_plot_data[[key]][[species_name]])) {
    return(NULL)
  }

  return(all_plot_data[[key]][[species_name]])
}

collect_y_values_for_chr_vec <- function(all_plot_data, species_name, chr_vec) {
  y_values <- c()

  for (chr_id in chr_vec) {
    one_chr <- get_chr_data(all_plot_data, species_name, chr_id)

    if (is.null(one_chr)) {
      next
    }

    y <- one_chr$sample_mean
    y <- y[is.finite(y)]

    if (length(y) > 0) {
      y_values <- c(y_values, y)
    }
  }

  return(y_values)
}

collect_y_values_for_group <- function(all_plot_data, group_info, species_names) {
  y_values <- c()

  for (species_name in species_names) {
    chr_vec <- group_info[[species_name]]

    if (is.null(chr_vec)) {
      next
    }

    y_values <- c(
      y_values,
      collect_y_values_for_chr_vec(all_plot_data, species_name, chr_vec)
    )
  }

  return(y_values)
}

make_concatenated_chr_curve <- function(all_plot_data,
                                        species_name,
                                        chr_vec,
                                        gap_ratio = 0.05) {
  out <- list()
  chr_lengths <- c()

  for (chr_id in chr_vec) {
    one_chr <- get_chr_data(all_plot_data, species_name, chr_id)

    if (is.null(one_chr)) {
      next
    }

    x <- one_chr$chromMid
    y <- one_chr$sample_mean
    is_flipped <- one_chr$is_flipped

    ok <- is.finite(x) & is.finite(y)

    if (sum(ok) < 2) {
      next
    }

    x <- x[ok]
    y <- y[ok]

    ## 转为每条染色体内部的相对坐标，方便多条染色体拼接
    x_local <- x - min(x, na.rm = TRUE)

    chr_len <- max(x_local, na.rm = TRUE)

    if (!is.finite(chr_len) || chr_len <= 0) {
      chr_len <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    }

    if (!is.finite(chr_len) || chr_len <= 0) {
      chr_len <- 1
    }

    ## 如果该物种该染色体被指定翻转，则反向绘制
    if (isTRUE(is_flipped)) {
      x_local <- chr_len - x_local
    }

    ## 排序，保证线条按画图坐标从左到右连接
    ord <- order(x_local)
    x_local <- x_local[ord]
    y <- y[ord]

    chr_lengths <- c(chr_lengths, chr_len)

    out[[length(out) + 1]] <- list(
      chr_id = chr_id,
      x_local = x_local,
      y = y,
      chr_len = chr_len,
      is_flipped = is_flipped
    )
  }

  if (length(out) == 0) {
    return(NULL)
  }

  max_len <- max(chr_lengths, na.rm = TRUE)
  gap <- max_len * gap_ratio

  offset <- 0
  boundary_positions <- c()
  label_positions <- c()

  for (i in seq_along(out)) {
    out[[i]]$x_plot <- out[[i]]$x_local + offset
    label_positions <- c(label_positions, offset + out[[i]]$chr_len / 2)

    offset <- offset + out[[i]]$chr_len

    if (i < length(out)) {
      boundary_positions <- c(boundary_positions, offset + gap / 2)
      offset <- offset + gap
    }
  }

  return(list(
    curves = out,
    boundary_positions = boundary_positions,
    label_positions = label_positions,
    xlim = c(0, offset)
  ))
}

draw_one_aligned_cell <- function(all_plot_data,
                                  species_name,
                                  chr_vec,
                                  species_config,
                                  group_label = "",
                                  show_group_title = FALSE,
                                  show_species_label = FALSE,
                                  y_lim = NULL,
                                  gap_ratio = 0.05) {
  species_color <- get_species_color(species_name, species_config)

  cell_data <- make_concatenated_chr_curve(
    all_plot_data = all_plot_data,
    species_name = species_name,
    chr_vec = chr_vec,
    gap_ratio = gap_ratio
  )

  if (is.null(cell_data)) {
    plot.new()
    box()

    if (show_group_title) {
      title(main = group_label, cex.main = 0.9)
    }

    if (show_species_label) {
      mtext(species_name, side = 2, line = 2.2, cex = 0.9, las = 0)
    }

    text(
      x = 0.5,
      y = 0.5,
      labels = paste0(species_name, "\nmissing"),
      cex = 0.7
    )

    return(invisible(NULL))
  }

  y_values <- unlist(lapply(cell_data$curves, function(z) z$y))

  if (is.null(y_lim)) {
    y_lim <- get_range_or_null(y_values)
  }

  if (is.null(y_lim)) {
    y_lim <- c(0, 1)
  }

  plot(
    NA,
    NA,
    xlim = cell_data$xlim,
    ylim = y_lim,
    xlab = "",
    ylab = "",
    xaxt = "n",
    yaxt = "n",
    main = ifelse(show_group_title, group_label, "")
  )

  box(col = "black", lty = 1, lwd = 1)

  ## y 轴只在每一行最左侧显示，避免太拥挤
  if (show_species_label) {
    axis(side = 2, las = 1, cex.axis = 0.55)
    mtext(species_name, side = 2, line = 2.3, cex = 0.8, las = 0)
  }

  ## 画每条染色体曲线
  for (i in seq_along(cell_data$curves)) {
    curve_i <- cell_data$curves[[i]]

    lines(
      curve_i$x_plot,
      curve_i$y,
      col = species_color,
      lwd = 0.8
    )
  }

  ##########################################################
  ## 修改点 1：
  ## 多条染色体之间的分割线改为黑色实线，
  ## 与小图边框风格一致。
  ##########################################################

  if (length(cell_data$boundary_positions) > 0) {
    abline(
      v = cell_data$boundary_positions,
      col = "black",
      lty = 1,
      lwd = 1
    )
  }

  ## 在每个小图底部标出染色体编号；
  ## 被翻转的染色体后面加 "_rev"
  chr_labels <- sapply(
    cell_data$curves,
    function(z) {
      paste0(
        "Chr",
        z$chr_id,
        ifelse(isTRUE(z$is_flipped), "_rev", "")
      )
    }
  )

  mtext(
    text = chr_labels,
    side = 1,
    at = cell_data$label_positions,
    line = 0.2,
    cex = 0.45
  )

  return(invisible(NULL))
}


############################################################
## 4. 主程序：逐物种、逐染色体读取数据并绘图
############################################################

pdf(file = pdf_file, width = 10, height = 5)
sink(file = r2_out_file)

for (sp in species_config) {

  species_name <- sp$name
  species_path <- sp$path
  species_nChr <- sp$nChr
  species_chrom_prefix <- sp$chrom_prefix
  species_color <- sp$color

  species_safe_name <- safe_file_name(species_name)
  rho_species_file <- file.path(out_dir, paste0("rho.", species_safe_name, ".txt"))
  rho_plot_species_file <- file.path(out_dir, paste0("rho.plot_window.", species_safe_name, ".txt"))

  ## 删除该物种旧的 rho 文件
  unlink(rho_species_file)
  unlink(rho_plot_species_file)

  cat("############################################################\n")
  cat("# Species:", species_name, "\n")
  cat("# Path:", species_path, "\n")
  cat("# Chromosome number:", species_nChr, "\n")
  cat("# Input window size:", input_window_size, "\n")
  cat("# Smooth curve:", smooth_curve, "\n")
  cat("# Plot window size:", plot_window_size, "\n")
  cat("# Flipped chromosomes:", paste(flip_chromosomes[[species_name]], collapse = ","), "\n")
  cat("############################################################\n")

  for (nChr in 1:species_nChr) {

    sF <- file.path(species_path, as.character(nChr), input_file_name)

    if (!file.exists(sF)) {
      cat("#Skip", species_name, "Chr", nChr, "file not found:", sF, "\n")
      next
    }

    dat20ind <- read.table(
      sF,
      header = TRUE,
      sep = "\t",
      fill = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    ## 检查必要列
    if (!"chromEnd" %in% colnames(dat20ind)) {
      cat("#Error", species_name, "Chr", nChr, "missing column: chromEnd in", sF, "\n")
      next
    }

    ## 根据原始输入窗口大小重新计算 chromStart
    dat20ind$chromStart <- dat20ind$chromEnd - (input_window_size - 1)

    ## 计算原始窗口中点
    dat20ind$chromMid <- (dat20ind$chromStart + dat20ind$chromEnd) / 2

    ## 删除空列名的列
    dat19ind <- dat20ind[, colnames(dat20ind) != "", drop = FALSE]

    ## 保留原脚本的样本列选择逻辑：
    ## 4:(ncol(dat19ind)-3)
    sample_col_end <- ncol(dat19ind) - n_tail_non_sample

    if (sample_col_end < sample_col_start) {
      cat(
        "#Error", species_name, "Chr", nChr,
        "cannot identify sample columns. ncol =", ncol(dat19ind),
        "sample_col_start =", sample_col_start,
        "sample_col_end =", sample_col_end, "\n"
      )
      next
    }

    sample_cols <- colnames(dat19ind)[sample_col_start:sample_col_end]

    ## 样本列转成 numeric，避免字符型导致 rowMeans 报错
    dat19ind[, sample_cols] <- lapply(
      dat19ind[, sample_cols, drop = FALSE],
      function(x) as.numeric(as.character(x))
    )

    ## 计算每个原始窗口跨样本平均 rho
    dat19ind$sample_mean <- rowMeans(
      dat19ind[, sample_cols, drop = FALSE],
      na.rm = TRUE
    )

    ## 设置输出用的染色体名
    dat19ind$chrom <- paste0(species_chrom_prefix, nChr)

    ## 判断当前物种当前染色体是否需要反向绘图
    is_flipped <- should_flip_chr(
      species_name = species_name,
      chr_id = nChr,
      flip_chromosomes = flip_chromosomes
    )

    ########################################################
    ## 新增：根据 plot_window_size 可选平滑，用于画图
    ########################################################

    dat_plot <- dat19ind[, c("chrom", "chromStart", "chromEnd", "chromMid", "sample_mean"), drop = FALSE]

    dat_plot <- smooth_data_by_window(
      dat = dat_plot,
      input_window_size = input_window_size,
      plot_window_size = plot_window_size,
      smooth_curve = smooth_curve
    )

    ## 绘图用坐标。
    ## chromMid 保持真实坐标；
    ## chromMid_plot 只用于画图。
    dat_plot$chromMid_plot <- dat_plot$chromMid

    if (isTRUE(is_flipped)) {
      dat_plot$chromMid_plot <- reverse_x_coord(dat_plot$chromMid)
    }

    ########################################################
    ## 4.1 绘图：每个物种每条染色体单独一页
    ########################################################

    plot(
      dat_plot$chromMid_plot,
      dat_plot$sample_mean,
      type = "l",
      col = species_color,
      xlab = ifelse(is_flipped, "pos, reversed", "pos"),
      ylab = "rho",
      main = paste0(
        species_name,
        "  Chr ",
        nChr,
        ifelse(is_flipped, "  reversed", ""),
        ifelse(
          isTRUE(smooth_curve) && plot_window_size > input_window_size,
          paste0("  plot_window=", plot_window_size / 1000, "kb"),
          paste0("  window=", input_window_size / 1000, "kb")
        )
      )
    )

    legend(
      "topright",
      legend = paste0(
        species_name,
        ifelse(is_flipped, " reversed", ""),
        ifelse(
          isTRUE(smooth_curve) && plot_window_size > input_window_size,
          paste0(" smoothed ", plot_window_size / 1000, "kb"),
          ""
        )
      ),
      col = species_color,
      lty = 1,
      bty = "n"
    )

    ########################################################
    ## 4.2 计算该染色体各样本与原始 sample_mean 的 R²
    ########################################################

    arrR2 <- calc_r2_for_samples(dat19ind, sample_cols)

    cat(
      "#Rsquare",
      species_name,
      "Chr", nChr,
      "mean_R2", mean(arrR2, na.rm = TRUE),
      "sd_R2", sd(arrR2, na.rm = TRUE),
      "n_samples", length(sample_cols),
      "n_valid_R2", sum(!is.na(arrR2)),
      "flipped", is_flipped,
      "input_window_size", input_window_size,
      "plot_window_size", plot_window_size,
      "smooth_curve", smooth_curve,
      "file", sF,
      "\n"
    )

    ########################################################
    ## 4.3 输出 rho 文件
    ########################################################

    ## 原始窗口输出：
    ## 每个物种单独输出一个 bedGraph 风格文件：
    ## chrom chromStart chromEnd sample_mean
    rho_species_dat <- dat19ind[, c("chrom", "chromStart", "chromEnd", "sample_mean")]

    append_table_no_header(
      rho_species_dat,
      rho_species_file
    )

    ## 所有物种合并输出一个带表头的原始窗口文件：
    ## species chrom chromStart chromEnd sample_mean
    rho_all_dat <- data.frame(
      species = species_name,
      chrom = dat19ind$chrom,
      chromStart = dat19ind$chromStart,
      chromEnd = dat19ind$chromEnd,
      sample_mean = dat19ind$sample_mean,
      stringsAsFactors = FALSE
    )

    append_table_with_header(
      rho_all_dat,
      rho_all_file
    )

    ## 新增：输出实际用于画图的窗口数据。
    ## 如果 smooth_curve=TRUE 且 plot_window_size > input_window_size，
    ## 这里就是平滑后的窗口数据。
    if (isTRUE(write_plot_window_rho_files)) {

      rho_plot_species_dat <- dat_plot[, c("chrom", "chromStart", "chromEnd", "sample_mean"), drop = FALSE]

      append_table_no_header(
        rho_plot_species_dat,
        rho_plot_species_file
      )

      rho_plot_all_dat <- data.frame(
        species = species_name,
        chrom = dat_plot$chrom,
        chromStart = dat_plot$chromStart,
        chromEnd = dat_plot$chromEnd,
        sample_mean = dat_plot$sample_mean,
        n_old_windows = dat_plot$n_old_windows,
        is_flipped = is_flipped,
        stringsAsFactors = FALSE
      )

      append_table_with_header(
        rho_plot_all_dat,
        rho_plot_all_file
      )
    }

    ########################################################
    ## 4.4 保存数据，供后续画图使用
    ########################################################

    key <- paste0("Chr_", nChr)

    if (is.null(all_plot_data[[key]])) {
      all_plot_data[[key]] <- list()
    }

    ## 后续所有图都使用 dat_plot。
    ## 也就是说，如果开启平滑，overlay 和 Tra-aligned 图也会使用平滑后的曲线。
    all_plot_data[[key]][[species_name]] <- list(
      species = species_name,
      chr = nChr,
      color = species_color,
      chromMid = dat_plot$chromMid,
      chromMid_plot = dat_plot$chromMid_plot,
      sample_mean = dat_plot$sample_mean,
      n_old_windows = dat_plot$n_old_windows,
      is_flipped = is_flipped
    )
  }

  cat("\n")
}

dev.off()
sink()


############################################################
## 5. 可选：按染色体编号叠加不同物种曲线
############################################################

if (make_overlay_plot) {

  pdf(file = overlay_pdf_file, width = 10, height = 5)

  for (chr_key in names(all_plot_data)) {

    chr_data <- all_plot_data[[chr_key]]

    if (length(chr_data) == 0) {
      next
    }

    ## 计算该页图的 x/y 范围
    x_range <- range(
      unlist(lapply(chr_data, function(x) x$chromMid_plot)),
      na.rm = TRUE
    )

    y_range <- range(
      unlist(lapply(chr_data, function(x) x$sample_mean)),
      na.rm = TRUE
    )

    first_sp <- chr_data[[1]]

    plot(
      first_sp$chromMid_plot,
      first_sp$sample_mean,
      type = "l",
      col = first_sp$color,
      xlab = "pos",
      ylab = "rho",
      main = paste0(
        "Overlay: ",
        chr_key,
        ifelse(
          isTRUE(smooth_curve) && plot_window_size > input_window_size,
          paste0("  smoothed ", plot_window_size / 1000, "kb"),
          ""
        )
      ),
      xlim = x_range,
      ylim = y_range
    )

    if (length(chr_data) > 1) {
      for (i in 2:length(chr_data)) {
        sp_dat <- chr_data[[i]]
        lines(
          sp_dat$chromMid_plot,
          sp_dat$sample_mean,
          col = sp_dat$color
        )
      }
    }

    legend(
      "topright",
      legend = sapply(
        chr_data,
        function(x) {
          paste0(
            x$species,
            ifelse(isTRUE(x$is_flipped), " reversed", "")
          )
        }
      ),
      col = sapply(chr_data, function(x) x$color),
      lty = 1,
      bty = "n"
    )
  }

  dev.off()
}


############################################################
## 6. 新增图：
##    以 Tra 染色体为锚点，
##    每一列是一组对应染色体，
##    每一行是一个物种，
##    所有对应关系画在同一页 PDF 中
############################################################

if (make_tra_aligned_one_page_plot) {

  n_group <- length(tra_chr_groups)
  n_species_row <- length(aligned_species_row_order)

  ## 根据列数自动调整 PDF 宽度。
  ## 如果列很多，PDF 会很宽，但仍然只有一页。
  pdf_width <- max(14, n_group * 1.7)
  pdf_height <- max(6, n_species_row * 2.0)

  pdf(
    file = tra_aligned_pdf_file,
    width = pdf_width,
    height = pdf_height,
    onefile = TRUE
  )

  old_par <- par(no.readonly = TRUE)

  par(
    mfrow = c(n_species_row, n_group),
    mar = c(1.4, 2.8, 2.0, 0.6),
    oma = c(3.5, 4.5, 3.5, 1.0),
    mgp = c(1.2, 0.3, 0),
    tcl = -0.2
  )

  ## 如果选择 row 或 global，需要预先计算 y 轴范围
  row_y_lim <- list()
  global_y_values <- c()

  if (aligned_y_scale %in% c("row", "global")) {
    for (species_name in aligned_species_row_order) {
      y_row <- c()

      for (group_info in tra_chr_groups) {
        chr_vec <- group_info[[species_name]]

        if (is.null(chr_vec)) {
          next
        }

        y_row <- c(
          y_row,
          collect_y_values_for_chr_vec(
            all_plot_data = all_plot_data,
            species_name = species_name,
            chr_vec = chr_vec
          )
        )
      }

      row_y_lim[[species_name]] <- get_range_or_null(y_row)
      global_y_values <- c(global_y_values, y_row)
    }
  }

  global_y_lim <- NULL

  if (aligned_y_scale == "global") {
    global_y_lim <- get_range_or_null(global_y_values)
  }

  ## 开始画 3 行 x N 列小图
  for (i_species in seq_along(aligned_species_row_order)) {

    species_name <- aligned_species_row_order[i_species]

    for (i_group in seq_along(tra_chr_groups)) {

      group_info <- tra_chr_groups[[i_group]]
      group_label <- group_info$group_label
      chr_vec <- group_info[[species_name]]

      if (is.null(chr_vec)) {
        chr_vec <- c()
      }

      show_group_title <- i_species == 1
      show_species_label <- i_group == 1

      y_lim_to_use <- NULL

      if (aligned_y_scale == "row") {
        y_lim_to_use <- row_y_lim[[species_name]]
      } else if (aligned_y_scale == "global") {
        y_lim_to_use <- global_y_lim
      }

      draw_one_aligned_cell(
        all_plot_data = all_plot_data,
        species_name = species_name,
        chr_vec = chr_vec,
        species_config = species_config,
        group_label = group_label,
        show_group_title = show_group_title,
        show_species_label = show_species_label,
        y_lim = y_lim_to_use,
        gap_ratio = aligned_chr_gap_ratio
      )
    }
  }

  mtext(
    "Tra-anchored corresponding chromosome groups",
    side = 3,
    outer = TRUE,
    line = 1.2,
    cex = 1.2,
    font = 2
  )

  mtext(
    paste0(
      "Each column is one Tra-anchored chromosome group; ",
      "chromosomes marked _rev are plotted in reversed x direction; ",
      ifelse(
        isTRUE(smooth_curve) && plot_window_size > input_window_size,
        paste0("curves are smoothed to ", plot_window_size / 1000, "kb windows"),
        paste0("curves use original ", input_window_size / 1000, "kb windows")
      )
    ),
    side = 1,
    outer = TRUE,
    line = 1.5,
    cex = 0.85
  )

  mtext(
    "rho",
    side = 2,
    outer = TRUE,
    line = 2.8,
    cex = 0.9
  )

  par(old_par)

  dev.off()
}


############################################################
## 7. 运行结束提示
############################################################

cat("Done!\n")
cat("Output directory:", out_dir, "\n")
cat("Main plot PDF:", pdf_file, "\n")
cat("All species raw rho file:", rho_all_file, "\n")
cat("R2 summary file:", r2_out_file, "\n")
cat("Input window size:", input_window_size, "\n")
cat("Smooth curve:", smooth_curve, "\n")
cat("Plot window size:", plot_window_size, "\n")

if (isTRUE(write_plot_window_rho_files)) {
  cat("All species plot-window rho file:", rho_plot_all_file, "\n")
}

if (make_overlay_plot) {
  cat("Overlay plot PDF:", overlay_pdf_file, "\n")
}

if (make_tra_aligned_one_page_plot) {
  cat("Tra-aligned one-page plot PDF:", tra_aligned_pdf_file, "\n")
}
