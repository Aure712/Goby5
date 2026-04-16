# ==========================
# Circos (multi-species, labels outside, configurable palettes + per-track spacing, robust save)
# ==========================

# ---- global params ----
win_size              <- 100000L   # window size (bp)
chrom_gap             <- 1         # uniform gap between all chromosomes (including last-to-first)
length_bar_height     <- 0.05      # outer length bar thickness
inner_track_height    <- 0.12      # thickness for all inner tracks (thin + uniform)
inner_track_alpha     <- 1.00      # 0..1, alpha for inner track colors (1 = opaque)

# spacing between rings (all tracks). The visible gap between two adjacent tracks ~= track_gap
track_gap             <- 0.003     # 0-1, increase to enlarge inter-track spacing

label_cex             <- 1.6       # chromosome label size (outside)
main_title_cex        <- 1.5       # species title size

# ---- palettes (customize here) ----
# at least 2 colors each; breaks will be spaced across [min, max]
palette_het           <- c("#c6dbef",  "#08519c")   # blue
palette_gene          <- c("#c7e9c0",  "#005a32")   # green
palette_cds           <- c("#fcae91",  "#a50f15")   # red
#length_bar_palette    <- c("white", "grey92")       # outer length bar alternating colors
length_bar_palette    <- c("#374151", "#D1D5DB")
#length_bar_palette    <- c("#1F2937", "#EAD7A3")

# ---- legend params (compact + vertical stack; customizable) ----
legend_cex_title      <- 0.9       # legend title text size
legend_cex_tick       <- 0.8       # legend tick text size
legend_bar_width_frac <- 0.15      # single bar width (0-1 of page width)
legend_bar_height     <- 0.06      # single bar height (0-1 of page height)
legend_vgap           <- 0.15      # vertical gap between bars (0-1 of page height)
legend_left           <- 0.22      # left margin for bars (0-1 of page width)
legend_ticks          <- 2         # number of ticks per legend (>=2)
legend_titles         <- c("Mean heterozygosity", "Gene Count", "CDS Density")

# ---- outputs ----
out_pdf               <- "circos_multi_species_side_by_side_heat_blocks.pdf"
out_png               <- "circos_multi_species_side_by_side_heat_blocks.png"
png_width_px          <- 2400
png_height_px         <- 1400
png_res               <- 220
pdf_width_in          <- 12
pdf_height_in         <- 8

set.seed(1)

# ---- species config ----
species_cfg <- list(
  list(
    id = "Tba",
    fai_file  = "Tba_sort.fa.fai",
    gff3_file = "Tba.prechange_longest_fixed.gff3",
    het_file  = "T-bar.hEst.bed",
    n_chromosomes = 22L,
    chr_prefix_to_strip = "TbaScf_"
  ),
  list(
    id = "Tbi",
    fai_file  = "Tbi_sort.fa.fai",
    gff3_file = "Tbi.prechange_longest_fixed.gff3",
    het_file  = "Tbifasciatus-muscle2.hEst.bed",
    n_chromosomes = 22L,
    chr_prefix_to_strip = "TbiScf_"
  ),
  list(
    id = "Tra",
    fai_file  = "Tra_sort.fa.fai",
    gff3_file = "Tra.prechange_longest_fixed.gff3",
    het_file  = "Tradiatus-muscle2.hEst.bed",
    n_chromosomes = 11L,
    chr_prefix_to_strip = "TraScf_"
  )
)

# ---- packages ----
pkgs <- c("data.table", "circlize")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
library(data.table)
library(circlize)

# ---- helpers ----
escape_regex <- function(s) gsub("([][{}()+*^$.|?\\\\])", "\\\\\\1", s)
simplify_chr_names <- function(x, prefix) {
  if (!nzchar(prefix)) return(as.character(x))
  sub(paste0("^", escape_regex(prefix)), "", x)
}

reduce_intervals <- function(dt) {
  if (nrow(dt) == 0L) return(dt[0])
  setorder(dt, chr, start, end)
  out <- dt[, {
    s <- start; e <- end
    cur_s <- s[1]; cur_e <- e[1]
    rs <- integer(); re <- integer()
    if (.N > 1L) {
      for (i in 2:.N) {
        if (s[i] <= cur_e + 1L) {
          if (e[i] > cur_e) cur_e <- e[i]
        } else {
          rs <- c(rs, cur_s); re <- c(re, cur_e)
          cur_s <- s[i]; cur_e <- e[i]
        }
      }
    }
    rs <- c(rs, cur_s); re <- c(re, cur_e)
    .(start = rs, end = re)
  }, by = chr]
  out[]
}

make_windows <- function(chr_len_dt, win_size) {
  rbindlist(lapply(seq_len(nrow(chr_len_dt)), function(i) {
    L   <- as.integer(chr_len_dt$length[i])
    chr <- chr_len_dt$chr[i]
    starts <- seq.int(1L, L, by = win_size)
    ends   <- pmin(starts + win_size - 1L, L)
    data.table(chr = chr, start = starts, end = ends,
               wid = seq_along(starts), wlen = ends - starts + 1L)
  }))
}

weighted_mean_in_windows <- function(x_dt, y_windows, value_col) {
  x <- copy(x_dt)
  x <- x[is.finite(x[[value_col]])]
  y <- copy(y_windows)
  if (nrow(x) == 0L) {
    z <- y[, .(chr, start, end, wid)]
    z[, wm := NA_real_]
    return(z)
  }
  setkey(x, chr, start, end)
  setkey(y, chr, start, end)
  ov <- foverlaps(x, y, nomatch = 0L)
  if (nrow(ov) == 0L) {
    z <- y[, .(chr, start, end, wid)]
    z[, wm := NA_real_]
    return(z)
  }
  ov[, ol_start := pmax(i.start, start)]
  ov[, ol_end   := pmin(i.end, end)]
  ov[, ol_len   := pmax(0L, ol_end - ol_start + 1L)]
  ov <- ov[ol_len > 0L]
  agg <- ov[, .(wm = sum(get(value_col) * ol_len) / sum(ol_len)),
            by = .(chr, start, end, wid)]
  merge(y[, .(chr, start, end, wid)], agg, by = c("chr","start","end","wid"), all.x = TRUE)[]
}

covered_bases_in_windows <- function(x_dt, y_windows) {
  y <- copy(y_windows)
  if (nrow(x_dt) == 0L) return(y[, .(chr, start, end, wid, wlen, cov_bp = 0L)])
  x <- copy(x_dt)
  setkey(x, chr, start, end)
  setkey(y, chr, start, end)
  ov <- foverlaps(x, y, nomatch = 0L)
  if (nrow(ov) == 0L) return(y[, .(chr, start, end, wid, wlen, cov_bp = 0L)])
  ov[, ol_start := pmax(i.start, start)]
  ov[, ol_end   := pmin(i.end, end)]
  ov[, ol_len   := pmax(0L, ol_end - ol_start + 1L)]
  cov_bp <- ov[ol_len > 0L, .(cov_bp = sum(ol_len)), by = .(chr, start, end, wid)]
  out <- merge(y[, .(chr, start, end, wid, wlen)], cov_bp,
               by = c("chr","start","end","wid"), all.x = TRUE)
  out[is.na(cov_bp), cov_bp := 0L][]
}

count_overlaps_in_windows <- function(x_dt, y_windows) {
  y <- copy(y_windows)
  if (nrow(x_dt) == 0L) return(y[, .(chr, start, end, wid, n = 0L)])
  x <- copy(x_dt)
  setkey(x, chr, start, end)
  setkey(y, chr, start, end)
  ov <- foverlaps(x, y, nomatch = 0L)
  if (nrow(ov) == 0L) return(y[, .(chr, start, end, wid, n = 0L)])
  cnt <- ov[, .N, by = .(chr, start, end, wid)]
  out <- merge(y[, .(chr, start, end, wid)], cnt,
               by = c("chr","start","end","wid"), all.x = TRUE)
  out[is.na(N), N := 0L]
  setnames(out, "N", "n")
  out[]
}

read_gff3_features <- function(gff3_file) {
  gff_try <- tryCatch(
    fread(gff3_file, sep = "\t", header = FALSE, fill = TRUE, quote = "", showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(gff_try) || ncol(gff_try) <= 1L) {
    gff_try <- as.data.table(read.delim(gff3_file, header = FALSE, sep = "\t",
                                        comment.char = "#", quote = "",
                                        stringsAsFactors = FALSE, fill = TRUE))
  }
  if (ncol(gff_try) == 1L) {
    lines <- readLines(gff3_file)
    lines <- lines[!grepl("^#", lines)]
    if (length(lines)) {
      parts <- strsplit(lines, "\t", fixed = TRUE)
      maxlen <- max(lengths(parts))
      mat <- t(vapply(parts, function(x) c(x, rep(NA_character_, maxlen - length(x))),
                      FUN.VALUE = character(maxlen)))
      gff_try <- as.data.table(mat)
    }
  }
  need_cols <- c("seqid","source","type","start","end","score","strand","phase","attr")
  if (ncol(gff_try) < 9L) for (j in (ncol(gff_try)+1L):9L) gff_try[[j]] <- NA_character_
  setnames(gff_try, 1:9, need_cols)
  gff <- gff_try[!(is.na(seqid) | seqid == "")]
  suppressWarnings(gff[, `:=`(start = as.integer(start), end = as.integer(end))])
  gff <- gff[!is.na(start) & !is.na(end)]
  list(
    gene = gff[type == "gene", .(chr = seqid, start, end)],
    cds  = gff[type == "CDS",  .(chr = seqid, start, end)]
  )
}

read_bed_het <- function(het_file) {
  het <- fread(het_file)
  nm <- names(het)
  chr_col  <- nm[grepl("CHROM", nm, ignore.case = TRUE)][1]
  if (is.na(chr_col)) chr_col <- nm[grepl("^chr$|^chrom$", nm, ignore.case = TRUE)][1]
  beg_col  <- nm[grepl("^BEGIN$|START|POS", nm, ignore.case = TRUE)][1]
  end_col  <- nm[grepl("^END$", nm, ignore.case = TRUE)][1]
  h_col    <- nm[grepl("^h$", nm, ignore.case = TRUE)][1]
  if (any(is.na(c(chr_col, beg_col, end_col, h_col)))) {
    stop("BED needs columns: CHROM, BEGIN/START, END, h")
  }
  dt <- het[, .(chr = as.character(get(chr_col)),
                start0 = as.integer(get(beg_col)),
                end    = as.integer(get(end_col)),
                h      = as.numeric(get(h_col)))]
  dt[, start := start0 + 1L][, start0 := NULL]
  dt[end >= start]
}

read_fai_topN <- function(fai_file, n) {
  fai <- fread(fai_file, header = FALSE,
               col.names = c("chr","length","offset","line_bases","line_width"))
  fai[, length := as.integer(length)]
  fai <- fai[order(-length)]
  fai[seq_len(min(n, .N)), .(chr, length)]
}

# Build color function from user palette (2+ colors)
make_col_fun <- function(vmin, vmax, palette) {
  if (length(palette) < 2) stop("palette must have >= 2 colors")
  if (!is.finite(vmin)) vmin <- 0
  if (!is.finite(vmax) || vmax <= vmin) vmax <- vmin + 1
  breaks <- seq(vmin, vmax, length.out = length(palette))
  circlize::colorRamp2(breaks, palette)
}

# Add alpha to hex colors
add_alpha <- function(cols, alpha = 1) {
  if (alpha >= 1) return(cols)
  rgb_mat <- col2rgb(cols, alpha = FALSE) / 255
  apply(rgb_mat, 2, function(x) rgb(x[1], x[2], x[3], alpha = alpha))
}

process_species <- function(sp, win_size) {
  fai_top <- read_fai_topN(sp$fai_file, sp$n_chromosomes)
  fai_top[, chr_s := simplify_chr_names(chr, sp$chr_prefix_to_strip)]
  
  gff <- read_gff3_features(sp$gff3_file)
  gene_dt <- gff$gene
  cds_dt  <- gff$cds
  het_dt  <- read_bed_het(sp$het_file)
  
  gene_dt[, chr := simplify_chr_names(chr, sp$chr_prefix_to_strip)]
  cds_dt[,  chr := simplify_chr_names(chr, sp$chr_prefix_to_strip)]
  het_dt[,  chr := simplify_chr_names(chr, sp$chr_prefix_to_strip)]
  
  gene_dt <- gene_dt[chr %in% fai_top$chr_s]
  cds_dt  <- cds_dt[ chr %in% fai_top$chr_s]
  het_dt  <- het_dt[ chr %in% fai_top$chr_s]
  
  kary <- fai_top[, .(chr = chr_s, length = as.integer(length))]
  windows <- make_windows(kary, win_size)
  
  het_wm <- weighted_mean_in_windows(het_dt[, .(chr, start, end, h)], windows, value_col = "h")
  setnames(het_wm, "wm", "het_mean")
  
  gene_cnt <- count_overlaps_in_windows(gene_dt, windows)
  setnames(gene_cnt, "n", "gene_count")
  
  cds_merged <- reduce_intervals(cds_dt)
  cds_cov <- covered_bases_in_windows(cds_merged, windows)
  cds_cov[, cds_density := pmin(1, cov_bp / wlen)]
  
  agg <- Reduce(function(a, b) merge(a, b, by = c("chr","start","end","wid"), all = TRUE),
                list(het_wm,
                     gene_cnt,
                     cds_cov[, .(chr, start, end, wid, cds_density)]))
  agg[is.na(het_mean), het_mean := NA_real_]
  agg[is.na(gene_count), gene_count := 0L]
  agg[is.na(cds_density), cds_density := 0]
  list(kary = kary, agg = agg, id = sp$id)
}

# ---- run per species ----
res_list <- lapply(species_cfg, process_species, win_size = win_size)

# ---- global min/max across species (for color mapping) ----
agg_all <- rbindlist(lapply(res_list, function(x) x$agg), fill = TRUE)
min_het  <- suppressWarnings(min(agg_all$het_mean,  na.rm = TRUE)); if (!is.finite(min_het)) min_het  <- 0
max_het  <- suppressWarnings(max(agg_all$het_mean,  na.rm = TRUE)); if (!is.finite(max_het)  || max_het  <= min_het) max_het  <- min_het + 1
min_gene <- suppressWarnings(min(agg_all$gene_count, na.rm = TRUE)); if (!is.finite(min_gene)) min_gene <- 0
max_gene <- suppressWarnings(max(agg_all$gene_count, na.rm = TRUE)); if (!is.finite(max_gene) || max_gene <= min_gene) max_gene <- min_gene + 1
min_cds  <- suppressWarnings(min(agg_all$cds_density,na.rm = TRUE)); if (!is.finite(min_cds)) min_cds  <- 0
max_cds  <- suppressWarnings(max(agg_all$cds_density,na.rm = TRUE)); if (!is.finite(max_cds)  || max_cds  <= min_cds) max_cds  <- min_cds + 1e-6

# Build color functions from user palettes
col_fun_het  <- make_col_fun(min_het,  max_het,  palette_het)
col_fun_gene <- make_col_fun(min_gene, max_gene, palette_gene)
col_fun_cds  <- make_col_fun(min_cds,  max_cds,  palette_cds)

# ---- compact stacked legends (base graphics) ----
draw_colorbar_h <- function(x0, y0, x1, y1, col_fun, min_v, max_v, title,
                            cex_title = 1, cex_tick = 0.9, n_breaks = 3) {
  op <- par(xpd = NA)
  on.exit(par(op))
  n <- 240
  vals <- seq(min_v, max_v, length.out = n)
  cols <- col_fun(vals)
  xs <- seq(x0, x1, length.out = n + 1L)
  for (i in 1:n) rect(xs[i], y0, xs[i+1], y1, col = cols[i], border = NA)
  rect(x0, y0, x1, y1, border = "grey30", lwd = 0.7)
  text(x0, y1 + 0.045, labels = title, adj = c(0, 0), cex = cex_title, font = 2)
  if (n_breaks < 2) n_breaks <- 2
  at_vals <- seq(min_v, max_v, length.out = n_breaks)
  at_pos  <- x0 + (at_vals - min_v) / (max_v - min_v) * (x1 - x0)
  segments(at_pos, y0, at_pos, y0 - 0.012, col = "grey20")
  text(at_pos, y0 - 0.032, labels = format(at_vals, digits = 3), cex = cex_tick)
}

render_all <- function() {
  nsp <- length(res_list)
  lay_mat <- rbind(1:nsp, rep(nsp + 1L, nsp))  # top: plots, bottom: legends
  layout(lay_mat, heights = c(7, 3))
  par(oma = c(1, 1, 1, 1))
  
  # top row: per-species circos
  for (i in seq_len(nsp)) {
    par(mar = c(1, 1, 3, 1))
    sp_id <- res_list[[i]]$id
    kary  <- res_list[[i]]$kary
    agg   <- res_list[[i]]$agg
    
    kary <- kary[order(-length)]
    gap_vec <- rep(chrom_gap, nrow(kary))  # uniform sector gaps (including last->first)
    
    circos.clear()
    circos.par(
      gap.after = gap_vec,
      start.degree = 90,
      track.height = 0.12,
      points.overflow.warning = FALSE
    )
    
    # labels outside
    karyo_df <- kary[, .(chr, start = 1L, end = as.integer(length))]
    circos.genomicInitialize(karyo_df, plotType = "labels", labels.cex = label_cex)
    
    # unify spacing per track via explicit track.margin
    tmargin <- c(track_gap/2, track_gap/2)
    
    # outer length bar (alternating colors) — remove border to avoid visual "narrower gap"
    length_bar_cols <- rep(length_bar_palette, length.out = nrow(kary))
    circos.track(
      ylim = c(0, 1),
      track.height = length_bar_height,
      bg.col = length_bar_cols,
      bg.border = NA,                 # was "#777777"; remove to make gaps visually equal
      track.margin = tmargin,
      panel.fun = function(x, y) { }
    )
    
    # Het (blocks)
    het_df <- agg[, .(chr, start, end, value = het_mean)]
    het_df[is.na(value), value := min_het]
    circos.genomicTrack(het_df, ylim = c(0, 1), bg.border = NA,
                        track.height = inner_track_height,
                        track.margin = tmargin,
                        panel.fun = function(region, value, ...) {
                          v <- as.numeric(value[[1]])
                          cols <- add_alpha(col_fun_het(v), inner_track_alpha)
                          circos.genomicRect(region, ybottom = 0, ytop = 1, col = cols, border = NA)
                        })
    
    # Gene Count (blocks)
    gene_df <- agg[, .(chr, start, end, value = gene_count)]
    gene_df[is.na(value), value := min_gene]
    circos.genomicTrack(gene_df, ylim = c(0, 1), bg.border = NA,
                        track.height = inner_track_height,
                        track.margin = tmargin,
                        panel.fun = function(region, value, ...) {
                          v <- as.numeric(value[[1]])
                          cols <- add_alpha(col_fun_gene(v), inner_track_alpha)
                          circos.genomicRect(region, ybottom = 0, ytop = 1, col = cols, border = NA)
                        })
    
    # CDS density (blocks)
    cds_df <- agg[, .(chr, start, end, value = cds_density)]
    cds_df[is.na(value), value := min_cds]
    circos.genomicTrack(cds_df, ylim = c(0, 1), bg.border = NA,
                        track.height = inner_track_height,
                        track.margin = tmargin,
                        panel.fun = function(region, value, ...) {
                          v <- as.numeric(value[[1]])
                          cols <- add_alpha(col_fun_cds(v), inner_track_alpha)
                          circos.genomicRect(region, ybottom = 0, ytop = 1, col = cols, border = NA)
                        })
    
    title(main = sp_id, cex.main = main_title_cex)
  }
  
  # bottom row: stacked compact legends
  par(mar = c(2, 2, 1, 2))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  
  left  <- legend_left
  right <- min(0.98, left + legend_bar_width_frac)
  h     <- legend_bar_height
  vgap  <- legend_vgap
  
  y_top    <- 0.78
  y_mid    <- y_top - (h + vgap)
  y_bottom <- y_mid - (h + vgap)
  
  draw_colorbar_h(left, y_top,    right, y_top    + h, col_fun_het,  min_het,  max_het,
                  title = legend_titles[1], cex_title = legend_cex_title, cex_tick = legend_cex_tick, n_breaks = legend_ticks)
  draw_colorbar_h(left, y_mid,    right, y_mid    + h, col_fun_gene, min_gene, max_gene,
                  title = legend_titles[2], cex_title = legend_cex_title, cex_tick = legend_cex_tick, n_breaks = legend_ticks)
  draw_colorbar_h(left, y_bottom, right, y_bottom + h, col_fun_cds,  min_cds,  max_cds,
                  title = legend_titles[3], cex_title = legend_cex_title, cex_tick = legend_cex_tick, n_breaks = legend_ticks)
}

# ---- robust save: ensure device closes even if rendering errors ----
safe_render <- function(dev_open, dev_close, file) {
  dev_open()
  ok <- FALSE
  tryCatch({
    render_all()
    ok <- TRUE
  }, error = function(e) {
    message("Error while rendering: ", conditionMessage(e))
  }, finally = {
    dev_close()
    if (!ok) message("The output file might be incomplete: ", file)
  })
}

# PDF
safe_render(
  dev_open  = function() pdf(out_pdf, width = pdf_width_in, height = pdf_height_in),
  dev_close = function() dev.off(),
  file = out_pdf
)

# PNG
safe_render(
  dev_open  = function() png(out_png, width = png_width_px, height = png_height_px, res = png_res),
  dev_close = function() dev.off(),
  file = out_png
)
