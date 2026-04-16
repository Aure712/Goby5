# Repro settings
Sys.setenv(OMP_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1")
if ("sample.kind" %in% names(formals(RNGkind))) {
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rounding")
} else {
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion")
}
options(stringsAsFactors = FALSE)

setwd("/public4/group_crf/home/g21shaoy23/cactus/phast_tra/")

library(rphast)
library(ape)

# --------------------- CLI args ---------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript rphastcon_24species_all_exceptT11.R <chr_id> [speciesSet|'all'|comma,list] [key=val ...]
  keys: seed=20250312 nrep=100000 split_len=50 site_q=0.95 site_min=1 site_nrep=50000 compat=0
        out_dir=./results_24species_exceptT2_tree22_test2bar2b cache=0 rebuild=0 cache_dir=<out_dir>/cache
        per_site=1 site_top_k=5 plot_tree=1 mark_rule=thr strict_present=off|ACGT|ACGTN strict_species=all enforce_strict=1
        msa_offset=0.03 msa_width=0.60 max_once_cols=1000000
        export_msa=1 msa_dir=<out_dir>/msa_exports msa_format=PHYLIP acc_msa_set=validated|sig
        win_fdr=0.05 win_q=NA win_p=NA")
}

nChr <- args[1]
sChr <- paste("TraScf_", nChr, sep = "")
sRefSp <- "Tridentiger_radiatus"
speciesArg <- if (length(args) >= 2) args[2] else "all"

# --------------------- Defaults ---------------------
seed_main <- 20250312L
nrep <- 100000L
splitLength <- 50L
site_q <- 0.95
site_min <- 1L
site_nrep <- 50000L
compat <- 0L

# 新增：统一输出目录（默认保持旧目录）
out_dir <- "./results_24species_final5_bar2_scalesame2_export2"
use_cache <- 0L
rebuild_cache <- 0L
cache_dir <- file.path(out_dir, "cache")
cache_dir_is_user <- FALSE

per_site <- 1L
site_top_k <- 5L
plot_tree <- 1L
mark_rule <- "thr"
strict_present <- "ACGT"
strict_species <- "all"
enforce_strict <- 1L
msa_offset <- 0.03
msa_width  <- 0.60
max_once_cols <- 1000000L

# 新增：MSA 导出参数
export_msa <- 1L
msa_dir <- file.path(out_dir, "msa_exports")
msa_dir_is_user <- FALSE
msa_format <- "PHYLIP"   # IQ-TREE 可直接用 "FASTA"
acc_msa_set <- "validated"  # "validated" or "sig"

# 新增：窗口级筛选阈值
win_fdr <- 0.05
win_q <- 0.95 #LRT > null 的分位 win_q <- 0.99
win_p <- 0.05 #右尾 p < 0.001 win_p <- 1e-3

if (length(args) >= 3) {
  for (kv in args[3:length(args)]) {
    if (!grepl("=", kv)) next
    key <- sub("=.*$", "", kv)
    val <- sub("^[^=]*=", "", kv)
    if (key == "seed")           { seed_main <- as.integer(val) }
    if (key == "nrep")           { nrep <- as.integer(val) }
    if (key == "split_len")      { splitLength <- as.integer(val) }
    if (key == "site_q")         { site_q <- as.numeric(val) }
    if (key == "site_min")       { site_min <- as.integer(val) }
    if (key == "site_nrep")      { site_nrep <- as.integer(val) }
    if (key == "compat")         { compat <- as.integer(val) }
    if (key == "out_dir")        {
      out_dir <- val
      if (!cache_dir_is_user) cache_dir <- file.path(out_dir, "cache")
      if (!msa_dir_is_user) msa_dir <- file.path(out_dir, "msa_exports")
    }
    if (key == "cache")          { use_cache <- as.integer(val) }
    if (key == "rebuild")        { rebuild_cache <- as.integer(val) }
    if (key == "cache_dir")      { cache_dir <- val; cache_dir_is_user <- TRUE }
    if (key == "per_site")       { per_site <- as.integer(val) }
    if (key == "site_top_k")     { site_top_k <- as.integer(val) }
    if (key == "plot_tree")      { plot_tree <- as.integer(val) }
    if (key == "mark_rule")      { mark_rule <- tolower(val) }
    if (key == "strict_present") { strict_present <- toupper(val) }
    if (key == "strict_species") { strict_species <- val }
    if (key == "enforce_strict") { enforce_strict <- as.integer(val) }
    if (key == "msa_offset")     { msa_offset <- as.numeric(val) }
    if (key == "msa_width")      { msa_width  <- as.numeric(val) }
    if (key == "max_once_cols")  { max_once_cols <- as.integer(val) }

    # 新增：MSA 导出
    if (key == "export_msa")     { export_msa <- as.integer(val) }
    if (key == "msa_dir")        { msa_dir <- val; msa_dir_is_user <- TRUE }
    if (key == "msa_format")     { msa_format <- toupper(val) }
    if (key == "acc_msa_set")    { acc_msa_set <- tolower(val) }

    # 新增：窗口级筛选
    if (key == "win_fdr")        { win_fdr <- as.numeric(val) }
    if (key == "win_q")          { win_q <- as.numeric(val) }
    if (key == "win_p")          { win_p <- as.numeric(val) }
  }
}

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed_main)
cat(sprintf("Params: speciesArg=%s, seed=%d, nrep=%d, split_len=%d, site_q=%.3f, site_min=%d, site_nrep=%d, compat=%d, out_dir=%s, cache=%d, rebuild=%d, cache_dir=%s, per_site=%d, site_top_k=%d, plot_tree=%d, mark_rule=%s, strict_present=%s, strict_species=%s, enforce_strict=%d, msa_offset=%.3f, msa_width=%.3f, max_once_cols=%d, export_msa=%d, msa_dir=%s, msa_format=%s, acc_msa_set=%s, win_fdr=%.4g, win_q=%.4g, win_p=%.4g\n",
            speciesArg, seed_main, nrep, splitLength, site_q, site_min, site_nrep, compat, out_dir, use_cache, rebuild_cache, cache_dir, per_site, site_top_k, plot_tree, mark_rule, strict_present, strict_species, enforce_strict, msa_offset, msa_width, max_once_cols, export_msa, msa_dir, msa_format, acc_msa_set, win_fdr, win_q, win_p))
if (use_cache == 1L && !dir.exists(cache_dir)) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
}
if (export_msa == 1L && !dir.exists(msa_dir)) {
  dir.create(msa_dir, recursive = TRUE, showWarnings = FALSE)
}

# --------------------- Constants & species groups ---------------------
win_null_chunk <- 5000L
site_null_chunk <- 20000L
site_obs_chunk  <- 10000L

tree <- "(((((Odontamblyopus_rebecca,(Taenioides_cirratus,Trypauchen_vagina)),Periophthalmus_magnuspinnatus),(Oligolepis_acutipinnis,(Stiphodon_pelewensis,Sicyopterus_lagocephalus))),((Eugnathogobius_siamensis,Redigobius_balteatus),(Chaeturichthys_stigmatias,((Rhinogobius_formosanus,Rhinogobius_similis),(Tridentiger_bifasciatus,(Tridentiger_radiatus,Tridentiger_barbatus)BAR))))),(Gobiodon_okinawae,(Neogobius_melanostomus,((Glossogobius_giuris,Bathygobius_meggitti),(Cryptocentrus_pavoninoides,(Exyrias_puntang,(Gobiopsis_macrostoma,(Amoya_brevirostris,Acentrogobius_pflaumii))))))));"

arrFreshs <- c('Tridentiger_bifasciatus','Rhinogobius_formosanus','Rhinogobius_similis','Redigobius_balteatus','Odontamblyopus_rebecca','Trypauchen_vagina','Periophthalmus_magnuspinnatus','Oligolepis_acutipinnis','Stiphodon_pelewensis','Sicyopterus_lagocephalus','Gobiodon_okinawae','Neogobius_melanostomus','Glossogobius_giuris','Bathygobius_meggitti','Cryptocentrus_pavoninoides','Exyrias_puntang','Amoya_brevirostris','Acentrogobius_pflaumii')
arrNonFreshs <- c('Tridentiger_radiatus','Tridentiger_barbatus')

# --------------------- Helpers (defined before use) ---------------------
empirical.pval <- function(x, dist) mean(dist >= x)

with_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(seed)
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old, envir = .GlobalEnv)
    }
  }, add = TRUE)
  eval.parent(substitute(expr))
}

sanitize <- function(x) gsub("[^A-Za-z0-9_+.-]", "_", x)

get_sub_msa_fun <- function() {
  if (!exists("sub.msa", where = asNamespace("rphast"), inherits = FALSE)) {
    stop("rphast::sub.msa not found. Please upgrade 'rphast'.")
  }
  get("sub.msa", envir = asNamespace("rphast"))
}

subset_msa_by_species <- function(aln_full, selSp) {
  subfun <- get_sub_msa_fun()
  subfun(aln_full, seqs = selSp)
}

write_gff10 <- function(feat_df, extra_vec, out_file, override_seqname = NULL) {
  has <- function(n) n %in% colnames(feat_df)
  seqname <- if (!is.null(override_seqname)) {
    rep(override_seqname, nrow(feat_df))
  } else if (has("seqname")) {
    feat_df$seqname
  } else {
    rep(".", nrow(feat_df))
  }
  source  <- if (has("source")) {
    feat_df$source
  } else if (has("src")) {
    feat_df$src
  } else {
    rep(".", nrow(feat_df))
  }
  feature <- if (has("feature")) {
    feat_df$feature
  } else if (has("feat")) {
    feat_df$feat
  } else {
    rep(".", nrow(feat_df))
  }
  start   <- if (has("start"))  feat_df$start else stop("feat_df missing start")
  end     <- if (has("end"))    feat_df$end   else stop("feat_df missing end")
  score   <- if (has("score"))  feat_df$score else rep(NA_real_, nrow(feat_df))
  strand  <- if (has("strand")) feat_df$strand else rep(".", nrow(feat_df))
  frame   <- if (has("frame"))  feat_df$frame  else rep(".", nrow(feat_df))
  group   <- if (has("group"))  feat_df$group  else rep(".", nrow(feat_df))
  df <- data.frame(seqname=seqname, source=source, feature=feature,
                   start=start, end=end, score=score, strand=strand,
                   frame=frame, group=group, info=extra_vec,
                   stringsAsFactors = FALSE)
  write.table(df, file = out_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

compute_split_parent_idx <- function(parent_feats, split_len) {
  if (!all(c("start","end") %in% colnames(parent_feats))) stop("parent_feats must have start/end")
  lens <- as.integer(parent_feats$end) - as.integer(parent_feats$start) + 1L
  nwin <- pmax(0L, lens %/% as.integer(split_len))
  rep(seq_len(nrow(parent_feats)), times = nwin)
}

read_fasta_to_named_vec <- function(path) {
  if (!file.exists(path)) return(character(0))
  lines <- readLines(path)
  hdr_idx <- which(substr(lines, 1, 1) == ">")
  if (length(hdr_idx) == 0) return(character(0))
  n <- length(hdr_idx)
  res_names <- character(n); res_seqs <- character(n)
  for (i in seq_len(n)) {
    start <- hdr_idx[i]
    end <- if (i < n) hdr_idx[i+1] - 1L else length(lines)
    name <- sub("^>", "", lines[start])
    seq_lines <- lines[(start+1L):end]
    seq_lines <- seq_lines[seq_lines != ""]
    res_names[i] <- name
    res_seqs[i]  <- paste(seq_lines, collapse = "")
  }
  names(res_seqs) <- res_names
  res_seqs
}

write_fasta_from_named_vec <- function(seqs_named, file, width = 60L) {
  con <- file(file, "w"); on.exit(close(con))
  nms <- names(seqs_named)
  if (length(nms) == 0) {
    writeLines("# empty window: no sequences extracted", con)
    return(invisible(NULL))
  }
  for (i in seq_along(seqs_named)) {
    writeLines(paste0(">", nms[i]), con)
    s <- seqs_named[i]
    L <- nchar(s, type = "chars")
    if (L == 0) { writeLines("", con); next }
    starts <- seq(1L, L, by = width)
    for (st in starts) {
      ed <- min(st + width - 1L, L)
      writeLines(substr(s, st, ed), con)
    }
  }
}

annotate_seq_by_positions <- function(seq_str, top_pos) {
  if (length(top_pos) == 0) return(tolower(seq_str))
  chars <- strsplit(seq_str, "", fixed = TRUE)[[1]]
  w <- length(chars)
  is_top <- rep(FALSE, w)
  top_pos <- top_pos[top_pos >= 1 & top_pos <= w]
  if (length(top_pos) > 0) is_top[top_pos] <- TRUE
  for (k in seq_len(w)) {
    ch <- chars[k]
    if (ch == "-" || ch == "") next
    if (is_top[k]) chars[k] <- toupper(ch) else chars[k] <- tolower(ch)
  }
  paste(chars, collapse = "")
}

find_seq_header_name <- function(headers, target_name) {
  if (length(headers) == 0) return(NA_character_)
  idx <- which(headers == target_name)
  if (length(idx) >= 1) return(headers[idx[1]])
  tnlen <- nchar(target_name, type = "chars")
  starts <- substr(headers, 1L, tnlen) == target_name
  if (any(starts)) {
    cand <- headers[starts]
    nextchar <- ifelse(nchar(cand) > tnlen, substring(cand, tnlen + 1L, tnlen + 1L), "")
    ok <- (nextchar == "") | (nextchar %in% c(".", ":", ";", "_", " ", "|"))
    cand_ok <- if (any(ok)) cand[ok] else cand
    cand_ok <- cand_ok[order(nchar(cand_ok))]
    return(cand_ok[1])
  }
  contains <- grepl(target_name, headers, fixed = TRUE)
  if (any(contains)) {
    cand <- headers[contains]
    cand <- cand[order(nchar(cand))]
    return(cand[1])
  }
  NA_character_
}

# -------- phyloP branch-only ACC wrapper --------
# Prefer 'branches='; else 'branch=' + subtree=FALSE; last resort falls back to 'branch=' (may default to subtree in older builds).
phyloP_acc_on_single_branch <- function(mod, msa, features, target_label) {
  f <- names(formals(phyloP))
  if ("branches" %in% f) {
    return(phyloP(mod, msa = msa, mode = "ACC", features = features, branches = target_label))
  } else if (all(c("branch","subtree") %in% f)) {
    return(phyloP(mod, msa = msa, mode = "ACC", features = features, branch = target_label, subtree = FALSE))
  } else {
    message("phyloP: this build cannot force single-branch; fallback to subtree for label=", target_label)
    return(phyloP(mod, msa = msa, mode = "ACC", features = features, branch = target_label))
  }
}

# -------- MSA 导出（保守/加速） --------
write_concat_msa <- function(aln_all, feats, out_file, fmt = "FASTA") {
  if (is.null(feats) || nrow(feats) == 0) {
    con <- file(out_file, "w"); writeLines("# empty MSA: no features", con); close(con)
    return(FALSE)
  }
  msa <- NULL
  try(msa <- extract.feature.msa(copy.msa(aln_all), feats, pointer.only = FALSE), silent = TRUE)
  if (is.null(msa) || ncol.msa(msa) <= 0 || length(names.msa(msa)) == 0) {
    con <- file(out_file, "w"); writeLines("# empty MSA: extract.feature.msa failed", con); close(con)
    return(FALSE)
  }
  write.msa(msa, file = out_file, format = fmt)
  TRUE
}

# -------- msaplot 支持（不叠加竖直红线） --------
has_TDbook <- requireNamespace("TDbook", quietly = TRUE)
has_ggtree <- requireNamespace("ggtree", quietly = TRUE)
has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)

build_msaplot_fasta <- function(phy, seqs_named, out_fa_path) {
  tips <- phy$tip.label
  hdrs <- names(seqs_named)
  out_names <- character(0); out_seqs <- character(0)
  for (sp in tips) {
    h <- find_seq_header_name(hdrs, sp)
    if (!is.na(h)) {
      out_names <- c(out_names, sp)
      out_seqs <- c(out_seqs, toupper(seqs_named[[h]]))
    }
  }
  if (length(out_names) == 0) return(FALSE)
  names(out_seqs) <- out_names
  write_fasta_from_named_vec(out_seqs, out_fa_path, width = 60L)
  TRUE
}

plot_tree_and_alignment_msaplot <- function(phy, msafasta, relpos_over_thr,
                                            outfile_png, title_txt = "",
                                            msa_offset = 0.03, msa_width = 0.60) {
  if (!(has_ggtree && has_ggplot2)) return(FALSE)
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggtree))
  seqs_named <- read_fasta_to_named_vec(msafasta)
  if (length(seqs_named) == 0) return(FALSE)
  Ls <- nchar(seqs_named, type = "chars")
  if (length(unique(Ls)) != 1) return(FALSE)
  p_tree <- ggtree(phy, size = 0.4) + geom_tiplab(size = 2.6, align = FALSE)
  msaplot_fn <- NULL
  if (has_TDbook && "msaplot" %in% getNamespaceExports("TDbook")) {
    msaplot_fn <- get("msaplot", asNamespace("TDbook"))
  } else if ("msaplot" %in% getNamespaceExports("ggtree")) {
    msaplot_fn <- get("msaplot", asNamespace("ggtree"))
  } else {
    return(FALSE)
  }
  p_msa <- NULL; ok <- FALSE
  try({ p_msa <- msaplot_fn(p_tree, fasta = msafasta, offset = msa_offset, width = msa_width); ok <- TRUE }, silent = TRUE)
  if (!ok || is.null(p_msa)) return(FALSE)
  p_msa <- p_msa + ggtitle(title_txt)
  ok2 <- TRUE
  tryCatch({
    ggsave(outfile_png, plot = p_msa, width = 10, height = 8, dpi = 150, limitsize = FALSE)
  }, error = function(e) { message("ggsave failed: ", conditionMessage(e)); ok2 <<- FALSE })
  ok2
}

# 两棵树并排绘制（统一横轴尺度）
plot_two_trees_side_by_side <- function(cons_tree, rar_tree, outfile_png) {
  to_phy <- function(obj) {
    if (inherits(obj, "phylo")) return(obj)
    if (is.character(obj)) return(ape::read.tree(text = obj))
    if (!is.null(obj$tree)) {
      if (inherits(obj$tree, "phylo")) return(obj$tree)
      if (is.character(obj$tree)) return(ape::read.tree(text = obj$tree))
    }
    stop("Invalid tree object.")
  }
  phy1 <- to_phy(cons_tree)
  phy2 <- to_phy(rar_tree)
  depth <- function(phy) max(ape::node.depth.edgelength(phy))
  xmax <- max(depth(phy1), depth(phy2))
  if (!is.finite(xmax) || xmax <= 0) xmax <- 1
  png(outfile_png, width = 1200, height = 600, res = 150)
  op <- par(no.readonly = TRUE); on.exit({ par(op); dev.off() }, add = TRUE)
  par(mfrow = c(1,2), mar = c(4,2,3,1))
  plot.phylo(phy1, x.lim = c(0, xmax * 1.05), main = "Conserved Elements", cex = 0.8)
  plot.phylo(phy2, x.lim = c(0, xmax * 1.05), main = "Accelerated Elements", cex = 0.8)
  invisible(TRUE)
}

# 提取窗口为文件（raw）——完全基于 feature
extract_window_to_files <- function(alnAll, s, e, out_raw_fa, ref_sp_for_feature) {
  nc <- ncol.msa(alnAll)
  s0 <- max(1L, min(as.integer(s), nc))
  e0 <- max(1L, min(as.integer(e), nc))
  if (e0 < s0) {
    con <- file(out_raw_fa, "w"); writeLines("# empty window: end<start", con); close(con)
    return(character(0))
  }
  anchor <- if (ref_sp_for_feature %in% names.msa(alnAll)) ref_sp_for_feature else names.msa(alnAll)[1]
  feat_win <- feat(seqname = anchor, src = "accwin", feat = ".", start = s0, end = e0)
  msa1 <- NULL
  try(msa1 <- extract.feature.msa(copy.msa(alnAll), feat_win, pointer.only = FALSE), silent = TRUE)
  wrote <- FALSE
  if (!is.null(msa1) && ncol.msa(msa1) > 0 && length(names.msa(msa1)) > 0) {
    try(write.msa(msa1, file = out_raw_fa, format = "FASTA"), silent = TRUE)
    wrote <- file.exists(out_raw_fa) && file.size(out_raw_fa) > 0
  }
  if (!wrote) {
    con <- file(out_raw_fa, "w"); writeLines("# empty window: no sequences extracted by extract.feature.msa", con); close(con)
    return(character(0))
  }
  read_fasta_to_named_vec(out_raw_fa)
}

# 提取窗口序列（不落盘）——完全基于 feature
extract_window_seqs_only <- function(alnAll, s, e, ref_sp_for_feature) {
  nc <- ncol.msa(alnAll)
  s0 <- max(1L, min(as.integer(s), nc))
  e0 <- max(1L, min(as.integer(e), nc))
  if (e0 < s0) return(character(0))
  anchor <- if (ref_sp_for_feature %in% names.msa(alnAll)) ref_sp_for_feature else names.msa(alnAll)[1]
  feat_win <- feat(seqname = anchor, src = "accwin", feat = ".", start = s0, end = e0)
  msa1 <- NULL
  try(msa1 <- extract.feature.msa(copy.msa(alnAll), feat_win, pointer.only = FALSE), silent = TRUE)
  if (is.null(msa1) || ncol.msa(msa1) <= 0 || length(names.msa(msa1)) == 0) return(character(0))
  tmpfa <- tempfile(fileext = ".fa"); on.exit(unlink(tmpfa), add = TRUE)
  write.msa(msa1, file = tmpfa, format = "FASTA")
  read_fasta_to_named_vec(tmpfa)
}

# 严格在场与正则
parse_strict_species <- function(arg, msa_names, arrFreshs) {
  al <- tolower(arg)
  if (al %in% c("off","","none")) return(character(0))
  if (al %in% c("arrfreshs","fresh","freshs")) return(arrFreshs)
  if (al %in% c("all","*")) return(msa_names)
  trimws(unlist(strsplit(arg, ",")))
}
allowed_regex_for <- function(present) {
  if (present == "ACGT")  return("[ACGT]")
  if (present == "ACGTN") return("[ACGTN]")
  return(NA_character_)
}
has_any_allowed_base <- function(seq_str, regex_pat) {
  if (is.na(regex_pat)) return(TRUE)
  grepl(regex_pat, toupper(seq_str))
}

# --------------------- Species set ---------------------
consSpecies <- NULL
la <- tolower(speciesArg)
if (la %in% c("all", "", "*")) {
  consSpecies <- NULL
} else if (la %in% c("arrfreshs", "fresh", "freshs")) {
  consSpecies <- arrFreshs
} else if (la %in% c("arrnonfreshs", "nonfresh", "marine")) {
  consSpecies <- arrNonFreshs
} else {
  consSpecies <- trimws(unlist(strsplit(speciesArg, ",")))
}

species_label <- sanitize(tolower(speciesArg))
make_cache_file <- function(kind, sChr, species_label, splitLength, nrep, cache_dir) {
  if (kind == "win") {
    file.path(cache_dir, sprintf("cache_null_win_%s_%s_split%d_nrep%d.rds", sChr, species_label, splitLength, nrep))
  } else if (kind == "site") {
    file.path(cache_dir, sprintf("cache_null_site_%s_%s_nrep%d.rds", sChr, species_label, nrep))
  } else {
    stop("unknown cache kind")
  }
}

# --------------------- Data load & neutral model ---------------------
feats <- read.feat(paste("./cds_gff/", sChr, ".cds.gff", sep = ""))
feats$seqname <- sRefSp

maf_file <- paste("./maf/", sChr, ".maf", sep = "")
alnAll <- read.msa(filename = maf_file, format = "MAF", ordered = TRUE)
aln4dAll <- get4d.msa(alnAll, feats)
neutralModAll <- phyloFit(aln4dAll, tree = tree, subst.mod = "REV")

msa_names <- names.msa(alnAll); print(msa_names)
tree_tips <- if (inherits(neutralModAll$tree, "phylo")) neutralModAll$tree$tip.label else ape::read.tree(text = neutralModAll$tree)$tip.label
print(tree_tips)
cat("In tree not in MSA: ", paste(setdiff(tree_tips, msa_names), collapse = ", "), "\n", sep = "")
cat("In MSA not in tree: ", paste(setdiff(msa_names, tree_tips), collapse = ", "), "\n", sep = "")

# --------------------- phastCons subset ---------------------
if (is.null(consSpecies)) {
  selSp <- msa_names
  alnCons <- alnAll
  cat(sprintf("Conservation species set: all (%d species)\n", length(selSp)))
} else {
  selSp <- intersect(msa_names, consSpecies)
  cat("Requested species: ", paste(consSpecies, collapse = ", "), "\n", sep = "")
  cat("Selected species in MSA: ", paste(selSp, collapse = ", "), "\n", sep = "")
  if (length(selSp) < 2) stop("Selected species not found sufficiently in the MSA; need >= 2 after intersection.")
  alnCons <- subset_msa_by_species(alnAll, selSp)
  cat("Conservation species set (subset): ", paste(names.msa(alnCons), collapse = ", "), "\n", sep = "")
}

tips_full <- if (inherits(neutralModAll$tree, "phylo")) neutralModAll$tree$tip.label else ape::read.tree(text = neutralModAll$tree)$tip.label
common_sp <- intersect(names.msa(alnCons), tips_full)
if (length(common_sp) < 2) stop("Model tree and alignment share fewer than 2 species; cannot run phastCons.")

pcEM_Cons <- phastCons(alnCons, neutralModAll,
                       viterbi = TRUE,
                       expected.length = 45,
                       rho = 0.3,
                       target.coverage = 0.3)

totSpCons <- length(names.msa(alnCons))
hasSelSpecies <- informative.regions.msa(alnCons, totSpCons)

informativeElements <- coverage.feat(pcEM_Cons[["most.conserved"]], hasSelSpecies, get.feats = TRUE)
write.feat(informativeElements, file.path(out_dir, paste("informativeElements_", sChr, ".gff", sep = "")))

# -------- 导出保守 MSA --------
if (export_msa == 1L) {
  out_cons_msa <- file.path(msa_dir, paste("msa_conserved_", sChr, ".", tolower(msa_format), sep = ""))
  ok_cons <- write_concat_msa(alnAll, informativeElements, out_cons_msa, fmt = msa_format)
  if (ok_cons) {
    cat("Wrote conserved MSA: ", out_cons_msa, "\n", sep = "")
  } else {
    cat("Failed to write conserved MSA: ", out_cons_msa, "\n", sep = "")
  }
}

# --------------------- Stats ---------------------
total_len_all <- ncol.msa(alnAll)
conserved_len_all <- coverage.feat(pcEM_Cons[["most.conserved"]])
prop_all <- conserved_len_all / total_len_all

total_len_inf <- coverage.feat(hasSelSpecies)
conserved_len_inf <- coverage.feat(informativeElements)
prop_inf <- conserved_len_inf / total_len_inf

cat(sprintf("[%-s] Conserved/All columns: %d/%d = %.6f\n", sChr, conserved_len_all, total_len_all, prop_all))
cat(sprintf("[%-s] Conserved/All-species-informative: %d/%d = %.6f\n", sChr, conserved_len_inf, total_len_inf, prop_inf))

stats_df <- data.frame(
  seq = sChr,
  measure = c("all_columns", "all_species_informative"),
  conserved_len = c(conserved_len_all, conserved_len_inf),
  total_len = c(total_len_all, total_len_inf),
  proportion = c(prop_all, prop_inf)
)
write.table(
  stats_df,
  file = file.path(out_dir, paste("conserved_proportion_", sChr, ".tsv", sep = "")),
  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
)

# --------------------- Nonparametric bootstrap + phyloP ---------------------
splitElements <- split.feat(informativeElements, f = splitLength, drop = TRUE)
split_parent_idx <- compute_split_parent_idx(informativeElements, splitLength)
if ("group" %in% colnames(informativeElements)) {
  splitElements$group <- informativeElements$group[split_parent_idx]
} else {
  splitElements$group <- paste0("ID=elem", split_parent_idx)
}

elementAlign <- extract.feature.msa(copy.msa(alnAll), informativeElements, pointer.only = TRUE)
target_label <- "BAR"

# 窗口级 null（缓存 + seed 校验 + 安全 compat）
make_cache_file <- function(kind, sChr, species_label, splitLength, nrep, cache_dir) {
  if (kind == "win") {
    file.path(cache_dir, sprintf("cache_null_win_%s_%s_split%d_nrep%d.rds", sChr, species_label, splitLength, nrep))
  } else if (kind == "site") {
    file.path(cache_dir, sprintf("cache_null_site_%s_%s_nrep%d.rds", sChr, species_label, nrep))
  } else {
    stop("unknown cache kind")
  }
}

win_cache_file <- make_cache_file("win", sChr, species_label, splitLength, nrep, cache_dir)
nonParaNullLRT <- NULL
if (use_cache == 1L && rebuild_cache == 0L && file.exists(win_cache_file)) {
  obj <- readRDS(win_cache_file)
  if (is.list(obj) && !is.null(obj$data) && length(obj$data) == nrep &&
      !is.null(obj$meta$seed) && as.integer(obj$meta$seed) == as.integer(seed_main)) {
    nonParaNullLRT <- obj$data
    cat("Loaded window-level null from cache (seed matched): ", win_cache_file, "\n", sep = "")
  } else {
    cat("Window-level cache incompatible or seed mismatched; will rebuild: ", win_cache_file, "\n", sep = "")
  }
}
if (is.null(nonParaNullLRT)) {
  nonParaNullLRT <- with_seed(seed_main + 1001L, {
    est_cols <- as.double(nrep) * as.double(splitLength)
    do_once <- (compat == 1L && est_cols <= max_once_cols)
    if (do_once) {
      cat(sprintf("Compat=1 & cols<=%d: window-level null once (nrep=%d,len=%d, cols=%d)\n", max_once_cols, nrep, splitLength, as.integer(est_cols)))
      simMsa <- sample.msa(elementAlign, nrep * splitLength, replace = TRUE)
      startIdx <- seq(from = 1L, by = splitLength, length.out = nrep)
      features <- feat(seqname = names.msa(simMsa)[1], src = "sim", feat = ".", start = startIdx, end = startIdx + splitLength - 1L)
      out <- phyloP_acc_on_single_branch(neutralModAll, msa = simMsa, features = features, target_label = target_label)$lnlratio
      rm(simMsa, features); gc()
      out
    } else {
      if (compat == 1L && est_cols > max_once_cols) {
        cat(sprintf("Compat=1 requested, cols=%d > max_once_cols=%d; using chunked mode.\n", as.integer(est_cols), as.integer(max_once_cols)))
      }
      n_chunks <- as.integer(ceiling(nrep / win_null_chunk))
      out <- numeric(0)
      for (cc in seq_len(n_chunks)) {
        this_n <- if (cc < n_chunks) win_null_chunk else (nrep - win_null_chunk * (n_chunks - 1L))
        simMsa <- sample.msa(elementAlign, this_n * splitLength, replace = TRUE)
        startIdx <- seq(from = 1L, by = splitLength, length.out = this_n)
        features <- feat(seqname = names.msa(simMsa)[1], src = "sim", feat = ".", start = startIdx, end = startIdx + splitLength - 1L)
        lrt_chunk <- phyloP_acc_on_single_branch(neutralModAll, msa = simMsa, features = features, target_label = target_label)$lnlratio
        out <- c(out, lrt_chunk)
        rm(simMsa, features, lrt_chunk); gc()
        if (cc %% 5 == 0 || cc == n_chunks) cat(sprintf("  null chunk %d/%d done\n", cc, n_chunks))
      }
      out
    }
  })
  if (use_cache == 1L) {
    meta <- list(kind="win", sChr=sChr, species_label=species_label, splitLength=splitLength, nrep=nrep, seed=seed_main, timestamp=Sys.time())
    saveRDS(list(meta=meta, data=nonParaNullLRT), file = win_cache_file)
    cat("Saved window-level null to cache: ", win_cache_file, "\n", sep = "")
  }
}

obsPhyloP <- phyloP_acc_on_single_branch(neutralModAll, msa = alnAll, features = splitElements, target_label = target_label)

# 右尾 p 值 + FDR
nonPara_ecdf <- ecdf(nonParaNullLRT)
eps_win <- 1.0 / max(1L, length(nonParaNullLRT))
nonParaPval <- 1 - nonPara_ecdf(obsPhyloP$lnlratio) + eps_win
nonParaFDR <- p.adjust(nonParaPval, method = "BH")

win_lrt_thr <- NA_real_
if (is.finite(win_q) && win_q > 0 && win_q < 1) {
  win_lrt_thr <- as.numeric(quantile(nonParaNullLRT, probs = win_q, names = FALSE, type = 7))
}

sig_idx <- which(nonParaFDR < win_fdr)
if (is.finite(win_lrt_thr)) {
  sig_idx <- sig_idx[obsPhyloP$lnlratio[sig_idx] >= win_lrt_thr]
}
if (is.finite(win_p) && win_p > 0 && win_p < 1) {
  sig_idx <- sig_idx[nonParaPval[sig_idx] < win_p]
}

out_gff <- file.path(out_dir, paste("Traref_", sChr, ".gff", sep = ""))

if (length(sig_idx) == 0) {
  cat(sprintf("No accelerated elements at thresholds: FDR<%.4g, win_q=%.4g, win_p=%.4g. Writing empty GFF: %s\n",
              win_fdr, win_q, win_p, out_gff))
  con <- file(out_gff, "w"); writeLines(paste0("# ", sChr, " no features at given thresholds"), con); close(con)
} else {
  sigWins <- splitElements[sig_idx,]
  nonParaSigFeats <- sigWins
  nonParaSigFeats$feature <- "acc"
  nonParaSigFeats$score <- obsPhyloP$lnlratio[sig_idx]
  extra_info <- sprintf("q=%.6g;p=%.6g", nonParaFDR[sig_idx], nonParaPval[sig_idx])
  write_gff10(nonParaSigFeats, extra_vec = extra_info, out_file = out_gff, override_seqname = sChr)
  cat("Wrote significant windows GFF: ", out_gff, " (n=", nrow(sigWins), ")\n", sep = "")
}

# --------------------- Site-level validation & export ---------------------
out_tsv <- file.path(out_dir, paste("validation_sites_", sChr, ".tsv", sep = ""))
out_gff_valid <- file.path(out_dir, paste("Traref_validated_", sChr, ".gff", sep = ""))
acc_align_dir <- file.path(out_dir, paste("acc_alignments_", sChr, sep = ""))
if (length(sig_idx) > 0 && !dir.exists(acc_align_dir)) {
  dir.create(acc_align_dir, recursive = TRUE, showWarnings = FALSE)
}

if (length(sig_idx) > 0) {
  # site-level null（缓存 + seed 校验）
  site_cache_file <- make_cache_file("site", sChr, species_label, splitLength, site_nrep, cache_dir)
  site_null_lrt <- NULL
  if (use_cache == 1L && rebuild_cache == 0L && file.exists(site_cache_file)) {
    obj <- readRDS(site_cache_file)
    if (is.list(obj) && !is.null(obj$data) && length(obj$data) == site_nrep &&
        !is.null(obj$meta$seed) && as.integer(obj$meta$seed) == as.integer(seed_main)) {
      site_null_lrt <- obj$data
      cat("Loaded site-level null from cache (seed matched): ", site_cache_file, "\n", sep = "")
    } else {
      cat("Site-level cache incompatible or seed mismatched; will rebuild: ", site_cache_file, "\n", sep = "")
    }
  }
  if (is.null(site_null_lrt)) {
    site_null_lrt <- with_seed(seed_main + 2001L, {
      do_once_site <- (compat == 1L && site_nrep <= site_null_chunk)
      if (do_once_site) {
        cat(sprintf("Compat=1: site-level null once (n=%d)\n", site_nrep))
        simMsaSite <- sample.msa(elementAlign, site_nrep, replace = TRUE)
        siteFeatNull <- feat(seqname = names.msa(simMsaSite)[1], src = "simsite", feat = ".", start = seq_len(site_nrep), end = seq_len(site_nrep))
        out <- phyloP_acc_on_single_branch(neutralModAll, msa = simMsaSite, features = siteFeatNull, target_label = target_label)$lnlratio
        rm(simMsaSite, siteFeatNull); gc()
        out
      } else {
        if (compat == 1L && site_nrep > site_null_chunk) {
          cat(sprintf("Compat=1 requested for site-level, n=%d > chunk=%d; chunked mode.\n", site_nrep, site_null_chunk))
        }
        n_chunks_site <- as.integer(ceiling(site_nrep / site_null_chunk))
        out <- numeric(0)
        for (cc in seq_len(n_chunks_site)) {
          this_n <- if (cc < n_chunks_site) site_null_chunk else (site_nrep - site_null_chunk * (n_chunks_site - 1L))
          simMsaSite <- sample.msa(elementAlign, this_n, replace = TRUE)
          siteFeatNull <- feat(seqname = names.msa(simMsaSite)[1], src = "simsite", feat = ".", start = seq_len(this_n), end = seq_len(this_n))
          lrt_chunk <- phyloP_acc_on_single_branch(neutralModAll, msa = simMsaSite, features = siteFeatNull, target_label = target_label)$lnlratio
          out <- c(out, lrt_chunk)
          rm(simMsaSite, siteFeatNull, lrt_chunk); gc()
          if (cc %% 5 == 0 || cc == n_chunks_site) cat(sprintf("  site-null chunk %d/%d done\n", cc, n_chunks_site))
        }
        out
      }
    })
    if (use_cache == 1L) {
      meta <- list(kind="site", sChr=sChr, species_label=species_label, nrep=site_nrep, seed=seed_main, timestamp=Sys.time())
      saveRDS(list(meta=meta, data=site_null_lrt), file = site_cache_file)
      cat("Saved site-level null to cache: ", site_cache_file, "\n", sep = "")
    }
  }

  site_lrt_thr <- as.numeric(quantile(site_null_lrt, probs = site_q, names = FALSE, type = 7))
  cat(sprintf("Site-level LRT threshold at q=%.3f is %.6f\n", site_q, site_lrt_thr))
  cdf_site_null <- ecdf(site_null_lrt)

  sigWins <- splitElements[sig_idx,]
  starts_vec <- integer(0); ends_vec <- integer(0); win_id <- integer(0)
  for (i in seq_len(nrow(sigWins))) {
    s <- sigWins$start[i]; e <- sigWins$end[i]
    if (is.na(s) || is.na(e) || e < s) next
    pos <- seq.int(s, e)
    starts_vec <- c(starts_vec, pos)
    ends_vec <- c(ends_vec, pos)
    win_id <- c(win_id, rep.int(i, length(pos)))
  }

  lrt_site_obs <- numeric(length(starts_vec))
  if (length(starts_vec) > 0) {
    n_sites <- length(starts_vec)
    n_chunks_obs <- as.integer(ceiling(n_sites / site_obs_chunk))
    cat(sprintf("Computing site-level observed LRT for %d sites (chunk=%d)...\n", n_sites, site_obs_chunk))
    for (cc in seq_len(n_chunks_obs)) {
      lo <- (cc - 1L) * site_obs_chunk + 1L
      hi <- min(cc * site_obs_chunk, n_sites)
      idx <- lo:hi
      siteFeatObs <- feat(seqname = rep(names.msa(alnAll)[1], length(idx)), src = "site", feat = ".", start = starts_vec[idx], end = ends_vec[idx])
      phy_site_obs <- phyloP_acc_on_single_branch(neutralModAll, msa = alnAll, features = siteFeatObs, target_label = target_label)
      lrt_site_obs[idx] <- phy_site_obs$lnlratio
      rm(siteFeatObs, phy_site_obs); gc()
      if (cc %% 10 == 0 || cc == n_chunks_obs) cat(sprintf("  site-obs chunk %d/%d done\n", cc, n_chunks_obs))
    }
  }

  nsig_by_win <- if (length(lrt_site_obs) > 0) tapply(lrt_site_obs > site_lrt_thr, INDEX = win_id, FUN = sum) else setNames(integer(nrow(sigWins)), as.character(seq_len(nrow(sigWins))))
  n_sig_sites <- as.integer(if (length(nsig_by_win) == 0) rep(0L, nrow(sigWins)) else nsig_by_win[as.character(seq_len(nrow(sigWins)))])
  top_k <- max(1L, min(site_top_k, splitLength))
  top_info <- vector("character", length = nrow(sigWins))
  if (length(lrt_site_obs) > 0) {
    lrt_split <- split(lrt_site_obs, win_id)
    pos_split <- split(starts_vec, win_id)
    names(top_info) <- as.character(seq_len(nrow(sigWins)))
    for (i in seq_len(nrow(sigWins))) {
      lrt_i <- lrt_split[[as.character(i)]]
      pos_i <- pos_split[[as.character(i)]]
      if (is.null(lrt_i) || length(lrt_i) == 0) {
        top_info[i] <- ""
      } else {
        ord <- order(lrt_i, decreasing = TRUE)
        k <- min(top_k, length(ord))
        picks <- ord[seq_len(k)]
        relpos <- pos_i[picks] - sigWins$start[i] + 1L
        top_info[i] <- paste(sprintf("%d:%.3f", relpos, lrt_i[picks]), collapse = ",")
      }
    }
  }

  qvals_sig <- nonParaFDR[sig_idx]
  pvals_sig <- nonParaPval[sig_idx]

  df_valid <- data.frame(
    species = sRefSp,
    seq = sChr,
    start = sigWins$start,
    end = sigWins$end,
    length = sigWins$end - sigWins$start + 1L,
    win_LRT = obsPhyloP$lnlratio[sig_idx],
    p = pvals_sig,
    q = qvals_sig,
    site_LRT_thr = rep(site_lrt_thr, length(sig_idx)),
    site_q = rep(site_q, length(sig_idx)),
    n_sig_sites = n_sig_sites,
    site_min = rep(site_min, length(sig_idx)),
    site_top_relpos_lrt = top_info,
    stringsAsFactors = FALSE
  )
  write.table(df_valid, file = out_tsv, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  cat("Wrote validation TSV: ", out_tsv, "\n", sep = "")

  pass <- which(n_sig_sites >= site_min)
  if (length(pass) > 0) {
    validatedFeats <- sigWins[pass,]
    validatedFeats$feature <- "acc"
    validatedFeats$score <- obsPhyloP$lnlratio[sig_idx[pass]]
    extra_valid <- sprintf("validate=1;nsig_sites=%d;q=%.6g;p=%.6g", n_sig_sites[pass], qvals_sig[pass], pvals_sig[pass])
    write_gff10(validatedFeats, extra_vec = extra_valid, out_file = out_gff_valid, override_seqname = sChr)
    cat("Validation summary: ", length(pass), " / ", nrow(sigWins), " windows passed (>= site_min). File: ", out_gff_valid, "\n", sep = "")
  } else {
    con <- file(out_gff_valid, "w")
    writeLines(paste0("# ", sChr, " no validated features at rule: site_min=", site_min, ", site_q=", site_q), con)
    close(con)
    cat("Validation summary: 0 / ", nrow(sigWins), " windows passed. Wrote empty file: ", out_gff_valid, "\n", sep = "")
  }

  # -------- 导出加速 MSA --------
  if (export_msa == 1L) {
    acc_feats <- NULL
    if (acc_msa_set == "validated" && exists("validatedFeats") && nrow(validatedFeats) > 0) {
      acc_feats <- validatedFeats
    } else if (nrow(sigWins) > 0) {
      acc_feats <- sigWins
    }
    out_acc_msa <- file.path(msa_dir, paste("msa_accelerated_", acc_msa_set, "_", sChr, ".", tolower(msa_format), sep = ""))
    ok_acc <- write_concat_msa(alnAll, acc_feats, out_acc_msa, fmt = msa_format)
    if (ok_acc) {
      cat("Wrote accelerated MSA: ", out_acc_msa, "\n", sep = "")
    } else {
      cat("Failed to write accelerated MSA: ", out_acc_msa, "\n", sep = "")
    }
  }

  # --------------------- Export（不叠加竖直红线） ---------------------
  eps <- 1.0 / max(1L, length(site_null_lrt))
  phy_obj <- if (inherits(neutralModAll$tree, "phylo")) neutralModAll$tree else ape::read.tree(text = neutralModAll$tree)

  present_regex <- allowed_regex_for(strict_present)
  strict_sp_list <- parse_strict_species(strict_species, msa_names, arrFreshs)

  align_index <- vector("list", 0)
  filtered_cnt <- 0L

  for (j in seq_len(nrow(sigWins))) {
    s <- sigWins$start[j]; e <- sigWins$end[j]
    idx_j <- which(win_id == j)

    seqs_named <- extract_window_seqs_only(alnAll, s, e, ref_sp_for_feature = sRefSp)

    strict_ok <- NA_integer_
    if (length(strict_sp_list) > 0 && !is.na(present_regex)) {
      if (length(seqs_named) == 0) {
        strict_ok <- 0L
      } else {
        ok <- TRUE
        hdrs <- names(seqs_named)
        for (sp in strict_sp_list) {
          h <- find_seq_header_name(hdrs, sp)
          if (is.na(h)) { ok <- FALSE; break }
          if (!has_any_allowed_base(seqs_named[[h]], present_regex)) { ok <- FALSE; break }
        }
        strict_ok <- as.integer(ok)
      }
    }

    if (enforce_strict == 1L && !is.na(strict_ok) && strict_ok == 0L) {
      filtered_cnt <- filtered_cnt + 1L
      cat(sprintf("SKIP window %d [%d-%d]: strict_present=%s failed for strict_species=%s\n",
                  j, s, e, strict_present, strict_species))
      next
    }

    # RAW FASTA
    out_raw_fa <- file.path(acc_align_dir, sprintf("win_%05d_%d_%d.raw.fasta", j, s, e))
    write_fasta_from_named_vec(seqs_named, out_raw_fa, width = 60L)

    # Annotated FASTA
    out_fa <- file.path(acc_align_dir, sprintf("win_%05d_%d_%d.fasta", j, s, e))
    top_pos <- integer(0)
    relpos_j <- integer(0); lrt_j <- numeric(0); over_thr_j <- logical(0); p_emp_j <- numeric(0)
    if (length(idx_j) > 0) {
      relpos_j <- starts_vec[idx_j] - s + 1L
      lrt_j <- lrt_site_obs[idx_j]
      over_thr_j <- lrt_j > site_lrt_thr
      p_emp_j <- 1 - cdf_site_null(lrt_j) + eps
      if (mark_rule == "thr") {
        top_pos <- relpos_j[over_thr_j]
      } else {
        ord <- order(lrt_j, decreasing = TRUE)
        k <- min(site_top_k, length(ord))
        top_pos <- relpos_j[ord[seq_len(k)]]
      }
    }
    seqs_ann <- vapply(seqs_named, annotate_seq_by_positions, character(1), top_pos = top_pos)
    write_fasta_from_named_vec(seqs_ann, out_fa, width = 60L)

    # per-site TSV
    out_site <- file.path(acc_align_dir, sprintf("win_%05d_%d_%d_sites.tsv", j, s, e))
    if (length(idx_j) == 0) {
      write.table(data.frame(), file = out_site, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    } else {
      nms <- names(seqs_named)
      rad_hdr <- find_seq_header_name(nms, "Tridentiger_radiatus")
      bar_hdr <- find_seq_header_name(nms, "Tridentiger_barbatus")
      rad_seq <- if (!is.na(rad_hdr)) seqs_named[[rad_hdr]] else NULL
      bar_seq <- if (!is.na(bar_hdr)) seqs_named[[bar_hdr]] else NULL
      base_rad <- if (!is.null(rad_seq)) substring(rad_seq, relpos_j, relpos_j) else rep(NA_character_, length(relpos_j))
      base_bar <- if (!is.null(bar_seq)) substring(bar_seq, relpos_j, relpos_j) else rep(NA_character_, length(relpos_j))
      has_rad <- if (!is.null(rad_seq)) as.integer(base_rad != "-") else rep(NA_integer_, length(relpos_j))
      has_bar <- if (!is.null(bar_seq)) as.integer(base_bar != "-") else rep(NA_integer_, length(relpos_j))
      df_site <- data.frame(
        seq = sChr, win_index = j, start = s, end = e,
        relpos = relpos_j, abs_col = starts_vec[idx_j],
        LRT = lrt_j, p_emp_right = p_emp_j, over_thr = as.integer(over_thr_j),
        base_rad = base_rad, base_bar = base_bar, has_rad = has_rad, has_bar = has_bar,
        stringsAsFactors = FALSE
      )
      write.table(df_site, file = out_site, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }

    # msaplot PNG（不叠加竖直红线）
    tree_png_path <- ""
    if (plot_tree == 1L) {
      msaplot_fa <- file.path(acc_align_dir, sprintf("win_%05d_%d_%d.msaplot.fasta", j, s, e))
      phy_obj <- if (inherits(neutralModAll$tree, "phylo")) neutralModAll$tree else ape::read.tree(text = neutralModAll$tree)
      ok_fa <- build_msaplot_fasta(phy_obj, seqs_named, msaplot_fa)
      if (ok_fa) {
        out_png <- file.path(acc_align_dir, sprintf("win_%05d_%d_%d_tree.png", j, s, e))
        main_t <- sprintf("%s | %s:%d-%d | LRT=%.3f q=%.3g p=%.3g",
                          speciesArg, sChr, s, e, obsPhyloP$lnlratio[sig_idx[j]],
                          nonParaFDR[sig_idx[j]], nonParaPval[sig_idx[j]])
        relpos_over <- if (length(idx_j) > 0) relpos_j[which(lrt_j > site_lrt_thr)] else integer(0)
        ok_png <- plot_tree_and_alignment_msaplot(
          phy_obj, msaplot_fa, relpos_over, out_png, title_txt = main_t,
          msa_offset = msa_offset, msa_width = msa_width
        )
        if (ok_png) {
          tree_png_path <- out_png
        } else {
          cat(sprintf("SKIP PNG window %d [%d-%d]: msaplot failed\n", j, s, e))
        }
      } else {
        cat(sprintf("SKIP PNG window %d [%d-%d]: no sequences for msaplot\n", j, s, e))
      }
    }

    grp_j <- if ("group" %in% colnames(sigWins)) sigWins$group[j] else "."
    align_index[[length(align_index)+1L]] <- data.frame(
      win_index = j, start = s, end = e, length = e - s + 1L,
      win_LRT = obsPhyloP$lnlratio[sig_idx[j]],
      p = nonParaPval[sig_idx[j]], q = nonParaFDR[sig_idx[j]],
      group = grp_j,
      file = file.path(acc_align_dir, sprintf("win_%05d_%d_%d.fasta", j, s, e)),
      raw_file = file.path(acc_align_dir, sprintf("win_%05d_%d_%d.raw.fasta", j, s, e)),
      sites = file.path(acc_align_dir, sprintf("win_%05d_%d_%d_sites.tsv", j, s, e)),
      msaplot_fa = if (plot_tree == 1L && tree_png_path != "") file.path(acc_align_dir, sprintf("win_%05d_%d_%d.msaplot.fasta", j, s, e)) else "",
      tree_png = tree_png_path,
      mark_rule = mark_rule,
      strict_present = strict_present,
      strict_species = strict_species,
      strict_ok = ifelse(is.na(strict_ok), NA_integer_, strict_ok),
      stringsAsFactors = FALSE
    )

    if (j %% 200 == 0 || j == nrow(sigWins)) {
      cat(sprintf("  export %d/%d\n", j, nrow(sigWins)))
    }
  }

  idx_file <- file.path(out_dir, paste("Traref_alignments_index_", sChr, ".tsv", sep = ""))
  if (length(align_index) > 0) {
    align_index_df <- do.call(rbind, align_index)
    write.table(align_index_df, file = idx_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  } else {
    write.table(data.frame(win_index=integer(), start=integer(), end=integer(), length=integer(),
                           win_LRT=numeric(), p=numeric(), q=numeric(), group=character(),
                           file=character(), raw_file=character(), sites=character(), msaplot_fa=character(), tree_png=character(),
                           mark_rule=character(), strict_present=character(), strict_species=character(), strict_ok=integer() ),
                file = idx_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
  cat("Raw/annotated FASTAs, per-site TSVs, and msaplot PNGs written to: ", acc_align_dir, " (index: ", idx_file, ")\n", sep = "")
  if (enforce_strict == 1L) {
    cat(sprintf("Strict-present enforcement filtered %d window(s).\n", filtered_cnt))
  }

  # --------------------- Extra diagnostics: signal vs model ---------------------
  acc_diag_file <- file.path(out_dir, paste("acc_signal_diagnostics_", sChr, ".tsv", sep = ""))

  sig_len <- if (exists("sigWins") && nrow(sigWins) > 0) sum(sigWins$end - sigWins$start + 1L) else 0L
  val_len <- if (exists("validatedFeats") && nrow(validatedFeats) > 0) sum(validatedFeats$end - validatedFeats$start + 1L) else 0L

  sig_win_LRT <- if (exists("sigWins") && nrow(sigWins) > 0) obsPhyloP$lnlratio[sig_idx] else numeric(0)
  val_win_LRT <- if (exists("validatedFeats") && nrow(validatedFeats) > 0 && exists("pass") && length(pass) > 0) obsPhyloP$lnlratio[sig_idx[pass]] else numeric(0)

  site_sig_ratio <- if (length(n_sig_sites) > 0 && sig_len > 0) {
    sum(n_sig_sites) / sig_len
  } else NA_real_

  site_lrt_stats <- if (length(lrt_site_obs) > 0) {
    as.numeric(quantile(lrt_site_obs, probs = c(0.5, 0.9, 0.99), names = FALSE))
  } else c(NA_real_, NA_real_, NA_real_)

  diag_df <- data.frame(
    seq = sChr,
    n_sig_windows = if (exists("sigWins")) nrow(sigWins) else 0L,
    n_valid_windows = if (exists("validatedFeats")) nrow(validatedFeats) else 0L,
    sig_total_len = sig_len,
    valid_total_len = val_len,
    sig_win_LRT_median = if (length(sig_win_LRT) > 0) median(sig_win_LRT) else NA_real_,
    sig_win_LRT_p90 = if (length(sig_win_LRT) > 0) quantile(sig_win_LRT, 0.9, names = FALSE) else NA_real_,
    valid_win_LRT_median = if (length(val_win_LRT) > 0) median(val_win_LRT) else NA_real_,
    site_sig_ratio = site_sig_ratio,
    site_LRT_q50 = site_lrt_stats[1],
    site_LRT_q90 = site_lrt_stats[2],
    site_LRT_q99 = site_lrt_stats[3],
    win_fdr = win_fdr,
    win_q = win_q,
    win_p = win_p,
    win_lrt_thr = win_lrt_thr,
    stringsAsFactors = FALSE
  )
  write.table(diag_df, file = acc_diag_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  cat("Wrote acceleration diagnostics TSV: ", acc_diag_file, "\n", sep = "")

  # --------------------- Two-tree comparison (Conserved vs Accelerated) ---------------------
  consEleModel <- NULL
  rarModel <- NULL
  try(consEleModel <- phyloFit(elementAlign, init.mod = neutralModAll, no.opt = c("backgd","ratematrix")), silent = TRUE)
  try({
    rarAlign <- extract.feature.msa(copy.msa(alnAll), sigWins, pointer.only = TRUE)
    rarModel <- phyloFit(rarAlign, init.mod = neutralModAll, no.opt = c("backgd","ratematrix"))
  }, silent = TRUE)
  cmp_png <- file.path(out_dir, paste("trees_compare_", sChr, ".png", sep = ""))
  if (!is.null(consEleModel) && !is.null(rarModel)) {
    ok_cmp <- FALSE
    try(ok_cmp <- plot_two_trees_side_by_side(consEleModel$tree, rarModel$tree, cmp_png), silent = TRUE)
    if (isTRUE(ok_cmp)) {
      cat("Wrote side-by-side tree comparison PNG: ", cmp_png, "\n", sep = "")
    } else {
      cat("Failed to write tree comparison PNG: ", cmp_png, "\n", sep = "")
    }
  } else {
    cat("Skip tree comparison: model fitting failed for conserved or accelerated set.\n")
  }

  # --------------------- BAR-branch-only scaling + highlight ---------------------
  to_phy <- function(obj) {
    if (inherits(obj, "phylo")) return(obj)
    if (is.character(obj)) return(ape::read.tree(text = obj))
    if (!is.null(obj$tree)) {
      if (inherits(obj$tree, "phylo")) return(obj$tree)
      if (is.character(obj$tree)) return(ape::read.tree(text = obj$tree))
    }
    stop("Invalid tree object.")
  }
  # Find MRCA node for BAR; fallback to internal node label "BAR"
  find_bar_node <- function(phy, rad_tip = "Tridentiger_radiatus", bar_tip = "Tridentiger_barbatus", label = "BAR") {
    tips <- phy$tip.label
    if (all(c(rad_tip, bar_tip) %in% tips)) {
      n <- ape::getMRCA(phy, c(rad_tip, bar_tip))
      if (is.finite(n)) return(n)
    }
    if (!is.null(phy$node.label)) {
      idx <- which(phy$node.label == label)
      if (length(idx) >= 1) return(ape::Ntip(phy) + idx[1])
    }
    NA_integer_
  }
  # Parent edge of a node
  parent_edge_of_node <- function(phy, node) {
    if (is.na(node)) return(NA_integer_)
    e <- phy$edge
    idx <- which(e[,2] == node)
    if (length(idx) >= 1) return(idx[1])
    NA_integer_
  }
  # Summaries with BAR single branch only
  bar_metrics <- function(phy, rad_tip = "Tridentiger_radiatus", bar_tip = "Tridentiger_barbatus") {
    n <- find_bar_node(phy, rad_tip, bar_tip, label = "BAR")
    eidx <- parent_edge_of_node(phy, n)
    Ltot <- sum(phy$edge.length)
    Lbar <- if (is.finite(eidx)) phy$edge.length[eidx] else NA_real_
    list(node = n, edge = eidx, L_total = Ltot, L_BAR = Lbar, L_nonBAR = Ltot - Lbar, frac_BAR = Lbar / Ltot)
  }
  # Plot two trees highlighting only the MRCA parent-edge (single branch)
  plot_two_trees_highlight_bar <- function(phy1, phy2, outfile_png,
                                          title1 = "Conserved (BAR scaling)",
                                          title2 = "Accelerated (BAR scaling)") {
    depth <- function(phy) max(ape::node.depth.edgelength(phy))
    xmax <- max(depth(phy1), depth(phy2))
    if (!is.finite(xmax) || xmax <= 0) xmax <- 1
    m1 <- bar_metrics(phy1); m2 <- bar_metrics(phy2)
    cols1 <- rep("gray40", nrow(phy1$edge)); cols2 <- rep("gray40", nrow(phy2$edge))
    w1 <- rep(1.5, nrow(phy1$edge));      w2 <- rep(1.5, nrow(phy2$edge))
    if (!is.na(m1$edge)) { cols1[m1$edge] <- "red"; w1[m1$edge] <- 3.5 }
    if (!is.na(m2$edge)) { cols2[m2$edge] <- "red"; w2[m2$edge] <- 3.5 }
    png(outfile_png, width = 1400, height = 700, res = 150)
    op <- par(no.readonly = TRUE); on.exit({ par(op); dev.off() }, add = TRUE)
    par(mfrow = c(1,2), mar = c(5,2,4,1))
    plot.phylo(phy1, x.lim = c(0, xmax * 1.05), main = title1, cex = 0.85,
               edge.color = cols1, edge.width = w1)
    mtext(sprintf("Total=%.3f, BAR=%.3f (%.2f%%)", m1$L_total, m1$L_BAR, 100*m1$frac_BAR),
          side = 1, line = 3, cex = 0.8)
    plot.phylo(phy2, x.lim = c(0, xmax * 1.05), main = title2, cex = 0.85,
               edge.color = cols2, edge.width = w2)
    mtext(sprintf("Total=%.3f, BAR=%.3f (%.2f%%)", m2$L_total, m2$L_BAR, 100*m2$frac_BAR),
          side = 1, line = 3, cex = 0.8)
    invisible(TRUE)
  }

  # Fit comparable models: conserved vs accelerated, both BAR scaling (fallback to subtree if needed)
  consScaled <- NULL
  rarBARScaled <- NULL
  ok_scale_args <- "scale.only" %in% names(formals(phyloFit))
  can_scale_branch <- "scale.branches" %in% names(formals(phyloFit))
  can_scale_subtree <- "scale.subtree" %in% names(formals(phyloFit))

  scale_arg_name <- NA_character_
  metBase <- NULL
  metC <- NULL
  metR <- NULL

  if (!ok_scale_args) {
    cat("Skip BAR scaling: rphast::phyloFit has no 'scale.only' parameter in this build.\n")
  } else if (!can_scale_branch && !can_scale_subtree) {
    cat("Skip BAR scaling: rphast::phyloFit has no 'scale.branches' or 'scale.subtree' parameter in this build.\n")
  } else {
    scale_arg_name <- if (can_scale_branch) "scale.branches" else "scale.subtree"
    if (!can_scale_branch && can_scale_subtree) {
      cat("phyloFit has no 'scale.branches'; fallback to 'scale.subtree' (model uses subtree).\n")
    }

    if (can_scale_branch) {
      try(consScaled <- phyloFit(elementAlign, init.mod = neutralModAll,
                                 scale.only = TRUE, scale.branches = target_label,
                                 no.opt = c("backgd","ratematrix")),
          silent = TRUE)
    } else {
      try(consScaled <- phyloFit(elementAlign, init.mod = neutralModAll,
                                 scale.only = TRUE, scale.subtree = target_label,
                                 no.opt = c("backgd","ratematrix")),
          silent = TRUE)
    }

    rarAlign <- extract.feature.msa(copy.msa(alnAll), sigWins, pointer.only = TRUE)
    if (can_scale_branch) {
      try(rarBARScaled <- phyloFit(rarAlign, init.mod = neutralModAll,
                                  scale.only = TRUE, scale.branches = target_label,
                                  no.opt = c("backgd","ratematrix")),
          silent = TRUE)
    } else {
      try(rarBARScaled <- phyloFit(rarAlign, init.mod = neutralModAll,
                                  scale.only = TRUE, scale.subtree = target_label,
                                  no.opt = c("backgd","ratematrix")),
          silent = TRUE)
    }

    if (!is.null(consScaled) && !is.null(rarBARScaled)) {
      phyBase <- to_phy(neutralModAll$tree)
      phyC <- to_phy(consScaled$tree)
      phyR <- to_phy(rarBARScaled$tree)

      metBase <- bar_metrics(phyBase)
      metC <- bar_metrics(phyC)
      metR <- bar_metrics(phyR)

      scaleC <- if (is.finite(metBase$L_BAR) && metBase$L_BAR > 0) metC$L_BAR / metBase$L_BAR else NA_real_
      scaleR <- if (is.finite(metBase$L_BAR) && metBase$L_BAR > 0) metR$L_BAR / metBase$L_BAR else NA_real_

      cat(sprintf("BAR scaling (scale.only + %s=%s)\n", scale_arg_name, target_label))
      cat(sprintf("  Base:        TotalLen=%.6f  BAR=%.6f  fracBAR=%.4f\n", metBase$L_total, metBase$L_BAR, metBase$frac_BAR))
      cat(sprintf("  Conserved:   TotalLen=%.6f  BAR=%.6f  fracBAR=%.4f  scaleBAR=%.4f\n",
                  metC$L_total, metC$L_BAR, metC$frac_BAR, scaleC))
      cat(sprintf("  Accelerated: TotalLen=%.6f  BAR=%.6f  fracBAR=%.4f  scaleBAR=%.4f\n",
                  metR$L_total, metR$L_BAR, metR$frac_BAR, scaleR))
      if (is.finite(metC$L_BAR) && metC$L_BAR > 0) {
        cat(sprintf("  BAR length fold-change (Acc/Cons): %.4f\n", metR$L_BAR / metC$L_BAR))
      }
      if (is.finite(scaleC) && is.finite(scaleR) && scaleC > 0) {
        cat(sprintf("  BAR scale fold-change (Acc/Cons): %.4f\n", scaleR / scaleC))
      }

      bar_stats_file <- file.path(out_dir, paste("trees_BAR_scaled_stats_", sChr, ".tsv", sep = ""))
      bar_stats_df <- data.frame(
        set = c("neutral_base", "conserved", "accelerated"),
        scale_method = c("none", scale_arg_name, scale_arg_name),
        total_len = c(metBase$L_total, metC$L_total, metR$L_total),
        bar_len = c(metBase$L_BAR, metC$L_BAR, metR$L_BAR),
        bar_frac = c(metBase$frac_BAR, metC$frac_BAR, metR$frac_BAR),
        bar_scale_vs_base = c(1, scaleC, scaleR),
        stringsAsFactors = FALSE
      )
      write.table(bar_stats_df, file = bar_stats_file, sep = "\t",
                  row.names = FALSE, col.names = TRUE, quote = FALSE)
      cat("Wrote BAR scaling stats TSV: ", bar_stats_file, "\n", sep = "")

      out_bar_png <- file.path(out_dir, paste("trees_BAR_scaled_highlight_", sChr, ".png", sep = ""))
      ok_bar <- FALSE
      try(ok_bar <- plot_two_trees_highlight_bar(phyC, phyR, out_bar_png), silent = TRUE)
      if (isTRUE(ok_bar)) {
        cat("Wrote BAR-branch-highlight comparison PNG: ", out_bar_png, "\n", sep = "")
      } else {
        cat("Failed to write BAR-branch-highlight PNG: ", out_bar_png, "\n", sep = "")
      }
    } else {
      cat("BAR scaling: model fitting failed for conserved or accelerated set.\n")
    }
  }

  # --------------------- BAR model diagnostics ---------------------
  bar_diag_file <- file.path(out_dir, paste("bar_model_diagnostics_", sChr, ".tsv", sep = ""))
  bar_diag <- data.frame(
    seq = sChr,
    scale_method = scale_arg_name,
    bar_edge_found_base = if (!is.null(metBase)) !is.na(metBase$edge) else NA,
    bar_edge_found_cons = if (!is.null(metC)) !is.na(metC$edge) else NA,
    bar_edge_found_acc = if (!is.null(metR)) !is.na(metR$edge) else NA,
    stringsAsFactors = FALSE
  )
  write.table(bar_diag, file = bar_diag_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  cat("Wrote BAR model diagnostics TSV: ", bar_diag_file, "\n", sep = "")

} else {
  con <- file(out_tsv, "w"); writeLines(paste0("species\tseq\tstart\tend\tlength\twin_LRT\tp\tq\tsite_LRT_thr\tsite_q\tn_sig_sites\tsite_min\tsite_top_relpos_lrt"), con); close(con)
  con <- file(out_gff_valid, "w"); writeLines(paste0("# ", sChr, " no validated features (no significant windows)"), con); close(con)
  cat("No significant windows; wrote empty validation TSV and GFF.\n")

  # 没有显著窗口时也给一个空的加速 MSA 文件（便于流水线）
  if (export_msa == 1L) {
    out_acc_msa <- file.path(msa_dir, paste("msa_accelerated_", acc_msa_set, "_", sChr, ".", tolower(msa_format), sep = ""))
    con <- file(out_acc_msa, "w"); writeLines("# empty MSA: no significant windows", con); close(con)
    cat("Wrote empty accelerated MSA: ", out_acc_msa, "\n", sep = "")
  }

  # 也输出空诊断文件（便于比较）
  acc_diag_file <- file.path(out_dir, paste("acc_signal_diagnostics_", sChr, ".tsv", sep = ""))
  diag_df <- data.frame(
    seq = sChr,
    n_sig_windows = 0L,
    n_valid_windows = 0L,
    sig_total_len = 0L,
    valid_total_len = 0L,
    sig_win_LRT_median = NA_real_,
    sig_win_LRT_p90 = NA_real_,
    valid_win_LRT_median = NA_real_,
    site_sig_ratio = NA_real_,
    site_LRT_q50 = NA_real_,
    site_LRT_q90 = NA_real_,
    site_LRT_q99 = NA_real_,
    win_fdr = win_fdr,
    win_q = win_q,
    win_p = win_p,
    win_lrt_thr = win_lrt_thr,
    stringsAsFactors = FALSE
  )
  write.table(diag_df, file = acc_diag_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  cat("Wrote acceleration diagnostics TSV: ", acc_diag_file, "\n", sep = "")
}
