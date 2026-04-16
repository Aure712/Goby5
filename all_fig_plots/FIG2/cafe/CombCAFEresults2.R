# ====================== 参数区（请按需修改） ======================
gene_list_files <- c(
  "TridentigerBarbatusExpanded_Contracted_summary.csv",
  "TridentigerBifasciatusExpanded_Contracted_summary.csv",
  "TridentigerRadiatusExpanded_Contracted_summary.csv"
)

suffixes <- c("tba", "tbi", "tra")   # 与上面文件一一对应

go_files <- c(
  "alltba_genes.fa.tsv",
  "alltbi_genes.fa.tsv",
  "alltra_genes.fa.tsv"
)

fa_files <- c(
  "alltba_genes.fa",
  "alltbi_genes.fa",
  "alltra_genes.fa"
)

go_out_dir   <- "go_summaries"
out_combined <- "all_species_gene_families_merged_with_genes.tsv"
# =================================================================

options(stringsAsFactors = FALSE)

# ====================== 工具函数 ======================
trim_ws <- function(x) gsub("^\\s+|\\s+$", "", x, perl = TRUE)

split_mrnas <- function(s) {
  if (length(s) == 0 || is.na(s) || !nzchar(s)) return(character(0))
  s <- gsub("[,|]", ";", s)
  parts <- trim_ws(unlist(strsplit(s, ";", fixed = TRUE)))
  unique(parts[nzchar(parts) & parts != "NA"])
}

extract_go_ids <- function(x) {
  if (is.null(x) || length(x) == 0) return(character(0))
  m <- regmatches(x, gregexpr("GO:[0-9]+", x, perl = TRUE))
  unique(trim_ws(unlist(m)))
}

# 从 fasta 提取该物种专属的 mrna -> gene 映射
# 头示例：">orig001309-T1 orig001309   TbaScf_1:26949829-26953361(+)"
extract_mrna_to_gene <- function(fa_path) {
  lines <- readLines(fa_path, warn = FALSE)
  hdr <- grep("^>", lines, value = TRUE)
  map <- character()
  for (h in hdr) {
    parts <- strsplit(trim_ws(sub("^>", "", h)), "\\s+")[[1]]
    if (length(parts) < 1) next
    mrna <- parts[1]
    gene <- if (length(parts) >= 2) parts[2] else sub("-T\\d+$", "", mrna)
    gene <- sub("-T\\d+$", "", gene)   # 去除 -T 后缀
    map[mrna] <- gene
  }
  cat("FA loaded:", basename(fa_path), "mrna count:", length(map), "\n")
  map
}

# 读取基因家族列表（兼容 BOM、空格、大小写、下划线等）
read_gene_list <- function(path) {
  df <- try(read.csv(path, header = TRUE, check.names = FALSE,
                     fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE),
            silent = TRUE)
  if (inherits(df, "try-error"))
    df <- read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  names(df) <- trim_ws(names(df))
  target <- c("ExpandedFam","ExpandedGenes","ContractedFam","ContractedGenes")
  canon  <- tolower(gsub("[^a-z]", "", names(df)))
  std    <- tolower(gsub("[^a-z]", "", target))
  rename <- character(length(target))
  for (i in seq_along(target)) {
    hit <- which(canon == std[i])[1]
    if (is.na(hit)) stop("Cannot find column for ", target[i], " in ", basename(path))
    rename[i] <- names(df)[hit]
  }
  df <- df[, rename, drop = FALSE]
  names(df) <- target
  for (col in target) df[[col]] <- ifelse(is.na(df[[col]]), "", as.character(df[[col]]))
  df
}

# GO 聚合（每物种独立）
summarize_go_file <- function(go_path, out_dir = NULL) {
  df <- try(read.table(go_path, sep = "\t", header = FALSE, quote = "",
                       fill = TRUE, comment.char = "", fileEncoding = "UTF-8-BOM"),
            silent = TRUE)
  if (inherits(df, "try-error"))
    df <- read.table(go_path, sep = "\t", header = FALSE, quote = "", fill = TRUE, comment.char = "")
  if (ncol(df) < 14) stop("GO file < 14 columns: ", basename(go_path))
  df$V1  <- as.character(df$V1)
  df$V13 <- as.character(df$V13); df$V13[is.na(df$V13)] <- "-"
  df$V14 <- as.character(df$V14); df$V14[is.na(df$V14)] <- "-"
  
  idx <- split(seq_len(nrow(df)), df$V1)
  res <- data.frame(mrna = names(idx),
                    mrna_function = character(length(idx)),
                    mrna_GO = character(length(idx)),
                    stringsAsFactors = FALSE)
  for (i in seq_along(idx)) {
    rows <- idx[[i]]
    funs <- unique(trim_ws(df$V13[rows])); funs <- funs[funs != "-" & nzchar(funs)]
    gos  <- extract_go_ids(df$V14[rows])
    res$mrna_function[i] <- paste(sort(unique(funs)), collapse = ";")
    res$mrna_GO[i]       <- paste(sort(unique(gos)),  collapse = ";")
  }
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    outf <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(go_path)),
                                      ".mrna_go_summary.tsv"))
    write.table(res, outf, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  res
}

# 合并同一文件内可能重复的家族行
collapse_family_map <- function(fam_col, gene_col) {
  fams  <- trim_ws(as.character(fam_col))
  genes <- as.character(gene_col)
  ok <- nzchar(fams)
  fams <- fams[ok]; genes <- genes[ok]
  if (!length(fams)) return(character(0))
  idx <- split(seq_along(fams), fams)
  vapply(names(idx), function(f) {
    gg <- character(0)
    for (j in idx[[f]]) gg <- c(gg, split_mrnas(genes[j]))
    paste(sort(unique(gg)), collapse = ";")
  }, FUN.VALUE = character(1))
}

# ====================== 主流程 ======================
stopifnot(length(gene_list_files) == length(suffixes),
          length(go_files) == length(suffixes),
          length(fa_files) == length(suffixes))

# 1. 每物种独立的 mrna -> gene 映射
cat("=== Loading mrna -> gene mappings ===\n")
mrna2gene_list <- setNames(lapply(fa_files, extract_mrna_to_gene), suffixes)

# 2. 每物种独立的 GO/功能映射
cat("\n=== Summarizing GO per species ===\n")
mrna2fun_list <- setNames(vector("list", length(suffixes)), suffixes)
mrna2go_list  <- setNames(vector("list", length(suffixes)), suffixes)
for (i in seq_along(suffixes)) {
  go_df <- summarize_go_file(go_files[i], out_dir = go_out_dir)
  mrna2fun_list[[i]] <- setNames(lapply(strsplit(go_df$mrna_function, ";"),
                                        function(x) trim_ws(x[nzchar(x)])), go_df$mrna)
  mrna2go_list[[i]]  <- setNames(lapply(strsplit(go_df$mrna_GO, ";"),
                                        function(x) trim_ws(x[nzchar(x)])), go_df$mrna)
}
names(mrna2fun_list) <- names(mrna2go_list) <- suffixes

# 3. 读取每物种的基因家族列表
cat("\n=== Reading gene family lists ===\n")
species_maps <- setNames(vector("list", length(suffixes)), suffixes)
for (i in seq_along(suffixes)) {
  df <- read_gene_list(gene_list_files[i])
  species_maps[[i]] <- list(
    expanded_map   = collapse_family_map(df$ExpandedFam,   df$ExpandedGenes),
    contracted_map = collapse_family_map(df$ContractedFam, df$ContractedGenes)
  )
}
names(species_maps) <- suffixes

# 所有基因家族并集（行索引）
all_families <- character(0)
for (suf in suffixes) {
  all_families <- c(all_families,
                    names(species_maps[[suf]]$expanded_map),
                    names(species_maps[[suf]]$contracted_map))
}
all_families <- sort(unique(trim_ws(all_families[nzchar(all_families)])))

# 4. 构建合并表
cat("\n=== Building final table ===\n")
combined <- data.frame(GeneFamily = all_families, stringsAsFactors = FALSE)

# 提前定义遗漏统计容器，避免找不到对象
missing_gene_by_sp <- setNames(vector("list", length(suffixes)), suffixes)
missing_go_by_sp   <- setNames(vector("list", length(suffixes)), suffixes)

for (suf in suffixes) {
  exp_flag <- paste0("ExpandedFam_", suf)
  exp_mrna <- paste0("ExpandedGenes_", suf)         # mrna 串
  exp_gene <- paste0("ExpandedGeneNames_", suf)     # 基因名串
  con_flag <- paste0("ContractedFam_", suf)
  con_mrna <- paste0("ContractedGenes_", suf)
  con_gene <- paste0("ContractedGeneNames_", suf)
  
  exp_map <- species_maps[[suf]]$expanded_map
  con_map <- species_maps[[suf]]$contracted_map
  m2g     <- mrna2gene_list[[suf]]
  
  combined[[exp_flag]] <- ifelse(combined$GeneFamily %in% names(exp_map), "Y", "N")
  combined[[con_flag]] <- ifelse(combined$GeneFamily %in% names(con_map), "Y", "N")
  
  combined[[exp_mrna]] <- vapply(
    combined$GeneFamily,
    function(f) { v <- unname(exp_map[f]); if (is.na(v)) "" else v },
    character(1)
  )
  combined[[con_mrna]] <- vapply(
    combined$GeneFamily,
    function(f) { v <- unname(con_map[f]); if (is.na(v)) "" else v },
    character(1)
  )
  
  # Expanded 基因名（在 vapply 闭包中使用 <<- 更新外层统计）
  combined[[exp_gene]] <- vapply(combined[[exp_mrna]], function(ms) {
    mrnas <- split_mrnas(ms)
    if (!length(mrnas)) return("")
    genes <- character(length(mrnas))
    for (j in seq_along(mrnas)) {
      g <- m2g[mrnas[j]]
      if (is.na(g) || length(g) == 0) {
        missing_gene_by_sp[[suf]] <<- c(missing_gene_by_sp[[suf]], mrnas[j])
        genes[j] <- mrnas[j]  # 找不到时用 mrna 兜底
      } else {
        genes[j] <- g
      }
    }
    paste(sort(unique(trim_ws(genes))), collapse = ";")
  }, character(1))
  
  # Contracted 基因名
  combined[[con_gene]] <- vapply(combined[[con_mrna]], function(ms) {
    mrnas <- split_mrnas(ms)
    if (!length(mrnas)) return("")
    genes <- character(length(mrnas))
    for (j in seq_along(mrnas)) {
      g <- m2g[mrnas[j]]
      if (is.na(g) || length(g) == 0) {
        missing_gene_by_sp[[suf]] <<- c(missing_gene_by_sp[[suf]], mrnas[j])
        genes[j] <- mrnas[j]
      } else {
        genes[j] <- g
      }
    }
    paste(sort(unique(trim_ws(genes))), collapse = ";")
  }, character(1))
}

# 5. 跨物种汇总功能与 GO（按物种各自查表再合并，避免跨物种同名混淆）
cat("\n=== Aggregating functions and GO across species ===\n")
all_funcs <- character(nrow(combined))
all_gos   <- character(nrow(combined))

for (i in seq_len(nrow(combined))) {
  fset <- character()
  gset <- character()
  for (suf in suffixes) {
    mrnas <- c(split_mrnas(combined[i, paste0("ExpandedGenes_", suf)]),
               split_mrnas(combined[i, paste0("ContractedGenes_", suf)]))
    if (!length(mrnas)) next
    fun_map <- mrna2fun_list[[suf]]
    go_map  <- mrna2go_list[[suf]]
    for (m in mrnas) {
      f <- fun_map[[m]]
      g <- go_map[[m]]
      if ((is.null(f) || length(f) == 0) && (is.null(g) || length(g) == 0)) {
        # 这里不在闭包中，无需 <<-
        missing_go_by_sp[[suf]] <- c(missing_go_by_sp[[suf]], m)
      }
      if (!is.null(f) && length(f)) fset <- c(fset, f)
      if (!is.null(g) && length(g)) gset <- c(gset, g)
    }
  }
  all_funcs[i] <- paste(sort(unique(trim_ws(fset[nzchar(fset)]))), collapse = ";")
  all_gos[i]   <- paste(sort(unique(trim_ws(gset[nzchar(gset)]))),   collapse = ";")
}
combined$All_mrna_functions <- all_funcs
combined$All_mrna_GOs       <- all_gos

# 6. 输出结果与遗漏统计
write.table(combined, out_combined, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n=== 完成 ===\n")
cat("合并表已写入:", normalizePath(out_combined), "\n")
cat("总基因家族数:", nrow(combined), "\n\n")

for (suf in suffixes) {
  miss_g <- unique(trim_ws(missing_gene_by_sp[[suf]]))
  miss_o <- unique(trim_ws(missing_go_by_sp[[suf]]))
  cat("[", suf, "] 未找到基因名的 mrna 数量:", length(miss_g), "\n", sep = "")
  if (length(miss_g)) cat("  例: ", paste(head(miss_g, 10), collapse = ", "), "\n", sep = "")
  cat("[", suf, "] 完全没有 GO/功能注释的 mrna 数量:", length(miss_o), "\n", sep = "")
  if (length(miss_o)) cat("  例: ", paste(head(miss_o, 10), collapse = ", "), "\n", sep = "")
}

