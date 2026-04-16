library(clusterProfiler)
library(ggplot2)
library(enrichplot)
library(GOSemSim)
library(DOSE)
library(stringr)
library(RColorBrewer)
library(dplyr)

# 设置全局绘图主题
theme_set(theme_bw())

# 读取数据
datGOMapall <- read.table("/fast3/group_crf/home/g21shaoy23/goby5rna/CEGA_ismc_snpeff_all/combinedGO.tsv", 
                          header=T, sep="\t")
datGOMap <- read.table("./all_combinedGO.tsv", header=T, sep="\t")
datRNA2GeneID <- read.table("./rna2geneid.txt", header=T, sep="\t")

# 数据处理
datGOMap <- merge(datGOMap, datRNA2GeneID, by.x="gene", by.y="rnaid")
datGOMap <- datGOMap[, c(2,3)]
datGOMap <- datGOMap[!duplicated(datGOMap),]
colnames(datGOMap) <- c('term', 'gene')

datGOMapall <- merge(datGOMapall, datRNA2GeneID, by.x="gene", by.y="rnaid")
datGOMapall <- datGOMapall[, c(2,3)]
datGOMapall <- datGOMapall[!duplicated(datGOMapall),]
colnames(datGOMapall) <- c('term', 'gene')

# GO富集分析函数
fnDEG <- function(arrExpandedGenes, arrControlGenes) {
  oEnrichRet <- enricher(
    gene = arrExpandedGenes,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    universe = arrControlGenes,
    minGSSize = 2,
    maxGSSize = 5000,
    qvalueCutoff = 0.2,
    TERM2GENE = datGOMapall,
    TERM2NAME = go2term(unique(datGOMap$term))
  )
  return(oEnrichRet)
}

# 执行富集分析
arrControlGenes <- unique(datGOMapall[, 'gene'])
arrExpandedGenes <- unique(datGOMap[, 'gene'])
enrichResult <- fnDEG(arrExpandedGenes, arrControlGenes)

# ========== 1. 结果筛选 ==========
# 提取富集结果数据框
datOut <- enrichResult@result

# 筛选显著富集的GO terms (p.adjust < 0.1)
datOut_sig <- datOut[datOut$p.adjust < 0.1, ]

# 按p.adjust值排序
datOut_sig <- datOut_sig[order(datOut_sig$p.adjust), ]

# 输出筛选后的结果
write.table(datOut_sig, file="GO_significant.tsv", sep="\t", 
            quote=F, row.names=F, col.names=T)

print(paste("Total GO terms:", nrow(datOut)))
print(paste("Significant GO terms (p.adjust < 0.1):", nrow(datOut_sig)))

# ========== 2. 可视化分析 ==========

# 如果有显著富集的结果，进行可视化
if(nrow(datOut_sig) > 0) {
  
  # 创建富集结果对象（仅包含显著结果）
  enrichResult_sig <- enrichResult
  enrichResult_sig@result <- datOut_sig
  
  # --- 图1: 条形图 (Bar Plot) ---
  p1 <- barplot(enrichResult_sig, showCategory = min(20, nrow(datOut_sig))) + 
    ggtitle("Top 20 Enriched GO Terms") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 12))
  
  ggsave("GO_barplot.pdf", p1, width = 10, height = 8, dpi = 300)
  ggsave("GO_barplot.png", p1, width = 10, height = 8, dpi = 300)
  
  # --- 图2: 点图 (Dot Plot) ---
  p2 <- dotplot(enrichResult_sig, showCategory = min(20, nrow(datOut_sig))) +
    ggtitle("GO Enrichment Dot Plot") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9))
  
  ggsave("GO_dotplot.pdf", p2, width = 10, height = 8, dpi = 300)
  ggsave("GO_dotplot.png", p2, width = 10, height = 8, dpi = 300)
  
  # --- 图3: 基因-概念网络图 (Gene-Concept Network) ---
  if(nrow(datOut_sig) >= 5) {
    p3 <- cnetplot(enrichResult_sig, 
                   showCategory = min(5, nrow(datOut_sig)),
                   foldChange = NULL,
                   circular = FALSE,
                   colorEdge = TRUE) +
      ggtitle("Gene-Concept Network") +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    
    ggsave("GO_cnetplot.pdf", p3, width = 12, height = 10, dpi = 300)
    ggsave("GO_cnetplot.png", p3, width = 12, height = 10, dpi = 300)
  }
  
  # --- 图4: 富集图谱 (Enrichment Map) ---
  if(nrow(datOut_sig) >= 10) {
    enrichResult_sig_em <- pairwise_termsim(enrichResult_sig)
    p4 <- emapplot(enrichResult_sig_em, 
                   showCategory = min(30, nrow(datOut_sig)),
                   layout = "nicely") +
      ggtitle("GO Enrichment Map") +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    
    ggsave("GO_emapplot.pdf", p4, width = 12, height = 10, dpi = 300)
    ggsave("GO_emapplot.png", p4, width = 12, height = 10, dpi = 300)
  }
  
  # --- 图5: 热图 (Heatmap) ---
  if(nrow(datOut_sig) >= 5) {
    p5 <- heatplot(enrichResult_sig, 
                   showCategory = min(20, nrow(datOut_sig)),
                   foldChange = NULL) +
      ggtitle("Gene-GO Term Heatmap") +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.text.y = element_text(size = 9),
            axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
    
    ggsave("GO_heatmap.pdf", p5, width = 12, height = 10, dpi = 300)
    ggsave("GO_heatmap.png", p5, width = 12, height = 10, dpi = 300)
  }
  
  # --- 图6: 自定义气泡图 (Custom Bubble Plot) ---
  # 选择top 20个显著的GO terms
  top_terms <- head(datOut_sig, min(20, nrow(datOut_sig)))
  
  # 计算富集倍数 (Fold Enrichment)
  top_terms$FoldEnrichment <- sapply(strsplit(top_terms$GeneRatio, "/"), function(x) {
    as.numeric(x[1])/as.numeric(x[2])
  }) / sapply(strsplit(top_terms$BgRatio, "/"), function(x) {
    as.numeric(x[1])/as.numeric(x[2])
  })
  
  # 截断过长的描述
  top_terms$Description <- str_wrap(top_terms$Description, width = 50)
  
  # 创建气泡图
  p6 <- ggplot(top_terms, aes(x = FoldEnrichment, 
                              y = reorder(Description, FoldEnrichment),
                              size = Count,
                              color = -log10(p.adjust))) +
    geom_point(alpha = 0.8) +
    scale_color_gradient(low = "blue", high = "red", 
                         name = "-log10(p.adjust)") +
    scale_size_continuous(range = c(3, 12), name = "Gene Count") +
    labs(x = "Fold Enrichment", 
         y = "GO Term",
         title = "GO Enrichment Analysis Results") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 11),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9),
          panel.grid.minor = element_blank())
  
  ggsave("GO_bubble_plot.pdf", p6, width = 10, height = 8, dpi = 300)
  ggsave("GO_bubble_plot.png", p6, width = 10, height = 8, dpi = 300)
  
  # --- 图7: GO三个本体分类图 ---
  # 提取GO本体信息（BP, CC, MF）
  top_terms$Ontology <- "Unknown"
  top_terms$Ontology[grep("^GO:000[0-9]", top_terms$ID)] <- "BP"  # Biological Process
  top_terms$Ontology[grep("^GO:000[5]", top_terms$ID)] <- "CC"    # Cellular Component  
  top_terms$Ontology[grep("^GO:000[3]", top_terms$ID)] <- "MF"    # Molecular Function
  
  # 按本体分类的条形图
  p7 <- ggplot(top_terms, aes(x = reorder(Description, -log10(p.adjust)), 
                              y = -log10(p.adjust),
                              fill = Ontology)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("BP" = "#E41A1C", "CC" = "#377EB8", 
                                 "MF" = "#4DAF4A", "Unknown" = "#999999")) +
    labs(x = "GO Term", 
         y = "-log10(p.adjust)",
         title = "GO Terms by Ontology") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.y = element_text(size = 9),
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 12),
          legend.position = "right")
  
  ggsave("GO_ontology_barplot.pdf", p7, width = 10, height = 8, dpi = 300)
  ggsave("GO_ontology_barplot.png", p7, width = 10, height = 8, dpi = 300)
  
  print("All visualizations have been saved successfully!")
  
} else {
  print("No significant GO terms found (p.adjust < 0.05)")
}

# ========== 3. 生成汇总报告 ==========
# 创建结果汇总
summary_stats <- data.frame(
  Metric = c("Total genes analyzed", 
             "Background genes",
             "Total GO terms tested",
             "Significant GO terms (p < 0.05)",
             "Significant GO terms (p.adjust < 0.1)",
             "Most significant GO term",
             "Most significant p.adjust value"),
  Value = c(length(arrExpandedGenes),
            length(arrControlGenes),
            nrow(datOut),
            sum(datOut$pvalue < 0.05),
            nrow(datOut_sig),
            ifelse(nrow(datOut_sig) > 0, datOut_sig$Description[1], "None"),
            ifelse(nrow(datOut_sig) > 0, format(datOut_sig$p.adjust[1], scientific = TRUE), "NA"))
)

write.table(summary_stats, file="GO_analysis_summary.txt", 
            sep="\t", quote=F, row.names=F, col.names=T)

print("Analysis complete! Check the output files:")
print("- GO_significant.tsv: Filtered significant results")
print("- GO_analysis_summary.txt: Analysis summary statistics")
print("- Various .pdf and .png files: Publication-quality figures")