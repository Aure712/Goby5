.libPaths("/data/projects/rcui/R/x86_64-pc-linux-gnu-library/4.1")
setwd("/data3/projects/yshao/cafe/cafe1/barbel_tridentiger/");
library(org.Dr.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)


sInDir <- "/data3/projects/yshao/cafe/cafe1/runresults_k6/1/";
#sRefSp <- sSp <- "MacropodusOpercularis";
sRefSp <- sSp <- "TridentigerBarbatus"
sSp <- "X.20."
nFDRCutoff <- 0.2
nPvalueCutoff <- 0.05
sCategory <- 'BP' ; #biological process
datPGID2Ref <- read.table("/data3/projects/yshao/Goby5RNA/go2/goby5/rundir2/orthofinder/Results_Jan13/Orthogroups/Orthogroups.tsv", header = T, stringsAsFactors = F, fill=T, sep="\t")
datFam2Gene <- read.table("/data3/projects/yshao/mcmctree/1_genespace/rundir/orthofinder/Results_Jan24/Phylogenetic_Hierarchical_Orthogroups/N0.tsv", header=T, sep="\t", stringsAsFactors = F);
datFam2Gene$OG <- gsub('N0.', '', datFam2Gene$OG , fixed = T)

datChange <- read.table(paste0(sInDir, "/Gamma_change.tab" ), sep="\t", header=T);
datChange2 <- datChange[, c(1, grep(sSp, colnames(datChange), fixed = T)) ];

datBranchProb <- read.table(paste0(sInDir, "/Gamma_branch_probabilities.tab" ), sep="\t", header=T, comment.char = '/', fill = T);

datChange3 <- merge(datChange2 , datBranchProb[, c(1,grep(sSp, colnames(datBranchProb)))], by.x = "FamilyID", by.y=1)

datChange4 <- datChange3[datChange3[,3] < nPvalueCutoff, ]
#datChange$fdr <- p.adjust(datChange[, 3])

arrExpandedFam <- datChange4[datChange4[,2] > 0  ,1]
arrContractedFam <- datChange4[datChange4[,2] < 0  ,1]
arrControlFam <- datChange4[ ,1]

fnFam2Gene <- function(arr) {
  arrGenes <- datFam2Gene[datFam2Gene$OG %in% arr, sRefSp];
  return(unique(str_trim(unlist(strsplit(x = arrGenes, split = ',')))));
}

arrExpandedGenes <- fnFam2Gene(arrExpandedFam);
arrContractedGenes <- fnFam2Gene(arrContractedFam);
arrControlGenes <- fnFam2Gene(arrControlFam)

