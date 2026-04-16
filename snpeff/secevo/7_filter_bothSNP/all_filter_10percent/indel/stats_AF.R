#setwd("/data2/projects/yshao/population/genetic_load/filter_bothSNP/indel/")

# 设置输出不使用科学计数法
options(scipen=999)

dat1 <- read.table("refTra.SV.w.indelmodifier.polarized.out.txt", header=F, sep="\t")
dat2 <- read.table("refTba.SV.w.indelmodifier.polarized.out.txt", header=F, sep="\t")

dat1_dedup <- dat1[!duplicated(dat1[, c(1,2)]), ]
dat2_dedup <- dat2[!duplicated(dat2[, c(1,2)]), ]

# 保存去重后的文件
write.table(dat1_dedup, "refTra.dedup.SV.w.indelmodifier.polarized.out.txt", 
            row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

write.table(dat2_dedup, "refTba.dedup.SV.w.indelmodifier.polarized.out.txt", 
            row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

dat1_dedup <- dat1_dedup[,c(1,2,2, 3:ncol(dat1_dedup))]
dat1_dedup[,2] <- dat1_dedup[,2]-1
#dat1_dedup <- dat1_dedup[complete.cases(dat1_dedup), ];
dat1_dedup$AF <- apply((dat1_dedup[, 10:ncol(dat1_dedup)]), 1, FUN= function(x) {non_na_x <- x[!is.na(x)]; nTotal = length(non_na_x)*2; p <- sum(non_na_x)/nTotal; return(p) }) # 去掉NA的个体再统计AF
dat1_dedup <- dat1_dedup[, c(1,2,3,31)]
# 过滤掉AF为NaN的行
#dat1_dedup_NNA <- dat1_dedup[!is.na(dat1_dedup$AF), ]

dat2_dedup <- dat2_dedup[,c(1,2,2, 3:ncol(dat2_dedup))]
dat2_dedup[,2] <- dat2_dedup[,2]-1
#dat2_dedup <- dat2_dedup[complete.cases(dat2_dedup), ];
dat2_dedup$AF <- apply((dat2_dedup[, 10:ncol(dat2_dedup)]), 1, FUN= function(x) {non_na_x <- x[!is.na(x)]; nTotal = length(non_na_x)*2; p <- sum(non_na_x)/nTotal; return(p) }) # 去掉NA的个体再统计AF
dat2_dedup <- dat2_dedup[, c(1,2,3,31)]
# 过滤掉AF为NaN的行
#dat2_dedup <- dat2_dedup[!is.na(dat2_dedup$AF), ]

# 保存为文本文件
write.table(dat1_dedup, file = "allpops.refTra.SV.withAF.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(dat2_dedup, file = "allpops.refTba.SV.withAF.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
