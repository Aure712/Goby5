#setwd("/data2/projects/yshao/population/genetic_load/filter_bothSNP/snp")
dat <- read.table("../allpops.cds.bothDP.bothAF.bed", header=F, sep="\t")
nrow(dat) 
quantile(dat$V8, probs = c(0.05, 0.95), na.rm = TRUE)
quantile(dat$V10, probs = c(0.05, 0.95), na.rm = TRUE)
quantile((dat$V8-dat$V10), probs = c(0.05, 0.95), na.rm = TRUE)
quantile(dat$V8, probs = c(0.1, 0.9), na.rm = TRUE)
quantile(dat$V10, probs = c(0.1, 0.9), na.rm = TRUE)
quantile((dat$V8-dat$V10), probs = c(0.1, 0.9), na.rm = TRUE)
#dat <- dat[dat$V8>532 & dat$V8<=847 & dat$V10>532 & dat$V10<=850 & (dat$V8-dat$V10)>-19 & (dat$V8-dat$V10)<=16,]
dat <- dat[dat$V8>567 & dat$V8<=808 & dat$V10>568 & dat$V10<=810 & (dat$V8-dat$V10)>-14 & (dat$V8-dat$V10)<=8,]
nrow(dat) 
dat$Fixed1 <- (dat$V14==0 | dat$V14==1); #单个数值用||，数组用|
dat$Fixed2 <- (dat$V18==0 | dat$V18==1);

# 筛选出Fixed1和Fixed2相等的行
dat_bothSNP <- dat[(dat$Fixed1 == dat$Fixed2) & !is.na(dat$Fixed1) & !is.na(dat$Fixed2), ]
nrow(dat_bothSNP)

# 计算dat的补集，包括NA行
dat_non_bothSNP <- dat[!(dat$Fixed1 == dat$Fixed2 & !is.na(dat$Fixed1) & !is.na(dat$Fixed2)), ]
nrow(dat_non_bothSNP)

dat1 <- dat_bothSNP[,c(1:3)]
dat2 <- dat_bothSNP[,c(4:6)]
dat3 <- dat_non_bothSNP[,c(1:3)]
dat4 <- dat_non_bothSNP[,c(4:6)]

write.table(dat1, file = "allpops.refTba.NONSYN.filtered.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(dat2, file = "allpops.refTra.NONSYN.filtered.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(dat3, file = "allpops.refTba.NONSYN.discarded.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(dat4, file = "allpops.refTra.NONSYN.discarded.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

