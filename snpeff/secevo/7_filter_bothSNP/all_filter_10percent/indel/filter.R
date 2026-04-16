#setwd("/data2/projects/yshao/population/genetic_load/filter_bothSNP/indel")
dat <- read.table("../allpops.cds.bothDP.bothAF.bed", header=F, sep="\t")
nrow(dat) 
quantile(dat$V8, probs = c(0.1, 0.9), na.rm = TRUE)
quantile(dat$V10, probs = c(0.1, 0.9), na.rm = TRUE)
quantile((dat$V8-dat$V10), probs = c(0.1, 0.9), na.rm = TRUE)
dat <- dat[dat$V8>568 & dat$V8<=805 & dat$V10>569 & dat$V10<=808 & (dat$V8-dat$V10)>-13 & (dat$V8-dat$V10)<=8,]
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

write.table(dat1, file = "allpops.refTba.SV.filtered.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(dat2, file = "allpops.refTra.SV.filtered.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(dat3, file = "allpops.refTba.SV.discarded.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(dat4, file = "allpops.refTra.SV.discarded.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

