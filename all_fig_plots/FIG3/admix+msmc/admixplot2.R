# 加载必要的包
library(RColorBrewer)

# 读取数据
ta1 <- read.table("Tbi.3.Q")

# 设置美观且色盲友好的颜色方案（使用ColorBrewer）
colors <- brewer.pal(3, "Set2")  # 使用Set2调色板

# 设置行名和排序
rownames(ta1)<-c("Tba9","Tba1","Tba10","Tba2","Tba3","Tba4","Tba5","Tba6","Tba7","Tba8","Tba11","Tba12","Tba13",
                 "Tbi1","Tbi10","Tbi2","Tbi3","Tbi4","Tbi5","Tbi6","Tbi7","Tbi8","Tbi9",
                 "Tra7","Tra1","Tra2","Tra3","Tra4","Tra5","Tra6")

ta1$order <- c(9,1,10,2,3,4,5,6,7,8,11,12,13,
               14,23,15,16,17,18,19,20,21,22,
               30,24,25,26,27,28,29)

ta1 <- ta1[order(ta1$order),]
ta2 <- ta1
ta2$order <- NULL

# 定义群体分组
groups <- list(
  Tra = 1:13,
  Tba = 14:23,
  Tbi = 24:30
)

# 创建分组间隔向量（群体内部间隔0，群体之间间隔1）
space_vec <- rep(0, nrow(ta2))
space_vec[head(groups$Tba,1) + 0.5] <- 1  # Tra-Tba之间
space_vec[head(groups$Tbi,1) + 0.5] <- 1  # Tba-Tbi之间

# 设置图形参数
png("admix-goby.png", width = 10, height = 6, units = "in", res = 300)
par(mar = c(8, 4, 3, 10) + 0.1,  # 增加边距用于标签和图例
    mgp = c(2.5, 0.7, 0),        # 调整坐标轴标签位置
    cex.lab = 1.2)                # 增大坐标轴标签字体

# 绘制admixture图
bp <- barplot(t(as.matrix(ta2)), 
              col = colors,
              space = space_vec,
              xlab = "", 
              ylab = "Ancestry Proportion",
              border = NA,
              xaxt = "n",  # 不显示x轴
              yaxt = "n")   # 不显示y轴

# 添加自定义y轴
axis(2, at = seq(0, 1, 0.2), las = 1, cex.axis = 0.9)

# 添加样本标签（旋转45度）
text(x = bp, 
     y = -0.05, 
     labels = rownames(ta2), 
     srt = 45, 
     adj = c(1, 0.5), 
     xpd = NA, 
     cex = 0.7)

# 添加分组标签
group_pos <- sapply(groups, function(idx) mean(range(bp[idx])))
text(x = group_pos,
     y = -0.15,
     labels = c("Tba", "Tbi", "Tra"),
     xpd = NA,
     font = 2,
     cex = 1.2)

# 添加图例
legend("right", 
       inset = c(-0.17, 0),  # 放置在图形外部
       legend = paste("Ancestry", 1:3), 
       fill = colors, 
       border = NA,
       bty = "n",
       xpd = NA,  # 允许在绘图区域外绘制
       cex = 0.9)

# 添加分隔线
abline(v = c(bp[head(groups$Tba,1)] - 0.5, bp[head(groups$Tbi,1)] - 0.5),
       lty = 2, col = "gray40", lwd = 1.5)

# 添加标题
title("Admixture Analysis (K=3)", cex.main = 1.3, line = 1)

dev.off()