# 设置输出PDF
pdf("multi_species2.pdf", width=7, height=6)

# 全局参数设置
arrYLim <- c(1e4, 1e9)
arrXLim <- c(1e3, 5e7)
nBSReps <- 30

# ============ 物种配置列表 ============
# 在这里配置所有需要绘制的物种
species_config <- list(
  # 物种1
  list(
    name = "Tra",                          # 物种名称（用于文件名）
    mainFolder = "tra/msmc2ret/",              # 主文件夹
    bsFolder = "tra/bootstrapped/",            # bootstrap文件夹
    mu = 2.45531e-9,                       # 突变率
    color = "yellow",                       # 绘图颜色
    plot_first = TRUE                      # 是否为第一个绘制（决定是plot还是lines）
  ),
  
  # 物种2（示例 - 取消注释并修改参数使用）
   list(
     name = "Tba",
     mainFolder = "tba/msmc2ret/",
     bsFolder = "tba/bootstrapped/",
     mu = 1.51533e-9,
     color = "blue",
     plot_first = FALSE
   ),
  
  # 物种3（示例 - 取消注释并修改参数使用）
   list(
     name = "Tbi",
     mainFolder = "tbi/msmc2ret/",
     bsFolder = "tbi/bootstrapped/",
     mu = 1.62435e-9,
     color = "green",
     plot_first = FALSE
   )
)

# ============ 辅助函数 ============
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

# ============ 主绘图函数 ============
fnPlotPop <- function(sPop, sMainFolder, sBSFolder, nMu, arrCol, 
                      bBSPlot = FALSE, bAppendPlot = TRUE, 
                      sThisBSFile1 = "") {
  
  # 构建文件路径
  sInFile1 <- paste(sMainFolder, "/", sPop, ".final.txt", sep="")
  arrCol <- add.alpha(arrCol, 0.8)
  
  # 如果是bootstrap绘图
  if (bBSPlot) {
    cat("Bootstrap: ", sThisBSFile1, "\n")
    sInFile1 <- sThisBSFile1
    arrCol <- add.alpha(arrCol, 0.2)
  }
  
  # 读取数据
  cat("Open ", sInFile1, "\n")
  if (!file.exists(sInFile1)) {
    warning("File not found: ", sInFile1)
    return()
  }
  
  datMSMC1 <- read.table(sInFile1, header=T, sep="\t")
  datMSMC1$gen <- datMSMC1$left_time_boundary / nMu
  datMSMC1$gen[1] <- 0.01
  datMSMC1$popsize <- (1 / datMSMC1$lambda) / (2 * nMu)
  
  # 绘图
  if (bBSPlot == FALSE && (!bAppendPlot)) {
    plot(datMSMC1$gen, datMSMC1$popsize, type = 's', log='xy', 
         xlim=arrXLim, ylim = arrYLim, 
         xlab="Generations", ylab = "Pop size", 
         lwd=3, col=arrCol[1])
  } else {
    lines(datMSMC1$gen, datMSMC1$popsize, type = 's', 
          lwd=3, col=arrCol[1])
  }
  
  # 绘制bootstrap结果
  if (bBSPlot == TRUE) {
    return()
  } else {
    cat("Plot bootstrap replicates...\n")
    for (nRep in 1:nBSReps) {
      sBSF1 <- paste(sBSFolder, '/', sPop, '/_', nRep, '/out.final.txt', sep="")
      if (file.exists(sBSF1)) {
        fnPlotPop(sPop, sMainFolder, sBSFolder, nMu, arrCol, 
                  bBSPlot = TRUE, sThisBSFile1=sBSF1)
      }
    }
  }
  
  return()
}

# ============ 主程序：循环绘制所有物种 ============
cat("Starting to plot multiple species...\n")
cat("Total species to plot:", length(species_config), "\n\n")

for (i in 1:length(species_config)) {
  species <- species_config[[i]]
  
  cat("=== Processing species", i, ":", species$name, "===\n")
  
  # 绘制该物种
  fnPlotPop(
    sPop = species$name,
    sMainFolder = species$mainFolder,
    sBSFolder = species$bsFolder,
    nMu = species$mu,
    arrCol = species$color,
    bAppendPlot = !species$plot_first
  )
  
  cat("\n")
}

# 添加图例（可选）
# 取消注释以添加图例
# legend("topright", 
#        legend = sapply(species_config, function(x) x$name),
#        col = sapply(species_config, function(x) x$color),
#        lwd = 3,
#        bty = "n")

cat("All species plotted successfully!\n")
dev.off()
