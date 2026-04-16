#setwd("/public4/group_crf/home/g21shaoy23/anavar/");

# 设置要分析的文件夹
arrSrcDirs <- c("../tra_out", "../tba_out")  # 两个文件夹
arrPlotNames <- c("Tra_Control", "Tba_Control")    # 对应的图表显示名称（可在此处自由修改）
arrRunTypes <- c('', '')  # 对应的运行类型
arrParamCounts <- c(59, 59)  # 对应的参数数量
nBestTolerance <- 0.05
nPerm <- 10000

options(digits=16)

# 频率分布图辅助函数
plot_hist_compare <- function(x1, x2, xlab, main, col1="blue", col2="red") {
  vals <- c(x1, x2)
  breaks <- pretty(range(vals, na.rm=TRUE), n=10)
  hist(x1, breaks=breaks, col=adjustcolor(col1, 0.4), border=col1,
       xlab=xlab, main=main, freq=TRUE)
  hist(x2, breaks=breaks, col=adjustcolor(col2, 0.4), border=col2,
       add=TRUE, freq=TRUE)
  legend("topright", legend=c(arrPlotNames[1], arrPlotNames[2]),
         col=c(col1, col2), lwd=2, bty="n", cex=0.9)
  grid()
}

plot_hist_multi <- function(series_list, labels, xlab, main, cols) {
  vals <- unlist(series_list)
  breaks <- pretty(range(vals, na.rm=TRUE), n=10)
  hist(series_list[[1]], breaks=breaks, col=adjustcolor(cols[1], 0.4),
       border=cols[1], xlab=xlab, main=main, freq=TRUE)
  if(length(series_list) > 1) {
    for(i in 2:length(series_list)) {
      hist(series_list[[i]], breaks=breaks, col=adjustcolor(cols[i], 0.4),
           border=cols[i], add=TRUE, freq=TRUE)
    }
  }
  legend("topright", legend=labels, col=cols, lwd=2, bty="n", cex=0.9)
  grid()
}

# 创建分析函数
analyzeFolder <- function(sSrcDir, sRunType, nParamCount, folderLabel) {
  cat("\n=================================\n")
  cat("Analyzing folder:", sSrcDir, "\n")
  cat("=================================\n")
  
  # 获取所有文件并提取control编号
  pattern <- paste(sSrcDir, '/', folderLabel, '_out.rep*', folderLabel, '_control_*.txt', sep="")
  arrFiles <- Sys.glob(pattern)
  
  if(length(arrFiles) == 0) {
    cat("No files found with pattern:", pattern, "\n")
    return(NULL)
  }
  
  cat("Found", length(arrFiles), "files\n")
  
  # 提取control编号
  getControlNumber <- function(filename) {
    match <- regmatches(filename, regexpr("control_[0-9]+", filename))
    if(length(match) > 0) {
      return(sub("control_", "", match))
    }
    return(NA)
  }
  
  # 提取rep编号
  getRepNumber <- function(filename) {
    match <- regmatches(filename, regexpr("rep[0-9]+", filename))
    if(length(match) > 0) {
      return(match)
    }
    return(NA)
  }
  
  # 将文件按control分组
  controlNumbers <- sapply(arrFiles, getControlNumber)
  uniqueControls <- unique(controlNumbers[!is.na(controlNumbers)])
  
  # 存储每个control的结果
  controlResults <- list()
  
  # 处理每个control组
  for(control in uniqueControls) {
    # 获取该control组的所有文件
    controlFiles <- arrFiles[controlNumbers == control]
    
    # 创建一个列表来存储每个文件的数据和文件名
    fileDataList <- list()
    
    for(sFile in controlFiles) {
      tempData <- read.table(sFile, header=T, skip = 6, sep="\t")
      tempData$sourceFile <- basename(sFile)
      fileDataList[[length(fileDataList) + 1]] <- tempData
    }
    
    # 合并所有数据
    datAnavar <- do.call(rbind, fileDataList)
    
    # 调整gamma类别的顺序
    arrColC1 <- grep('sel_.*_1' , colnames(datAnavar))
    arrColC2 <- grep('sel_.*_2' , colnames(datAnavar))
    arrColC3 <- grep('sel_.*_3' , colnames(datAnavar))
    datCols <- rbind(arrColC1, arrColC2, arrColC3)
    
    for(nRow in 1:nrow(datAnavar)) {
      datVals <- data.frame()
      for(nR2 in 1:nrow(datCols) ) {
        datVals <- rbind(datVals, as.numeric(datAnavar[nRow ,datCols[nR2,]]) )
      }
      colnames(datVals) <- c('sel_theta' , 'sel_gamma' , 'sel_e' )
      datVals <- datVals[order(datVals$sel_gamma) , ]
      
      for(nR2 in 1:nrow(datVals) ) {
        datAnavar[nRow, datCols[nR2,] ] <- datVals[nR2, ]
      }
    }
    
    # 按似然值排序并获取最佳结果
    datAnavar <- datAnavar[order(-datAnavar$lnL) , ]
    datBest <- datAnavar[1,]
    
    # 计算theta比例
    theta_props <- datBest[, grep('sel_theta_', colnames(datBest))]
    theta_props <- theta_props / sum(theta_props)
    
    # 保存该control的结果
    controlResults[[control]] <- list(
      control = control,
      lnL = datBest$lnL,
      AIC = 2*nParamCount - 2*datBest$lnL,
      theta_prop_1 = theta_props[[1]],
      theta_prop_2 = theta_props[[2]], 
      theta_prop_3 = theta_props[[3]],
      gamma_1 = datBest$sel_gamma_1,
      gamma_2 = datBest$sel_gamma_2,
      gamma_3 = datBest$sel_gamma_3,
      e_1 = datBest$sel_e_1,
      e_2 = datBest$sel_e_2,
      e_3 = datBest$sel_e_3,
      bestFileName = datBest$sourceFile
    )
  }
  
  # 转换为数据框
  plotData <- do.call(rbind, lapply(controlResults, function(x) {
    data.frame(
      control = as.numeric(x$control),
      gamma_1 = x$gamma_1,
      gamma_2 = x$gamma_2,
      gamma_3 = x$gamma_3,
      theta_prop_1 = x$theta_prop_1,
      theta_prop_2 = x$theta_prop_2,
      theta_prop_3 = x$theta_prop_3,
      lnL = x$lnL,
      AIC = x$AIC,
      bestFileName = x$bestFileName,
      folder = sSrcDir  # 添加文件夹标识
    )
  }))
  
  # 排序
  plotData <- plotData[order(plotData$control), ]
  
  return(plotData)
}

# 分析两个文件夹
results_tra <- analyzeFolder(arrSrcDirs[1], arrRunTypes[1], arrParamCounts[1], "tra_out")
results_tba <- analyzeFolder(arrSrcDirs[2], arrRunTypes[2], arrParamCounts[2], "tba_out")

# 确保两个结果都存在
if(is.null(results_tra) || is.null(results_tba)) {
  stop("One or both folders did not yield results. Please check file paths and patterns.")
}

# 经验法检验
cat("\n=================================\n")
cat("Empirical Paired Test\n")
cat("=================================\n")

# 获取共同的control编号
common_controls <- intersect(results_tra$control, results_tba$control)
cat("Common controls found:", length(common_controls), "\n")
cat("Controls:", paste(common_controls, collapse=", "), "\n\n")

if(length(common_controls) > 0) {
  # 确保顺序一致
  tra_gamma1_props <- results_tra$theta_prop_1[match(common_controls, results_tra$control)]
  tba_gamma1_props <- results_tba$theta_prop_1[match(common_controls, results_tba$control)]
  
  # 配对差异
  differences <- tra_gamma1_props - tba_gamma1_props
  
  # 经验置换法（随机符号翻转）
  obs_stat <- mean(differences)
  perm_stats <- replicate(nPerm, {
    mean(differences * sample(c(-1, 1), length(differences), replace=TRUE))
  })
  empirical_p <- (sum(abs(perm_stats) >= abs(obs_stat)) + 1) / (nPerm + 1)
  
  empirical_test <- list(
    statistic = obs_stat,
    p.value = empirical_p,
    method = "Empirical paired permutation test (mean difference)",
    nPerm = nPerm
  )
  
  cat("Testing Gamma 1 proportions between", arrSrcDirs[1], "and", arrSrcDirs[2], "\n")
  cat("Number of paired samples:", length(common_controls), "\n")
  cat("Test statistic (mean difference) =", empirical_test$statistic, "\n")
  cat("Empirical P-value =", empirical_test$p.value, "\n")
  cat("Method:", empirical_test$method, "\n")
  cat("Permutations:", empirical_test$nPerm, "\n\n")
  
  # 计算差异统计
  cat("Summary of differences (", arrPlotNames[1], " - ", arrPlotNames[2], "):\n", sep="")
  cat("  Mean difference:", mean(differences), "\n")
  cat("  Median difference:", median(differences), "\n")
  cat("  SD of differences:", sd(differences), "\n")
  cat("  Range:", min(differences), "to", max(differences), "\n\n")
  
  # 显示配对数据
  paired_data <- data.frame(
    Control = common_controls,
    tra_in_gamma1_prop = round(tra_gamma1_props * 100, 1),
    tba_in_gamma1_prop = round(tba_gamma1_props * 100, 1),
    Difference = round(differences * 100, 1)
  )
  print(paired_data)
} else {
  cat("No common controls found between folders. Cannot perform empirical test.\n")
  empirical_test <- NULL
  paired_data <- NULL
}

# 输出文本结果
sOutStem <- "combined_analysis"
sink(file=paste(sOutStem, ".txt", sep=""))

cat("COMBINED ANALYSIS RESULTS\n")
cat("========================\n\n")

# tra_in 结果
cat("FOLDER:", arrSrcDirs[1], "\n")
cat("========================\n")
for(i in 1:nrow(results_tra)) {
  cat("Control", results_tra$control[i], ":\n")
  cat("  Best file:", results_tra$bestFileName[i], "\n")
  cat("  AIC =", results_tra$AIC[i], "\n")
  cat("  lnL =", results_tra$lnL[i], "\n")
  cat("  Gamma values:", results_tra$gamma_1[i], results_tra$gamma_2[i], results_tra$gamma_3[i], "\n")
  cat("  Theta proportions:", 
      round(results_tra$theta_prop_1[i]*100), "%,",
      round(results_tra$theta_prop_2[i]*100), "%,",
      round(results_tra$theta_prop_3[i]*100), "%\n")
  cat("\n")
}

# tba_in 结果
cat("\nFOLDER:", arrSrcDirs[2], "\n")
cat("========================\n")
for(i in 1:nrow(results_tba)) {
  cat("Control", results_tba$control[i], ":\n")
  cat("  Best file:", results_tba$bestFileName[i], "\n")
  cat("  AIC =", results_tba$AIC[i], "\n")
  cat("  lnL =", results_tba$lnL[i], "\n")
  cat("  Gamma values:", results_tba$gamma_1[i], results_tba$gamma_2[i], results_tba$gamma_3[i], "\n")
  cat("  Theta proportions:", 
      round(results_tba$theta_prop_1[i]*100), "%,",
      round(results_tba$theta_prop_2[i]*100), "%,",
      round(results_tba$theta_prop_3[i]*100), "%\n")
  cat("\n")
}

# 经验法检验结果
cat("\nEMPIRICAL PAIRED TEST RESULTS\n")
cat("==================================\n")
if(!is.null(empirical_test)) {
  cat("Comparing Gamma 1 proportions between", arrSrcDirs[1], "and", arrSrcDirs[2], "\n")
  cat("Number of paired samples:", length(common_controls), "\n")
  cat("Test statistic (mean difference) =", empirical_test$statistic, "\n")
  cat("Empirical P-value =", empirical_test$p.value, "\n")
  cat("Method:", empirical_test$method, "\n")
  cat("Permutations:", empirical_test$nPerm, "\n")
  cat("Interpretation: ")
  if(empirical_test$p.value < 0.05) {
    cat("Significant difference at α = 0.05\n")
  } else {
    cat("No significant difference at α = 0.05\n")
  }
  cat("\nPaired Data:\n")
  print(paired_data)
} else {
  cat("Test could not be performed due to lack of common controls.\n")
}

sink()

# 创建比较图
pdf(file=paste(sOutStem, "_comparison.pdf", sep=""), width=14, height=10)
par(mfrow=c(2,3))

# 图1: 两个文件夹的Gamma 1分布
plot_hist_compare(results_tra$gamma_1, results_tba$gamma_1,
                  xlab="Gamma 1", main="Gamma 1 Distribution")

# 图2: 两个文件夹的Theta Proportion 1分布
plot_hist_compare(results_tra$theta_prop_1, results_tba$theta_prop_1,
                  xlab="Theta Proportion 1", main="Theta Proportion 1 Distribution")

# 图3: 配对差异分布（仅对共同controls）
if(length(common_controls) > 0) {
  hist(differences, breaks=10, col='plum', border="white",
       xlab="Difference in Theta Prop 1", 
       main=paste("Paired Differences (", arrPlotNames[1], " - ", arrPlotNames[2], ")", sep=""),
       freq=TRUE)
  abline(v=0, lty=2, col='gray')
  abline(v=mean(differences), lty=3, col='blue')
  mtext(paste("Mean =", round(mean(differences), 4)), side=3, line=0, cex=0.8)
  grid()
} else {
  plot.new()
  text(0.5, 0.5, "No common controls for paired comparison", cex=1.2)
}

# 图4: AIC分布比较
plot_hist_compare(results_tra$AIC, results_tba$AIC,
                  xlab="AIC", main="AIC Distribution")

# 图5: 箱型图比较Gamma 1 proportions
if(length(common_controls) > 0) {
  boxplot(list(tra_in = results_tra$theta_prop_1,
               tba_in = results_tba$theta_prop_1),
          col = c("lightblue", "lightcoral"),
          main = "Distribution of Gamma 1 Proportions",
          ylab = "Theta Proportion 1",
          names = arrPlotNames)
  
  # 添加经验法检验的p值
  if(!is.null(empirical_test)) {
    mtext(paste("Empirical p-value =", round(empirical_test$p.value, 4)), 
          side = 3, line = 0, cex = 0.9)
  }
  grid()
} else {
  plot.new()
  text(0.5, 0.5, "No data for boxplot comparison", cex=1.2)
}

# 图6: 所有Gamma值的综合分布比较
tra_gamma_all <- c(results_tra$gamma_1, results_tra$gamma_2, results_tra$gamma_3)
tba_gamma_all <- c(results_tba$gamma_1, results_tba$gamma_2, results_tba$gamma_3)
plot_hist_compare(tra_gamma_all, tba_gamma_all,
                  xlab="Gamma", main="All Gamma Values Distribution")

dev.off()

# 创建详细的配对分析图（如果有共同controls）
if(length(common_controls) > 0) {
  pdf(file=paste(sOutStem, "_paired_analysis.pdf", sep=""), width=12, height=10)
  par(mfrow=c(2,2))
  
  # 图1: 两组Gamma 1比例分布
  plot_hist_compare(tra_gamma1_props, tba_gamma1_props,
                    xlab="Gamma 1 Proportion", main="Gamma 1 Proportion Distribution")
  
  # 图2: 平均值分布
  mean_props <- (tra_gamma1_props + tba_gamma1_props) / 2
  hist(mean_props, breaks=10, col='lightgreen', border="white",
       xlab="Mean of Two Measurements",
       main="Mean Proportion Distribution",
       freq=TRUE)
  abline(v=mean(mean_props), lty=2, col='blue')
  grid()
  
  # 图3: 差异的直方图
  diff_props <- tra_gamma1_props - tba_gamma1_props
  hist(diff_props, breaks=10, col='lightblue',
       main="Distribution of Differences",
       xlab="Difference in Gamma 1 Proportion",
       ylab="Frequency",
       freq=TRUE)
  abline(v=0, lty=2, col='red')
  abline(v=mean(diff_props), lty=2, col='blue')
  
  # 添加正态性检验
  shapiro_test <- shapiro.test(diff_props)
  text(min(diff_props), par("usr")[4]*0.9,
       paste("Shapiro-Wilk p =", round(shapiro_test$p.value, 4)),
       pos=4, cex=0.9)
  grid()
  
  # 图4: 配对差异的条形图
  barplot(diff_props, names.arg=common_controls,
          col=ifelse(diff_props > 0, "lightgreen", "lightcoral"),
          main="Individual Differences by Control",
          xlab="Control", ylab=paste("Difference (", arrPlotNames[1], " - ", arrPlotNames[2], ")", sep=""),
          las=2)
  abline(h=0, lty=1, col='black')
  
  # 添加经验法检验结果
  if(!is.null(empirical_test)) {
    mtext(paste("Empirical test: p =", round(empirical_test$p.value, 4),
                ifelse(empirical_test$p.value < 0.05, "(Significant)", "(Not significant)")),
          side=3, line=-2, cex=0.9, col='darkblue')
  }
  
  dev.off()
}

# 为每个文件夹单独生成详细图表
for(i in 1:2) {
  if(i == 1) {
    plotData <- results_tra
    folderName <- arrSrcDirs[1]
    plotName <- arrPlotNames[1]
  } else {
    plotData <- results_tba
    folderName <- arrSrcDirs[2]
    plotName <- arrPlotNames[2]
  }
  
  if(is.null(plotData)) next
  
  sOutStem_folder <- paste(folderName, ".out", sep="")
  
  # 基本分布图
  pdf(file=paste(sOutStem_folder, ".pdf", sep=""), width=10, height=8)
  par(mfrow=c(2,2))
  
  # Gamma值分布
  plot_hist_multi(
    list(plotData$gamma_1, plotData$gamma_2, plotData$gamma_3),
    c("Gamma 1", "Gamma 2", "Gamma 3"),
    xlab="Gamma",
    main=paste("Gamma Values Distribution -", plotName),
    cols=c("blue", "darkgreen", "red")
  )
  
  # Theta比例分布
  plot_hist_multi(
    list(plotData$theta_prop_1, plotData$theta_prop_2, plotData$theta_prop_3),
    c("Theta 1", "Theta 2", "Theta 3"),
    xlab="Theta Proportion",
    main=paste("Theta Proportions Distribution -", plotName),
    cols=c("blue", "darkgreen", "red")
  )
  
  # AIC值分布
  hist(plotData$AIC, breaks=10, col='purple', border="white",
       xlab="AIC", main=paste("AIC Distribution -", plotName), freq=TRUE)
  grid()
  
  # 对数似然值分布
  hist(plotData$lnL, breaks=10, col='orange', border="white",
       xlab="Log Likelihood", main=paste("Log Likelihood Distribution -", plotName), freq=TRUE)
  grid()
  
  dev.off()
}

cat("\nAnalysis complete. Results saved to:\n")
cat("  ", paste(sOutStem, ".txt", sep=""), " - Combined text output with empirical test\n")
cat("  ", paste(sOutStem, "_comparison.pdf", sep=""), " - Comparison plots\n")
if(length(common_controls) > 0) {
  cat("  ", paste(sOutStem, "_paired_analysis.pdf", sep=""), " - Detailed paired analysis\n")
}
cat("  ", paste(arrSrcDirs[1], ".out.pdf", sep=""), " - Individual plots for", arrSrcDirs[1], "\n")
cat("  ", paste(arrSrcDirs[2], ".out.pdf", sep=""), " - Individual plots for", arrSrcDirs[2], "\n")

