#setwd("/public4/group_crf/home/g21shaoy23/anavar/");
sSrcDir <- "tra_out";
arrRunTypes <- c('', 'sametheta_', 'sametheta_nopolerr_');
arrParamCounts <- c(59, 58, 54);
nBestTolerance <- 0.05;
lsDat <- list();
lsBest <- list();
lsPlotFiles <- list();

options(digits=16)

# ============ 文件信息报告 ============
cat(paste(rep("=", 70), collapse=""), "\n")
cat("文件信息报告\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

cat("工作目录:", getwd(), "\n")
cat("输入目录:", sSrcDir, "\n\n")

allInputFiles <- list()

for( i in 1:length(arrRunTypes) ) {
  sRunType <- arrRunTypes[i];
  if (sRunType == 'full') {
    sRunType <- '';
  }
  nParamCount <- arrParamCounts[i];
  
  arrFiles <- Sys.glob(paste(sSrcDir,'/tra_out_',sRunType, 'out.rep*.txt', sep=""));
  
  cat("--------------------------------------------------\n")
  cat("模型类型:", ifelse(arrRunTypes[i] == '', "完整模型", arrRunTypes[i]), "\n")
  cat("参数数量:", nParamCount, "\n")
  cat("匹配的文件数:", length(arrFiles), "\n")
  if(length(arrFiles) > 0) {
    cat("输入文件列表:\n")
    for(f in arrFiles) {
      cat("  -", f, "\n")
    }
  } else {
    cat("  警告: 未找到匹配的文件!\n")
  }
  cat("\n")
  
  allInputFiles[[i]] <- arrFiles
  
  datAnavar <- NULL;
  for(sFile in arrFiles) {
    datTemp <- read.table(sFile, header=T, skip = 6, sep="\t");
    datTemp$source_file <- sFile;
    datAnavar <- rbind(datAnavar, datTemp);
  }
  
  # 调整gamma类别的顺序
  arrColC1 <- grep('sel_.*_1' , colnames(datAnavar));
  arrColC2 <- grep('sel_.*_2' , colnames(datAnavar));
  arrColC3 <- grep('sel_.*_3' , colnames(datAnavar));
  datCols <- rbind(arrColC1, arrColC2, arrColC3);
  
  for(nRow in 1:nrow(datAnavar)) {
    datVals <- data.frame(); 
    for(nR2 in 1:nrow(datCols) ) {
      datVals <- rbind(datVals, as.numeric(datAnavar[nRow ,datCols[nR2,]]) );
    }
    colnames(datVals) <- c('sel_theta' , 'sel_gamma' , 'sel_e' );
    datVals <- datVals[order(datVals$sel_gamma) , ];
    
    for(nR2 in 1:nrow(datVals) ) {
      datAnavar[nRow, datCols[nR2,] ] <- datVals[nR2, ];
    }
  }
  
  datAnavar <- datAnavar[order(-datAnavar$lnL) , ];
  datBest <- datAnavar[1,];
  lsBest[[i]] <- datBest;
  lsDat[[i]] <- datAnavar;
}

# ============ 输出文件信息 ============
sOutStem <- paste(basename(sSrcDir), ".out", sep="");
sOutTxt <- paste(sOutStem, ".txt", sep="")
sOutPdf <- paste(sOutStem, ".pdf", sep="")

cat(paste(rep("=", 70), collapse=""), "\n")
cat("输出文件:\n")
cat("  文本结果:", file.path(getwd(), sOutTxt), "\n")
cat("  PDF图表:", file.path(getwd(), sOutPdf), "\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

sink(file=sOutTxt);

cat(paste(rep("=", 70), collapse=""), "\n")
cat("ANAVAR结果分析报告\n")
cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

cat("输入文件摘要:\n")
for(i in 1:length(arrRunTypes)) {
  cat("\n模型", i, ":", ifelse(arrRunTypes[i] == '', "完整模型", arrRunTypes[i]), "\n")
  cat("  文件数:", length(allInputFiles[[i]]), "\n")
  if(length(allInputFiles[[i]]) > 0) {
    for(f in allInputFiles[[i]]) {
      cat("  ", f, "\n")
    }
  }
}
cat("\n", paste(rep("=", 70), collapse=""), "\n\n")

cat("最佳拟合结果来源文件:\n")
cat(paste(rep("-", 70), collapse=""), "\n")
for(i in 1:length(arrRunTypes)) {
  cat("\n模型", i, ":", ifelse(arrRunTypes[i] == '', "完整模型", arrRunTypes[i]), "\n")
  cat("  最佳lnL:", lsBest[[i]]$lnL, "\n")
  cat("  来源文件:", as.character(lsBest[[i]]$source_file), "\n")
}
cat("\n", paste(rep("=", 70), collapse=""), "\n\n")

# ============ 增大图形尺寸，改善布局 ============
pdf(file=sOutPdf, width=10, height=7);

# 统计分析
for( i in 1:length(arrRunTypes) ) {
  sRunType1 <- arrRunTypes[i];
  nParamCount1 <- arrParamCounts[i];
  nLnL1 <- lsBest[[i]]$lnL;
  nAIC <- 2*nParamCount1 - 2*nLnL1;
  cat(sRunType1, "AIC=",nAIC, "\n");
}

for( i in 1:length(arrRunTypes) ) {
  sRunType1 <- arrRunTypes[i];
  nParamCount1 <- arrParamCounts[i];
  nLnL1 <- lsBest[[i]]$lnL;

  for( j in 1:length(arrRunTypes) ) {
    if (i == j) {
      next;
    }
    
    sRunType2 <- arrRunTypes[j];
    nParamCount2 <- arrParamCounts[j];
    nLnL2 <- lsBest[[j]]$lnL;
    
    if (nParamCount2 > nParamCount1) {
      next;
    }
    
    nStat <- 2 * (nLnL1 - nLnL2);
    nDF <- nParamCount1 - nParamCount2;
    nP <- 1 - pchisq(nStat,nDF);
    cat(sRunType1, "vs.", sRunType2, "chistat=", nStat , "LRT=", nP, "d.f.=",nDF, "\n");
  }
}

# ============ 优化的画图部分 ============
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("!!! 重要:画图使用的数据来源文件 !!!\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

for( i in 1:length(arrRunTypes) ) {
  sRunType1 <- arrRunTypes[i];
  nParamCount1 <- arrParamCounts[i];
  datAnavar <- lsDat[[i]][ (lsBest[[i]]$lnL - lsDat[[i]]$lnL) < nBestTolerance , ];
  
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat(">>> 图表 ", i, " <<<\n", sep="")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("模型类型:", ifelse(arrRunTypes[i] == '', "完整模型", arrRunTypes[i]), "\n")
  cat("参数数量:", nParamCount1, "\n\n")
  
  cat("【画图使用的主要数据】\n")
  cat("  数据行索引: 第1行 (最佳拟合)\n")
  cat("  lnL值:", datAnavar[1, 'lnL'], "\n")
  cat("  *** 数据来源文件: ", as.character(datAnavar[1, 'source_file']), " ***\n\n", sep="")
  
  cat("【该模型筛选后的数据统计】\n")
  cat("  符合tolerance条件的结果总数:", nrow(datAnavar), "\n")
  cat("  筛选条件: lnL >= ", lsBest[[i]]$lnL - nBestTolerance, "\n", sep="")
  
  if(nrow(datAnavar) > 1) {
    cat("\n  注意: 有", nrow(datAnavar), "个结果符合条件,但画图只使用第1行数据\n")
    cat("  这", nrow(datAnavar), "个结果来自以下文件:\n")
    unique_files <- unique(as.character(datAnavar$source_file))
    for(uf in unique_files) {
      count <- sum(datAnavar$source_file == uf)
      cat("    - ", uf, " (", count, "个结果)\n", sep="")
    }
  }
  
  cat("\n【画图使用的参数值】\n")
  cat("  Gamma 1:", datAnavar$sel_gamma_1[1], "\n")
  cat("  Gamma 2:", datAnavar$sel_gamma_2[1], "\n")
  cat("  Gamma 3:", datAnavar$sel_gamma_3[1], "\n")
  cat("\n")
  
  lsPlotFiles[[i]] <- list(
    model = ifelse(arrRunTypes[i] == '', "完整模型", arrRunTypes[i]),
    file = as.character(datAnavar[1, 'source_file']),
    lnL = datAnavar[1, 'lnL']
  )
  
  theta_props <- datAnavar[, grep('sel_theta_', colnames(datAnavar))];
  theta_props <- theta_props / rowSums(theta_props);
  datAnavar[, c('sel_theta_prop_1','sel_theta_prop_2', 'sel_theta_prop_3')] <- theta_props;
  
  nTotal <- 10000;
  
  n1 <- round(datAnavar[1, 'sel_theta_prop_1'] * nTotal)
  n2 <- round(datAnavar[1, 'sel_theta_prop_2'] * nTotal)
  n3 <- round(datAnavar[1, 'sel_theta_prop_3'] * nTotal)
  
  total_n <- n1 + n2 + n3
  if (total_n != nTotal) {
    max_idx <- which.max(c(n1, n2, n3))
    if (max_idx == 1) n1 <- n1 + (nTotal - total_n)
    if (max_idx == 2) n2 <- n2 + (nTotal - total_n)
    if (max_idx == 3) n3 <- n3 + (nTotal - total_n)
  }
  
  arrG1 <- if (n1 > 0) rep(datAnavar$sel_gamma_1[1], n1) else numeric(0)
  arrG2 <- if (n2 > 0) rep(datAnavar$sel_gamma_2[1], n2) else numeric(0)
  arrG3 <- if (n3 > 0) rep(datAnavar$sel_gamma_3[1], n3) else numeric(0)
  
  all_gammas <- c(
    if (length(arrG1) > 0) arrG1 else numeric(0),
    if (length(arrG2) > 0) arrG2 else numeric(0),
    if (length(arrG3) > 0) arrG3 else numeric(0)
  )
  
  # ============ 优化:更智能的坐标轴范围计算 ============
  if (length(all_gammas) > 0) {
    gamma_range <- max(all_gammas) - min(all_gammas)
    margin <- max(gamma_range * 0.3, 20)  # 至少留20单位的边距
    xlim_range <- c(min(all_gammas) - margin, max(all_gammas) + margin)
  } else {
    xlim_range <- c(0, 1)
  }
  
  # ============ 优化:增加图形边距 ============
  par(mar=c(5, 4, 4, 2) + 0.5)  # 增加边距
  
  # ============ 优化:更清晰的标题 ============
  model_name <- ifelse(sRunType1 == '', "Full Model", sRunType1)
  plot_title <- paste(model_name, " (lnL = ", round(datAnavar[1, 'lnL'], 2), ")", sep="")
  
  # ============ 绘制直方图 ============
  if (length(arrG1) > 0) {
    hist(arrG1, col=rgb(0, 0, 1, 0.6), border="blue", xlim = xlim_range, 
         main=plot_title, xlab="Selection Coefficient (Gamma)", 
         ylab="Frequency", freq = TRUE, cex.main=1.3, cex.lab=1.1)
  } else {
    plot(1, type="n", xlim=xlim_range, ylim=c(0, nTotal/2), 
         main=plot_title, xlab="Selection Coefficient (Gamma)", 
         ylab="Frequency", cex.main=1.3, cex.lab=1.1)
  }
  
  if (length(arrG2) > 0) {
    hist(arrG2, col=rgb(0, 0.5, 0, 0.6), border="darkgreen", add=TRUE, freq = TRUE)
  }
  
  if (length(arrG3) > 0) {
    hist(arrG3, col=rgb(1, 0, 0, 0.6), border="red", add=TRUE, freq = TRUE)
  }
  
  # ============ 修正:使用标准ASCII字符 ============
  # 确定最佳文本位置(避开主要数据区域)
  if (length(all_gammas) > 0) {
    gamma_center <- mean(all_gammas)
    # 如果数据集中在右侧,文本放左侧;反之放右侧
    if (gamma_center > mean(xlim_range)) {
      text_x <- xlim_range[1] + (xlim_range[2] - xlim_range[1]) * 0.05
      text_adj <- c(0, 1)  # 左对齐,顶部对齐
    } else {
      text_x <- xlim_range[2] - (xlim_range[2] - xlim_range[1]) * 0.05
      text_adj <- c(1, 1)  # 右对齐,顶部对齐
    }
  } else {
    text_x <- xlim_range[1] + (xlim_range[2] - xlim_range[1]) * 0.05
    text_adj <- c(0, 1)
  }
  
  # 获取当前y轴范围
  y_range <- par("usr")[3:4]
  text_y <- y_range[2] * 0.95  # 从顶部开始
  
  # ============ 修正:使用标准字符表示 ============
  source_filename <- basename(as.character(datAnavar[1, 'source_file']))
  
  # 构建文本内容(分开主要信息和文件信息)
  main_text <- paste(
    "Gamma Distribution:\n",
    sprintf("%.1f%% : Gamma_1 = %.2f", datAnavar[1,'sel_theta_prop_1'] * 100, datAnavar$sel_gamma_1[1]), "\n",
    sprintf("%.1f%% : Gamma_2 = %.2f", datAnavar[1,'sel_theta_prop_2'] * 100, datAnavar$sel_gamma_2[1]), "\n",
    sprintf("%.1f%% : Gamma_3 = %.2f", datAnavar[1,'sel_theta_prop_3'] * 100, datAnavar$sel_gamma_3[1]),
    sep=""
  )
  
  file_text <- paste("\nSource:\n", source_filename, sep="")
  
  # 绘制主要信息文本框
  text(text_x, text_y, main_text, 
       cex=1.0, adj=text_adj, font=2, family="mono")
  
  # 绘制文件信息(稍小一些)
  file_y <- text_y - (y_range[2] - y_range[1]) * 0.25
  text(text_x, file_y, file_text, 
       cex=0.75, adj=text_adj, font=1, col="gray30")
  
  # ============ 修正:图例也使用标准字符 ============
  # 只为有数据的gamma添加图例
  legend_colors <- c()
  legend_labels <- c()
  legend_border <- c()
  
  if(n1 > 0) {
    legend_colors <- c(legend_colors, rgb(0, 0, 1, 0.6))
    legend_labels <- c(legend_labels, sprintf("Gamma_1 = %.2f", datAnavar$sel_gamma_1[1]))
    legend_border <- c(legend_border, "blue")
  }
  if(n2 > 0) {
    legend_colors <- c(legend_colors, rgb(0, 0.5, 0, 0.6))
    legend_labels <- c(legend_labels, sprintf("Gamma_2 = %.2f", datAnavar$sel_gamma_2[1]))
    legend_border <- c(legend_border, "darkgreen")
  }
  if(n3 > 0) {
    legend_colors <- c(legend_colors, rgb(1, 0, 0, 0.6))
    legend_labels <- c(legend_labels, sprintf("Gamma_3 = %.2f", datAnavar$sel_gamma_3[1]))
    legend_border <- c(legend_border, "red")
  }
  
  # 智能放置图例(与文本相对)
  if (text_adj[1] == 0) {
    legend_pos <- "topright"
  } else {
    legend_pos <- "topleft"
  }
  
  if(length(legend_colors) > 0) {
    legend(legend_pos, legend=legend_labels, 
           fill=legend_colors, border=legend_border,
           bty="n", cex=0.9)
  }
}

cat("\n", paste(rep("=", 70), collapse=""), "\n\n")

sink();
dev.off();

# ============ 画图数据来源摘要 ============
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("【画图数据来源摘要】\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")
for(i in 1:length(lsPlotFiles)) {
  cat("图", i, ":", lsPlotFiles[[i]]$model, "\n")
  cat("  使用文件:", lsPlotFiles[[i]]$file, "\n")
  cat("  lnL值:", lsPlotFiles[[i]]$lnL, "\n\n")
}

cat(paste(rep("=", 70), collapse=""), "\n")
cat("分析完成!\n\n")
cat("请查看以下文件:\n")
cat("  [1] 文本结果: ", sOutTxt, "\n", sep="")
cat("      (包含详细的画图数据来源信息)\n\n")
cat("  [2] PDF图表: ", sOutPdf, "\n", sep="")
cat("      (每个图中标注了数据来源文件)\n\n")
cat(paste(rep("=", 70), collapse=""), "\n")
