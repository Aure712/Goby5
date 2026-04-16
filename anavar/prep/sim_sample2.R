# ========== 参数设置区 ==========
# 原始数组
original_array <- c(0,176,13,8,5,3,0,2,0,3,3,0,1,0,3,0,1,1,1,0,0,1,0,1,0,8,82)

# 抽样次数
sample_size <- 312

# 循环次数
iterations <- 100

# 输出文件名
output_file <- "tba_pureFSR_sel_m_sampling_results.txt"

# 随机种子
# 设为固定整数可保证结果可复现；如果不想固定结果，可改为 NULL
random_seed <- 20260401

# ========== 主程序 ==========

# 设置随机种子以保证可重复性
if (!is.null(random_seed)) {
  set.seed(random_seed)
  cat("已设置随机种子：", random_seed, "\n")
} else {
  cat("未设置随机种子，本次结果不可复现\n")
}

# 步骤1: 创建扩展数组
# 根据原数组的值和序号创建新数组
cat("正在创建扩展数组...\n")
expanded_array <- c()

# 对原数组的每个元素进行处理
for (i in 1:length(original_array)) {
  # 填充original_array[i]个i值
  expanded_array <- c(expanded_array, rep(i, original_array[i]))
}

# 检查扩展数组的长度
cat("扩展数组长度：", length(expanded_array), "\n")
cat("原数组元素总和：", sum(original_array), "\n\n")

# 步骤2: 执行多次抽样并记录结果
cat("开始执行", iterations, "次抽样...\n")

# 创建结果矩阵，用于存储每次的统计结果
results_matrix <- matrix(0, nrow = iterations, ncol = length(original_array))

# 进度显示设置
progress_interval <- max(1, floor(iterations / 10))

# 执行循环抽样
for (iter in 1:iterations) {
  # 显示进度
  if (iter %% progress_interval == 0 || iter == iterations) {
    cat("进度：", iter, "/", iterations, " (", round(iter / iterations * 100), "%)\n", sep = "")
  }
  
  # 步骤3: 使用sample函数进行抽样
  # 参数replace=TRUE表示有放回抽样
  sampled_data <- sample(expanded_array, size = sample_size, replace = TRUE)
  
  # 统计各个数字出现的次数
  # 初始化计数数组
  count_array <- rep(0, length(original_array))
  
  # 统计每个数字的出现次数
  for (j in 1:length(original_array)) {
    count_array[j] <- sum(sampled_data == j)
  }
  
  # 将结果存储到结果矩阵中
  results_matrix[iter, ] <- count_array
}

cat("抽样完成！\n\n")

# 步骤4: 输出结果到文本文件（只输出第2个到倒数第2个元素）
cat("正在写入结果文件：", output_file, "\n")
cat("注意：输出仅包含第2个到倒数第2个元素\n")

# 创建输出字符串向量
output_lines <- c()

# 可选：把随机种子信息写入文件头部
if (!is.null(random_seed)) {
  output_lines <- c(output_lines, paste0("# random_seed = ", random_seed))
} else {
  output_lines <- c(output_lines, "# random_seed = NULL")
}

# 确定输出的列范围（第2列到倒数第2列）
start_col <- 2
end_col <- ncol(results_matrix) - 1

# 格式化每行结果
for (i in 1:iterations) {
  # 只选择第2个到倒数第2个元素，用逗号和空格连接
  line <- paste(results_matrix[i, start_col:end_col], collapse = ", ")
  output_lines <- c(output_lines, line)
}

# 写入文件
writeLines(output_lines, output_file)

cat("结果已保存到文件：", output_file, "\n\n")

# ========== 统计分析（可选） ==========
cat("=============== 统计摘要 ===============\n")
cat("原数组长度：", length(original_array), "\n")
cat("扩展数组总长度：", length(expanded_array), "\n")
cat("每次抽样数量：", sample_size, "\n")
cat("总循环次数：", iterations, "\n")
cat("输出元素范围：第", start_col, "个到第", end_col, "个元素\n")
if (!is.null(random_seed)) {
  cat("随机种子：", random_seed, "\n\n")
} else {
  cat("随机种子：NULL（结果不可复现）\n\n")
}

# 计算每个位置的平均抽样次数（全部位置）
mean_counts <- colMeans(results_matrix)
cat("各位置平均抽样次数（全部", length(mean_counts), "个位置）：\n")
for (i in 1:length(mean_counts)) {
  cat("位置", i, ": ", round(mean_counts[i], 2),
      " (理论值: ", round(original_array[i] * sample_size / sum(original_array), 2), ")\n", sep = "")
}

# 显示前5行结果作为示例（只显示输出的部分）
cat("\n前5行输出结果示例（第", start_col, "到第", end_col, "个元素）：\n")
for (i in 1:min(5, iterations)) {
  cat("第", i, "次: ", paste(results_matrix[i, start_col:end_col], collapse = ", "), "\n", sep = "")
}

cat("\n程序执行完成！\n")
cat("输出文件包含", end_col - start_col + 1, "个元素（原数组去掉首尾）\n")
