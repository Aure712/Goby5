# 数据
a_del <- 130
a_syn <- 86250
b_del <- 64
b_syn <- 48674

# 比例
a_prop <- a_del / (a_del + a_syn)
b_prop <- b_del / (b_del + b_syn)

cat("a proportion =", a_prop, "\n")
cat("b proportion =", b_prop, "\n")

cat("a deleterious/synonymous =", a_del / a_syn, "\n")
cat("b deleterious/synonymous =", b_del / b_syn, "\n")

# 2×2表
mat <- matrix(c(a_del, a_syn,
                b_del, b_syn),
              nrow = 2, byrow = TRUE)

rownames(mat) <- c("species_a", "species_b")
colnames(mat) <- c("deleterious", "synonymous")

print(mat)

# Fisher精确检验
fisher_res <- fisher.test(mat)
print(fisher_res)

# 两样本比例检验
prop_res <- prop.test(x = c(a_del, b_del),
                      n = c(a_del + a_syn, b_del + b_syn))
print(prop_res)

# 卡方检验
chisq_res <- chisq.test(mat)
print(chisq_res)
