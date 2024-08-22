setwd("D:\\生信分析数据\\12 免疫浸润")
rm(list=ls())

# 加载需要的R包 -----------------------------------------------------------------
library(ggcorrplot)
library(tidyr)


#IS
# 导入数据 --------------------------------------------------------------------
DEG_expr <- read.csv("IS array.csv",row.names = 1)#行名基因名，列名样本名
# cibersort运算结果
ssgsea <- read.csv("ssgsea_IS_std.csv",row.names = 1)#行名样本名，列名免疫细胞浸润评分
#数据格式整理
data <- as.data.frame(t(DEG_expr))#数据转置，行名是样本名，列名基因名
ssgsea <- as.data.frame(t(ssgsea))#数据转置，行名是样本名，列名基因名
#check一下
identical(rownames(data),rownames(ssgsea))
data_to_calculate <- data[, c("TM9SF2", "CCL8")]

# 一种简单的方法 （感受R包的魅力）----------------------------------------------
# 计算相关系数矩阵
corr_mat <- cor(data_to_calculate, ssgsea, method = "spearman")
# 绘制相关性矩阵图
plot<-ggcorrplot(t(corr_mat), 
                 # 显示颜色图例
                 show.legend = T, 
                 # 设置相关系数矩阵的颜色，其中第一个为最小值对应的颜色，
                 # 中间的白色为0对应的颜色，最后一个为最大值对应的颜色
                 colors = c("#2166AC", "white", "#B2182B"), 
                 # 设置数字显示的位数
                 digits = 2, 
                 # 显示变量标签（默认为TRUE）
                 lab = T) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))#修改横轴标签方向
ggsave("基因相关性_IS.pdf", plot = plot, width = 15, height = 3)


#OSA
DEG_expr <- read.csv("OSA array.csv",row.names = 1)#行名基因名，列名样本名
# cibersort运算结果
ssgsea <- read.csv("ssgsea_OSA_std.csv",row.names = 1)#行名样本名，列名免疫细胞浸润评分
#数据格式整理
data <- as.data.frame(t(DEG_expr))#数据转置，行名是样本名，列名基因名
ssgsea <- as.data.frame(t(ssgsea))#数据转置，行名是样本名，列名基因名
#check一下
identical(rownames(data),rownames(ssgsea))
data_to_calculate <- data[, c("TM9SF2", "CCL8")]

# 一种简单的方法 （感受R包的魅力）----------------------------------------------
# 计算相关系数矩阵
corr_mat <- cor(data_to_calculate, ssgsea, method = "spearman")
# 绘制相关性矩阵图
plot<-ggcorrplot(t(corr_mat), 
                 # 显示颜色图例
                 show.legend = T, 
                 # 设置相关系数矩阵的颜色，其中第一个为最小值对应的颜色，
                 # 中间的白色为0对应的颜色，最后一个为最大值对应的颜色
                 colors = c("#2166AC", "white", "#B2182B"), 
                 # 设置数字显示的位数
                 digits = 2, 
                 # 显示变量标签（默认为TRUE）
                 lab = T) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))#修改横轴标签方向
ggsave("基因相关性_OSA.pdf", plot = plot, width = 15, height = 3)


#另一种复杂方法
#数据简单处理
# Shapiro-Wilk检验
shapiro.test(data$TM9SF2)#适用于小样本数据，N≤50
# Kolmogorov-Smirnov检验(柯尔莫哥洛夫检验)
ks.test(data$TM9SF2, pnorm, mean = mean(data$TM9SF2), sd = sd(data$TM9SF2))#适用于大样本数据，N＞50

# 计算相关性 -------------------------------------------------------------------
colnames(data)
colnames(ssgsea)
target_Gene <- "TM9SF2"#确定目标基因
expr <- as.numeric(data[,target_Gene])
cor.test(ssgsea$'Activated CD8 T cell', expr, use = "everything", method = "spearman")

# 封装一个函数来快速计算 -------------------------------------------------------------
# spearman相关性分析(非正态)
calculate_correlation <- function(Gene_expr, ssgsea_result) {
  cor_matrix <- data.frame("Gene" = character(),"im_cell" = character(),"Cor" = numeric(),"p-value" = numeric(), stringsAsFactors = FALSE)
  
  for (i in 1:ncol(ssgsea_result)) {
    result <- cor.test(Gene_expr, ssgsea_result[, i], method = "spearman")
    new_row <- data.frame("Gene" = "Gene", "im_cell" = colnames(ssgsea_result)[i], "Cor" = result$estimate, "p-value" = result$p.value)
    cor_matrix <- rbind(cor_matrix, new_row)
  }
  return(cor_matrix)
}
# pearson相关性分析(正态)
calculate_correlation <- function(Gene_expr, ssgsea_result) {
  cor_matrix <- data.frame("Gene" = character(),"im_cell" = character(),"Cor" = numeric(),"p-value" = numeric(), stringsAsFactors = FALSE)
  
  for (i in 1:ncol(ssgsea_result)) {
    result <- cor.test(Gene_expr, ssgsea_result[, i], method = "pearson")
    new_row <- data.frame("Gene" = "Gene", "im_cell" = colnames(ssgsea_result)[i], "Cor" = result$estimate, "p-value" = result$p.value)
    cor_matrix <- rbind(cor_matrix, new_row)
  }
  return(cor_matrix)
}

# 计算5个基因与免疫细胞浸润分数的相关性 -----------------------------------------------------
# 选择要计算相关性的列
data_to_calculate <- data[, c("TM9SF2", "CCL8")]
# 新建一个空的数据框保存结果
results <- data.frame("Gene" = character(),"im_cell" = character(),"Cor" = numeric(),"p-value" = numeric(), stringsAsFactors = FALSE)
# 使用for循环遍历数据框中的每一列，并计算相关性
for (i in 1:ncol(data_to_calculate)) {
  print(i)
  gene_expr <- data_to_calculate[, i]
  corr_result <- calculate_correlation(gene_expr, ssgsea)
  # 将每次计算的结果添加到新的数据框中
  results <- rbind(results, corr_result)
}
# 查看列名，下方修改的依据
colnames(data_to_calculate)
# 手动修改每一个基因名称
results$Gene <- c(rep("TM9SF2", 28),
                  rep("CCL8", 28))

results <- results[,c(1,2,3)]#转换数据，这里不需要显著性p值，也可以通过显著性P值进行筛选
colnames(results)
# 使用spread函数将长数据转换为宽格式
results_wide <- spread(results, key = im_cell, value = Cor)
# 将基因名设置为行名
rownames(results_wide) <- results_wide$Gene
# 删除原始数据框中的基因名列
results_wide$Gene <- NULL

# 绘制相关性矩阵图
plot<-ggcorrplot(t(results_wide), 
           # 显示颜色图例
           show.legend = T, 
           # 设置相关系数矩阵的颜色，其中第一个为最小值对应的颜色，
           # 中间的白色为0对应的颜色，最后一个为最大值对应的颜色
           colors = c("#2166AC", "white", "#B2182B"), 
           # 设置数字显示的位数
           digits = 2, 
           # 显示变量标签（默认为TRUE）
           lab = T)
ggsave("基因相关性_IS.pdf", plot = plot, width = 15, height = 6)
