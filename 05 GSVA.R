setwd("D:\\生信分析数据\\05 GSVA")
rm(list=ls())

library(msigdbr)
library(dplyr)
library(data.table)
## 
## Attaching package: 'data.table'
## The following objects are masked from 'package:dplyr':
## 
##     between, first, last
library(GSVA)
library(limma)
library(stringr)
library(ggplot2)
#查看msigdbr包里自带的物种
msigdbr_show_species()
h <- msigdbr(species = "Homo sapiens", # 物种拉丁名
             category = "H") #此处以hallmark为例，你也可以选择MSigDB的其他注释

# 示例数据表达矩阵的基因名是gene symbol，这里就选gene_symbol。
# 如果你的表达矩阵以ENTREZ ID作为基因名，就把下面这段的gene_symbol换成entrez_gene
h <- select(h, gs_name, gene_symbol) %>% #或entrez_gene
  as.data.frame %>% 
  split(., .$gs_name) %>% 
  lapply(., function(x)(x$gene_symbol)) #或entrez_gene

# 在每个geneset里面去掉重复的基因
gs <- lapply(h, unique)

# 接下来去掉那些在两个或更多个pathways里出现过的genes
count <- table(unlist(gs))
keep <- names(which(table(unlist(gs)) < 2))
gs <- lapply(gs, function(x) intersect(keep, x))

# 过滤之后，很多pathway一个gene都不剩了，去掉这些
gs <- gs[lapply(gs, length) > 0]

# 预览过滤后的结果
head(gs)


#IS
# 读入基因表达矩阵
a<-read.csv("IS array.csv",row.names = 1)
# 这一句就完成了GSVA分析
gsva_es <- gsva(as.matrix(a), gs)
# 把通路的表达量保存到文件
write.csv(gsva_es, "gsva_output_IS.csv", quote = F)
# 分组
group_list <- data.frame(sample = colnames(gsva_es), group = c(rep("Control", 23), rep("IS", 23)))
head(group_list)
# 设置对比
design <- model.matrix(~ 0 + factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva_es)
design
# 构建差异比较矩阵
contrast.matrix <- makeContrasts(IS-Control, levels = design)

# 差异分析，b vs. a
fit <- lmFit(gsva_es, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(x)

#把通路的limma分析结果保存到文件
write.csv(x, "gsva_limma_IS.csv", quote = F)

#输出t值，用做FigureYa39bar的输入数据
pathway <- str_replace(row.names(x), "HALLMARK_", "")
df <- data.frame(ID = pathway, score = x$t, adj.P.Val= x$adj.P.Val)
write.csv(df, "easy_input2_for39bar_IS.csv", quote = F, row.names = F)
head(df)

#按照score的值分组
cutoff <- 1
df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf),labels = c(1,2,3))
df$group2 <- cut(df$adj.P.Val, breaks = c(0.05, 1), labels = 2)
df$group[!is.na(df$group2)] <- df$group2[!is.na(df$group2)]

#按照score排序
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)

ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c('palegreen3', 'snow3', 'dodgerblue4'), guide = FALSE) + 
  
  #画2条虚线
  geom_hline(yintercept = c(-cutoff,cutoff), 
             color="white",
             linetype = 2, #画虚线
             size = 0.3) + #线的粗细
  
  #写label
  geom_text(data = subset(df, score < 0),
            aes(x=ID, y= 0, label= paste0(" ", ID) , color = group),#bar跟坐标轴间留出间隙
            size = 3, #字的大小
            hjust = "inward" ) +  #字的对齐方式
  geom_text(data = subset(df, score > 0),
            aes(x=ID, y= -0.1,label=ID, color = group),
            size = 3, hjust = "outward") +  
  scale_colour_manual(values = c("black","snow3","black"), guide = FALSE) +
  
  xlab("") +ylab("t value of GSVA score, IS \n versus Control")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_rect(size = 0.6)) + #边框粗细
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) #去除y轴

ggsave("gsva_IS.pdf", width = 9, height = 8)



#OSA
# 读入基因表达矩阵
a<-read.csv("OSA array.csv",row.names = 1)
# 这一句就完成了GSVA分析
gsva_es <- gsva(as.matrix(a), gs)
# 把通路的表达量保存到文件
write.csv(gsva_es, "gsva_output_OSA.csv", quote = F)
# 分组
group_list <- data.frame(sample = colnames(gsva_es), group = c(rep("control", 8), rep("stroke", 34)))
head(group_list)
# 设置对比
design <- model.matrix(~ 0 + factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva_es)
design
# 构建差异比较矩阵
contrast.matrix <- makeContrasts(stroke-control, levels = design)

# 差异分析，b vs. a
fit <- lmFit(gsva_es, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(x)

#把通路的limma分析结果保存到文件
write.csv(x, "gsva_limma_OSA.csv", quote = F)

#输出t值，用做FigureYa39bar的输入数据
pathway <- str_replace(row.names(x), "HALLMARK_", "")
df <- data.frame(ID = pathway, score = x$t, adj.P.Val= x$adj.P.Val)
write.csv(df, "easy_input2_for39bar_OSA.csv", quote = F, row.names = F)
head(df)

#按照score的值分组
cutoff <- 1
df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf),labels = c(1,2,3))
df$group2 <- cut(df$adj.P.Val, breaks = c(0.05, 1), labels = 2)
df$group[!is.na(df$group2)] <- df$group2[!is.na(df$group2)]

#按照score排序
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)

ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c('palegreen3', 'snow3', 'dodgerblue4'), guide = FALSE) + 
  
  #画2条虚线
  geom_hline(yintercept = c(-cutoff,cutoff), 
             color="white",
             linetype = 2, #画虚线
             size = 0.3) + #线的粗细
  
  #写label
  geom_text(data = subset(df, score < 0),
            aes(x=ID, y= 0, label= paste0(" ", ID) , color = group),#bar跟坐标轴间留出间隙
            size = 3, #字的大小
            hjust = "inward" ) +  #字的对齐方式
  geom_text(data = subset(df, score > 0),
            aes(x=ID, y= -0.1,label=ID, color = group),
            size = 3, hjust = "outward") +  
  scale_colour_manual(values = c("black","snow3","black"), guide = FALSE) +
  
  xlab("") +ylab("t value of GSVA score, OSA \n versus Control")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_rect(size = 0.6)) + #边框粗细
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) #去除y轴

ggsave("gsva_OSA.pdf", width = 9, height = 8)



#取重复上调下调
a1<-read.csv("easy_input2_for39bar_IS.csv",row.names = 1)
a2<-read.csv("easy_input2_for39bar_OSA.csv",row.names = 1)
up1<-rownames(a1[a1[,1]>1,,drop=FALSE])
up2<-rownames(a2[a2[,1]>1,,drop=FALSE])         
up<-intersect(up1,up2)
print(up)
down1<-rownames(a1[a1[,1] < -1,,drop=FALSE])
down2<-rownames(a2[a2[,1] < -1,,drop=FALSE]) 
down<-intersect(down1,down2)
print(down)
not1<-rownames(a1[a1[,1]<1 & a1[,1]>-1,,drop=FALSE])
not2<-rownames(a2[a2[,1]<1 & a2[,1]>-1,,drop=FALSE]) 
not<-intersect(not1,not2)
write.csv(up, "up.csv", quote = F, row.names = F)
write.csv(down, "down.csv", quote = F, row.names = F)
