setwd("D:\\生信分析数据\\08 ROC")
rm(list = ls())

#建模
hubgene <- read.csv("lasso.csv",row.names = 1)
genename <- rownames(hubgene)
genecoef <- hubgene[,1]
myfun = function(x){crossprod(as.numeric(x), genecoef)}


#ROC
library(pROC); library(ggplot2)

#GSE58294 IS
IS <- read.csv("IS array.csv", row.names = 1)
data <- t(IS)
data = data[,genename]
data <- as.data.frame(data)
riskscore = apply(data,1,myfun)
group <- read.csv("fen_GSE58294.csv", row.names = 1)
model = cbind(group$Groups, data, riskscore=as.vector(riskscore))

roc <- roc(group$Groups ~ riskscore, data = model)
pdf("GSE58294.pdf", width = 8, height = 8)
plot(roc,
     print.auc = TRUE,
     smooth = F,
     col = "red",
     legacy.axes = TRUE,
     main = "GSE58294")
dev.off()


#GSE135917 OSA
OSA <- read.csv("OSA array.csv", row.names = 1)
data <- t(OSA)
data = data[,genename]
data <- as.data.frame(data)
riskscore = apply(data, 1, myfun)
group <- read.csv("fen_GSE135917.csv", row.names = 1)
model = cbind(group$Groups, data, riskscore=as.vector(riskscore))

roc <- roc(group$Groups ~ riskscore, data = model)
pdf("GSE135917.pdf", width = 8, height = 8)
plot(roc,
     print.auc = TRUE,
     smooth = F,
     col = "red",
     legacy.axes = TRUE,
     main = "GSE135917")
dev.off()


#IS验证集GSE16561
stroke <- read.csv("IS_GSE16561.csv", row.names = 1)
data <- t(stroke)
data = data[,genename]
data <- as.data.frame(data)
riskscore = apply(data, 1, myfun)
group <- read.csv("fen_GSE16561.csv", row.names = 1)
model = cbind(group$Groups, data, riskscore=as.vector(riskscore))

roc <- roc(group$Groups ~ riskscore, data = model)
pdf("GSE16561.pdf", width = 8, height = 8)
plot(roc,
     print.auc = TRUE,
     smooth = F,
     col = "red",
     legacy.axes = TRUE,
     main = "GSE16561")
dev.off()


#OSA验证集GSE38792
OSA <- read.csv("OSA_GSE38792.csv", row.names = 1)
data <- t(OSA)
data = data[,genename]
data <- as.data.frame(data)
riskscore = apply(data,1,myfun)
group <- read.csv("fen_GSE38792.csv", row.names = 1)
model = cbind(group$Groups, data, riskscore=as.vector(riskscore))

roc <- roc(group$Groups ~ riskscore, data = model)
pdf("GSE38792.pdf", width = 8, height = 8)
plot(roc,
     print.auc = TRUE,
     smooth = F,
     col = "red",
     legacy.axes = TRUE,
     main = "GSE38792")
dev.off()
