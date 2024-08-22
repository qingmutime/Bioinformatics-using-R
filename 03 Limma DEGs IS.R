setwd("D:\\生信分析数据\\03 IS DEGs")
rm(list = ls())
library(limma)
library(dplyr)
#以整理好的GSE46960为例
a<-read.csv("IS array.csv",row.names = 1)#64个患者和14个疾病对照


list <- c(rep("Control", 23), rep("IS",23)) %>% factor(., levels = c("Control", "IS"), ordered = F)
list1 <- model.matrix(~factor(list)+0)  #把group设置成一个model matrix
colnames(list1) <- c("Control", "IS")
df.fit <- lmFit(a, list1)  ## 数据与list进行匹配
df.matrix <- makeContrasts(IS - Control, levels = list1)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
tempOutput <- topTable(fit,coef=1,n = Inf, adjust = "BH")
nrDEG = na.omit(tempOutput) ## 去掉数据中有NA的行或列
diffsig <- nrDEG  
write.csv(diffsig, "all.limmaOut.csv")
## 我们使用|logFC| > 1，padj < 0.05（矫正后P值）
foldChange = 0.585
padj = 0.05
## 筛选出所有差异基因的结果
All_diffSig <- diffsig[(diffsig$adj.P.Val < padj & (diffsig$logFC>foldChange | diffsig$logFC < (-foldChange))),]
write.csv(All_diffSig, "all.diffsig.csv")  ##输出差异基因数据集
diffup <-  All_diffSig[(All_diffSig$adj.P.Val < padj & (All_diffSig$logFC > foldChange)),]
write.csv(diffup, "diffup.csv")
#
diffdown <- All_diffSig[(All_diffSig$adj.P.Val < padj & (All_diffSig$logFC < -foldChange)),]
write.csv(diffdown, "diffdown.csv")


#热图
Groups<-list
Groups<-as.data.frame(Groups)
rownames(Groups)<-colnames(a)
heat<-a[rownames(a) %in% c(head(rownames(subset(All_diffSig,All_diffSig$logFC>0)),25),head(rownames(subset(All_diffSig,All_diffSig$logFC<0)),25)),]
library(pheatmap)
x <- t(scale(t(heat)))
gg<-pheatmap(x,annotation_col=Groups,fontsize = 4)

#火山图
png(file="火山图_IS.png",width = 600,height = 600)
yMax=max(-log10(diffsig$adj.P.Val))
xMax=max(abs(diffsig$logFC))
plot(diffsig$logFC,-log10(diffsig$adj.P.Val), xlab="logFC",ylab="-log10(adj.P.Val)",main="Volcano", xlim=c(-xMax,xMax),ylim=c(0,yMax),yaxs="i",pch=19, cex=1.2)
diffSub=subset(diffsig, adj.P.Val<padj & logFC>foldChange)
points(diffSub$logFC,-log10(diffSub$adj.P.Val), pch=19, col="red",cex=1.2)
diffSub=subset(diffsig, adj.P.Val<padj & logFC<(-foldChange))
points(diffSub$logFC,-log10(diffSub$adj.P.Val),  pch=19, col="green",cex=1.2)
abline(h=-log10(padj),lty=2,lwd=2)
abline(v=c(-0.585,0.585),lty=2,lwd=2)
dev.off()
