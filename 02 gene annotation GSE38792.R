setwd("D:\\生信分析数据\\15 OSA验证")
rm(list = ls())
gene<-read.csv("GSE38792_series_matrix.csv",row.names = 1)
zhu<-read.csv("GPL6244-17930.csv")

#log2
ex <- gene
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NA
gene <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}

#如果有多个基因，只取前面第一个
k<-zhu$gene_assignment
k<-strsplit(k," // ")
o<-data.frame()
for(i in 1:length(k)){
  k1<-k[[i]][2]
  o<-rbind(o,k1)
}
zhu$gene<-o$NA_character_.


gene$ID<-rownames(gene)
m<-merge(zhu,gene,"ID")
m<-m[,-c(1:2)]
k=aggregate(.~gene,mean,data=m)
rownames(k)<-k$gene
k<-k[,-1]
write.csv(k,"k.csv")
