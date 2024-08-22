setwd("D:\\生信分析数据\\02 OSA ID转换")
rm(list = ls())
gene<-read.csv("OSA matrix.csv",row.names = 1)
zhu<-read.csv("GPL6244-17930.csv")

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
