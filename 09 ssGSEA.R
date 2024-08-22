setwd("D:\\生信分析数据\\12 免疫浸润")
rm(list=ls())

#IS
library(GSVA)
library(limma)
exp<-read.csv("IS array.csv",row.names = 1)
list<- read.csv("CellReports.csv",row.names = 1,header=F)
list<-as.data.frame(t(list))
list<-as.list(list)
ssgsea.res <-gsva(
  as.matrix(exp),
  list,
  method = "ssgsea",
  kcdf = "Gaussian",
  abs.ranking = T
)


library(pheatmap)
library(limma)
library(dplyr)
library(vegan)
library(ggplot2)
annotation_col<-read.csv("fen_GSE58294.csv",row.names = 1)
data.1 <- decostand(ssgsea.res,"standardize",MARGIN = 1)
color.key <- c("#3300CC", "#3399FF", "white", "#FF3333", "#CC0000")
gg<-pheatmap(
  data.1,
  color = colorRampPalette(color.key)(50),
  border_color = NA,
  annotation_col = annotation_col,
  labels_row = NULL,
  clustering_method = "ward.D2",
  show_rownames = T,
  show_colnames = F,
  fontsize_col = 5,
  cluster_cols = F
)
ggsave(gg,filename = '热图_IS.pdf', height = 6, width = 8)
dev.off()

write.csv(ssgsea.res,"ssgsea_IS.csv")
write.csv(data.1,"ssgsea_IS_std.csv")


#OSA
exp<-read.csv("OSA array.csv",row.names = 1)
list<- read.csv("CellReports.csv",row.names = 1,header=F)
list<-as.data.frame(t(list))
list<-as.list(list)
ssgsea.res <-gsva(
  as.matrix(exp),
  list,
  method = "ssgsea",
  kcdf = "Gaussian",
  abs.ranking = T
)


annotation_col<-read.csv("fen_GSE135917.csv",row.names = 1)
data.1 <- decostand(ssgsea.res,"standardize",MARGIN = 1)
color.key <- c("#3300CC", "#3399FF", "white", "#FF3333", "#CC0000")
gg<-pheatmap(
  data.1,
  color = colorRampPalette(color.key)(50),
  border_color = NA,
  annotation_col = annotation_col,
  labels_row = NULL,
  clustering_method = "ward.D2",
  show_rownames = T,
  show_colnames = F,
  fontsize_col = 5,
  cluster_cols = F
)
ggsave(gg,filename = '热图_OSA.pdf', height = 6, width = 8)
dev.off()

write.csv(ssgsea.res,"ssgsea_OSA.csv")
write.csv(data.1,"ssgsea_OSA_std.csv")
