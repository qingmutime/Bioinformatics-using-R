setwd("D:\\生信分析数据\\08 交集venn")
rm(list = ls())

a1<-read.csv("diffup_IS.csv",row.names = 1)
up_IS<-row.names(a1)
a2<-read.csv("diffdown_IS.csv",row.names = 1)
down_IS<-row.names(a2)
a3<-read.csv("turquoise_IS.csv",row.names = 1)
IS_WGCNA<-row.names(a3)

b1<-read.csv("diffup_OSA.csv",row.names = 1)
up_OSA<-row.names(b1)
b2<-read.csv("diffdown_OSA.csv",row.names = 1)
down_OSA<-row.names(b2)
b3<-read.csv("turquoise_OSA.csv",row.names = 1)
OSA_WGCNA<-row.names(b3)


#venn图绘制
library(VennDiagram)
venn.diagram(list( up_IS=up_IS ,
                   IS_WGCNA=IS_WGCNA ,
                   up_OSA=up_OSA ,
                   OSA_WGCNA=OSA_WGCNA ),
             fill=c("red","green","yellow","blue"),
             filename="up_VennDiagram.tiff")

venn.diagram(list( down_IS=down_IS ,
                   IS_WGCNA=IS_WGCNA ,
                   down_OSA=down_OSA ,
                   OSA_WGCNA=OSA_WGCNA ),
             fill=c("red","green","yellow","blue"),
             filename="down_VennDiagram.tiff")
