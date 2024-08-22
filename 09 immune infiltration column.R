setwd("D:\\生信分析数据\\12 免疫浸润")
rm(list=ls())

#IS
library(tidyr)
library(reshape2)

a <- read.csv("ssgsea_IS_std.csv", row.names = 1, header = T, as.is = F)
a<-as.data.frame(t(a))
a$X<-rownames(a)
fen<-read.csv("fen_GSE58294.csv")
a<-merge(fen,a,"X")

mydata1<-melt(
  a,
  id.vars=c("X","Groups"),
  variable.name="immunecell",
  value.name="tpm"
)

library(ggpubr)
library(ggplot2)

ylabname <- paste("immunecell", "expression")
colnames(mydata1) <- c("Sample", "Groups", "immunecell","tpm")
# 计算p value
pvalues <- sapply(mydata1$immunecell, function(x) {
  res <- wilcox.test(as.numeric(tpm) ~ Groups, data = subset(mydata1, immunecell == x)) #两组，wilcox.test或t.test；多组，kruskal.test或aov(one-way ANOVA test)
  res$p.value
})
pv <- data.frame(gene = mydata1$immunecell, pvalue = pvalues)
pv$sigcode <- cut(pv$pvalue, c(0,0.0001, 0.001, 0.01, 0.05, 1), 
                  labels=c('****','***', '**', '*', 'ns'))
mydata1<-mydata1[,-1]
# 画box plot
p.box <- ggplot(mydata1, aes(x=immunecell, y=tpm, color=Groups, fill=Groups)) +
  geom_boxplot(alpha = .5) + #半透明
  theme_classic() + #或theme_bw()
  scale_fill_brewer(palette = "Set1") + #按类填充颜色
  scale_color_brewer(palette = "Set1") + #按类给边框着色
  
  theme(axis.text.x = element_text(colour="black", size = 11,
                                   #名太挤，旋转45度
                                   angle = 45, hjust = .5, vjust = .5)) +
  geom_text(aes(x=gene, y=max(mydata1$tpm) * 1.1,
                label = pv$sigcode),
            data=pv, 
            inherit.aes=F) +
  ylab(ylabname)
p.box
ggsave("immunecellbox_IS.pdf", width = 14, height = 5)
# 画带散点的box plot
p.box.dot <- p.box + geom_point(shape = 21, size=.5, # 点的形状和大小
                                position = position_jitterdodge(), # 让点散开
                                alpha = .5) #半透明
p.box.dot
ggsave("immunecellsanbox_OSA.pdf", width = 14, height = 5)


#OSA
a <- read.csv("ssgsea_OSA_std.csv", row.names = 1, header = T, as.is = F)
a<-as.data.frame(t(a))
a$X<-rownames(a)
fen<-read.csv("fen_GSE135917.csv")
a<-merge(fen,a,"X")

mydata1<-melt(
  a,
  id.vars=c("X","Groups"),
  variable.name="immunecell",
  value.name="tpm"
)


ylabname <- paste("immunecell", "expression")
colnames(mydata1) <- c("Sample", "Groups", "immunecell","tpm")
# 计算p value
pvalues <- sapply(mydata1$immunecell, function(x) {
  res <- wilcox.test(as.numeric(tpm) ~ Groups, data = subset(mydata1, immunecell == x)) #两组，wilcox.test或t.test；多组，kruskal.test或aov(one-way ANOVA test)
  res$p.value
})
pv <- data.frame(gene = mydata1$immunecell, pvalue = pvalues)
pv$sigcode <- cut(pv$pvalue, c(0,0.0001, 0.001, 0.01, 0.05, 1), 
                  labels=c('****','***', '**', '*', 'ns'))
mydata1<-mydata1[,-1]
# 画box plot
p.box <- ggplot(mydata1, aes(x=immunecell, y=tpm, color=Groups, fill=Groups)) +
  geom_boxplot(alpha = .5) + #半透明
  theme_classic() + #或theme_bw()
  scale_fill_brewer(palette = "Set1") + #按类填充颜色
  scale_color_brewer(palette = "Set1") + #按类给边框着色
  
  theme(axis.text.x = element_text(colour="black", size = 11,
                                   #名太挤，旋转45度
                                   angle = 45, hjust = .5, vjust = .5)) +
  geom_text(aes(x=gene, y=max(mydata1$tpm) * 1.1,
                label = pv$sigcode),
            data=pv, 
            inherit.aes=F) +
  ylab(ylabname)
p.box
ggsave("immunecellbox_OSA.pdf", width = 14, height = 5)
# 画带散点的box plot
p.box.dot <- p.box + geom_point(shape = 21, size=.5, # 点的形状和大小
                                position = position_jitterdodge(), # 让点散开
                                alpha = .5) #半透明
p.box.dot
ggsave("immunecellsanbox_OSA.pdf", width = 14, height = 5)
