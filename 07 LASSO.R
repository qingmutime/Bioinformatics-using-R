#f1 = glmnet(x, y, family="binomial", nlambda=100, alpha=1) #这里alpha=1为LASSO回归，如果等于0就是岭回归
#参数 family 规定了回归模型的类型：
#family="gaussian" 适用于一维连续因变量（univariate）
#family="mgaussian" 适用于多维连续因变量（multivariate）
#family="poisson" 适用于非负次数因变量（count）
#family="binomial" 适用于二元离散因变量（binary）
#family="multinomial" 适用于多元离散因变量（category）
#我们这里结局指标是2分类变量，所以使用binomial
#print(f1)#把f1结果输出

setwd("D:\\生信分析数据\\07 LASSO")
rm(list = ls())
options(stringsAsFactors=F)

library(tidyverse)
library(broom)
library(glmnet)
k<-read.csv("IS array.csv",row.names = 1)
gene<-read.csv("交集基因.csv")
k<-k[rownames(k)%in%gene$gene,]
k<-as.data.frame(t(k))
fen<-read.csv("fen.csv",row.names = 1)
k<-k[rownames(fen),]
k<-cbind(fen$lasso,k)

states<-as.matrix(k)
x<-states[,-1]
y<-states[,1]
set.seed(1234567)
cvfit=cv.glmnet(x,y,type.measure = "mse",nfolds = 5,alpha=1)
plot(cvfit)
c(cvfit$lambda.min, cvfit$lambda.1se)
coef(cvfit$glmnet.fit, s=c(cvfit$lambda.min ,cvfit$lambda.1se))

print(cvfit$glmnet.fit) #结果中依次特征数，偏差解释比，参数λ，通常看最后一
plot(cvfit$glmnet.fit,label=T, cex = 0) #L1范数
plot(cvfit$glmnet.fit,xvar="lambda",label=T) #lambda参数
lasso.coef<-predict(cvfit$glmnet.fit,s=cvfit$lambda.min ,type="coefficients") #回归系数
plot(cvfit$glmnet.fit,xvar="dev",label=T) #解释偏差和回归系数的关系
lasso.y<-predict(cvfit$glmnet.fit,newx=x,type="response",s=cvfit$lambda.min ) #拟合
plot(lasso.y,y,xlab="Predicted",ylab="Actual",main="Lasso Regression")

#导出获取最优 lambda 的系数文件
coef_min=coef(cvfit$glmnet.fit,s=cvfit$lambda.min)
coef_min <- coef_min[-1,]
coef_min<-as.data.frame(coef_min)
coef_min<-coef_min[coef_min[, 1] != 0, , drop = FALSE]
write.csv(coef_min, "lasso.csv")
