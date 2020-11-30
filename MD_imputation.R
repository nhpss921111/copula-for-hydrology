# Missing data impution
# 開始撰寫日期：2020/11/12
# 完成撰寫日期：2020/11/30
# 讀取MD
# 選擇Q與Qs的邊際分布及其參數
# 選擇聯結函數及其參數
rm(list=ls())
library(copula)
library(CoImp)
library(stats) # 機率分布 
library(actuar) # 機率分布 
library(dplyr) # 資料整理 
library(FAdist) # 水文上常用機率分布
library(fitdistrplus) # 估計參數
library(EnvStats) # Environmental Statistics, Including US EPA Guidance
library(reliaR) # AIC of gumbel
library(lestat) # inverse CDF
library(goftest) # 適合度檢定
library(goft) # 適合度檢定
library(gumbel) # for gumbel copula
library(ggplot2) # 繪圖用
library(VineCopula) # for high dimension copula function
library(copula) # 聯結函數
library(scatterplot3d) #畫3D圖
library(locfit) # local polynomial estimation

year <- c(2012:2019)
MD.input <- c(paste0(year,"QandQs.csv"))
output <- c(paste0(year,"imp.csv"))
MD <- c()
# 主要迴圈(以年份為底)
for( y in 1:length(year)){
  # 讀取有缺失的資料
  setwd("F:/R_reading/CHIA-YUANG/missingdata")
  data <- read.csv(file.path(getwd(),MD.input[y]))
  data <- data[,-1]
  MD <- rbind(MD,data)
}
#
# ------------------------------- Imputation ----------------------------------------
n.marg <- 2
samp.miss <- as.matrix(MD[,4:5])
imp <- CoImp(samp.miss, n.marg=n.marg, q.up=rep(0.5,n.marg), q.lo=rep(0.1,n.marg), 
             smoothing = rep(0.7,n.marg),type.data="continuous", 
             model=list(gumbelCopula(),frankCopula(),claytonCopula()))
show(imp)
plot(imp)
imp@Imputed.data.matrix[,2]
imp.data <- cbind(MD,imp@Imputed.data.matrix[,2])
imp.data <- imp.data[,-5]
colnames(imp.data) <- c("Year","Month","Day","Discharge","Imp.SSL")
ggplot(data=imp.data)+
  geom_point(aes(x=Imp.SSL,y=Discharge))
# 存檔
file <- paste("F:/R_output/CHIA-YUANG/imputation/", output[y], sep="") #存檔路徑  
write.csv(imp.data,file)

