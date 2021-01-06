# Missing data impution 驗證
# 開始撰寫日期：2020/12/24
# 完成撰寫日期：2020/12/29
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
library(Metrics) # for RMSE
library(nls2)
library(Metrics) # 評估指標
# ===========
# 家源橋CHIA-YUANG year <- c(1974:2009,2012:2019)
# 彰雲橋CHUNYUN BRIDGE year <- c(1987:2019)
# ===========

year <- c(1974:2009,2012:2019)
MD.input <- c(paste0(year,"QandQs.csv"))
output <- c(paste0(year,"imp.csv"))
ob.data <- c()
set.seed(101)
# set.seed(101)
# set.seed(102)
# set.seed(103)
# set.seed(104)
# set.seed(105)
# set.seed(106)
# set.seed(107)
# set.seed(108)
# set.seed(109)
#
# 主要迴圈(以年份為底)
for( y in 1:length(year)){
  # 讀取有缺失的資料
  setwd("F:/R_reading/CHIA-YUANG/missingdata")
  data <- read.csv(file.path(getwd(),MD.input[y]))
  data <- data[,-1]
  ob.data <- rbind(ob.data,data)
}
ob.data <- subset(ob.data,ob.data[,4]>0) # 把流量為0(無觀測資料)刪除
ob.data <- ob.data[complete.cases(ob.data), ] # 移除全部NA 
#ob.data <- subset(ob.data,ob.data[,4]>20 & ob.data[,5]>1000)
#
# ---- 不分組資料畫成 scatter plot ----
#
ggplot(data=ob.data,aes(x=Discharge,y=Suspended.Load))+
  geom_point()+
  scale_color_discrete(name="年",labels=c(""))+
  labs(x="流量(cms)",y="輸砂量Qs (公噸)") + # 座標軸名稱
  #geom_mark_ellipse(expand = 0,aes(fill=group))+ # 橢圓形圈
  #geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=group))+ # 多邊形
  #geom_mark_hull(expand=0.01,aes(fill=group))+ # 凹凸多邊形
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=30))  # 字體大小
#
# ----------- 將60%資料總數建模，剩下40%當成驗證 -------------
#
perc.mis    <- 0.3 # 40%的資料當成NA 
x.samp <- as.matrix(ob.data[,4:5])
miss.row    <- sample(1:length(ob.data$Discharge), perc.mis*length(ob.data$Discharge), replace=FALSE)
miss.col    <- rep(2,perc.mis*length(ob.data$Discharge))
miss        <- cbind(miss.row,miss.col)
samp.miss <- replace(x.samp,miss,NA) # NA的欄位座標
MD <- cbind(ob.data[,1:3],samp.miss) # 將40%觀測資料轉換成NA

MD60 <- MD[complete.cases(MD), ] # 移除全部NA (60%剩餘資料)



#  -------------- 決定流量分組範圍 (來源：samp.miss) ---------------
# group 1： 0% ~  20%
# group 2：20% ~  40%
# group 3：40% ~  60%
# group 4：60% ~  80%
# group 5：80% ~  90%
# group 6：90% ~  95%
# group 8：95% ~  98%
# group 8：98% ~  99%
# group 9：99% ~ 100%

rank.data <- cbind(MD60,rank(MD60$Discharge))
persent <- (rank.data$`rank(MD60$Discharge)`) / length(rank.data$Discharge)
per.data <- cbind(rank.data,persent)

data.1 <- data.frame(subset(per.data, persent<=0.2),group="group1")
data.2 <- data.frame(subset(per.data, persent>0.2 & persent<=0.4),group="group2")
data.3 <- data.frame(subset(per.data, persent>0.4 & persent<=0.6),group="group3")
data.4 <- data.frame(subset(per.data, persent>0.6 & persent<=0.8),group="group4")
data.5 <- data.frame(subset(per.data, persent>0.8 & persent<=0.9),group="group5")
data.6 <- data.frame(subset(per.data, persent>0.9 & persent<=0.96),group="group6")
data.7 <- data.frame(subset(per.data, persent>0.96 & persent<=0.98),group="group7")
data.8 <- data.frame(subset(per.data, persent>0.98 & persent<=0.99),group="group8")
data.9 <- data.frame(subset(per.data, persent>0.99),group="group9")

data.group <- rbind(data.1, data.2, data.3,
                    data.4, data.5, data.6,
                    data.7, data.8, data.9)

# ------------------------------- imputation_validation ----------------------------------------
# =========== 分一組 (建立 rating curve)===============
rating <- nls(Suspended.Load ~ a*Discharge^b,  algorithm="port",
              control = list(maxiter = 200,minFactor = 1/2000000,warnOnly = TRUE),
              start=list(a=0.01,b=1), data=MD60,trace=T) # 乘冪迴歸式
summary(rating) # 觀察
a <- environment(rating[["m"]][["resid"]])[["env"]][["a"]] # 迴歸係數a
b <- environment(rating[["m"]][["resid"]])[["env"]][["b"]] # 迴歸係數b
ratingSSL <- a*(MD$Discharge)^b # 計算率定曲線推估的輸砂量
rating.1 <- cbind(subset(ob.data),subset(MD[,4:5])[,2],ratingSSL)

rmse <- function(actual, predicted) { # 均方根誤差公式
  sqrt(mean((actual - predicted) ^ 2))
}
# rmse (root mean square error)

rmse.rating <- rmse(rating.1$Suspended.Load, rating.1$ratingSSL) #原本移除的觀測值與率定值

# mae (mean absolute error)

mae.rating <- mae(rating.1$Suspended.Load, rating.1$ratingSSL) #原本移除的觀測值與率定值

# mape (mean absolute persentage error)

mape.rating <- mape(rating.1$Suspended.Load, rating.1$ratingSSL) #原本移除的觀測值與率定值

rating.table <- cbind(rmse.rating,mae.rating,mape.rating)

# =========== 分九組 ===============
#分組
# ---- group1：0%~20% ----
# 將一開始移除的40% 加回來
maxBC1 <- max(data.1$Discharge)
MD.1 <- as.matrix(subset(MD[,4:5],Discharge<maxBC1)) #提取Q與Qs出來，並限制流量範圍
n.marg <- 2 # 兩個變數(Q、Qs)
imp1 <- CoImp(MD.1, n.marg=n.marg, q.up=c(0.5,0.5), q.lo=c(0.01,0.01), 
             smoothing = c(1,1),type.data="continuous", 
             model=list(gumbelCopula(),frankCopula(),claytonCopula())) # 補遺計算
plot(imp1)
imp.1 <- cbind(subset(ob.data,Discharge<maxBC1),
               subset(MD[,4:5],Discharge<maxBC1)[,2],
               imp1@Imputed.data.matrix[,2]) # 提取補遺值，並合併成結果表格
colnames(imp.1) <- c("Year","Month","Day","Discharge",
                     "Suspended.Load","asNA","Imp.SSL") # 為框架的行命名
rating <- nls(Suspended.Load ~ a*Discharge^b,  algorithm="port",
              control = list(maxiter = 200,minFactor = 1/200,warnOnly = TRUE),
              start=list(a=1,b=1), data=data.1[,4:5],trace=T) # 乘冪迴歸式
summary(rating) # 觀察
a <- environment(rating[["m"]][["resid"]])[["env"]][["a"]] # 迴歸係數a
b <- environment(rating[["m"]][["resid"]])[["env"]][["b"]] # 迴歸係數b
ratingSSL <- a*(imp.1$Discharge)^b # 計算率定曲線推估的輸砂量
imp.1 <- cbind(imp.1,ratingSSL)

rmse <- function(actual, predicted) { # 均方根誤差公式
  sqrt(mean((actual - predicted) ^ 2))
}
# rmse (root mean square error)
rmse.imp <- rmse(imp.1$Suspended.Load, imp.1$Imp.SSL) #原本移除的觀測值與補遺值
rmse.rating <- rmse(imp.1$Suspended.Load, imp.1$ratingSSL) #原本移除的觀測值與率定值
rmse1 <- cbind(rmse.imp, rmse.rating)
# mae (mean absolute error)
mae.imp <- mae(imp.1$Suspended.Load, imp.1$Imp.SSL) #原本移除的觀測值與補遺值
mae.rating <- mae(imp.1$Suspended.Load, imp.1$ratingSSL) #原本移除的觀測值與率定值
mae1 <- cbind(mae.imp, mae.rating)
# mape (mean absolute persentage error)
mape.imp <- mape(imp.1$Suspended.Load, imp.1$Imp.SSL) #原本移除的觀測值與補遺值
mape.rating <- mape(imp.1$Suspended.Load, imp.1$ratingSSL) #原本移除的觀測值與率定值
mape1 <- cbind(mape.imp, mape.rating)
# 補遺資料出圖
setwd("F:/R_output/CHIA-YUANG/imputation_validation") # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(0~20)3.png"),width = 1250, height = 700, units = "px", pointsize = 12)
Imp.1 <- ggplot(data=imp.1)+
  geom_point(aes(x=Discharge,y=Imp.SSL,color="補遺後")) +
  geom_point(aes(x=Discharge,y=Suspended.Load,color="補遺前")) +
  geom_line(aes(x=Discharge,y=ratingSSL,color="率定曲線")) +
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(0%~20%)")) +
  theme(text=element_text(size=20))  # 字體大小
plot(Imp.1)
dev.off()
  
# ---- group2：20%~40% ----

# 將一開始移除的40% 加回來
maxBC2 <- max(data.2$Discharge)
MD.2 <- as.matrix(subset(MD[,4:5],Discharge>=maxBC1 & Discharge<maxBC2)) #提取Q與Qs出來，並限制流量範圍
n.marg <- 2 # 兩個變數(Q、Qs)
imp2 <- CoImp(MD.2, n.marg=n.marg, q.up=c(0.5,0.5), q.lo=c(0.01,0.01), 
              smoothing = c(1,1),type.data="continuous", 
              model=list(gumbelCopula(),frankCopula(),claytonCopula())) # 補遺計算
imp.2 <- cbind(subset(ob.data,Discharge>=maxBC1 & Discharge<maxBC2),
               subset(MD[,4:5],Discharge>=maxBC1 & Discharge<maxBC2)[,2],
               imp2@Imputed.data.matrix[,2]) # 提取補遺值，並合併成結果表格
colnames(imp.2) <- c("Year","Month","Day","Discharge",
                     "Suspended.Load","asNA","Imp.SSL") # 為框架的行命名
rating <- nls(Suspended.Load ~ a*Discharge^b, algorithm="port",
              control = list(maxiter = 200,minFactor = 1/200,warnOnly = TRUE),
              start=list(a=0.01,b=1), data=data.2[,4:5],trace=T) # 乘冪迴歸式
#save <- data.2[,4:5]
#file <- paste("F:/R_output/CHIA-YUANG/imputation_validation_validation_validation/ttt.csv", sep="") # 存檔路徑
#write.csv(save, file)

summary(rating) # 觀察
a <- environment(rating[["m"]][["resid"]])[["env"]][["a"]] # 迴歸係數a
b <- environment(rating[["m"]][["resid"]])[["env"]][["b"]] # 迴歸係數b
ratingSSL <- a*(imp.2$Discharge)^b # 計算率定曲線推估的輸砂量
imp.2 <- cbind(imp.2,ratingSSL)

rmse <- function(actual, predicted) { # 均方根誤差公式
  sqrt(mean((actual - predicted) ^ 2))
}
# rmse (root mean square error)
rmse.imp <- rmse(imp.2$Suspended.Load, imp.2$Imp.SSL) #原本移除的觀測值與補遺值
rmse.rating <- rmse(imp.2$Suspended.Load, imp.2$ratingSSL) #原本移除的觀測值與率定值
rmse2 <- cbind(rmse.imp,rmse.rating)
# mae (mean absolute error)
mae.imp <- mae(imp.2$Suspended.Load, imp.2$Imp.SSL) #原本移除的觀測值與補遺值
mae.rating <- mae(imp.2$Suspended.Load, imp.2$ratingSSL) #原本移除的觀測值與率定值
mae2 <- cbind(mae.imp, mae.rating)
# mape (mean absolute persentage error)
mape.imp <- mape(imp.2$Suspended.Load, imp.2$Imp.SSL) #原本移除的觀測值與補遺值
mape.rating <- mape(imp.2$Suspended.Load, imp.2$ratingSSL) #原本移除的觀測值與率定值
mape2 <- cbind(mape.imp, mape.rating)
# 補遺資料出圖
setwd("F:/R_output/CHIA-YUANG/imputation_validation") # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(20~40)3.png"),width = 1250, height = 700, units = "px", pointsize = 12)
Imp.2 <- ggplot(data=imp.2)+
  geom_point(aes(x=Discharge,y=Imp.SSL,color="補遺後"))+
  geom_point(aes(x=Discharge,y=Suspended.Load,color="補遺前"))+
  geom_line(aes(x=Discharge,y=ratingSSL,color="率定曲線")) + #
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(20%~40%)"))+
  theme(text=element_text(size=20))  # 字體大小
plot(Imp.2)
dev.off()

# ---- group3：40% ~ 60% ----

# 將一開始移除的40% 加回來
maxBC3 <- max(data.3$Discharge)
MD.3 <- as.matrix(subset(MD[,4:5],Discharge>=maxBC2 & Discharge<maxBC3)) #提取Q與Qs出來，並限制流量範圍
n.marg <- 2 # 兩個變數(Q、Qs)
imp3 <- CoImp(MD.3, n.marg=n.marg, q.up=c(0.5,0.5), q.lo=c(0.01,0.01), 
              smoothing = c(1,1),type.data="continuous", 
              model=list(gumbelCopula(),frankCopula(),claytonCopula())) # 補遺計算
imp.3 <- cbind(subset(ob.data,Discharge>=maxBC2 & Discharge<maxBC3),
               subset(MD[,4:5],Discharge>=maxBC2 & Discharge<maxBC3)[,2],
               imp3@Imputed.data.matrix[,2]) # 提取補遺值，並合併成結果表格
colnames(imp.3) <- c("Year","Month","Day","Discharge",
                     "Suspended.Load","asNA","Imp.SSL") # 為框架的行命名
rating <- nls(Suspended.Load ~ a*Discharge^b,algorithm="port",
              control = list(maxiter = 200,minFactor = 1/200,warnOnly = TRUE),
              start=list(a=1,b=1), data=data.3[,4:5],trace=T) # 乘冪迴歸式
#save <- data.3[,4:5]
#file <- paste("F:/R_output/CHIA-YUANG/imputation_validation/ttt.csv", sep="") # 存檔路徑
#write.csv(save, file)

summary(rating) # 觀察
a <- environment(rating[["m"]][["resid"]])[["env"]][["a"]] # 迴歸係數a
b <- environment(rating[["m"]][["resid"]])[["env"]][["b"]] # 迴歸係數b
ratingSSL <- a*(imp.3$Discharge)^b # 計算率定曲線推估的輸砂量
imp.3 <- cbind(imp.3,ratingSSL)

rmse <- function(actual, predicted) { # 均方根誤差公式
  sqrt(mean((actual - predicted) ^ 2))
}
# rmse (root mean square error)
rmse.imp <- rmse(imp.3$Suspended.Load, imp.3$Imp.SSL) #原本移除的觀測值與補遺值
rmse.rating <- rmse(imp.3$Suspended.Load, imp.3$ratingSSL) #原本移除的觀測值與率定值
rmse3 <- cbind(rmse.imp,rmse.rating)
# mae (mean absolute error)
mae.imp <- mae(imp.3$Suspended.Load, imp.3$Imp.SSL) #原本移除的觀測值與補遺值
mae.rating <- mae(imp.3$Suspended.Load, imp.3$ratingSSL) #原本移除的觀測值與率定值
mae3 <- cbind(mae.imp, mae.rating)
# mape (mean absolute persentage error)
mape.imp <- mape(imp.3$Suspended.Load, imp.3$Imp.SSL) #原本移除的觀測值與補遺值
mape.rating <- mape(imp.3$Suspended.Load, imp.3$ratingSSL) #原本移除的觀測值與率定值
mape3 <- cbind(mape.imp, mape.rating)
# 補遺資料出圖
setwd("F:/R_output/CHIA-YUANG/imputation_validation") # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(40~60)3.png"),width = 1250, height = 700, units = "px", pointsize = 12)
Imp.3 <- ggplot(data=imp.3)+
  geom_point(aes(x=Discharge,y=Imp.SSL,color="補遺後"))+
  geom_point(aes(x=Discharge,y=Suspended.Load,color="補遺前"))+
  geom_line(aes(x=Discharge,y=ratingSSL,color="率定曲線")) + #
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(40%~60%)"))+
  theme(text=element_text(size=20))  # 字體大小
plot(Imp.3)
dev.off()

# ---- group4：60% ~ 80% ----

# 將一開始移除的40% 加回來
maxBC4 <- max(data.4$Discharge)
MD.4 <- as.matrix(subset(MD[,4:5],Discharge>=maxBC3 & Discharge<maxBC4)) #提取Q與Qs出來，並限制流量範圍
n.marg <- 2 # 兩個變數(Q、Qs)
imp4 <- CoImp(MD.4, n.marg=n.marg, q.up=c(0.5,0.5), q.lo=c(0.01,0.01), 
              smoothing = c(1,1),type.data="continuous", 
              model=list(gumbelCopula(),frankCopula(),claytonCopula())) # 補遺計算
imp.4 <- cbind(subset(ob.data,Discharge>=maxBC3 & Discharge<maxBC4),
               subset(MD[,4:5],Discharge>=maxBC3 & Discharge<maxBC4)[,2],
               imp4@Imputed.data.matrix[,2]) # 提取補遺值，並合併成結果表格
colnames(imp.4) <- c("Year","Month","Day","Discharge",
                     "Suspended.Load","asNA","Imp.SSL") # 為框架的行命名
rating <- nls(Suspended.Load ~ a*Discharge^b,algorithm="port",
              control = list(maxiter = 200,minFactor = 1/200,warnOnly = TRUE),
              start=list(a=1,b=1), data=data.4[,4:5],trace=T) # 乘冪迴歸式
coef(rating)
#save <- data.4[,4:5]
#file <- paste("F:/R_output/CHIA-YUANG/imputation_validation/ttt.csv", sep="") # 存檔路徑
#write.csv(save, file)

summary(rating) # 觀察
a <- environment(rating[["m"]][["resid"]])[["env"]][["a"]] # 迴歸係數a
b <- environment(rating[["m"]][["resid"]])[["env"]][["b"]] # 迴歸係數b
ratingSSL <- a*(imp.4$Discharge)^b # 計算率定曲線推估的輸砂量
imp.4 <- cbind(imp.4,ratingSSL)

rmse <- function(actual, predicted) { # 均方根誤差公式
  sqrt(mean((actual - predicted) ^ 2))
}
# rmse (root mean square error)
rmse.imp <- rmse(imp.4$Suspended.Load, imp.4$Imp.SSL) #原本移除的觀測值與補遺值
rmse.rating <- rmse(imp.4$Suspended.Load, imp.4$ratingSSL) #原本移除的觀測值與率定值
rmse4 <- cbind(rmse.imp,rmse.rating)
# mae (mean absolute error)
mae.imp <- mae(imp.4$Suspended.Load, imp.4$Imp.SSL) #原本移除的觀測值與補遺值
mae.rating <- mae(imp.4$Suspended.Load, imp.4$ratingSSL) #原本移除的觀測值與率定值
mae4 <- cbind(mae.imp, mae.rating)
# mape (mean absolute persentage error)
mape.imp <- mape(imp.4$Suspended.Load, imp.4$Imp.SSL) #原本移除的觀測值與補遺值
mape.rating <- mape(imp.4$Suspended.Load, imp.4$ratingSSL) #原本移除的觀測值與率定值
mape4 <- cbind(mape.imp, mape.rating)
# 補遺資料出圖
setwd("F:/R_output/CHIA-YUANG/imputation_validation") # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(60~80)3.png"),width = 1250, height = 700, units = "px", pointsize = 12)
Imp.4 <- ggplot(data=imp.4)+
  geom_point(aes(x=Discharge,y=Imp.SSL,color="補遺後"))+
  geom_point(aes(x=Discharge,y=Suspended.Load,color="補遺前"))+
  geom_line(aes(x=Discharge,y=ratingSSL,color="率定曲線")) + #
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(60%~80%)"))+
  theme(text=element_text(size=20))  # 字體大小
plot(Imp.4)
dev.off()

# ---- group5：80% ~ 90% ----

# 將一開始移除的40% 加回來
maxBC5 <- max(data.5$Discharge)
MD.5 <- as.matrix(subset(MD[,4:5],Discharge>=maxBC4 & Discharge<maxBC5)) #提取Q與Qs出來，並限制流量範圍
n.marg <- 2 # 兩個變數(Q、Qs)
imp5 <- CoImp(MD.5, n.marg=n.marg, q.up=c(0.5,0.5), q.lo=c(0.01,0.01), 
              smoothing = c(1,1),type.data="continuous", 
              model=list(gumbelCopula(),frankCopula(),claytonCopula())) # 補遺計算
imp.5 <- cbind(subset(ob.data,Discharge>=maxBC4 & Discharge<maxBC5),
               subset(MD[,4:5],Discharge>=maxBC4 & Discharge<maxBC5)[,2],
               imp5@Imputed.data.matrix[,2]) # 提取補遺值，並合併成結果表格
colnames(imp.5) <- c("Year","Month","Day","Discharge",
                     "Suspended.Load","asNA","Imp.SSL") # 為框架的行命名
rating <- nls(Suspended.Load ~ a*Discharge^b,algorithm="port",
              control = list(maxiter = 200,minFactor = 1/200,warnOnly = TRUE),
              start=list(a=1,b=1), data=data.5[,4:5],trace=T) # 乘冪迴歸式
coef(rating)
#save <- data.5[,4:5]
#file <- paste("F:/R_output/CHIA-YUANG/imputation_validation/ttt.csv", sep="") # 存檔路徑
#write.csv(save, file)

summary(rating) # 觀察
a <- environment(rating[["m"]][["resid"]])[["env"]][["a"]] # 迴歸係數a
b <- environment(rating[["m"]][["resid"]])[["env"]][["b"]] # 迴歸係數b
ratingSSL <- a*(imp.5$Discharge)^b # 計算率定曲線推估的輸砂量
imp.5 <- cbind(imp.5,ratingSSL)

rmse <- function(actual, predicted) { # 方均根誤差公式
  sqrt(mean((actual - predicted) ^ 2))
}
# rmse (root mean square error)
rmse.imp <- rmse(imp.5$Suspended.Load, imp.5$Imp.SSL) #原本移除的觀測值與補遺值
rmse.rating <- rmse(imp.5$Suspended.Load, imp.5$ratingSSL) #原本移除的觀測值與率定值
rmse5 <- cbind(rmse.imp,rmse.rating)
# mae (mean absolute error)
mae.imp <- mae(imp.5$Suspended.Load, imp.5$Imp.SSL) #原本移除的觀測值與補遺值
mae.rating <- mae(imp.5$Suspended.Load, imp.5$ratingSSL) #原本移除的觀測值與率定值
mae5 <- cbind(mae.imp, mae.rating)
# mape (mean absolute persentage error)
mape.imp <- mape(imp.5$Suspended.Load, imp.5$Imp.SSL) #原本移除的觀測值與補遺值
mape.rating <- mape(imp.5$Suspended.Load, imp.5$ratingSSL) #原本移除的觀測值與率定值
mape5 <- cbind(mape.imp, mape.rating)
# 補遺資料出圖
setwd("F:/R_output/CHIA-YUANG/imputation_validation") # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(80~90)3.png"),width = 1250, height = 700, units = "px", pointsize = 12)
Imp.5 <- ggplot(data=imp.5)+
  geom_point(aes(x=Discharge,y=Imp.SSL,color="補遺後"))+
  geom_point(aes(x=Discharge,y=Suspended.Load,color="補遺前"))+
  geom_line(aes(x=Discharge,y=ratingSSL,color="率定曲線")) + #
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(80%~90%)"))+
  theme(text=element_text(size=20))  # 字體大小
plot(Imp.5)
dev.off()

# ---- group6：90% ~ 95% ----

# 將一開始移除的40% 加回來
maxBC6 <- max(data.6$Discharge)
MD.6 <- as.matrix(subset(MD[,4:5],Discharge>=maxBC5 & Discharge<maxBC6)) #提取Q與Qs出來，並限制流量範圍
n.marg <- 2 # 兩個變數(Q、Qs)
imp6 <- CoImp(MD.6, n.marg=n.marg, q.up=c(0.5,0.5), q.lo=c(0.01,0.01), 
              smoothing = c(1,1),type.data="continuous", 
              model=list(gumbelCopula(),frankCopula(),claytonCopula())) # 補遺計算
imp.6 <- cbind(subset(ob.data,Discharge>=maxBC5 & Discharge<maxBC6),
               subset(MD[,4:5],Discharge>=maxBC5 & Discharge<maxBC6)[,2],
               imp6@Imputed.data.matrix[,2]) # 提取補遺值，並合併成結果表格
colnames(imp.6) <- c("Year","Month","Day","Discharge",
                     "Suspended.Load","asNA","Imp.SSL") # 為框架的行命名
rating <- nls(Suspended.Load ~ a*Discharge^b,algorithm="port",
              control = list(maxiter = 200,minFactor = 1/200,warnOnly = TRUE),
              start=list(a=1,b=1), data=data.6[,4:5],trace=T) # 乘冪迴歸式
coef(rating)
#save <- data.6[,4:5]
#file <- paste("F:/R_output/CHIA-YUANG/imputation_validation_validation/ttt.csv", sep="") # 存檔路徑
#write.csv(save, file)

summary(rating) # 觀察
a <- environment(rating[["m"]][["resid"]])[["env"]][["a"]] # 迴歸係數a
b <- environment(rating[["m"]][["resid"]])[["env"]][["b"]] # 迴歸係數b
ratingSSL <- a*(imp.6$Discharge)^b # 計算率定曲線推估的輸砂量
imp.6 <- cbind(imp.6,ratingSSL)

rmse <- function(actual, predicted) { # 方均根誤差公式
  sqrt(mean((actual - predicted) ^ 2))
}
# rmse (root mean square error)
rmse.imp <- rmse(imp.6$Suspended.Load, imp.6$Imp.SSL) #原本移除的觀測值與補遺值
rmse.rating <- rmse(imp.6$Suspended.Load, imp.6$ratingSSL) #原本移除的觀測值與率定值
rmse6 <- cbind(rmse.imp,rmse.rating)
# mae (mean absolute error)
mae.imp <- mae(imp.6$Suspended.Load, imp.6$Imp.SSL) #原本移除的觀測值與補遺值
mae.rating <- mae(imp.6$Suspended.Load, imp.6$ratingSSL) #原本移除的觀測值與率定值
mae6 <- cbind(mae.imp, mae.rating)
# mape (mean absolute persentage error)
mape.imp <- mape(imp.6$Suspended.Load, imp.6$Imp.SSL) #原本移除的觀測值與補遺值
mape.rating <- mape(imp.6$Suspended.Load, imp.6$ratingSSL) #原本移除的觀測值與率定值
mape6 <- cbind(mape.imp, mape.rating)
# 補遺資料出圖
setwd("F:/R_output/CHIA-YUANG/imputation_validation") # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(90~95)3.png"),width = 1250, height = 700, units = "px", pointsize = 12)
Imp.6 <- ggplot(data=imp.6)+
  geom_point(aes(x=Discharge,y=Imp.SSL,color="補遺後"))+
  geom_point(aes(x=Discharge,y=Suspended.Load,color="補遺前"))+
  geom_line(aes(x=Discharge,y=ratingSSL,color="率定曲線")) + #
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(90%~95%)"))+
  theme(text=element_text(size=20))  # 字體大小
plot(Imp.6)
dev.off()

# ---- group7：95% ~ 98% ----

# 將一開始移除的40% 加回來
maxBC7 <- max(data.7$Discharge)
MD.7 <- as.matrix(subset(MD[,4:5],Discharge>=maxBC6 & Discharge<maxBC7)) #提取Q與Qs出來，並限制流量範圍
n.marg <- 2 # 兩個變數(Q、Qs)
imp7 <- CoImp(MD.7, n.marg=n.marg, q.up=c(0.5,0.5), q.lo=c(0.01,0.01), 
              smoothing = c(1,1),type.data="continuous", 
              model=list(gumbelCopula(),frankCopula(),claytonCopula())) # 補遺計算
imp.7 <- cbind(subset(ob.data,Discharge>=maxBC6 & Discharge<maxBC7),
               subset(MD[,4:5],Discharge>=maxBC6 & Discharge<maxBC7)[,2],
               imp7@Imputed.data.matrix[,2]) # 提取補遺值，並合併成結果表格
colnames(imp.7) <- c("Year","Month","Day","Discharge",
                     "Suspended.Load","asNA","Imp.SSL") # 為框架的行命名
rating <- nls(Suspended.Load ~ a*Discharge^b,algorithm="port",
              control = list(maxiter = 200,minFactor = 1/200,warnOnly = TRUE),
              start=list(a=1,b=1), data=data.7[,4:5],trace=T) # 乘冪迴歸式
coef(rating)
#save <- data.7[,4:5]
#file <- paste("F:/R_output/CHIA-YUANG/imputation_validation_validation/ttt.csv", sep="") # 存檔路徑
#write.csv(save, file)

summary(rating) # 觀察
a <- environment(rating[["m"]][["resid"]])[["env"]][["a"]] # 迴歸係數a
b <- environment(rating[["m"]][["resid"]])[["env"]][["b"]] # 迴歸係數b
ratingSSL <- a*(imp.7$Discharge)^b # 計算率定曲線推估的輸砂量
imp.7 <- cbind(imp.7,ratingSSL)

rmse <- function(actual, predicted) { # 方均根誤差公式
  sqrt(mean((actual - predicted) ^ 2))
}
# rmse (root mean square error)
rmse.imp <- rmse(imp.7$Suspended.Load, imp.7$Imp.SSL) #原本移除的觀測值與補遺值
rmse.rating <- rmse(imp.7$Suspended.Load, imp.7$ratingSSL) #原本移除的觀測值與率定值
rmse7 <- cbind(rmse.imp,rmse.rating)
# mae (mean absolute error)
mae.imp <- mae(imp.7$Suspended.Load, imp.7$Imp.SSL) #原本移除的觀測值與補遺值
mae.rating <- mae(imp.7$Suspended.Load, imp.7$ratingSSL) #原本移除的觀測值與率定值
mae7 <- cbind(mae.imp, mae.rating)
# mape (mean absolute persentage error)
mape.imp <- mape(imp.7$Suspended.Load, imp.7$Imp.SSL) #原本移除的觀測值與補遺值
mape.rating <- mape(imp.7$Suspended.Load, imp.7$ratingSSL) #原本移除的觀測值與率定值
mape7 <- cbind(mape.imp, mape.rating)
# 補遺資料出圖
setwd("F:/R_output/CHIA-YUANG/imputation_validation") # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(95~98)3.png"),width = 1250, height = 700, units = "px", pointsize = 12)
Imp.7 <- ggplot(data=imp.7)+
  geom_point(aes(x=Discharge,y=Imp.SSL,color="補遺後"))+
  geom_point(aes(x=Discharge,y=Suspended.Load,color="補遺前"))+
  geom_line(aes(x=Discharge,y=ratingSSL,color="率定曲線")) + #
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(95%~98%)"))+
  theme(text=element_text(size=20))  # 字體大小
plot(Imp.7)
dev.off()

# ---- group8：98% ~99% ----

# 將一開始移除的40% 加回來
maxBC8 <- max(data.8$Discharge)
MD.8 <- as.matrix(subset(MD[,4:5],Discharge>=maxBC7 & Discharge<maxBC8)) #提取Q與Qs出來，並限制流量範圍
n.marg <- 2 # 兩個變數(Q、Qs)
imp8 <- CoImp(MD.8, n.marg=n.marg, q.up=c(0.5,0.5), q.lo=c(0.01,0.01), 
              smoothing = c(1,1),type.data="continuous", 
              model=list(gumbelCopula(),frankCopula(),claytonCopula())) # 補遺計算
imp.8 <- cbind(subset(ob.data,Discharge>=maxBC7 & Discharge<maxBC8),
               subset(MD[,4:5],Discharge>=maxBC7 & Discharge<maxBC8)[,2],
               imp8@Imputed.data.matrix[,2]) # 提取補遺值，並合併成結果表格
colnames(imp.8) <- c("Year","Month","Day","Discharge",
                     "Suspended.Load","asNA","Imp.SSL") # 為框架的行命名
rating <- nls(Suspended.Load ~ a*Discharge^b,algorithm="port",
              control = list(maxiter = 200,minFactor = 1/200,warnOnly = TRUE),
              start=list(a=1,b=1), data=data.8[,4:5],trace=T) # 乘冪迴歸式
coef(rating)
#save <- data.8[,4:5]
#file <- paste("F:/R_output/CHIA-YUANG/imputation_validation_validation/ttt.csv", sep="") # 存檔路徑
#write.csv(save, file)

summary(rating) # 觀察
a <- environment(rating[["m"]][["resid"]])[["env"]][["a"]] # 迴歸係數a
b <- environment(rating[["m"]][["resid"]])[["env"]][["b"]] # 迴歸係數b
ratingSSL <- a*(imp.8$Discharge)^b # 計算率定曲線推估的輸砂量
imp.8 <- cbind(imp.8,ratingSSL)

rmse <- function(actual, predicted) { # 方均根誤差公式
  sqrt(mean((actual - predicted) ^ 2))
}
# rmse (root mean square error)
rmse.imp <- rmse(imp.8$Suspended.Load, imp.8$Imp.SSL) #原本移除的觀測值與補遺值
rmse.rating <- rmse(imp.8$Suspended.Load, imp.8$ratingSSL) #原本移除的觀測值與率定值
rmse8 <- cbind(rmse.imp,rmse.rating)
# mae (mean absolute error)
mae.imp <- mae(imp.8$Suspended.Load, imp.8$Imp.SSL) #原本移除的觀測值與補遺值
mae.rating <- mae(imp.8$Suspended.Load, imp.8$ratingSSL) #原本移除的觀測值與率定值
mae8 <- cbind(mae.imp, mae.rating)
# mape (mean absolute persentage error)
mape.imp <- mape(imp.8$Suspended.Load, imp.8$Imp.SSL) #原本移除的觀測值與補遺值
mape.rating <- mape(imp.8$Suspended.Load, imp.8$ratingSSL) #原本移除的觀測值與率定值
mape8 <- cbind(mape.imp, mape.rating)
# 補遺資料出圖
setwd("F:/R_output/CHIA-YUANG/imputation_validation") # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(98~99)3.png"),width = 1250, height = 700, units = "px", pointsize = 12)
Imp.8 <- ggplot(data=imp.8)+
  geom_point(aes(x=Discharge,y=Imp.SSL,color="補遺後"))+
  geom_point(aes(x=Discharge,y=Suspended.Load,color="補遺前"))+
  geom_line(aes(x=Discharge,y=ratingSSL,color="率定曲線")) + #
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(98%~99%)"))+
  theme(text=element_text(size=20))  # 字體大小
plot(Imp.8)
dev.off()

# ---- group9：99% ~100% ----

# 將一開始移除的40% 加回來
maxBC9 <- max(data.9$Discharge)
MD.9 <- as.matrix(subset(MD[,4:5],Discharge>=maxBC8 & Discharge<maxBC9)) #提取Q與Qs出來，並限制流量範圍
n.marg <- 2 # 兩個變數(Q、Qs)
imp9 <- CoImp(MD.9, n.marg=n.marg, q.up=c(0.5,0.5), q.lo=c(0.01,0.01), 
              smoothing = c(1,1),type.data="continuous", 
              model=list(gumbelCopula(),frankCopula(),claytonCopula())) # 補遺計算
imp.9 <- cbind(subset(ob.data,Discharge>=maxBC8 & Discharge<maxBC9),
               subset(MD[,4:5],Discharge>=maxBC8 & Discharge<maxBC9)[,2],
               imp9@Imputed.data.matrix[,2]) # 提取補遺值，並合併成結果表格
colnames(imp.9) <- c("Year","Month","Day","Discharge",
                     "Suspended.Load","asNA","Imp.SSL") # 為框架的行命名
rating <- nls(Suspended.Load ~ a*Discharge^b,algorithm="port",
              control = list(maxiter = 200,minFactor = 1/200,warnOnly = TRUE),
              start=list(a=1,b=1), data=data.9[,4:5],trace=T) # 乘冪迴歸式
coef(rating)
#save <- data.9[,4:5]
#file <- paste("F:/R_output/CHIA-YUANG/imputation_validation_validation/ttt.csv", sep="") # 存檔路徑
#write.csv(save, file)

summary(rating) # 觀察
a <- environment(rating[["m"]][["resid"]])[["env"]][["a"]] # 迴歸係數a
b <- environment(rating[["m"]][["resid"]])[["env"]][["b"]] # 迴歸係數b
ratingSSL <- a*(imp.9$Discharge)^b # 計算率定曲線推估的輸砂量
imp.9 <- cbind(imp.9,ratingSSL)

rmse <- function(actual, predicted) { # 方均根誤差公式
  sqrt(mean((actual - predicted) ^ 2))
}
# rmse (root mean square error)
rmse.imp <- rmse(imp.9$Suspended.Load, imp.9$Imp.SSL) #原本移除的觀測值與補遺值
rmse.rating <- rmse(imp.9$Suspended.Load, imp.9$ratingSSL) #原本移除的觀測值與率定值
rmse9 <- cbind(rmse.imp,rmse.rating)
# mae (mean absolute error)
mae.imp <- mae(imp.9$Suspended.Load, imp.9$Imp.SSL) #原本移除的觀測值與補遺值
mae.rating <- mae(imp.9$Suspended.Load, imp.9$ratingSSL) #原本移除的觀測值與率定值
mae9 <- cbind(mae.imp, mae.rating)
# mape (mean absolute persentage error)
mape.imp <- mape(imp.9$Suspended.Load, imp.9$Imp.SSL) #原本移除的觀測值與補遺值
mape.rating <- mape(imp.9$Suspended.Load, imp.9$ratingSSL) #原本移除的觀測值與率定值
mape9 <- cbind(mape.imp, mape.rating)
# 補遺資料出圖
setwd("F:/R_output/CHIA-YUANG/imputation_validation") # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(99~100)3.png"),width = 1250, height = 700, units = "px", pointsize = 12)
Imp.9 <- ggplot(data=imp.9)+
  geom_point(aes(x=Discharge,y=Imp.SSL,color="補遺後"))+
  geom_point(aes(x=Discharge,y=Suspended.Load,color="補遺前"))+
  geom_line(aes(x=Discharge,y=ratingSSL,color="率定曲線")) + #
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(99%~100%)"))+
  theme(text=element_text(size=20))  # 字體大小
plot(Imp.9)
dev.off()

# ----- all groups together --------

imp.sum <- rbind(imp.1,imp.2,imp.3,imp.4,imp.5,imp.6,imp.7,imp.8,imp.9)
file <- paste("F:/R_output/CHIA-YUANG/imputation_validation/imp_sum3.csv", sep="") # 存檔路徑
write.csv(imp.sum, file)
# rmse (root mean square error)
rmse.imp <- rmse(imp.sum$Suspended.Load, imp.sum$Imp.SSL) #原本移除的觀測值與補遺值
rmse.rating <- rmse(imp.sum$Suspended.Load, imp.sum$ratingSSL) #原本移除的觀測值與率定值
rmse.sum <- cbind(rmse.imp,rmse.rating)
# mae (mean absolute error)
mae.imp <- mae(imp.sum$Suspended.Load, imp.sum$Imp.SSL) #原本移除的觀測值與補遺值
mae.rating <- mae(imp.sum$Suspended.Load, imp.sum$ratingSSL) #原本移除的觀測值與率定值
mae.sum <- cbind(mae.imp, mae.rating)
# mape (mean absolute persentage error)
mape.imp <- mape(imp.sum$Suspended.Load, imp.sum$Imp.SSL) #原本移除的觀測值與補遺值
mape.rating <- mape(imp.sum$Suspended.Load, imp.sum$ratingSSL) #原本移除的觀測值與率定值
mape.sum <- cbind(mape.imp, mape.rating)

# 補遺資料出圖
setwd("F:/R_output/CHIA-YUANG/imputation_validation") # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(all group)3.png"),width = 1250, height = 700, units = "px", pointsize = 12)
Imp.sum <- ggplot(data=imp.sum)+
  geom_point(aes(x=Discharge,y=Imp.SSL,color="補遺值"))+
  geom_point(aes(x=Discharge,y=Suspended.Load,color="觀測值"))+
  geom_line(aes(x=Discharge,y=ratingSSL,color="率定曲線")) + #
  geom_vline(xintercept=c(maxBC1,maxBC2,maxBC3,maxBC4,maxBC5,maxBC6,maxBC7,maxBC8)) +
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(all group)"))+
  theme(text=element_text(size=20))  # 字體大小
plot(Imp.sum)
dev.off()

# ----- 表格整理 -----

rmse.table <- rbind(rmse1,rmse2,rmse3,rmse4,rmse5,rmse6,rmse7,rmse8,rmse9,rmse.sum)
rownames(rmse.table) <-c("group1","group2","group3","group4","group5",
                         "group6","group7","group8","group9","total")
file <- paste("F:/R_output/CHIA-YUANG/imputation_validation/imp_rmse3.csv", sep="") # 存檔路徑
write.csv(rmse.table, file)

mae.table <- rbind(mae1,mae2,mae3,mae4,mae5,mae6,mae7,mae8,mae9,mae.sum)
rownames(mae.table) <-c("group1","group2","group3","group4","group5",
                         "group6","group7","group8","group9","total")
file <- paste("F:/R_output/CHIA-YUANG/imputation_validation/imp_mae3.csv", sep="") # 存檔路徑
write.csv(mae.table, file)

mape.table <- rbind(mape1,mape2,mape3,mape4,mape5,mape6,mape7,mape8,mape9,mape.sum)
rownames(mape.table) <-c("group1","group2","group3","group4","group5",
                         "group6","group7","group8","group9","total")
file <- paste("F:/R_output/CHIA-YUANG/imputation_validation/imp_mape3.csv", sep="") # 存檔路徑
write.csv(mape.table, file)
