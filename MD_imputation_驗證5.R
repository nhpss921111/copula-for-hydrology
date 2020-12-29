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

# ===========
# 家源橋CHIA-YUANG year <- c(1974:2009,2012:2019)
# 彰雲橋CHUNYUN BRIDGE year <- c(1987:2019)
# ===========

year <- c(1974:2009,2012:2019)
MD.input <- c(paste0(year,"QandQs.csv"))
output <- c(paste0(year,"imp.csv"))
ob.data <- c()
set.seed(100)
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
perc.mis    <- 0.4 # 40%的資料當成NA 
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

# ------------------------------- Imputation ----------------------------------------

# ---- group1：2yr ----
# 將一開始移除的40% 加回來
maxBC1 <- max(data.1$Discharge)
MD.1 <- as.matrix(subset(MD[,4:5],Discharge<maxBC1)) #提取Q與Qs出來，並限制流量範圍
n.marg <- 2 # 兩個變數(Q、Qs)
imp1 <- CoImp(MD.1, n.marg=n.marg, q.up=c(0.5,0.5), q.lo=c(0.01,0.01), 
             smoothing = c(1,1),type.data="continuous", 
             model=list(gumbelCopula(),frankCopula(),claytonCopula())) # 補遺計算
imp.1 <- cbind(subset(ob.data,Discharge<maxBC1),
               subset(MD[,4:5],Discharge<maxBC1)[,2],
               imp1@Imputed.data.matrix[,2]) # 提取補遺值，並合併成結果表格
colnames(imp.1) <- c("Year","Month","Day","Discharge",
                     "Suspended.Load","asNA","Imp.SSL") # 為框架的行命名
rating <- nls(Suspended.Load ~ a*Discharge^b, 
              start=list(a=1,b=1), data=data.1[,4:5]) # 乘冪迴歸式
summary(rating) # 觀察
a <- environment(rating[["m"]][["resid"]])[["env"]][["a"]] # 迴歸係數a
b <- environment(rating[["m"]][["resid"]])[["env"]][["b"]] # 迴歸係數b
ratingSSL <- a*(imp.1$Discharge)^b # 計算率定曲線推估的輸砂量
imp.1 <- cbind(imp.1,ratingSSL)

rmse <- function(actual, predicted) { # 均方根誤差公式
  sqrt(mean((actual - predicted) ^ 2))
}
rmse.imp <- rmse(imp.1$Suspended.Load, imp.1$Imp.SSL) #原本移除的觀測值與補遺值
rmse.rating <- rmse(imp.1$Suspended.Load, imp.1$ratingSSL) #原本移除的觀測值與率定值
plot(rmse.imp,rmse.rating)
# 補遺資料出圖
#setwd("F:/R_output/CHUNYUN BRIDGE/imputation") # 請修改儲存路徑：
#png(paste0(year[1],"到",year[y],"驗證-1-no1.png"),width = 1250, height = 700, units = "px", pointsize = 12)
Imp.1 <- ggplot(data=imp.1)+
  geom_point(aes(x=Discharge,y=Imp.SSL,color="補遺後"))+
  geom_point(aes(x=Discharge,y=Suspended.Load,color="補遺前"))+
  stat_smooth(data=data.1[,4:5],
              aes(x=Discharge,y=Suspended.Load,color="率定曲線"),
              method = 'nls', formula = 'y~a*x^b', 
              start=c(a = 0.1,b=0.1),se=FALSE) +  #
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年家源橋測站補遺"))+
  theme(text=element_text(size=20))  # 字體大小
plot(Imp.1)
#dev.off()
  
# ---- group2：5yr ----

MD.5yr <- as.matrix(data.5yr[,4:5])
n.marg <- 2
imp5yr <- CoImp(MD.5yr, n.marg=n.marg, q.up=c(0.5,0.5), q.lo=c(0.01,0.01), 
             smoothing = c(1,1),type.data="continuous", 
             model=list(gumbelCopula(),frankCopula(),claytonCopula()))

imp.5yr <- cbind(data.5yr,imp5yr@Imputed.data.matrix[,2])
imp.5yr <- imp.5yr[,-6:-8]
colnames(imp.5yr) <- c("Year","Month","Day","Discharge","Suspended.Load","Imp.SSL")


# 補遺資料出圖
#setwd("F:/R_output/CHUNYUN BRIDGE/imputation") # 請修改儲存路徑：
#png(paste0(year[1],"到",year[y],"驗證-1-no1.png"),width = 1250, height = 700, units = "px", pointsize = 12)
Imp.5yr <- ggplot()+
  geom_point(data=imp.5yr,aes(x=Discharge,y=Imp.SSL,color="補遺值"))+
  geom_point(data=imp.5yr,aes(x=Discharge,y=Suspended.Load,color="觀測值"))+
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年家源橋測站補遺"))+
  theme(text=element_text(size=20))  # 字體大小
plot(Imp.5yr)
#dev.off()

# ---- group3：10yr ----

MD.10yr <- as.matrix(data.10yr[,4:5])
n.marg <- 2
imp10yr <- CoImp(MD.10yr, n.marg=n.marg, q.up=c(0.5,0.5), q.lo=c(0.01,0.01), 
                smoothing = c(1,1),type.data="continuous", 
                model=list(gumbelCopula(),frankCopula(),claytonCopula()))

imp.10yr <- cbind(data.10yr,imp10yr@Imputed.data.matrix[,2])
imp.10yr <- imp.10yr[,-6:-8]
colnames(imp.10yr) <- c("Year","Month","Day","Discharge","Suspended.Load","Imp.SSL")


# 補遺資料出圖
#setwd("F:/R_output/CHUNYUN BRIDGE/imputation") # 請修改儲存路徑：
#png(paste0(year[1],"到",year[y],"驗證-1-no1.png"),width = 1250, height = 700, units = "px", pointsize = 12)
Imp.10yr <- ggplot()+
  geom_point(data=imp.10yr,aes(x=Discharge,y=Imp.SSL,color="補遺值"))+
  geom_point(data=imp.10yr,aes(x=Discharge,y=Suspended.Load,color="觀測值"))+
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年家源橋測站補遺"))+
  theme(text=element_text(size=20))  # 字體大小
plot(Imp.10yr)
#dev.off()

# ---- group4：25yr ----

MD.25yr <- as.matrix(data.25yr[,4:5])
n.marg <- 2
imp25yr <- CoImp(MD.25yr, n.marg=n.marg, q.up=c(0.5,0.5), q.lo=c(0.01,0.01), 
                 smoothing = c(1,1),type.data="continuous", 
                 model=list(gumbelCopula(),frankCopula(),claytonCopula()))

imp.25yr <- cbind(data.25yr,imp25yr@Imputed.data.matrix[,2])
imp.25yr <- imp.25yr[,-6:-8]
colnames(imp.25yr) <- c("Year","Month","Day","Discharge","Suspended.Load","Imp.SSL")


# 補遺資料出圖
#setwd("F:/R_output/CHUNYUN BRIDGE/imputation") # 請修改儲存路徑：
#png(paste0(year[1],"到",year[y],"驗證-1-no1.png"),width = 1250, height = 700, units = "px", pointsize = 12)
Imp.25yr <- ggplot()+
  geom_point(data=imp.25yr,aes(x=Discharge,y=Imp.SSL,color="補遺值"))+
  geom_point(data=imp.25yr,aes(x=Discharge,y=Suspended.Load,color="觀測值"))+
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年家源橋測站補遺"))+
  theme(text=element_text(size=20))  # 字體大小
plot(Imp.25yr)
#dev.off()

# ---- group5：50yr ----

MD.50yr <- as.matrix(data.50yr[,4:5])
n.marg <- 2
imp50yr <- CoImp(MD.50yr, n.marg=n.marg, q.up=c(0.5,0.5), q.lo=c(0.01,0.01), 
                 smoothing = c(1,1),type.data="continuous", 
                 model=list(gumbelCopula(),frankCopula(),claytonCopula()))

imp.50yr <- cbind(data.50yr,imp50yr@Imputed.data.matrix[,2])
imp.50yr <- imp.50yr[,-6:-8]
colnames(imp.50yr) <- c("Year","Month","Day","Discharge","Suspended.Load","Imp.SSL")


# 補遺資料出圖
#setwd("F:/R_output/CHUNYUN BRIDGE/imputation") # 請修改儲存路徑：
#png(paste0(year[1],"到",year[y],"驗證-1-no1.png"),width = 1250, height = 700, units = "px", pointsize = 12)
Imp.50yr <- ggplot()+
  geom_point(data=imp.50yr,aes(x=Discharge,y=Imp.SSL,color="補遺值"))+
  geom_point(data=imp.50yr,aes(x=Discharge,y=Suspended.Load,color="觀測值"))+
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年家源橋測站補遺"))+
  theme(text=element_text(size=20))  # 字體大小
plot(Imp.50yr)
#dev.off()

# ---- group6：100yr ----

MD.100yr <- as.matrix(data.100yr[,4:5])
n.marg <- 2
imp100yr <- CoImp(MD.100yr, n.marg=n.marg, q.up=c(0.5,0.5), q.lo=c(0.01,0.01), 
                 smoothing = c(1,1),type.data="continuous", 
                 model=list(gumbelCopula(),frankCopula(),claytonCopula()))

imp.100yr <- cbind(data.100yr,imp100yr@Imputed.data.matrix[,2])
imp.100yr <- imp.100yr[,-6:-8]
colnames(imp.100yr) <- c("Year","Month","Day","Discharge","Suspended.Load","Imp.SSL")


# 補遺資料出圖
#setwd("F:/R_output/CHUNYUN BRIDGE/imputation") # 請修改儲存路徑：
#png(paste0(year[1],"到",year[y],"驗證-1-no1.png"),width = 1250, height = 700, units = "px", pointsize = 12)
Imp.100yr <- ggplot()+
  geom_point(data=imp.100yr,aes(x=Discharge,y=Imp.SSL,color="補遺值"))+
  geom_point(data=imp.100yr,aes(x=Discharge,y=Suspended.Load,color="觀測值"))+
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年家源橋測站補遺"))+
  theme(text=element_text(size=20))  # 字體大小
plot(Imp.100yr)
#dev.off()


# ---- group7：200yr ----

MD.200yr <- as.matrix(data.200yr[,4:5])
n.marg <- 2
imp200yr <- CoImp(MD.200yr, n.marg=n.marg, q.up=c(0.5,0.5), q.lo=c(0.01,0.01), 
                  smoothing = c(1,1),type.data="continuous", 
                  model=list(gumbelCopula(),frankCopula(),claytonCopula()))

imp.200yr <- cbind(data.200yr,imp200yr@Imputed.data.matrix[,2])
imp.200yr <- imp.200yr[,-6:-8]
colnames(imp.200yr) <- c("Year","Month","Day","Discharge","Suspended.Load","Imp.SSL")


# 補遺資料出圖
#setwd("F:/R_output/CHUNYUN BRIDGE/imputation") # 請修改儲存路徑：
#png(paste0(year[1],"到",year[y],"驗證-1-no1.png"),width = 1250, height = 700, units = "px", pointsize = 12)
Imp.200yr <- ggplot()+
  geom_point(data=imp.200yr,aes(x=Discharge,y=Imp.SSL,color="補遺值"))+
  geom_point(data=imp.200yr,aes(x=Discharge,y=Suspended.Load,color="觀測值"))+
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年家源橋測站補遺"))+
  theme(text=element_text(size=20))  # 字體大小
plot(Imp.200yr)
#dev.off()

imp.sum <- rbind(imp.2yr,imp.5yr,imp.10yr,imp.25yr,imp.50yr,imp.100yr,imp.200yr)

Imp.sum <- ggplot()+
  geom_point(data=imp.sum,aes(x=Discharge,y=Imp.SSL,color="補遺值"))+
  geom_point(data=imp.sum,aes(x=Discharge,y=Suspended.Load,color="觀測值"))+
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年家源橋測站補遺"))+
  theme(text=element_text(size=20))  # 字體大小
plot(Imp.sum)


file <- paste("F:/R_output/CHIA-YUANG/imputation/",year[1],"到", year[y],"imp成果-1-no2.csv", sep="") # 存檔路徑
write.csv(imp.sum, file)

summary(imp.sum$Discharge)
