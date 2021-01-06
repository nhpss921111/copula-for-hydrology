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
# =========== 不分組 (建立 rating curve)===============
rating <- nls(Suspended.Load ~ a*Discharge^b,  algorithm="port",
              control = list(maxiter = 200,minFactor = 1/2000000,warnOnly = TRUE),
              start=list(a=0.1,b=0.1), data=MD60,trace=T) # 乘冪迴歸式
summary(rating) # 觀察
a <- environment(rating[["m"]][["resid"]])[["env"]][["a"]] # 迴歸係數a
b <- environment(rating[["m"]][["resid"]])[["env"]][["b"]] # 迴歸係數b
ratingSSL <- a*(MD$Discharge)^b # 計算率定曲線推估的輸砂量
rating.1 <- cbind(subset(ob.data),subset(MD[,4:5])[,2],ratingSSL)
colnames(rating.1) <- c("Year","Month","Day","Discharge",
                         "Suspended.Load","asNA","ratingSSL1")
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
# 補遺資料出圖
setwd("F:/R_output/CHIA-YUANG/imputation_validation") # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(不分組).png"),width = 1250, height = 700, units = "px", pointsize = 12)
Rating.1 <- ggplot(data=rating.1)+
  geom_point(aes(x=Discharge,y=asNA,color="觀測資料")) +
  geom_line(aes(x=Discharge,y=ratingSSL,color="率定曲線")) +
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(不分組)")) +
  theme(text=element_text(size=20))  # 字體大小
plot(Rating.1)
dev.off()

# =========== 分九組 ===============
# 分組數
group.BC <- c(0,max(data.1$Discharge),max(data.2$Discharge),max(data.3$Discharge),
              max(data.4$Discharge),max(data.5$Discharge),max(data.6$Discharge),
              max(data.7$Discharge),max(data.8$Discharge),max(data.9$Discharge))
persent.BC <- c(0,20,40,60,80,90,95,98,99,100)
imp.table <- c() # 儲存資料
copula.parameter.table <- c() # 儲存資料
loss.table <- c() # 儲存資料
# ---- 分組迴圈  ----
# 將一開始移除的40% 加回來
for (g in 1:(length(group.BC)-1)){
  #maxBC1 <- max(data.1$Discharge)
  MD.bygroup <- as.matrix(subset(MD[,4:5],Discharge>group.BC[g] & Discharge<=group.BC[g+1])) #提取Q與Qs出來，並限制流量範圍
  colnames(MD.bygroup) <- c("Discharge","Suspended.Load")
  n.marg <- 2 # 兩個變數(Q、Qs)
  imp <- CoImp(MD.bygroup, n.marg=n.marg, q.up=c(0.5,0.5), q.lo=c(0.01,0.01), 
                smoothing = c(1,1),type.data="continuous", 
                model=list(gumbelCopula(),frankCopula(),claytonCopula())) # 補遺計算
  plot(imp)
  copula.func <- imp@Estimated.Model.Imp[["model"]]
  copula.para <- imp@Estimated.Model.Imp[["parameter"]]
  copula.list <- cbind(copula.func, copula.para)
  rownames(copula.list) <- paste0("group",g)
  copula.parameter.table <- rbind(copula.parameter.table,copula.list) # 傳承
  
  imp.group <- cbind(subset(ob.data,Discharge>group.BC[g] & Discharge<=group.BC[g+1]),
                 subset(MD[,4:5],Discharge>group.BC[g] & Discharge<=group.BC[g+1])[,2],
                 imp@Imputed.data.matrix[,2]) # 提取補遺值，並合併成結果表格
  
  colnames(imp.group) <- c("Year","Month","Day","Discharge",
                           "Suspended.Load","asNA","Imp.SSL") # 為框架的行命名
  rating <- nls(Suspended.Load ~ a*Discharge^b, algorithm="port",
                control = list(maxiter = 200, minFactor = 1/2^30, warnOnly = TRUE),
                start=list(a=0.1,b=0.1), data=as.data.frame(MD.bygroup[complete.cases(MD.bygroup), ]),
                trace=T) # 乘冪迴歸式
  summary(rating) # 觀察
  a <- environment(rating[["m"]][["resid"]])[["env"]][["a"]] # 迴歸係數a
  b <- environment(rating[["m"]][["resid"]])[["env"]][["b"]] # 迴歸係數b
  ratingSSL <- a*(imp.group$Discharge)^b # 計算率定曲線推估的輸砂量
  imp.group <- cbind(imp.group,ratingSSL)
  imp.table <- rbind(imp.table,imp.group) # 傳承
  
  rmse <- function(actual, predicted) { # 均方根誤差公式
    sqrt(mean((actual - predicted) ^ 2))
  }
  # rmse (root mean square error)
  rmse.imp <- rmse(imp.group$Suspended.Load, imp.group$Imp.SSL) #原本移除的觀測值與補遺值
  rmse.rating <- rmse(imp.group$Suspended.Load, imp.group$ratingSSL) #原本移除的觀測值與率定值
  
  # mae (mean absolute error)
  mae.imp <- mae(imp.group$Suspended.Load, imp.group$Imp.SSL) #原本移除的觀測值與補遺值
  mae.rating <- mae(imp.group$Suspended.Load, imp.group$ratingSSL) #原本移除的觀測值與率定值
  
  # mape (mean absolute persentage error)
  mape.imp <- mape(imp.group$Suspended.Load, imp.group$Imp.SSL) #原本移除的觀測值與補遺值
  mape.rating <- mape(imp.group$Suspended.Load, imp.group$ratingSSL) #原本移除的觀測值與率定值
  
  # loss function
  loss <- cbind(mae.imp, mae.rating,mape.imp, mape.rating,rmse.imp, rmse.rating)
  rownames(loss) <- paste0("group",g)
  loss.table <- rbind(loss.table,loss)
  # 補遺資料出圖
  setwd("F:/R_output/CHIA-YUANG/imputation_validation") # 請修改儲存路徑：
  png(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證group",g,".png"),width = 1250, height = 700, units = "px", pointsize = 12)
  Imp.group <- ggplot(data=imp.group)+
    geom_point(aes(x=Discharge,y=Imp.SSL,color="補遺後")) +
    geom_point(aes(x=Discharge,y=Suspended.Load,color="補遺前")) +
    geom_line(aes(x=Discharge,y=ratingSSL,color="率定曲線")) +
    scale_color_discrete(name="圖例") + #圖例名稱
    ggtitle(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證",persent.BC[g],"% ~ ",persent.BC[g+1],"%")) +
    theme(text=element_text(size=20))  # 字體大小
  plot(Imp.group)
  dev.off()
}

# ------------- 把分組合併 -------------------
# rmse (root mean square error)

rmse.rating <- rmse(imp.table$Suspended.Load, imp.table$ratingSSL) #原本移除的觀測值與率定值

# mae (mean absolute error)

mae.rating <- mae(imp.table$Suspended.Load, imp.table$ratingSSL) #原本移除的觀測值與率定值

# mape (mean absolute persentage error)

mape.rating <- mape(imp.table$Suspended.Load, imp.table$ratingSSL) #原本移除的觀測值與率定值

rating.bygroup <- cbind(rmse.rating,mae.rating,mape.rating)
rating.cf <- rbind(rating.table,rating.bygroup)
rownames(rating.cf) <- c("不分組","有分組")

# 組數合併資料出圖
setwd("F:/R_output/CHIA-YUANG/imputation_validation") # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(all group).png"),width = 1250, height = 700, units = "px", pointsize = 12)
Imp.sum <- ggplot(data=imp.table)+
  #geom_point(aes(x=Discharge,y=Imp.SSL,color="補遺值"))+
  geom_point(aes(x=Discharge,y=asNA,color="觀測值"))+
  geom_line(aes(x=Discharge,y=ratingSSL,color="率定曲線")) + #
  geom_vline(xintercept=group.BC) +
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(all group)"))+
  theme(text=element_text(size=20))  # 字體大小
plot(Imp.sum)
dev.off()

# 率定曲線不分組 vs 率定曲線有分組
rating.sum <- full_join(imp.table,rating.1)
#rating.sum <- cbind(,$ratingSSL)
colnames(rating.sum) <- c("Year","Month","Day","Discharge",
                        "Suspended.Load","asNA","ImpSSL","ratingSSL1","ratingSSL2")
setwd("F:/R_output/CHIA-YUANG/imputation_validation") # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(no group vs all group ).png"),width = 1250, height = 700, units = "px", pointsize = 12)
Rating.cf <- ggplot(data=rating.sum)+
  #geom_point(aes(x=Discharge,y=Imp.SSL,color="補遺值"))+
  geom_point(aes(x=Discharge,y=asNA,color="觀測值"))+
  geom_line(aes(x=Discharge,y=ratingSSL1,color="率定曲線(有分組)")) + #
  geom_line(aes(x=Discharge,y=ratingSSL2,color="率定曲線(不分組)")) + #
  #geom_vline(xintercept=group.BC) +
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年家源橋測站補遺驗證(all group)"))+
  theme(text=element_text(size=20))  # 字體大小
plot(Rating.cf)
dev.off()
