# 邊際分布：parametric method
# copula函數：IFM method
# 開始撰寫日期：2020/02/16
# 完成撰寫日期：2021/06/09
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
library(nls2)
library(kdensity)
library(spatstat)
library(Metrics) # 評估指標
library(ggforce)
library(numDeriv)
library(latex2exp) # labtex語法
library(Hmisc)
library(ggprism)
library(patchwork)
library(ggrepel)
library(dynatopmodel) #Nash-Sutcliffe Efficiency
library(hydroGOF)
library(ggpmisc) # ggplot peak point
# =======================================================================
# 家源橋("CHIA-YUANG")：       year <- c(1974:2019)            都沒過
# 彰雲橋("CHUNYUN BRIDGE")：   year <- c(1987:2019)            copula沒過CVM
# 內茅埔("NEI-MAO-PU")：       year <- c(1972:2001,2003:2019)  Q沒過copula沒過
# 仁壽橋("JEN-SHOU BRIDGE")：  year <- c(1960:2019)            Q沒過
# 六龜("LIU-KWEI")：           year <- c(1982:2009,2011:2019)  Q沒過copula沒過
# 內灣("NEI-WAN")：            year <- c(1971:2005,2009:2019)  Qs沒過copula沒過
# 經國橋("JEIN-KUO BRIDGE")：  year <- c(1990:2005,2009:2019)  Qs沒過copula沒過
# 蘭陽大橋("LAN-YANG BRIDGE")：year <- c(1949:1957,1959:2017,2019) Q、Qs沒過copula沒過CVM
# 橫溪("HENG CHI")：           year <- c(1974:2004,2006:2019)  都沒過
# 秀朗("HSIU-LUNG")：          year <- c(1970:2004,2006:2019)  copula沒過
# 義里("I-LI")：               year <- c(1966:2003,2006:2010,2013:2019) Qs沒過copula沒過
# 荖濃(新發大橋)("LAO-NUNG")： year <- c(1956:2009)            Qs沒過copula沒過
# 里嶺大橋("LI-LIN BRIDGE")：  year <- c(1991:2004,2007:2019)  copula沒過CVM
# ======================================================================
station <- c("JEN-SHOU BRIDGE") # 測站名稱
station_ch <-c("仁壽橋")
year <- c(1960:2019) # 年分

group.number <- c(9) # 分組的組數
log.group.number <- c(9) # 分組的組數
set.seed(100)
perc.mis <- 0.3 # 多少%的資料當成NA
build <- ("deterministic.1") # deterministic / stochastic
MD.input <- c(paste0(year,"QandQs.csv"))
output <- c(paste0(year,"imp.csv"))
ob.data <- c()
# ---- 誤差指標公式 ----
mse <- function(actual, predicted) { # 均方誤差
  mean((actual - predicted) ^ 2)
}
rmse <- function(actual, predicted) { # 均方根誤差
  sqrt(mean((actual - predicted) ^ 2))
}
nmse <- function(actual, predicted) { # 正歸化均方誤差
  mean((actual - predicted) ^ 2)/mean((actual - mean(actual)) ^ 2)
}
mape <- function(actual, predicted) { # 平均絕對百分誤差
  100*mean(abs((actual - predicted)/actual))
}
nse <- function(sim, obs){ # Nash-Sutcliffe Efficiency
  1-sum((obs-sim)^2)/sum((obs-mean(obs))^2)
}
nse.abs <- function(sim, obs){ # Nash-Sutcliffe Efficiency
  1-sum(abs(obs-sim))/sum(abs(obs-mean(obs)))
}
#
# --------------
# 主要迴圈(以年份為底)
for( y in 1:length(year)){
  # 讀取有缺失的資料
  setwd(paste0("F:/R_reading/",station,"/missingdata"))
  data <- read.csv(file.path(getwd(),MD.input[y]))
  data <- data[,-1]
  ob.data <- rbind(ob.data,data)
}
ob.data <- unique(ob.data) # 移除重複的觀測資料
ob.data <- subset(ob.data,ob.data$Discharge>0) # 把流量為0(無觀測資料)刪除
ob.data[ob.data==0] <- NA #將輸砂量0的資料當成NA
#ob.data$Suspended.Load <- ob.data$Suspended.Load/1000 #將公噸換算成千噸
log.data <- cbind(ob.data[,1:3],log10(ob.data[,4:5]))
rm.ob.data <- ob.data[complete.cases(ob.data), ] # 移除原始觀測資料中全部NA
file <- paste("F:/R_output/",station,"/vinecopula/",
              year[1],"到", year[y],station_ch,"流量&輸砂量(100%觀測資料).csv", sep="") #存檔路徑
write.csv(rm.ob.data,file)
rm.log.data <- log.data[complete.cases(log.data), ] # 移除觀測資料取對數後全部NA
#ob.data <- subset(ob.data,ob.data[,4]>20 & ob.data[,5]>1000)
summary(rm.ob.data)
# ----------- 將同時有Q與Qs的資料分兩組 (80%資料總數建模，剩下20%當成驗證) -------------
if (build == "stochastic"){
  x.samp <- as.matrix(rm.ob.data[,4:5])
  x.samp.log <- as.matrix(rm.log.data[,4:5])
  miss.row    <- sample(1:length(rm.ob.data$Discharge), perc.mis*length(rm.ob.data$Discharge), replace=FALSE)
  miss.col    <- rep(2,perc.mis*length(rm.ob.data$Discharge))
  miss        <- cbind(miss.row,miss.col) # NA的欄位座標
  samp.miss <- replace(x.samp,miss,NA) # 將 ?% 觀測資料轉換成NA
  samp.miss.log <- replace(x.samp.log,miss,NA) # 將 ?% 觀測資料轉換成NA
  MD <- cbind(rm.ob.data[,1:3],samp.miss) 
  MD.log <- cbind(rm.log.data[,1:3],samp.miss.log) 
  MD.withNA <- MD[,4:5]
  MD.log.withNA <- MD.log[,4:5]
  MD.rmNA <- MD[complete.cases(MD), ] # 移除全部NA (剩餘資料)
  MD.log.rmNA <- MD.log[complete.cases(MD.log), ] # 移除全部NA (剩餘資料)
  file <- paste("F:/R_output/",station,"/vinecopula/",
                year[1],"到", year[y],station_ch,"率定模式(70%觀測資料).csv", sep="") #存檔路徑
  write.csv(MD.rmNA,file)
}
if (build == "deterministic"){
  x.samp <- as.matrix(rm.ob.data[,4:5])
  miss.year.start <- round(max(rm.ob.data$Year)-(max(rm.ob.data$Year)-min(rm.ob.data$Year)+1)*perc.mis+1)
  # --- 原始尺度
  miss <- subset(rm.ob.data,Year>=miss.year.start)
  miss.data <- cbind(miss[,-5],Suspended.Load=NA)
  MD <- rbind(subset(rm.ob.data,Year<miss.year.start),miss.data)
  MD.withNA <- MD[,4:5]
  MD.rmNA <- MD[complete.cases(MD), ] # 移除全部NA (剩餘資料)
  # --- 對數尺度
  miss.log <- subset(rm.log.data,Year>=miss.year.start)
  miss.data.log <- cbind(miss.log[,-5],Suspended.Load=NA)
  MD.log <- rbind(subset(rm.log.data,Year<miss.year.start),miss.data.log)
  MD.log.withNA <- MD.log[,4:5]
  MD.log.rmNA <- MD.log[complete.cases(MD.log), ] # 移除全部NA (剩餘資料)
  # --- 原始尺度存檔
  file <- paste("F:/R_output/",station,"/vinecopula/",
                year[1],"到", year[y],station_ch,"率定模式(70%觀測資料).csv", sep="") #存檔路徑
  write.csv(MD.rmNA,file)
  
}
if (build == "deterministic.1"){
  x.samp <- as.matrix(rm.ob.data[,4:5])
  miss.start <- round(length(rm.ob.data$Discharge)*(1-perc.mis))
  miss.year.start <- rm.ob.data[miss.start,]$Year
  miss.month.start <- rm.ob.data[miss.start,]$Month
  miss.day.start <- rm.ob.data[miss.start,]$Day
  year.data <- subset(rm.ob.data, Year==miss.year.start)
  month.data <- subset(year.data, Month<=miss.month.start)
  i <- (length(year.data$Discharge)-length(month.data$Discharge))/length(year.data$Discharge)
  if(i<=0.5){
    miss.year.start <- miss.year.start+1
  }
  # --- 原始尺度
  miss <- subset(rm.ob.data,Year>=miss.year.start)
  build.perc <- length(miss$Discharge)/length(rm.ob.data$Discharge)
  miss.data <- cbind(miss[,-5],Suspended.Load=NA)
  MD <- rbind(subset(rm.ob.data,Year<miss.year.start),miss.data)
  MD.withNA <- MD[,4:5]
  MD.rmNA <- MD[complete.cases(MD), ] # 移除全部NA (剩餘資料)
  # --- 對數尺度
  miss.log <- subset(rm.log.data,Year>=miss.year.start)
  miss.data.log <- cbind(miss.log[,-5],Suspended.Load=NA)
  MD.log <- rbind(subset(rm.log.data,Year<miss.year.start),miss.data.log)
  MD.log.withNA <- MD.log[,4:5]
  MD.log.rmNA <- MD.log[complete.cases(MD.log), ] # 移除全部NA (剩餘資料)
  # --- 原始尺度存檔
  file <- paste("F:/R_output/",station,"/vinecopula/",
                year[1],"到", year[y],station_ch,"率定模式(70%觀測資料).csv", sep="") #存檔路徑
  write.csv(MD.rmNA,file)
  
}
# 分組的目的：讓CoImp計算用
if (group.number==9){
  #  -------------- 決定流量分組範圍 (來源：log.data) ---------------
  # group 1： 0% ~  20%
  # group 2：20% ~  40%
  # group 3：40% ~  60%
  # group 4：60% ~  80%
  # group 5：80% ~  90%
  # group 6：90% ~  95%
  # group 8：95% ~  98%
  # group 8：98% ~  100%
  
  rank.data <- cbind(rm.ob.data,rank(rm.ob.data$Discharge))
  persent <- (rank.data$`rank(rm.ob.data$Discharge)`) / length(rank.data$Discharge)
  per.data <- cbind(rank.data,persent)
  
  data.1 <- data.frame(subset(per.data, persent<=0.2),group="group1")
  data.2 <- data.frame(subset(per.data, persent>0.2 & persent<=0.4),group="group2")
  data.3 <- data.frame(subset(per.data, persent>0.4 & persent<=0.6),group="group3")
  data.4 <- data.frame(subset(per.data, persent>0.6 & persent<=0.8),group="group4")
  data.5 <- data.frame(subset(per.data, persent>0.8 & persent<=0.9),group="group5")
  data.6 <- data.frame(subset(per.data, persent>0.9 & persent<=0.96),group="group6")
  data.7 <- data.frame(subset(per.data, persent>0.96 & persent<=0.98),group="group7")
  data.8 <- data.frame(subset(per.data, persent>0.98),group="group8")
  
  data.group <- rbind(data.1, data.2, data.3,data.4, 
                      data.5, data.6,data.7, data.8)
  
  group.BC <- c(0,max(data.1$Discharge),max(data.2$Discharge),max(data.3$Discharge),
                max(data.4$Discharge),max(data.5$Discharge),max(data.6$Discharge),
                max(data.7$Discharge),max(data.8$Discharge))
  persent.BC <- c(0,20,40,60,80,90,95,98,100)
}

ggplot(data=MD.rmNA)+
  geom_point(aes(x=Discharge,y=Suspended.Load))+
  scale_color_discrete(name="年",labels=c(""))+
  labs(x="日流量Q (cms)",y="日懸浮載輸砂量L (公噸/日)") + # 座標軸名稱
  #geom_mark_ellipse(expand = 0,aes(fill=group))+ # 橢圓形圈
  #geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=group))+ # 多邊形
  #geom_mark_hull(expand=0.01,aes(fill=group))+ # 凹凸多邊形
  ggtitle(paste0(station_ch,"測站建立模型")) +
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=40))  # 字體大小

# =========== 不分組 (建立 rating curve)===============
# 從"率定曲線初始值.csv"找到a,b的初始值再帶入
file <- paste("F:/R_output/",station,"/vinecopula/",
              year[1],"到", year[y],"年不分組率定曲線初始值查找.csv", sep="") #存檔路徑
write.csv(MD.rmNA,file)
# 從線性模型計算log10a,b 初始值
rating <- lm(Suspended.Load ~ Discharge,data=MD.log.rmNA)
summary(rating) # 觀察

# 將線性計算的係數當成初始值，再帶入非線性模型中
log10.a <- rating$coefficients[1] # 迴歸係數log10(a)
a.start <- 10^log10.a # 線性迴歸係數log10(a) -> 非線性迴歸係數a
b.start <- rating$coefficients[2] # 線性迴歸係數b -> 非線性迴歸係數b

rating <- nls(Suspended.Load ~ a*Discharge^b, algorithm="default",
              control = list(maxiter = 200,minFactor = 1/2^300,warnOnly = TRUE),
              start=list(a=a.start,b=b.start), data=MD.rmNA,trace=T) # 初始值要給好!!!

summary(rating)
a <- environment(rating[["m"]][["resid"]])[["env"]][["a"]]
b <- environment(rating[["m"]][["resid"]])[["env"]][["b"]]
rating.par <- cbind(a,b)
file <- paste("F:/R_output/",station,"/vinecopula/",
              year[1],"到", year[y],"年不分組率定曲線係數.csv", sep="") #存檔路徑
write.csv(rating.par,file)
ratingSSL <- a *(MD$Discharge)^b # 計算率定曲線推估的輸砂量
rating.all <- cbind(rm.ob.data,MD[,5],ratingSSL)
colnames(rating.all) <- c("Year","Month","Day","Discharge",
                          "Suspended.Load","asNA","ratingSSL_all")


# 計算損失函數
rating.all_vs_observation <- rating.all[(is.na(rating.all$asNA)),] #保留 na
# rmse (root mean square error)

rmse.rating <- rmse(rating.all_vs_observation$Suspended.Load, rating.all_vs_observation$ratingSSL_all) #原本移除的觀測值與率定值

# nmse (normalize mean squared error)

nmse.rating <- nmse(rating.all_vs_observation$Suspended.Load, rating.all_vs_observation$ratingSSL_all) #原本移除的觀測值與率定值
# mse (mean absolute error)

mse.rating <- mse(rating.all_vs_observation$Suspended.Load, rating.all_vs_observation$ratingSSL_all) #原本移除的觀測值與率定值

# mape (mean absolute persentage error)

mape.rating <- mape(rating.all_vs_observation$Suspended.Load, rating.all_vs_observation$ratingSSL_all) #原本移除的觀測值與率定值

rating.table <- cbind(rmse.rating,
                      nmse.rating,
                      mse.rating,
                      mape.rating) #沒分組的率定曲線vs觀測資料
# 率定曲線(不分組)資料出圖
setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站率定曲線.png"),width = 1000, height = 700, units = "px", pointsize = 12)
ggplot(data=rating.all)+
  geom_point(aes(x=Discharge,y=Suspended.Load/10000),size=4) +
  geom_line(aes(x=Discharge,y=ratingSSL_all/10000),size=2) +
  scale_x_continuous(guide = "prism_minor", #x軸副刻度
                     limits = c(0, max(rating.all$Discharge)),
                     minor_breaks = seq(0, max(rating.all$Discharge), 100))+
  scale_y_continuous(guide = "prism_minor", #y軸副刻度
                     limits = c(0, max(rating.all$Suspended.Load)/10000),
                     minor_breaks = seq(0, max(rating.all$Suspended.Load)/10000, 20))+
  labs(x="日流量Q(立方米/秒)",y=TeX("$日懸浮載輸砂量L(10^4公頓/日)$")) + # 座標軸名稱
  theme_bw()+ #白底
  theme(legend.position = c(0.2,0.9))+ #圖例位置座標
  theme(legend.text = element_text("1","2"))+
  theme(legend.title=element_blank())+ #隱藏圖例標題
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(prism.ticks.length=unit(.7,"lines"))+
  theme(axis.ticks.length=unit(1.2,"lines"))+
  theme(text=element_text(size=40,color="black"))+  # 字體大小
  theme(axis.text = element_text(colour = "black"))
dev.off()

#
# =============== Determine margianl distribution (parametric) ==============
#

print("邊際分布之參數估計")
candidate <- c("norm","lnorm","gumbel","weibull","gamma")
dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b)) # a:location, b:scale
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))
#
MD.anal <- as.data.frame(MD.rmNA[,4:5])
# 邊際分布參數表格
margin.par <- matrix(nrow=length(candidate),ncol=2*2)
rownames(margin.par) <- c(candidate)
colnames(margin.par) <- c("par1","par2","par1","par2") 
# margin.ks 放置p-value
margin.ks <- matrix(nrow=length(candidate)+1,ncol=dim(MD.anal)[2]) # +1 是為了要放最佳分布的欄位
rownames(margin.ks) <- c(candidate,"good dist")
colnames(margin.ks) <- c(colnames(MD.anal))
# margin.aic  放置AIC value
margin.aic <- matrix(nrow=length(candidate)+1,ncol=dim(MD.anal)[2]) # +1 是為了要放最小AIC的欄位
rownames(margin.aic) <- c(candidate,"good dist")
colnames(margin.aic) <- c(colnames(MD.anal))
# 邊際分布編號
margin.num <- matrix(ncol=2)

#
for(i in 1:dim(MD.anal)[2]){ # Q ,QS
  var <- MD.anal[,i]
  print(paste0("第",i,"個變數：",colnames(MD.anal)[i])) #顯示第幾個及變數名稱
  # By Maximun Likelihood Estimate Method
  for(dist in c(1:length(candidate))){
    # -------------------------- parameter estimate ----------------------------------
    print(candidate[dist])
    dist.char <- c(candidate[dist],
                   paste0("d", candidate[dist]), 
                   paste0("p", candidate[dist]),
                   paste0("q", candidate[dist]))
    #md <- fitdistr(x, distribution, start = list(parameter1 = 1, parameter2 = 1))
    if(candidate[dist] == "norm"){
      md <- fitdist(var, dist = dist.char[1])}
    if(candidate[dist] == "lnorm"){
      md <- fitdist(var, dist = dist.char[1], start = list(meanlog=1, sdlog=1))}
    if(candidate[dist] == "gumbel"){
      fitgumbel <- eevd(var,method = "mle")  # 先計算初始值
      md <- fitdist(var, dist = dist.char[1], 
                    start = list(a=as.numeric(fitgumbel$parameters[1]),
                                 b=as.numeric(fitgumbel$parameters[2])))}
    if(candidate[dist] == "weibull"){
      md <- fitdist(var, dist = dist.char[1])}
    if(candidate[dist] == "gamma"){
      md <- fitdist(var, dist = dist.char[1],lower=0, upper=Inf)}
    par1 <- md$estimate[1] #fitting參數1
    par2 <- md$estimate[2] #參數2
    print(c(par1, par2))
    # parameter 擺放設定
    if(i==1){
      margin.par[dist,1] <- as.numeric(par1) 
      margin.par[dist,2] <- as.numeric(par2)}
    if(i==2){
      margin.par[dist,3] <- as.numeric(par1) 
      margin.par[dist,4] <- as.numeric(par2)}
    #
    # ------------------------------- K-S test ------------------------
    #
    print("KS test")
    result <- ks.test(var, dist.char[3], par1, par2)
    print(paste0(candidate[dist], "的KS檢定P-value: ", result$p.value))
    #
    # 將P-value整理成表格
    margin.ks[dist,i] <- result$p.value
    #margin.ks.D[dist,i] <- result$p.value
    #
    # --------------------------------- AIC -----------------------------
    print("AIC")
    if(candidate[dist] != "gumbel"){
      aic <- gofstat(md)
      margin.aic[dist,i] <- aic$aic # 將AIC-value整理成表格
      print(paste0(candidate[dist], "AIC值: ", aic$aic))}
    if(candidate[dist] == "gumbel"){
      aic <- abic.gumbel(var,par1,par2)
      margin.aic[dist,i] <- aic$AIC # 將AIC-value整理成表格
      print(paste0(candidate[dist], "AIC值: ", aic$AIC))}
  }
  # 每個延時P-value排序，變成數值再排序
  ks.choice <- rank(as.numeric(margin.ks[(1:dist),i])) 
  margin.ks[length(candidate)+1,i] <- candidate[which.max(ks.choice)] #最大的P-value對應的機率分布
  # 每個延時AIC-value排序，變成數值再排序
  aic.choice <- rank(as.numeric(margin.aic[(1:dist),i])) 
  margin.aic[length(candidate)+1,i] <- candidate[which.min(aic.choice)] #最小的AIC-value對應的機率分布
  # 邊際分布代號(norm=1, lnorm=2, gumbel=3, weibull=4, gamma=5)
  margin.num[,i] <- c(which.min(aic.choice)) 
}

margin.dist <-c(margin.aic[length(candidate)+1,1],margin.aic[length(candidate)+1,2]) # 請輸入邊際分布：

# 畫最佳邊際分布分析圖(PDFs.CDFs,Q-Q plot,P-P plot)
# 注意：因尺度不同，需要手動調整，各自出圖
# for(i in 1:dim(MD.anal)[2]){ # Q ,QS
#   i <- 1
#   var <- MD.anal[,i]
#   print(paste0("第",i,"個變數：",colnames(MD.anal)[i])) #顯示第幾個及變數名稱
#   if(margin.dist[i] != "gumbel"){
#     md.plot <- fitdist(var, dist = margin.dist[i])}
# 
#   if(margin.dist[i] == "gumbel"){
#     fitgumbel <- eevd(var,method = "mle")  # 先計算初始值
#     md.plot <- fitdist(var, dist = margin.dist[i],
#                        start = list(a=as.numeric(fitgumbel$parameters[1]),
#                                     b=as.numeric(fitgumbel$parameters[2])))}
#   
#   setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
#   png(paste0(year[1],"到",year[y],"年",station_ch,"測站",colnames(MD.anal)[i],"邊際分布分析.png"),width = 800, height = 600, units = "px", pointsize = 12)
#   # ---- PDF、CDF、Q-Q plot、P-P plot ----
#   par(mfcol=c(2,2))
#   options(scipen = 200)
#   denscomp(md.plot,xlab=TeX("$日流量(m^3/s)$"),main="PDF",breaks=200,xlim=c(0,300),cex=2)
#   cdfcomp(md.plot,xlab=TeX("$日流量(m^3/s)$"),main="CDF",xlim=c(0,300),horizontals=F,lines01=T,xlegend=NULL)
#   qqcomp(md.plot,xlab=TeX("$理論分布之分位數(m^3/s)$"),ylab=TeX("$觀測值之分位數(m^3/s)$"),
#          xlim=c(0,300), ylim=c(0,300),addlegend = F)
#   ppcomp(md.plot,xlab=TeX("$理論值累積機率$"),ylab=TeX("$觀測值累積機率$"),
#          addlegend = F)
#   
#   # par(mfcol=c(2,2))
#   # options(scipen = 200)
#   # denscomp(md.plot,xlab=TeX("$日輸砂量(ton/d)$"),main="PDF",breaks=2000,xlim=c(0,100000),cex=2)
#   # cdfcomp(md.plot,xlab=TeX("$日輸砂量(ton/d)$"),main="CDF",xlim=c(0,100000),horizontals=F,lines01=T,xlegend=NULL)
#   # qqcomp(md.plot,xlab=TeX("$理論分布之分位數(ton/d)$"),ylab=TeX("$觀測值之分位數(ton/d)$"),
#   #        xlim=c(0,100000),ylim=c(0,100000),addlegend = F)
#   # ppcomp(md.plot,xlab=TeX("$理論值累積機率$"),ylab=TeX("$觀測值累積機率$"),
#   #        addlegend = F)
#   
#   # --- CDF plot ----
#   #cdfcomp(md.plot,xlim=c(0,800),xlab=TeX("$Discharge(m^3/s)$"),
#   #        legendtext = NULL,datapch=16,horizontals=F,verticals=F,cex.lab=1.4,cex.axis=1.5)
#   #legend("bottomright", legend = c("obsreved", "fitted"),
#   #       pch = c(19, NA), lty = c(NA, 1),col=c(1,2),cex=1.5)
#   
#   # options(scipen = 200)
#   # cdfcomp(md.plot,xlim=c(0,1000000),xlab=TeX("Suspended sediment $(ton/d)$"),
#   #         legendtext = NULL,datapch=16,horizontals=F,verticals=F,cex.lab=1.4,cex.axis=1.5)
#   # legend("bottomright", legend = c("obsreved", "fitted"),
#   #        pch = c(19, NA), lty = c(NA, 1),col=c(1,2),cex=1.5)
#   
#   dev.off()
# }
f <- get(paste0("p",margin.dist[1]))(MD.anal$Discharge, margin.par[margin.num[1],1], margin.par[margin.num[1],2])
ggplot(data.frame(MD.anal$Discharge,f),aes(x=MD.anal$Discharge,y=f))+
  geom_line()+
  xlim(0,500)+
  xlab("流量 Q (cms)")+
  ylab("累積機率")+
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=20))  # 字體大小

f <- get(paste0("p",margin.dist[2]))(MD.anal$Suspended.Load, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
options(scipen = 200)
ggplot(data.frame(MD.anal$Suspended.Load,f),aes(x=MD.anal$Suspended.Load,y=f))+
  geom_line()+
  xlim(0,100000)+
  xlab("輸砂量 L (ton/d)")+
  ylab("累積機率")+
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=15))  # 字體大小

# 儲存 margins 相關資料
file <- paste("F:/R_output/",station,"/vinecopula/",
              year[1],"到", year[y],station_ch," margins parameter.csv", sep="") #存檔路徑
write.csv(margin.par,file)
file <- paste("F:/R_output/",station,"/vinecopula/",
              year[1],"到", year[y],station_ch," margins ks_test.csv", sep="") #存檔路徑
write.csv(margin.ks,file)
file <- paste("F:/R_output/",station,"/vinecopula/",
              year[1],"到", year[y],station_ch," margins AIC.csv", sep="") #存檔路徑
write.csv(margin.aic,file)
#
# ================== Determine copula function ==========================
#
#
# 建立copula參數估計表格(參數估計方法：itau, irho, mpl, ml)
fitcopula.par <- matrix(nrow=1,ncol=20)
colnames(fitcopula.par) <- c("gumbel.itau","gumbel.irho","gumbel.mpl","gumbel.ml",
                             "frank.itau","frank.irho","frank.mpl","frank.ml",
                             "clayton.itau","clayton.irho","clayton.mpl","clayton.ml",
                             "amh.itau","amh.irho","amh.mpl","amh.ml",
                             "joe.itau","joe.irho","joe.mpl","joe.ml")
# 新增最後一行放選擇的參數值(預設值為0)
fitcopula.par <- data.frame(fitcopula.par, copula.parameter=0)

# 建立p-value表格(參數估計方法：itau, irho, mpl, ml)
gof.pvalue <- matrix(nrow=1,ncol=20)
colnames(gof.pvalue) <- c("gumbel.itau","gumbel.irho","gumbel.mpl","gumbel.ml",
                          "frank.itau","frank.irho","frank.mpl","frank.ml",
                          "clayton.itau","clayton.irho","clayton.mpl","clayton.ml",
                          "amh.itau","amh.irho","amh.mpl","amh.ml",
                          "joe.itau","joe.irho","joe.mpl","joe.ml")
# 新增最後一行放選擇的聯結函數(預設值為0)
gof.pvalue <- data.frame(gof.pvalue, resultcopula=0)
# ----- fitCopula (Inference For Margin estimation) ------ 
candidate.copula <- as.data.frame(cbind(number=c(3,4,5),copula=c("clayton","gumbel","frank")))
#Q與Qs的邊際CDF
F_Q <- get(paste0("p",margin.dist[1]))(MD.anal$Discharge, margin.par[margin.num[1],1], margin.par[margin.num[1],2])
F_L <- get(paste0("p",margin.dist[2]))(MD.anal$Suspended.Load, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
margins.cdf <- cbind(F_Q,F_L)
# select the bivariate copula family and estimate the parameter(s)
cop <- BiCopSelect(F_Q, F_L, familyset = c(3,4,5), 
                    indeptest = FALSE, level = 0.05, method="mle")
                     # (3：clayton) , (4：gumbel), (5：frank)
family.num <- cop$family  # AIC值最小之聯結函數代碼
cop.gof <- BiCopGofTest(F_Q, F_L, family = family.num,par=cop$par, method = "kendall")

if(cop.gof$p.value.KS>0.05){
  if(family.num==3){
    print("Clayton copula function 通過 K-S test")
  }
  if(family.num==4){
    print("Gumbel copula function 通過 K-S test")
  }
  if(family.num==5){
    print("Frank copula function 通過 K-S test")
  }
  copula.func <- get(paste0(subset(candidate.copula,number==family.num)[,2],"Copula"))(cop$par,dim=2)
  pvalue <- cop.gof$p.value.KS
  Dvalue <- cop.gof$statistic.KS
  copula.gof.ks <- cbind(pvalue,Dvalue)
  file <- paste("F:/R_output/",station,"/vinecopula/",
                year[1],"到", year[y],station_ch," copula_gof(ks).csv", sep="") #存檔路徑
  write.csv(copula.gof.ks,file)
}
if(cop.gof$p.value.CvM>0.05){
  if(family.num==3){
    print("clayton copula function 通過 Cramer-Von Mises test")
  }
  if(family.num==4){
    print("gumbel copula function 通過 Cramer-Von Mises test")
  }
  if(family.num==5){
    print("frank copula function 通過 Cramer-Von Mises test")
  }
  copula.func <- get(paste0(subset(candidate.copula,number==family.num)[,2],"Copula"))(cop$par,dim=2)
  pvalue <- cop.gof$p.value.CvM
  snvalue <- cop.gof$statistic.CvM
  copula.gof.cvm <- cbind(pvalue,snvalue)
  file <- paste("F:/R_output/",station,"/vinecopula/",
                year[1],"到", year[y],station_ch," copula_gof(cvm).csv", sep="") #存檔路徑
  write.csv(copula.gof.cvm,file)
}
copula.func

# 儲存copula function相關資料
file <- paste("F:/R_output/",station,"/vinecopula/",
              year[1],"到", year[y],station_ch," copula parameter.csv", sep="") #存檔路徑
write.csv(fitcopula.par,file)
file <- paste("F:/R_output/",station,"/vinecopula/",
              year[1],"到", year[y],station_ch," copula gof.csv", sep="") #存檔路徑
write.csv(gof.pvalue,file)

# ============= Build Mvdc (parametric margins + parametric copula) ==============
# ---- 1. variable----
# L：random variable of Suspended sediment load
# Q：random variable od streamflow discharge

# ---- 2.probability distribution -----
# f_L(l)： PDF of L (Suspended.Load)
# f_Q(q)： PDF of Q (Discharge)
# F_L(l)： CDF of L (Suspended.Load)
# F_Q(q)： CDF of Q (Discharge)

# ---- 3. joint probability distribution ----
# joint CDF of L and Q：
# F_{L,Q}(l,q) = C(F_L(l),F_Q(q)), C：copula
# joint PDF of L and Q：
# f_{L,Q}(l,q) = c(F_L(l),F_Q(q))*f_L(l)*f_Q(q) , c：copula density

# ---- 4. conditional probability distribution ----
# conditional PDF of L given the observation discharge q_0：
# f_{L|q_0}(l) = f_L(l)*c(F_L(l),F_Q(q))
# conditional CDF of L given the observation discharge q_0：
# F_{L|q_0}(l) = C_{L|q_0}(F_L(l)|F_Q(q_0))
# ============================================================

# ---- 2.probability distribution -----
# f_Q(q)： PDF of Q (Discharge)
# f_L(l)： PDF of L (Suspended.Load)
# F_Q(q)： CDF of Q (Discharge)
# F_L(l)： CDF of L (Suspended.Load)
f_Q <- get(paste0("d",margin.dist[1]))(MD.anal$Discharge, margin.par[margin.num[1],1], margin.par[margin.num[1],2])
f_L <- get(paste0("d",margin.dist[2]))(MD.anal$Suspended.Load, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
F_Q <- get(paste0("p",margin.dist[1]))(MD.anal$Discharge, margin.par[margin.num[1],1], margin.par[margin.num[1],2])
F_L <- get(paste0("p",margin.dist[2]))(MD.anal$Suspended.Load, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
# f_Q
setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站流量之機率密度函數.png"),width = 1250, height = 700, units = "px", pointsize = 12)
ggplot()+
  geom_line(aes(x=MD.anal$Discharge,y=f_Q))+
  labs(x="日流量Q (cms)",y="PDF函數值") + # 座標軸名稱
  ggtitle(paste0(station_ch,"測站流量之機率密度函數")) +
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=40))  # 字體大小
dev.off()
# f_L
setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站輸砂量之機率密度函數.png"),width = 1250, height = 700, units = "px", pointsize = 12)
ggplot()+
  geom_line(aes(x=MD.anal$Suspended.Load,y=f_L))+
  labs(x="日懸浮載輸砂量L(ton)",y="PDF函數值") + # 座標軸名稱
  xlim(0,10000) +
  ggtitle(paste0(station_ch,"測站輸砂量之機率密度函數")) +
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=40))  # 字體大小
dev.off()
# F_Q
setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站流量之累積分布函數.png"),width = 1250, height = 700, units = "px", pointsize = 12)
ggplot()+
  geom_line(aes(x=MD.anal$Discharge,y=F_Q))+
  labs(x="日流量Q(cms)",y="累積機率") + # 座標軸名稱
  ggtitle(paste0(station_ch,"測站流量之累積分布函數")) +
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=40))  # 字體大小
dev.off()
# F_L
setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站輸砂量之累積分布函數.png"),width = 1250, height = 700, units = "px", pointsize = 12)
ggplot()+
  geom_line(aes(x=MD.anal$Suspended.Load,y=F_L))+
  labs(x="日懸浮載輸砂量L(ton)",y="累積機率") + # 座標軸名稱
  #xlim(0,2000) +
  ggtitle(paste0(station_ch,"測站輸砂量之累積分布函數")) +
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=40))  # 字體大小
dev.off()
# ---- 3. joint probability distribution ----
# joint CDF of L and Q：
# F_{L,Q}(l,q) = C(F_L(l),F_Q(q)), C：copula
# joint PDF of L and Q：
# f_{L,Q}(l,q) = c(F_L(l),F_Q(q))*f_L(l)*f_Q(q) , c：copula density
#
# plot copula function
simcopula <- pCopula(margins.cdf,copula.func)
df <- data.frame(FQ=margins.cdf[,1],FQs=margins.cdf[,2],simcopula)
setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站 copula CDF 等高線圖.png"),width = 500, height = 350, units = "px", pointsize = 12)
contour(copula.func, pCopula, xlim = c(0,1), ylim=c(0,1),xlab = TeX('$\\F_Q(q)$'), 
        col=c(2),ylab=TeX('$\\F_L(l)$'),labcex = 1.3,cex.lab=1.15,cex.axis=1.2)
points(df$FQ,df$FQs,pch=16)
dev.off()

ggplot()+
  geom_point(data=as.data.frame(margins.cdf),aes(x=F_Q,y=F_L))+
  geom_contour_filled(data=df,aes(x=FQ,y=FQs,z=simcopula),breaks=simcopula)
dev.off()

#建立Mvdc
mymvdc <- mvdc(copula.func,margin.dist,
               paramMargins =list(list(margin.par[margin.num[1],1], margin.par[margin.num[1],2]),
                                  list(margin.par[margin.num[2],3], margin.par[margin.num[2],4])))
set.seed(100)
#n <- 5000
# v <- rMvdc(n, mymvdc) # 隨機生成5000點的模型
v <- as.matrix(MD.anal) # 率定資料

# Compute the density
mvd.pdf <- dMvdc(v, mymvdc)

# Compute the CDF
mvd.cdf <- pMvdc(v, mymvdc)

# plot joint PDF and joint CDF

setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站 mvd 散佈圖.png"),width = 1250, height = 700, units = "px", pointsize = 12)
par(mfrow = c(1, 2))
scatterplot3d(v[,1],v[,2], mvd.pdf, color="red", main= paste0("joint PDF 散佈圖"), xlab = "Q(cms)", ylab="L(ton)", zlab="pMvdc",pch=".")
scatterplot3d(v[,1],v[,2], mvd.cdf, color="red", main=paste0("joint CDF 散佈圖"), xlab = "Q(cms)", ylab="L(ton)", zlab="pMvdc",pch=".")
dev.off()

setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站joint PDF.png"),width = 500, height = 350, units = "px", pointsize = 12)
contour(mymvdc, dMvdc, xlim = c(0,50), ylim=c(0, 5000), xlab =  TeX("$流量Q \\,(m^3/s)$"), ylab=TeX("$輸砂量L\\,(ton/d)$"),col=c(2),labcex = 1.2,cex.lab=1.2,cex.axis=1.2)
points(MD.anal$Discharge,MD.anal$Suspended.Load,pch=16)
minor.tick(nx=4, ny=4, tick.ratio=.5)
dev.off()

setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站joint CDF.png"),width = 500, height = 350, units = "px", pointsize = 12)
contour(mymvdc, pMvdc, xlim = c(0,300), ylim=c(0, 300000),xlab =  TeX("$流量Q \\,(m^3/s)$"), ylab=TeX("$輸砂量L\\,(ton/d)$"),col=c(2),labcex = 1.2,cex.lab=1.2,cex.axis=1.2)
points(MD.anal$Discharge,MD.anal$Suspended.Load,pch=16)
minor.tick(nx=4, ny=4, tick.ratio=.5)
dev.off()

# ---- 4. conditional probability distribution ----

## conditional pdf of L given q0
## conditional PDF of L given the observation discharge q_0：
## f_{L|q_0}(l) = f_L(l)*c(F_L(l),F_Q(q))
# copula cdf and pdf function

# ---- 老師推導gumbelcopula pdf function ----
# gumbelcopula.pdf <- function(u,v,theta){
#   exp(-((-log(u))^theta+(-log(v))^theta)^(1/theta)) *
#     (((-log(u))*(-log(v)))^(theta-1))/(u*v) *
#     ((-log(u))^theta+(-log(v))^theta)^(2/theta-2) *
#     ((theta-1)*((-log(u))^theta+(-log(v))^theta)^(-1/theta)+1)
# }
# ------- 程式從CDF 對u,v微分成PDF -----------
gumbelcopula.cdf <- expression(exp(-((-log(u))^theta+(-log(v))^theta)^(1/theta)))
D(D(gumbelcopula.cdf,"v"),"u")
gumbelcopula.pdf <- function(u,v,theta){
  exp(-((-log(u))^theta + (-log(v))^theta)^(1/theta)) * (((-log(u))^theta +
  (-log(v))^theta)^((1/theta) - 1) * ((1/theta) * ((-log(u))^(theta -1) *
  (theta * (1/u))))) * (((-log(u))^theta + (-log(v))^theta)^((1/theta) -1) *
  ((1/theta) * ((-log(v))^(theta - 1) * (theta * (1/v))))) -
  exp(-((-log(u))^theta + (-log(v))^theta)^(1/theta)) * (((-log(u))^theta +
  (-log(v))^theta)^(((1/theta) - 1) - 1) * (((1/theta) -1) *
  ((-log(u))^(theta - 1) * (theta * (1/u)))) * ((1/theta) *
  ((-log(v))^(theta - 1) * (theta * (1/v)))))
}
gumbel.pdf <- expression( exp(-((-log(u))^theta + (-log(v))^theta)^(1/theta)) * (((-log(u))^theta +
              (-log(v))^theta)^((1/theta) - 1) * ((1/theta) * ((-log(u))^(theta -1) *
              (theta * (1/u))))) * (((-log(u))^theta + (-log(v))^theta)^((1/theta) -1) *
              ((1/theta) * ((-log(v))^(theta - 1) * (theta * (1/v))))) -
               exp(-((-log(u))^theta + (-log(v))^theta)^(1/theta)) * (((-log(u))^theta +
               (-log(v))^theta)^(((1/theta) - 1) - 1) * (((1/theta) -1) *
               ((-log(u))^(theta - 1) * (theta * (1/u)))) * ((1/theta) *
                ((-log(v))^(theta - 1) * (theta * (1/v))))))
D(gumbel.pdf,"u")

claytoncopula.cdf <- expression((u^(-theta)+v^(-theta)-1)^(-1/theta))
D(D(claytoncopula.cdf,"v"),"u")
claytoncopula.pdf <- function(u,v,theta){
  (u^(-theta) + v^(-theta) - 1)^(((-1/theta) - 1) - 1) * (((-1/theta) - 1) *
  (u^((-theta) - 1) * (-theta))) * ((-1/theta) * (v^((-theta) - 1) * (-theta)))
}
frankcopula.cdf <- expression((-1/theta) * log(1+(exp(-theta*u)-1)*(exp(-theta*v)-1)/(exp(-theta)-1)))
D(D(frankcopula.cdf,"v"),"u")
frankcopula.pdf <- function(u,v,theta){
  (-1/theta) * (exp(-theta * u) * theta * (exp(-theta * v) * theta)/
  (exp(-theta) - 1)/(1 + (exp(-theta * u) - 1) * (exp(-theta * v) - 1)/
  (exp(-theta) - 1)) - (exp(-theta * u) - 1) * (exp(-theta * v) * theta)/
  (exp(-theta) - 1) * (exp(-theta * u) * theta * (exp(-theta * v) - 1)/
  (exp(-theta) - 1))/(1 + (exp(-theta * u) - 1) * (exp(-theta * v) - 1)/
  (exp(-theta) - 1))^2)
}

amhcopula.cdf <- expression((u*v)/(1-theta*(1-u)*(1-v)))
D(D(amhcopula.cdf,"v"),"u")
amhcopula.pdf <- function(u,v,theta){
  1/(1 - theta * (1 - u) * (1 - v)) - u * (theta * (1 - v))/
    (1 - theta * (1 - u) * (1 - v))^2 - ((v * (theta * (1 - u)) -
    (u * v) * theta)/(1 - theta * (1 - u) * (1 - v))^2 -
    (u * v) * (theta * (1 - u)) *
    (2 * (theta * (1 - v) * (1 - theta * (1 - u) * (1 - v))))/
    ((1 - theta * (1 - u) * (1 - v))^2)^2)
}

joecopula.cdf <- expression(1-((1-u)^theta+(1-v)^theta-(1-u)^theta*(1-v)^theta)^(1/theta))
D(D(joecopula.cdf,"v"),"u")
joecopula.pdf <- function(u,v,theta){
  ((1 - u)^theta + (1 - v)^theta - (1 - u)^theta * (1 - v)^theta)^
  ((1/theta) -1) * ((1/theta) * ((1 - u)^(theta - 1) * theta *
  ((1 - v)^(theta - 1) * theta))) - ((1 - u)^theta + (1 - v)^theta -
  (1 - u)^theta * (1 - v)^theta)^(((1/theta) - 1) - 1) * (((1/theta) - 1) *
  ((1 - u)^(theta - 1) * theta - (1 - u)^(theta - 1) * theta *
  (1 - v)^theta)) * ((1/theta) * ((1 - v)^(theta - 1) * theta -
  (1 - u)^theta * ((1 - v)^(theta - 1) * theta)))
}

# ---- 單一流量組 ----
givenQ <- c(5) # 設定條件流量大小：
small.Qs <- seq(from=1,to=99999,by=1)
middle.Qs <- append(small.Qs,seq(from=100000,to=999999,by=10))
big.Qs <- append(middle.Qs,seq(from=1000000,to=9999999,by=100))
all.Qs <- append(big.Qs,seq(from=10000000,to=100000000,by=1000))

## caculate conditional PDF of L given the q0
F_Q <- get(paste0("p",margin.dist[1]))(givenQ, margin.par[margin.num[1],1], margin.par[margin.num[1],2])
F_L <- get(paste0("p",margin.dist[2]))(all.Qs, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
f_L <- get(paste0("d",margin.dist[2]))(all.Qs, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
con.pdf <-get(paste0(subset(candidate.copula,number==family.num)[,2],"copula.pdf"))(F_Q,F_L,cop$par)*f_L

# plot conditional PDF of L given the q0

setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站Q=",givenQ,"cms下conditional PDF.png"),width = 1250, height = 700, units = "px", pointsize = 12)
Con.pdf <- ggplot()+
  geom_line(aes(x=all.Qs,y=con.pdf),size=1.5)+
  labs(x="輸砂量Qs (公噸)",y="PDF函數值") + # 座標軸名稱
  xlim(0,10000) +
  ggtitle(paste0(station_ch,"測站Q=",givenQ,"cms下條件機率密度函數")) +
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=40))  # 字體大小
print(Con.pdf)
dev.off()

## caculate conditional CDF of L given the q0
F_Q <- get(paste0("p",margin.dist[1]))(givenQ, margin.par[margin.num[1],1], margin.par[margin.num[1],2])
F_L <- get(paste0("p",margin.dist[2]))(all.Qs, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
con.cdf <- cCopula(cbind(F_Q,F_L), copula = copula.func,indices = dim(copula.func), inverse = F)

# plot conditional CDF of L given the q0
setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站Q=",givenQ,"cms下conditional CDF.png"),width = 1250, height = 700, units = "px", pointsize = 12)
ggplot()+
  geom_line(aes(x=all.Qs,y=con.cdf),size=1.5) +
  labs(x="輸砂量Qs (公噸)",y="CDF累積機率值") + # 座標軸名稱
  xlim(0,200000)+
  ggtitle(paste0(station_ch,"測站Q=",givenQ,"cms下條件累積分布函數")) +
  scale_color_discrete(name="圖例")+
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=40))  # 字體大小
print(Con.cdf)
dev.off()

# ------ 比較多個流量組 ------
givenQ <- c(20,30,50) # 設定條件流量大小：

small.Qs <- seq(from=1,to=99999,by=1)
middle.Qs <- append(small.Qs,seq(from=100000,to=999999,by=10))
big.Qs <- append(middle.Qs,seq(from=1000000,to=9999999,by=100))
all.Qs <- append(big.Qs,seq(from=10000000,to=100000000,by=1000))
con.pdf.table <- c()
con.cdf.table <- c()
for ( i in 1:length(givenQ)){
  ## caculate conditional PDF of L given the q0
  F_Q <- get(paste0("p",margin.dist[1]))(givenQ[i], margin.par[margin.num[1],1], margin.par[margin.num[1],2])
  F_L <- get(paste0("p",margin.dist[2]))(all.Qs, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  f_L <- get(paste0("d",margin.dist[2]))(all.Qs, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  con.pdf <-get(paste0(subset(candidate.copula,number==family.num)[,2],"copula.pdf"))(F_Q,F_L,cop$par)*f_L
  each.pdf <- data.frame(paste0(givenQ[i],"cms"),all.Qs,con.pdf)
  colnames(each.pdf) <- c("Discharge","SSL","PDF")
  con.pdf.table <- rbind(con.pdf.table,each.pdf)

  ## caculate conditional CDF of L given the q0
  F_Q <- get(paste0("p",margin.dist[1]))(givenQ[i], margin.par[margin.num[1],1], margin.par[margin.num[1],2])
  F_L <- get(paste0("p",margin.dist[2]))(all.Qs, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  con.cdf <- cCopula(cbind(F_Q,F_L), copula = copula.func,indices = dim(copula.func), inverse = F)
  each.cdf <- data.frame(paste0(givenQ[i],"cms"),all.Qs,con.cdf)
  colnames(each.cdf) <- c("Discharge","SSL","CDF")
  con.cdf.table <- rbind(con.cdf.table,each.cdf)
}

# 描述CDF之用
which(subset(con.cdf.table,Discharge=="50cms")$CDF>0.75)
(subset(con.cdf.table,Discharge=="50cms")[5000,]$CDF)-subset(con.cdf.table,Discharge=="50cms")[1000,]$CDF

# plot conditional PDF of L given the q0
setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站不同流量情況下conditional PDF.png"),width = 1000, height = 600, units = "px", pointsize = 12)
ggplot(data=con.pdf.table,aes(x=SSL,y=PDF,group=Discharge,color=Discharge))+
  geom_line(aes(size=Discharge))+
  stat_peaks(span = NULL,
             geom = "text_repel",
             mapping = aes(label = paste(..y.label.., ..x.label..)),
             x.label.fmt = "at %.0d ton",
             y.label.fmt = "max PDF value = %.5f",
             segment.colour = "black",
             arrow = grid::arrow(angle=10,length = unit(1, "cm"),ends="last",type="closed"),
             nudge_x = 5,size = 8, vjust =.25,hjust=-.8)+
  scale_linetype_manual(values=c("twodash", "dotted","solid")) +
  scale_size_manual(values=c(1.3,1.3,1.3))+
  labs(x=TeX("$日懸浮載輸砂量 \\,(ton/d)$"),y="PDF") + # 座標軸名稱
  xlim(0,10000)+
  scale_x_continuous(guide = "prism_minor", #x軸副刻度
                     limits = c(0, 10000),
                     minor_breaks = seq(0, 10000, 500))+
  scale_y_continuous(guide = "prism_minor", #x軸副刻度
                     limits = c(0, max(con.pdf.table$PDF)),
                     minor_breaks = seq(0, max(con.pdf.table$PDF), .00005))+
  theme_bw()+ #白底
  theme(legend.position = c(.8,.8))+ #圖例位置座標
  theme(legend.title=element_blank())+ #隱藏圖例標題
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(prism.ticks.length=unit(.7,"lines"))+
  theme(axis.ticks.length=unit(1.2,"lines"))+
  theme(text=element_text(size=40,color="black"))+  # 字體大小
  theme(axis.text = element_text(colour = "black"))
dev.off()

# plot conditional CDF of L given the q0
setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站不同流量情況下conditional CDF.png"),width = 950, height = 600, units = "px", pointsize = 12)
ggplot(data=con.cdf.table,aes(x=SSL,y=CDF,color=Discharge))+
  geom_line(aes(size=Discharge))+
  scale_linetype_manual(values=c("twodash", "dotted","solid")) +
  scale_size_manual(values=c(1.3,1.3,1.3))+
  labs(x=TeX("$日懸浮載輸砂量 \\,(ton/d)$"),y="CDF") + # 座標軸名稱
  xlim(0,50000)+
  scale_x_continuous(guide = "prism_minor", #x軸副刻度
                     limits = c(0, 50000),
                     minor_breaks = seq(0, 50000, 1000))+
  scale_y_continuous(guide = "prism_minor", #x軸副刻度
                     limits = c(0, max(con.cdf.table$CDF)),
                     minor_breaks = seq(0, max(con.cdf.table$CDF), .05))+
  theme_bw()+ #白底
  theme(legend.position = c(.85,.2))+ #圖例位置座標
  theme(legend.title=element_blank())+ #隱藏圖例標題
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(prism.ticks.length=unit(.7,"lines"))+
  theme(axis.ticks.length=unit(1.2,"lines"))+
  theme(text=element_text(size=40,color="black"))+  # 字體大小
  theme(axis.text = element_text(colour = "black"))
dev.off()

#  =========== 輸砂量推估 ============
# 1. mode(中位數)
# 2. (Peng et al. 2020)
# 3. (Bezak et al. 2017)
# 4. Suspended load rating curve(非線性最小平方法)
# 5. (Di Lascio et al. 2015)原始尺度

# ====================================
print("開始率定")
givenQ.table <- MD.rmNA$Discharge
MD.rmNA.rmSSL <- MD.rmNA[,-5]
rat.SSL0 <- data.frame(MD.rmNA,group="觀測資料")
estimate.SSL <- c()

# 1. Suspended load rating curve 率定曲線
# by rating curve coefficient a, b
Suspended.Load <- c()
for (q in 1:length(givenQ.table)){
  each.m1 <- a * (givenQ.table[q])^b
  Suspended.Load <- rbind(Suspended.Load,each.m1)
}
rat.SSL1 <- data.frame(MD.rmNA.rmSSL,Suspended.Load,group="率定曲線")
colnames(rat.SSL1) <- c("Year","Month","Day","Discharge","Suspended.Load","group")
print("率定曲線率定完畢")

# 2.方法1 (Peng et al. 2020) 從CDF隨機取一個點
Suspended.Load <- c()
for (q in 1:length(givenQ.table)){
  r2 <- runif(1,min=0.001,max=0.995) #隨機產生1個均勻分布
  F_Q <- get(paste0("p",margin.dist[1]))(givenQ.table[q], margin.par[margin.num[1],1], margin.par[margin.num[1],2])
  con.cdf2 <- cCopula(cbind(F_Q,r2), copula = copula.func,indices = 2, inverse = T)
  if(con.cdf2>0.995){
    con.cdf2 <- 0.995
  }
  each.m2 <- get(paste0("q",margin.dist[2]))(con.cdf2,margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  Suspended.Load <- rbind(Suspended.Load,each.m2)
}  
rat.SSL2 <- data.frame(MD.rmNA.rmSSL,Suspended.Load,group="方法 1")
colnames(rat.SSL2) <- c("Year","Month","Day","Discharge","Suspended.Load","group")
print("方法1率定完畢")

# 3.方法2 (Bezak et al. 2017) 從CDF隨機取10000個點，再取中位數
Suspended.Load <- c()
for (q in 1:length(givenQ.table)){
  r3 <- runif(10000,min=0.001,max=0.995) #隨機產生10000個均勻分布
  r3.median <- median(r3) # 隨機亂數中先取中位數
  F_Q <- get(paste0("p",margin.dist[1]))(givenQ.table[q], margin.par[margin.num[1],1], margin.par[margin.num[1],2])
  con.cdf3 <- cCopula(cbind(F_Q,r3.median), copula = copula.func,indices = 2, inverse = T)
  if(con.cdf3>0.995){
    con.cdf3 <- 0.995
  }
  each.m3 <- get(paste0("q",margin.dist[2]))(con.cdf3,margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  Suspended.Load <- rbind(Suspended.Load,each.m3)
}
rat.SSL3 <- data.frame(MD.rmNA.rmSSL,Suspended.Load,group="方法 2")
colnames(rat.SSL3) <- c("Year","Month","Day","Discharge","Suspended.Load","group")
print("方法2率定完畢")

# 4. 方法3 mode(眾數)
Suspended.Load <- c()
for (q in 1:length(givenQ.table)){
  F_Q <- get(paste0("p",margin.dist[1]))(givenQ.table[q], margin.par[margin.num[1],1], margin.par[margin.num[1],2])
  min.Qs <- get(paste0("q",margin.dist[2]))(0.001, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  max.Qs <- get(paste0("q",margin.dist[2]))(0.995, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  small.Qs <- seq(from=round(min.Qs),to=99999,by=1)
  middle.Qs <- append(small.Qs,seq(from=100000,to=999999,by=10))
  all.Qs <- append(middle.Qs,seq(from=1000000,to=max.Qs,by=100))
  F_L <- get(paste0("p",margin.dist[2]))(all.Qs, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  f_L <- get(paste0("d",margin.dist[2]))(all.Qs, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  con.pdf <-get(paste0(subset(candidate.copula,number==family.num)[,2],"copula.pdf"))(F_Q,F_L,cop$par)*f_L
  #max(con.pdf) # 最大PDF值
  each.m4 <- all.Qs[which.max(con.pdf)] # PDF最大值所推估的輸砂量
  Suspended.Load <- rbind(Suspended.Load,each.m4)
}
rat.SSL4 <- data.frame(MD.rmNA.rmSSL,Suspended.Load,group="方法 3")
colnames(rat.SSL4) <- c("Year","Month","Day","Discharge","Suspended.Load","group")
print("方法3率定完畢")

# 5. 方法4 Hit or Miss
Suspended.Load <- c()
for (q in 1:length(givenQ.table)){
  F_Q <- get(paste0("p",margin.dist[1]))(givenQ.table[q], margin.par[margin.num[1],1], margin.par[margin.num[1],2])
  r0.001 <- c(0.001)
  min.Qs <- get(paste0("q",margin.dist[2]))(r0.001,margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  r0.995 <- c(0.995)
  max.Qs <- get(paste0("q",margin.dist[2]))(r0.995,margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  small.Qs <- seq(from=round(min.Qs),to=99999,by=1)
  middle.Qs <- append(small.Qs,seq(from=100000,to=999999,by=10))
  all.Qs <- append(middle.Qs,seq(from=1000000,to=max.Qs,by=100))
  F_L <- get(paste0("p",margin.dist[2]))(all.Qs, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  f_L <- get(paste0("d",margin.dist[2]))(all.Qs, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  con.pdf <-get(paste0(subset(candidate.copula,number==family.num)[,2],"copula.pdf"))(F_Q,F_L,cop$par)*f_L
  mode.pdf <- max(con.pdf)
  r <- c(0,1)
  con.pdf1 <- 0
  t <- 1
  while(r[2]*mode.pdf > con.pdf1){
    r <- runif(2,max=0.995)
    u <- min.Qs+r[1]*(max.Qs-min.Qs)
    F_L1 <- get(paste0("p",margin.dist[2]))(u, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
    f_L1 <- get(paste0("d",margin.dist[2]))(u, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
    con.pdf1 <-get(paste0(subset(candidate.copula,number==family.num)[,2],"copula.pdf"))(F_Q,F_L1,cop$par)*f_L1
    t <- t+1
  }
  each.m5 <- u
  Suspended.Load <- rbind(Suspended.Load,each.m5)
  print(paste0("第",q,"次完成"))
}
rat.SSL5 <- data.frame(MD.rmNA.rmSSL,Suspended.Load,group="方法 4")
colnames(rat.SSL5) <- c("Year","Month","Day","Discharge","Suspended.Load","group")
print("方法4率定完畢")

# 合併5種率定推估值
rating.table <- rbind(rat.SSL0,rat.SSL1,rat.SSL2,rat.SSL3,rat.SSL4,rat.SSL5)
rating.table[is.na(rating.table)] <-0

file <- paste("F:/R_output/",station,"/vinecopula/",
              year[1],"到", year[y],station_ch,"率定總表(70%觀測資料).csv", sep="") #存檔路徑
write.csv(rating.table,file)

# 率定資料分組
group.obser  <- subset(rating.table,group=="觀測資料")
group.rating <- subset(rating.table,group=="率定曲線")
group.m1     <- subset(rating.table,group=="方法 1")
group.m2     <- subset(rating.table,group=="方法 2")
group.m3     <- subset(rating.table,group=="方法 3")
group.m4     <- subset(rating.table,group=="方法 4")

# 1. MSE
mse1 <- mse(group.obser$Suspended.Load,group.rating$Suspended.Load)
mse2 <- mse(group.obser$Suspended.Load,group.m1$Suspended.Load)
mse3 <- mse(group.obser$Suspended.Load,group.m2$Suspended.Load)
mse4 <- mse(group.obser$Suspended.Load,group.m3$Suspended.Load)
mse5 <- mse(group.obser$Suspended.Load,group.m4$Suspended.Load)

mse.table <- cbind(mse1,mse2,mse3,mse4,mse5)

# 2. RMSE
rmse1 <- rmse(group.obser$Suspended.Load,group.rating$Suspended.Load)
rmse2 <- rmse(group.obser$Suspended.Load,group.m1$Suspended.Load)
rmse3 <- rmse(group.obser$Suspended.Load,group.m2$Suspended.Load)
rmse4 <- rmse(group.obser$Suspended.Load,group.m3$Suspended.Load)
rmse5 <- rmse(group.obser$Suspended.Load,group.m4$Suspended.Load)

rmse.table <- cbind(rmse1,rmse2,rmse3,rmse4,rmse5)

# 3. NMSE
nmse1 <- nmse(group.obser$Suspended.Load,group.rating$Suspended.Load)
nmse2 <- nmse(group.obser$Suspended.Load,group.m1$Suspended.Load)
nmse3 <- nmse(group.obser$Suspended.Load,group.m2$Suspended.Load)
nmse4 <- nmse(group.obser$Suspended.Load,group.m3$Suspended.Load)
nmse5 <- nmse(group.obser$Suspended.Load,group.m4$Suspended.Load)
nmse.table <- cbind(nmse1,nmse2,nmse3,nmse4,nmse5)

# 4. MAPE
mape1 <- mape(group.obser$Suspended.Load, group.rating$Suspended.Load)
mape2 <- mape(group.obser$Suspended.Load, group.m1$Suspended.Load)
mape3 <- mape(group.obser$Suspended.Load, group.m2$Suspended.Load)
mape4 <- mape(group.obser$Suspended.Load, group.m3$Suspended.Load)
mape5 <- mape(group.obser$Suspended.Load, group.m4$Suspended.Load)
mape.table <- cbind(mape1,mape2,mape3,mape4,mape5)

# 5. NSE
nse1 <- nse(group.rating$Suspended.Load,group.obser$Suspended.Load)
nse2 <- nse(group.m1$Suspended.Load,group.obser$Suspended.Load)
nse3 <- nse(group.m2$Suspended.Load,group.obser$Suspended.Load)
nse4 <- nse(group.m3$Suspended.Load,group.obser$Suspended.Load)
nse5 <- nse(group.m4$Suspended.Load,group.obser$Suspended.Load)
nse.table <- cbind(nse1,nse2,nse3,nse4,nse5)

# 6.NSE.abs
nse.abs1 <- nse.abs(group.rating$Suspended.Load,group.obser$Suspended.Load)
nse.abs2 <- nse.abs(group.m1$Suspended.Load,group.obser$Suspended.Load)
nse.abs3 <- nse.abs(group.m2$Suspended.Load,group.obser$Suspended.Load)
nse.abs4 <- nse.abs(group.m3$Suspended.Load,group.obser$Suspended.Load)
nse.abs5 <- nse.abs(group.m4$Suspended.Load,group.obser$Suspended.Load)
nse.abs.table <- cbind(nse.abs1,nse.abs2,nse.abs3,nse.abs4,nse.abs5)

# 7.KGE
kge1 <- KGE(group.rating$Suspended.Load,group.obser$Suspended.Load)
kge2 <- KGE(group.m1$Suspended.Load,group.obser$Suspended.Load)
kge3 <- KGE(group.m2$Suspended.Load,group.obser$Suspended.Load)
kge4 <- KGE(group.m3$Suspended.Load,group.obser$Suspended.Load)
kge5 <- KGE(group.m4$Suspended.Load,group.obser$Suspended.Load)
kge.table <- cbind(kge1,kge2,kge3,kge4,kge5)

error.table <- rbind(mse.table,rmse.table,nmse.table,mape.table,nse.table,nse.abs.table,kge.table)
colnames(error.table) <- c("率定曲線","est.SSL1","est.SSL2","est.SSL3","est.SSL4")
rownames(error.table) <- c("mse","rmse","nmse","mape","nse","nse.abs","kge")
file <- paste("F:/R_output/",station,"/vinecopula/",
              year[1],"到", year[y],station_ch,"率定_誤差指標(70%觀測資料).csv", sep="") #存檔路徑
write.csv(error.table,file)

# 率定輸砂量推估值
setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站率定(5種推估方法).png"),width = 800, height = 600, units = "px", pointsize = 12)
ggplot(rating.table,aes(x=Discharge,y=Suspended.Load/10000,shape=group))+
  geom_point(aes(shape=group,color=group,size=group))+
  geom_line(aes(color=group))+
  scale_x_continuous(guide = "prism_minor", #x軸副刻度
                     limits = c(0, max(rating.table$Discharge)),
                     minor_breaks = seq(0, max(rating.table$Discharge), 100))+
  scale_y_continuous(guide = "prism_minor", #y軸副刻度
                     limits = c(0, max(rating.table$Suspended.Load)/10000),
                     minor_breaks = seq(0, max(rating.table$Suspended.Load)/10000, 20))+
  labs(x=TeX("$日流量Q \\,(m^3/s)$"),y=TeX("$日懸浮載輸砂量L \\,(10^4 \\,ton/d)$")) + # 座標軸名稱
  scale_color_manual(values=c("black","red","orange","darkgreen","blue","purple"))+
  scale_shape_manual(values=c(19,0,1,2,4,6))+
  scale_size_manual(values=c(4,3,3,3,3,3))+
  theme_bw()+ #白底
  theme(legend.position = c(.15,.75))+ #圖例位置座標
  theme(legend.title=element_blank())+ #隱藏圖例標題
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(prism.ticks.length=unit(.7,"lines"))+
  theme(axis.ticks.length=unit(1.2,"lines"))+
  theme(text=element_text(size=30,color="black"))+  # 字體大小
  theme(axis.text = element_text(colour = "black"))
dev.off()

#  =========== 輸砂量推估 ============
# 1. mode(中位數)
# 2. (Peng et al. 2020)
# 3. (Bezak et al. 2017)
# 4. Suspended load rating curve(非線性最小平方法)
# 5. (Di Lascio et al. 2015)原始尺度

# ====================================
print("開始驗證")
MD.onlyNA <- MD[complete.cases(MD)==F, ]
MD.onlyNA.rmNASSL <- MD.onlyNA[,-5]
ori.data.vad <- left_join(MD.onlyNA.rmNASSL,rm.ob.data)
vad.SSL0 <- data.frame(ori.data.vad,group="觀測資料")
givenQ.table <- MD.onlyNA$Discharge
estimate.SSL <- c()

# 1. Suspended load rating curve 率定曲線
# by rating curve coefficient a, b
Suspended.Load <- c()
for (q in 1:length(givenQ.table)){
  each.L1 <- a * (givenQ.table[q])^b
  Suspended.Load <- rbind(Suspended.Load,each.L1)
}
vad.SSL1 <- data.frame(MD.onlyNA.rmNASSL,Suspended.Load,group="率定曲線")
colnames(vad.SSL1) <- c("Year","Month","Day","Discharge","Suspended.Load","group")
print("率定曲線驗證完畢")

# 2.方法1 (Peng et al. 2020) 從CDF隨機取一個點
Suspended.Load <- c()
for (q in 1:length(givenQ.table)){
  r2 <- runif(1,min=0.001,max=0.995) #隨機產生1個均勻分布
  F_Q <- get(paste0("p",margin.dist[1]))(givenQ.table[q], margin.par[margin.num[1],1], margin.par[margin.num[1],2])
  con.cdf2 <- cCopula(cbind(F_Q,r2), copula = copula.func,indices = 2, inverse = T)
  if(con.cdf2>0.995){
    con.cdf2 <- 0.995
  }
  each.L2 <- get(paste0("q",margin.dist[2]))(con.cdf2,margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  get(paste0("q",margin.dist[2]))(0.995,margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  Suspended.Load <- rbind(Suspended.Load,each.L2)
}  
vad.SSL2 <- data.frame(MD.onlyNA.rmNASSL,Suspended.Load,group="方法 1")
colnames(vad.SSL2) <- c("Year","Month","Day","Discharge","Suspended.Load","group")
print("方法1驗證完畢")

# 3.方法2 (Bezak et al. 2017) 從CDF隨機取10000個點，再取中位數
Suspended.Load <- c()
for (q in 1:length(givenQ.table)){
  r3 <- runif(10000,min=0.001,max=0.995) #隨機產生10000個均勻分布
  r3.median <- median(r3) # 隨機亂數中先取中位數
  F_Q <- get(paste0("p",margin.dist[1]))(givenQ.table[q], margin.par[margin.num[1],1], margin.par[margin.num[1],2])
  con.cdf3 <- cCopula(cbind(F_Q,r3.median), copula = copula.func,indices = 2, inverse = T)
  if(con.cdf3>0.995){
    con.cdf3 <- 0.995
  }
  each.L3 <- get(paste0("q",margin.dist[2]))(con.cdf3,margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  Suspended.Load <- rbind(Suspended.Load,each.L3)
}
vad.SSL3 <- data.frame(MD.onlyNA.rmNASSL,Suspended.Load,group="方法 2")
colnames(vad.SSL3) <- c("Year","Month","Day","Discharge","Suspended.Load","group")
print("方法2驗證完畢")

# 4. 方法3 mode(眾數)
Suspended.Load <- c()
for (q in 1:length(givenQ.table)){
  F_Q <- get(paste0("p",margin.dist[1]))(givenQ.table[q], margin.par[margin.num[1],1], margin.par[margin.num[1],2])
  min.Qs <- get(paste0("q",margin.dist[2]))(0.001, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  max.Qs <- get(paste0("q",margin.dist[2]))(0.995, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  small.Qs <- seq(from=round(min.Qs),to=99999,by=1)
  middle.Qs <- append(small.Qs,seq(from=100000,to=999999,by=10))
  all.Qs <- append(middle.Qs,seq(from=1000000,to=max.Qs,by=100))
  F_L <- get(paste0("p",margin.dist[2]))(all.Qs, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  f_L <- get(paste0("d",margin.dist[2]))(all.Qs, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  con.pdf <-get(paste0(subset(candidate.copula,number==family.num)[,2],"copula.pdf"))(F_Q,F_L,cop$par)*f_L
  #max(con.pdf) # 最大PDF值
  each.L4 <- all.Qs[which.max(con.pdf)] # PDF最大值所推估的輸砂量
  Suspended.Load <- rbind(Suspended.Load,each.L4)
}
vad.SSL4 <- data.frame(MD.onlyNA.rmNASSL,Suspended.Load,group="方法 3")
colnames(vad.SSL4) <- c("Year","Month","Day","Discharge","Suspended.Load","group")
print("方法3驗證完畢")
  
# 5. 方法4 Hit or Miss
Suspended.Load <- c()
for (q in 1:length(givenQ.table)){
  r0.001 <- c(0.001)
  F_Q <- get(paste0("p",margin.dist[1]))(givenQ.table[q], margin.par[margin.num[1],1], margin.par[margin.num[1],2])
  min.Qs <- get(paste0("q",margin.dist[2]))(r0.001,margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  r0.995 <- c(0.995)
  max.Qs <- get(paste0("q",margin.dist[2]))(r0.995,margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  small.Qs <- seq(from=round(min.Qs),to=99999,by=1)
  middle.Qs <- append(small.Qs,seq(from=100000,to=999999,by=10))
  all.Qs <- append(middle.Qs,seq(from=1000000,to=max.Qs,by=100))
  F_L <- get(paste0("p",margin.dist[2]))(all.Qs, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  f_L <- get(paste0("d",margin.dist[2]))(all.Qs, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  con.pdf <-get(paste0(subset(candidate.copula,number==family.num)[,2],"copula.pdf"))(F_Q,F_L,cop$par)*f_L
  mode.pdf <- max(con.pdf)
  r <- c(0,1)
  con.pdf1 <- 0
  t <- 1
  while(r[2]*mode.pdf > con.pdf1){
    r <- runif(2,max=0.995)
    u <- min.Qs+r[1]*(max.Qs-min.Qs)
    F_L1 <- get(paste0("p",margin.dist[2]))(u, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
    f_L1 <- get(paste0("d",margin.dist[2]))(u, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
    con.pdf1 <-get(paste0(subset(candidate.copula,number==family.num)[,2],"copula.pdf"))(F_Q,F_L1,cop$par)*f_L1
    t <- t+1
  }
  print(paste0("第",q,"次完成"))
  each.L5 <- u
  Suspended.Load <- rbind(Suspended.Load,each.L5)
}
vad.SSL5 <- data.frame(MD.onlyNA.rmNASSL,Suspended.Load,group="方法 4")
colnames(vad.SSL5) <- c("Year","Month","Day","Discharge","Suspended.Load","group")
print("方法4驗證完畢")

# 合併5種方法的驗證推估值

validation.table <- rbind(vad.SSL0,vad.SSL1,vad.SSL2,vad.SSL3,vad.SSL4,vad.SSL5)
validation.table[is.na(validation.table)] <-0

file <- paste("F:/R_output/",station,"/vinecopula/",
              year[1],"到", year[y],station_ch,"驗證總表(30%觀測資料).csv", sep="") #存檔路徑
write.csv(validation.table,file)

# 驗證資料分組
group.obser <- subset(validation.table,group=="觀測資料")
group.rating <- subset(validation.table,group=="率定曲線")
group.m1 <- subset(validation.table,group=="方法 1")
group.m2 <- subset(validation.table,group=="方法 2")
group.m3 <- subset(validation.table,group=="方法 3")
group.m4 <- subset(validation.table,group=="方法 4")

# 1. MSE
mse1 <- mse(group.obser$Suspended.Load,group.rating$Suspended.Load)
mse2 <- mse(group.obser$Suspended.Load,group.m1$Suspended.Load)
mse3 <- mse(group.obser$Suspended.Load,group.m2$Suspended.Load)
mse4 <- mse(group.obser$Suspended.Load,group.m3$Suspended.Load)
mse5 <- mse(group.obser$Suspended.Load,group.m4$Suspended.Load)

mse.table <- cbind(mse1,mse2,mse3,mse4,mse5)

# 2. RMSE
rmse1 <- rmse(group.obser$Suspended.Load,group.rating$Suspended.Load)
rmse2 <- rmse(group.obser$Suspended.Load,group.m1$Suspended.Load)
rmse3 <- rmse(group.obser$Suspended.Load,group.m2$Suspended.Load)
rmse4 <- rmse(group.obser$Suspended.Load,group.m3$Suspended.Load)
rmse5 <- rmse(group.obser$Suspended.Load,group.m4$Suspended.Load)

rmse.table <- cbind(rmse1,rmse2,rmse3,rmse4,rmse5)

# 3. NMSE
nmse1 <- nmse(group.obser$Suspended.Load,group.rating$Suspended.Load)
nmse2 <- nmse(group.obser$Suspended.Load,group.m1$Suspended.Load)
nmse3 <- nmse(group.obser$Suspended.Load,group.m2$Suspended.Load)
nmse4 <- nmse(group.obser$Suspended.Load,group.m3$Suspended.Load)
nmse5 <- nmse(group.obser$Suspended.Load,group.m4$Suspended.Load)
nmse.table <- cbind(nmse1,nmse2,nmse3,nmse4,nmse5)

# 4. MAPE
mape1 <- mape(group.obser$Suspended.Load, group.rating$Suspended.Load)
mape2 <- mape(group.obser$Suspended.Load, group.m1$Suspended.Load)
mape3 <- mape(group.obser$Suspended.Load, group.m2$Suspended.Load)
mape4 <- mape(group.obser$Suspended.Load, group.m3$Suspended.Load)
mape5 <- mape(group.obser$Suspended.Load, group.m4$Suspended.Load)
mape.table <- cbind(mape1,mape2,mape3,mape4,mape5)

# 5. NSE
nse1 <- nse(group.rating$Suspended.Load,group.obser$Suspended.Load)
nse2 <- nse(group.m1$Suspended.Load,group.obser$Suspended.Load)
nse3 <- nse(group.m2$Suspended.Load,group.obser$Suspended.Load)
nse4 <- nse(group.m3$Suspended.Load,group.obser$Suspended.Load)
nse5 <- nse(group.m4$Suspended.Load,group.obser$Suspended.Load)
nse.table <- cbind(nse1,nse2,nse3,nse4,nse5)

# 6. NSE.abs
nse.abs1 <- nse.abs(group.rating$Suspended.Load,group.obser$Suspended.Load)
nse.abs2 <- nse.abs(group.m1$Suspended.Load,group.obser$Suspended.Load)
nse.abs3 <- nse.abs(group.m2$Suspended.Load,group.obser$Suspended.Load)
nse.abs4 <- nse.abs(group.m3$Suspended.Load,group.obser$Suspended.Load)
nse.abs5 <- nse.abs(group.m4$Suspended.Load,group.obser$Suspended.Load)
nse.abs.table <- cbind(nse.abs1,nse.abs2,nse.abs3,nse.abs4,nse.abs5)

# 7.KGE
kge1 <- KGE(group.rating$Suspended.Load,group.obser$Suspended.Load)
kge2 <- KGE(group.m1$Suspended.Load,group.obser$Suspended.Load)
kge3 <- KGE(group.m2$Suspended.Load,group.obser$Suspended.Load)
kge4 <- KGE(group.m3$Suspended.Load,group.obser$Suspended.Load)
kge5 <- KGE(group.m4$Suspended.Load,group.obser$Suspended.Load)
kge.table <- cbind(kge1,kge2,kge3,kge4,kge5)

error.table <- rbind(mse.table,rmse.table,nmse.table,mape.table,nse.table,nse.abs.table,kge.table)
colnames(error.table) <- c("率定曲線","est.SSL1","est.SSL2","est.SSL3","est.SSL4")
rownames(error.table) <- c("mse","rmse","nmse","mape","nse","nse.abs","kge")
file <- paste("F:/R_output/",station,"/vinecopula/",
              year[1],"到", year[y],station_ch,"驗證_誤差指標(30%觀測資料).csv", sep="") #存檔路徑
write.csv(error.table,file)

table <- rbind(data.frame(MD.rmNA,type="率定時期 (1960-2000)"),
               data.frame(ori.data.vad,type="驗證時期 (2001-2019)"))

setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站模式率定&驗證.png"),width = 800, height = 600, units = "px", pointsize = 12)
ggplot(data=table,aes(x=Discharge,y=Suspended.Load/10000,shape=type,color=type))+
  geom_point(size=3)+
  scale_x_continuous(guide = "prism_minor",
                     limits = c(0, max(table$Discharge)),
                     minor_breaks = seq(0, max(table$Discharge), 100))+
  scale_y_continuous(guide = "prism_minor",
                     limits = c(0, max(table$Suspended.Load)/10000),
                     minor_breaks = seq(0, max(table$Suspended.Load), 20))+
  labs(x=TeX("日流量Q (cms)"),y=TeX("$日懸浮載輸砂量L \\,(10^4 \\,ton/d )$")) + # 座標軸名稱
  scale_shape_manual(values=c(1,4))+  #圖例
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=30))+# 字體大小
  theme(axis.text = element_text(colour = "black"))+
  theme(legend.position=c(0.25,0.8))+
  theme(legend.title=element_blank())+
  theme(prism.ticks.length=unit(.5,"lines"))+
  theme(axis.ticks.length=unit(1,"lines"))
dev.off()

rating.curve <- data.frame(rm.ob.data[,-5],Suspended.Load=ratingSSL,type="率定曲線")
rating.curve.table <- rbind(table,rating.curve)
setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站率定曲線&觀測值.png"),width = 800, height = 600, units = "px", pointsize = 12)
ggplot(data=rating.curve.table,aes(x=Discharge,y=Suspended.Load/10000))+
  geom_point(data=subset(rating.curve.table,type!="率定曲線"),aes(shape=type),size=3)+
  geom_line(data=subset(rating.curve.table,type=="率定曲線"),aes(color=type),size=1.2)+
  scale_x_continuous(guide = "prism_minor",
                     limits = c(0, max(rating.curve.table$Discharge)),
                     minor_breaks = seq(0, max(rating.curve.table$Discharge), 100))+
  scale_y_continuous(guide = "prism_minor",
                     limits = c(0, max(rating.curve.table$Suspended.Load)/10000),
                     minor_breaks = seq(0, max(rating.curve.table$Suspended.Load), 20))+
  labs(x=TeX("$日流量Q \\,(m^3/s)$"),y=TeX("$日懸浮載輸砂量L \\, (10^4 \\,ton/d)$")) + # 座標軸名稱
  scale_shape_manual(values=c(1,4))+  #圖例
  theme_bw() + # 白底
  #theme(legend.position = "none")+ #圖例位置座標
  #theme(legend.title=element_blank())+ #隱藏圖例標題
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=30))+# 字體大小
  theme(legend.position=c(0.25,0.8))+
  theme(legend.title=element_blank())+
  theme(prism.ticks.length=unit(.5,"lines"))+
  theme(axis.ticks.length=unit(1,"lines"))
dev.off()

# 驗證輸砂量推估值
setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站驗證(5種推估方法).png"),width = 800, height = 600, units = "px", pointsize = 12)
ggplot(data=validation.table,aes(x=Discharge,y=Suspended.Load/10000))+
  geom_point(aes(shape=group,color=group,size=group))+
  geom_line(aes(color=group))+
  scale_x_continuous(guide = "prism_minor", #x軸副刻度
                     limits = c(0, max(validation.table$Discharge)),
                     minor_breaks = seq(0, max(validation.table$Discharge), 100))+
  scale_y_continuous(guide = "prism_minor", #y軸副刻度
                     limits = c(0, max(validation.table$Suspended.Load)/10000),
                     minor_breaks = seq(0, max(validation.table$Suspended.Load)/10000, 20))+
  scale_color_manual(values=c("black","red","orange","darkgreen","blue","purple"))+
  scale_shape_manual(values=c(19,0,1,2,4,6))+
  scale_size_manual(values=c(4,3,3,3,3,3))+
  labs(x=TeX("$日流量Q \\,(m^3/s)$"),y=TeX("$日懸浮載輸砂量L \\,(10^4 \\,ton/d)$")) + # 座標軸名稱
  theme_bw()+ #白底
  theme(legend.position = c(.15,.75))+ #圖例位置座標
  theme(legend.title=element_blank())+ #隱藏圖例標題
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(prism.ticks.length=unit(.7,"lines"))+
  theme(axis.ticks.length=unit(1.2,"lines"))+
  theme(text=element_text(size=30,color="black"))+  # 字體大小
  theme(axis.text = element_text(colour = "black"))
dev.off()

#
# =========================== 討論 ==============================================
# 率定+驗證一起執行誤差指標
all.table <- rbind(rating.table,validation.table)

# 驗證資料分組
group.obser <- subset(all.table,group=="觀測資料")
group.rating <- subset(all.table,group=="率定曲線")
group.m1 <- subset(all.table,group=="方法 1")
group.m2 <- subset(all.table,group=="方法 2")
group.m3 <- subset(all.table,group=="方法 3")
group.m4 <- subset(all.table,group=="方法 4")

# 1. MSE
mse1 <- mse(group.obser$Suspended.Load,group.rating$Suspended.Load)
mse2 <- mse(group.obser$Suspended.Load,group.m1$Suspended.Load)
mse3 <- mse(group.obser$Suspended.Load,group.m2$Suspended.Load)
mse4 <- mse(group.obser$Suspended.Load,group.m3$Suspended.Load)
mse5 <- mse(group.obser$Suspended.Load,group.m4$Suspended.Load)

mse.table <- cbind(mse1,mse2,mse3,mse4,mse5)

# 2. RMSE
rmse1 <- rmse(group.obser$Suspended.Load,group.rating$Suspended.Load)
rmse2 <- rmse(group.obser$Suspended.Load,group.m1$Suspended.Load)
rmse3 <- rmse(group.obser$Suspended.Load,group.m2$Suspended.Load)
rmse4 <- rmse(group.obser$Suspended.Load,group.m3$Suspended.Load)
rmse5 <- rmse(group.obser$Suspended.Load,group.m4$Suspended.Load)

rmse.table <- cbind(rmse1,rmse2,rmse3,rmse4,rmse5)

# 3. NMSE
nmse1 <- nmse(group.obser$Suspended.Load,group.rating$Suspended.Load)
nmse2 <- nmse(group.obser$Suspended.Load,group.m1$Suspended.Load)
nmse3 <- nmse(group.obser$Suspended.Load,group.m2$Suspended.Load)
nmse4 <- nmse(group.obser$Suspended.Load,group.m3$Suspended.Load)
nmse5 <- nmse(group.obser$Suspended.Load,group.m4$Suspended.Load)
nmse.table <- cbind(nmse1,nmse2,nmse3,nmse4,nmse5)

# 4. MAPE
mape1 <- mape(group.obser$Suspended.Load, group.rating$Suspended.Load)
mape2 <- mape(group.obser$Suspended.Load, group.m1$Suspended.Load)
mape3 <- mape(group.obser$Suspended.Load, group.m2$Suspended.Load)
mape4 <- mape(group.obser$Suspended.Load, group.m3$Suspended.Load)
mape5 <- mape(group.obser$Suspended.Load, group.m4$Suspended.Load)
mape.table <- cbind(mape1,mape2,mape3,mape4,mape5)

# 5. NSE
nse1 <- nse(group.rating$Suspended.Load,group.obser$Suspended.Load)
nse2 <- nse(group.m1$Suspended.Load,group.obser$Suspended.Load)
nse3 <- nse(group.m2$Suspended.Load,group.obser$Suspended.Load)
nse4 <- nse(group.m3$Suspended.Load,group.obser$Suspended.Load)
nse5 <- nse(group.m4$Suspended.Load,group.obser$Suspended.Load)
nse.table <- cbind(nse1,nse2,nse3,nse4,nse5)

# 6. NSE.abs
nse.abs1 <- nse.abs(group.rating$Suspended.Load,group.obser$Suspended.Load)
nse.abs2 <- nse.abs(group.m1$Suspended.Load,group.obser$Suspended.Load)
nse.abs3 <- nse.abs(group.m2$Suspended.Load,group.obser$Suspended.Load)
nse.abs4 <- nse.abs(group.m3$Suspended.Load,group.obser$Suspended.Load)
nse.abs5 <- nse.abs(group.m4$Suspended.Load,group.obser$Suspended.Load)
nse.abs.table <- cbind(nse.abs1,nse.abs2,nse.abs3,nse.abs4,nse.abs5)

# 7.KGE
kge1 <- KGE(group.rating$Suspended.Load,group.obser$Suspended.Load)
kge2 <- KGE(group.m1$Suspended.Load,group.obser$Suspended.Load)
kge3 <- KGE(group.m2$Suspended.Load,group.obser$Suspended.Load)
kge4 <- KGE(group.m3$Suspended.Load,group.obser$Suspended.Load)
kge5 <- KGE(group.m4$Suspended.Load,group.obser$Suspended.Load)
kge.table <- cbind(kge1,kge2,kge3,kge4,kge5)

error.table <- rbind(mse.table,rmse.table,nmse.table,mape.table,nse.table,nse.abs.table,kge.table)
colnames(error.table) <- c("率定曲線","est.SSL1","est.SSL2","est.SSL3","est.SSL4")
rownames(error.table) <- c("mse","rmse","nmse","mape","nse","nse.abs","kge")
file <- paste("F:/R_output/",station,"/vinecopula/",
              year[1],"到", year[y],station_ch,"率定+驗證_誤差指標(100%觀測資料).csv", sep="") #存檔路徑
write.csv(error.table,file)

#

# 1. 以不同流量分組比較誤差指標
# ---- 1.1率定+驗證時期 -------
Q0.3 <- sort(MD$Discharge)[round(length(MD$Discharge)*0.3)]
Q0.7 <- sort(MD$Discharge)[round(length(MD$Discharge)*0.7)]
group0 <- all.table
group1 <- subset(all.table,Discharge<Q0.3)
group2 <- subset(all.table,Discharge>=Q0.3&Discharge<Q0.7)
group3 <- subset(all.table,Discharge>Q0.7)

for(i in 0:3){
  group.obser <- subset(get(paste0("group",i)),group=="觀測資料")
  group.rating <- subset(get(paste0("group",i)),group=="率定曲線")
  group.m1 <- subset(get(paste0("group",i)),group=="方法 1")
  group.m2 <- subset(get(paste0("group",i)),group=="方法 2")
  group.m3 <- subset(get(paste0("group",i)),group=="方法 3")
  group.m4 <- subset(get(paste0("group",i)),group=="方法 4")
  
  # 1. MSE
  mse1 <- mse(group.obser$Suspended.Load,group.rating$Suspended.Load)
  mse2 <- mse(group.obser$Suspended.Load,group.m1$Suspended.Load)
  mse3 <- mse(group.obser$Suspended.Load,group.m2$Suspended.Load)
  mse4 <- mse(group.obser$Suspended.Load,group.m3$Suspended.Load)
  mse5 <- mse(group.obser$Suspended.Load,group.m4$Suspended.Load)
  
  mse.table <- cbind(mse1,mse2,mse3,mse4,mse5)
  
  # 2. RMSE
  rmse1 <- rmse(group.obser$Suspended.Load,group.rating$Suspended.Load)
  rmse2 <- rmse(group.obser$Suspended.Load,group.m1$Suspended.Load)
  rmse3 <- rmse(group.obser$Suspended.Load,group.m2$Suspended.Load)
  rmse4 <- rmse(group.obser$Suspended.Load,group.m3$Suspended.Load)
  rmse5 <- rmse(group.obser$Suspended.Load,group.m4$Suspended.Load)
  
  rmse.table <- cbind(rmse1,rmse2,rmse3,rmse4,rmse5)
  
  # 3. NMSE
  nmse1 <- nmse(group.obser$Suspended.Load,group.rating$Suspended.Load)
  nmse2 <- nmse(group.obser$Suspended.Load,group.m1$Suspended.Load)
  nmse3 <- nmse(group.obser$Suspended.Load,group.m2$Suspended.Load)
  nmse4 <- nmse(group.obser$Suspended.Load,group.m3$Suspended.Load)
  nmse5 <- nmse(group.obser$Suspended.Load,group.m4$Suspended.Load)
  nmse.table <- cbind(nmse1,nmse2,nmse3,nmse4,nmse5)
  
  # 4. MAPE
  mape1 <- mape(group.obser$Suspended.Load, group.rating$Suspended.Load)
  mape2 <- mape(group.obser$Suspended.Load, group.m1$Suspended.Load)
  mape3 <- mape(group.obser$Suspended.Load, group.m2$Suspended.Load)
  mape4 <- mape(group.obser$Suspended.Load, group.m3$Suspended.Load)
  mape5 <- mape(group.obser$Suspended.Load, group.m4$Suspended.Load)
  mape.table <- cbind(mape1,mape2,mape3,mape4,mape5)
  
  # 5. NSE
  nse1 <- nse(group.rating$Suspended.Load,group.obser$Suspended.Load)
  nse2 <- nse(group.m1$Suspended.Load,group.obser$Suspended.Load)
  nse3 <- nse(group.m2$Suspended.Load,group.obser$Suspended.Load)
  nse4 <- nse(group.m3$Suspended.Load,group.obser$Suspended.Load)
  nse5 <- nse(group.m4$Suspended.Load,group.obser$Suspended.Load)
  nse.table <- cbind(nse1,nse2,nse3,nse4,nse5)
  
  # 6. NSE.abs
  nse.abs1 <- nse.abs(group.rating$Suspended.Load,group.obser$Suspended.Load)
  nse.abs2 <- nse.abs(group.m1$Suspended.Load,group.obser$Suspended.Load)
  nse.abs3 <- nse.abs(group.m2$Suspended.Load,group.obser$Suspended.Load)
  nse.abs4 <- nse.abs(group.m3$Suspended.Load,group.obser$Suspended.Load)
  nse.abs5 <- nse.abs(group.m4$Suspended.Load,group.obser$Suspended.Load)
  nse.abs.table <- cbind(nse.abs1,nse.abs2,nse.abs3,nse.abs4,nse.abs5)
  
  # 7.KGE
  kge1 <- KGE(group.rating$Suspended.Load,group.obser$Suspended.Load)
  kge2 <- KGE(group.m1$Suspended.Load,group.obser$Suspended.Load)
  kge3 <- KGE(group.m2$Suspended.Load,group.obser$Suspended.Load)
  kge4 <- KGE(group.m3$Suspended.Load,group.obser$Suspended.Load)
  kge5 <- KGE(group.m4$Suspended.Load,group.obser$Suspended.Load)
  kge.table <- cbind(kge1,kge2,kge3,kge4,kge5)
  
  error.table <- rbind(mse.table,rmse.table,nmse.table,mape.table,nse.table,nse.abs.table,kge.table)
  colnames(error.table) <- c("率定曲線","est.SSL1","est.SSL2","est.SSL3","est.SSL4")
  rownames(error.table) <- c("mse","rmse","nmse","mape","nse","nse.abs","kge")
  file <- paste("F:/R_output/",station,"/vinecopula/",
                year[1],"到", year[y],station_ch,"率定+驗證_誤差指標(第",i,"組).csv", sep="") #存檔路徑
  write.csv(error.table,file)
}

# ---- 1.2.1 率定時期 -------
# # (all-30%, 70%)
# Q0.3 <- sort(MD$Discharge)[round(length(MD$Discharge)*0.3)]
# Q0.7 <- sort(MD$Discharge)[round(length(MD$Discharge)*0.7)]

#(rating-30%, 70%)
Q0.3 <- sort(MD.rmNA$Discharge)[round(length(MD.rmNA$Discharge)*0.3)]
Q0.7 <- sort(MD.rmNA$Discharge)[round(length(MD.rmNA$Discharge)*0.7)]

group0 <- rating.table
group1 <- subset(rating.table,Discharge<Q0.3)
group2 <- subset(rating.table,Discharge>=Q0.3&Discharge<Q0.7)
group3 <- subset(rating.table,Discharge>Q0.7)

for(i in 0:3){
  group.obser <- subset(get(paste0("group",i)),group=="觀測資料")
  group.rating <- subset(get(paste0("group",i)),group=="率定曲線")
  group.m1 <- subset(get(paste0("group",i)),group=="方法 1")
  group.m2 <- subset(get(paste0("group",i)),group=="方法 2")
  group.m3 <- subset(get(paste0("group",i)),group=="方法 3")
  group.m4 <- subset(get(paste0("group",i)),group=="方法 4")
  
  # 1. MSE
  mse1 <- mse(group.obser$Suspended.Load,group.rating$Suspended.Load)
  mse2 <- mse(group.obser$Suspended.Load,group.m1$Suspended.Load)
  mse3 <- mse(group.obser$Suspended.Load,group.m2$Suspended.Load)
  mse4 <- mse(group.obser$Suspended.Load,group.m3$Suspended.Load)
  mse5 <- mse(group.obser$Suspended.Load,group.m4$Suspended.Load)
  
  mse.table <- cbind(mse1,mse2,mse3,mse4,mse5)
  
  # 2. RMSE
  rmse1 <- rmse(group.obser$Suspended.Load,group.rating$Suspended.Load)
  rmse2 <- rmse(group.obser$Suspended.Load,group.m1$Suspended.Load)
  rmse3 <- rmse(group.obser$Suspended.Load,group.m2$Suspended.Load)
  rmse4 <- rmse(group.obser$Suspended.Load,group.m3$Suspended.Load)
  rmse5 <- rmse(group.obser$Suspended.Load,group.m4$Suspended.Load)
  
  rmse.table <- cbind(rmse1,rmse2,rmse3,rmse4,rmse5)
  
  # 3. NMSE
  nmse1 <- nmse(group.obser$Suspended.Load,group.rating$Suspended.Load)
  nmse2 <- nmse(group.obser$Suspended.Load,group.m1$Suspended.Load)
  nmse3 <- nmse(group.obser$Suspended.Load,group.m2$Suspended.Load)
  nmse4 <- nmse(group.obser$Suspended.Load,group.m3$Suspended.Load)
  nmse5 <- nmse(group.obser$Suspended.Load,group.m4$Suspended.Load)
  nmse.table <- cbind(nmse1,nmse2,nmse3,nmse4,nmse5)
  
  # 4. MAPE
  mape1 <- mape(group.obser$Suspended.Load, group.rating$Suspended.Load)
  mape2 <- mape(group.obser$Suspended.Load, group.m1$Suspended.Load)
  mape3 <- mape(group.obser$Suspended.Load, group.m2$Suspended.Load)
  mape4 <- mape(group.obser$Suspended.Load, group.m3$Suspended.Load)
  mape5 <- mape(group.obser$Suspended.Load, group.m4$Suspended.Load)
  mape.table <- cbind(mape1,mape2,mape3,mape4,mape5)
  
  # 5. NSE
  nse1 <- nse(group.rating$Suspended.Load,group.obser$Suspended.Load)
  nse2 <- nse(group.m1$Suspended.Load,group.obser$Suspended.Load)
  nse3 <- nse(group.m2$Suspended.Load,group.obser$Suspended.Load)
  nse4 <- nse(group.m3$Suspended.Load,group.obser$Suspended.Load)
  nse5 <- nse(group.m4$Suspended.Load,group.obser$Suspended.Load)
  nse.table <- cbind(nse1,nse2,nse3,nse4,nse5)
  
  # 6. NSE.abs
  nse.abs1 <- nse.abs(group.rating$Suspended.Load,group.obser$Suspended.Load)
  nse.abs2 <- nse.abs(group.m1$Suspended.Load,group.obser$Suspended.Load)
  nse.abs3 <- nse.abs(group.m2$Suspended.Load,group.obser$Suspended.Load)
  nse.abs4 <- nse.abs(group.m3$Suspended.Load,group.obser$Suspended.Load)
  nse.abs5 <- nse.abs(group.m4$Suspended.Load,group.obser$Suspended.Load)
  nse.abs.table <- cbind(nse.abs1,nse.abs2,nse.abs3,nse.abs4,nse.abs5)
  
  # 7.KGE
  kge1 <- KGE(group.rating$Suspended.Load,group.obser$Suspended.Load)
  kge2 <- KGE(group.m1$Suspended.Load,group.obser$Suspended.Load)
  kge3 <- KGE(group.m2$Suspended.Load,group.obser$Suspended.Load)
  kge4 <- KGE(group.m3$Suspended.Load,group.obser$Suspended.Load)
  kge5 <- KGE(group.m4$Suspended.Load,group.obser$Suspended.Load)
  kge.table <- cbind(kge1,kge2,kge3,kge4,kge5)
  
  error.table <- rbind(mse.table,rmse.table,nmse.table,mape.table,nse.table,nse.abs.table,kge.table)
  colnames(error.table) <- c("率定曲線","est.SSL1","est.SSL2","est.SSL3","est.SSL4")
  rownames(error.table) <- c("mse","rmse","nmse","mape","nse","nse.abs","kge")
  file <- paste("F:/R_output/",station,"/vinecopula/",
                year[1],"到", year[y],station_ch,"率定+驗證_誤差指標(第",i,"組).csv", sep="") #存檔路徑
  write.csv(error.table,file)
}

# ---- 1.3 驗證時期 -------

# # (all-30%, 70%)
# Q0.3 <- sort(MD$Discharge)[round(length(MD$Discharge)*0.3)]
# Q0.7 <- sort(MD$Discharge)[round(length(MD$Discharge)*0.7)]

#(validation-30%, 70%)
Q0.3 <- sort(MD.onlyNA$Discharge)[round(length(MD.onlyNA$Discharge)*0.3)]
Q0.7 <- sort(MD.onlyNA$Discharge)[round(length(MD.onlyNA$Discharge)*0.7)]

group0 <- validation.table
group1 <- subset(validation.table,Discharge<Q0.3)
group2 <- subset(validation.table,Discharge>=Q0.3&Discharge<Q0.7)
group3 <- subset(validation.table,Discharge>Q0.7)

for(i in 0:3){
  group.obser <- subset(get(paste0("group",i)),group=="觀測資料")
  group.rating <- subset(get(paste0("group",i)),group=="率定曲線")
  group.m1 <- subset(get(paste0("group",i)),group=="方法 1")
  group.m2 <- subset(get(paste0("group",i)),group=="方法 2")
  group.m3 <- subset(get(paste0("group",i)),group=="方法 3")
  group.m4 <- subset(get(paste0("group",i)),group=="方法 4")
  
  # 1. MSE
  mse1 <- mse(group.obser$Suspended.Load,group.rating$Suspended.Load)
  mse2 <- mse(group.obser$Suspended.Load,group.m1$Suspended.Load)
  mse3 <- mse(group.obser$Suspended.Load,group.m2$Suspended.Load)
  mse4 <- mse(group.obser$Suspended.Load,group.m3$Suspended.Load)
  mse5 <- mse(group.obser$Suspended.Load,group.m4$Suspended.Load)
  
  mse.table <- cbind(mse1,mse2,mse3,mse4,mse5)
  
  # 2. RMSE
  rmse1 <- rmse(group.obser$Suspended.Load,group.rating$Suspended.Load)
  rmse2 <- rmse(group.obser$Suspended.Load,group.m1$Suspended.Load)
  rmse3 <- rmse(group.obser$Suspended.Load,group.m2$Suspended.Load)
  rmse4 <- rmse(group.obser$Suspended.Load,group.m3$Suspended.Load)
  rmse5 <- rmse(group.obser$Suspended.Load,group.m4$Suspended.Load)
  
  rmse.table <- cbind(rmse1,rmse2,rmse3,rmse4,rmse5)
  
  # 3. NMSE
  nmse1 <- nmse(group.obser$Suspended.Load,group.rating$Suspended.Load)
  nmse2 <- nmse(group.obser$Suspended.Load,group.m1$Suspended.Load)
  nmse3 <- nmse(group.obser$Suspended.Load,group.m2$Suspended.Load)
  nmse4 <- nmse(group.obser$Suspended.Load,group.m3$Suspended.Load)
  nmse5 <- nmse(group.obser$Suspended.Load,group.m4$Suspended.Load)
  nmse.table <- cbind(nmse1,nmse2,nmse3,nmse4,nmse5)
  
  # 4. MAPE
  mape1 <- mape(group.obser$Suspended.Load, group.rating$Suspended.Load)
  mape2 <- mape(group.obser$Suspended.Load, group.m1$Suspended.Load)
  mape3 <- mape(group.obser$Suspended.Load, group.m2$Suspended.Load)
  mape4 <- mape(group.obser$Suspended.Load, group.m3$Suspended.Load)
  mape5 <- mape(group.obser$Suspended.Load, group.m4$Suspended.Load)
  mape.table <- cbind(mape1,mape2,mape3,mape4,mape5)
  
  # 5. NSE
  nse1 <- nse(group.rating$Suspended.Load,group.obser$Suspended.Load)
  nse2 <- nse(group.m1$Suspended.Load,group.obser$Suspended.Load)
  nse3 <- nse(group.m2$Suspended.Load,group.obser$Suspended.Load)
  nse4 <- nse(group.m3$Suspended.Load,group.obser$Suspended.Load)
  nse5 <- nse(group.m4$Suspended.Load,group.obser$Suspended.Load)
  nse.table <- cbind(nse1,nse2,nse3,nse4,nse5)
  
  # 6. NSE.abs
  nse.abs1 <- nse.abs(group.rating$Suspended.Load,group.obser$Suspended.Load)
  nse.abs2 <- nse.abs(group.m1$Suspended.Load,group.obser$Suspended.Load)
  nse.abs3 <- nse.abs(group.m2$Suspended.Load,group.obser$Suspended.Load)
  nse.abs4 <- nse.abs(group.m3$Suspended.Load,group.obser$Suspended.Load)
  nse.abs5 <- nse.abs(group.m4$Suspended.Load,group.obser$Suspended.Load)
  nse.abs.table <- cbind(nse.abs1,nse.abs2,nse.abs3,nse.abs4,nse.abs5)
  
  # 7.KGE
  kge1 <- KGE(group.rating$Suspended.Load,group.obser$Suspended.Load)
  kge2 <- KGE(group.m1$Suspended.Load,group.obser$Suspended.Load)
  kge3 <- KGE(group.m2$Suspended.Load,group.obser$Suspended.Load)
  kge4 <- KGE(group.m3$Suspended.Load,group.obser$Suspended.Load)
  kge5 <- KGE(group.m4$Suspended.Load,group.obser$Suspended.Load)
  kge.table <- cbind(kge1,kge2,kge3,kge4,kge5)
  
  error.table <- rbind(mse.table,rmse.table,nmse.table,mape.table,nse.table,nse.abs.table,kge.table)
  colnames(error.table) <- c("率定曲線","est.SSL1","est.SSL2","est.SSL3","est.SSL4")
  rownames(error.table) <- c("mse","rmse","nmse","mape","nse","nse.abs","kge")
  file <- paste("F:/R_output/",station,"/vinecopula/",
                year[1],"到", year[y],station_ch,"率定+驗證_誤差指標(第",i,"組).csv", sep="") #存檔路徑
  write.csv(error.table,file)
}
# 驗證輸砂量推估值
setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站率定+驗證(5種推估方法).png"),width = 800, height = 600, units = "px", pointsize = 12)
ggplot(all.table,aes(x=Discharge,y=Suspended.Load/10000,shape=group,color=group))+
  geom_point(size=2.5)+
  scale_x_continuous(guide = "prism_minor", #x軸副刻度
                     limits = c(0, max(validation.table$Discharge)),
                     minor_breaks = seq(0, max(validation.table$Discharge), 100))+
  scale_y_continuous(guide = "prism_minor", #y軸副刻度
                     limits = c(0, max(validation.table$Suspended.Load)/10000),
                     minor_breaks = seq(0, max(validation.table$Suspended.Load)/10000, 20))+
  labs(x=TeX("$日流量Q \\,(m^3/s)$"),y=TeX("$日懸浮載輸砂量L \\,(10^4 \\,ton/d)$")) + # 座標軸名稱
  theme_bw()+ #白底
  theme(legend.position = c(.15,.75))+ #圖例位置座標
  theme(legend.title=element_blank())+ #隱藏圖例標題
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(prism.ticks.length=unit(.7,"lines"))+
  theme(axis.ticks.length=unit(1.2,"lines"))+
  theme(text=element_text(size=30,color="black"))+  # 字體大小
  theme(axis.text = element_text(colour = "black"))
dev.off()

  #推估值與觀測值
i <- 2 # 1 or2
setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站討論(5種推估方法)第",i,"組.png"),width = 800, height = 600, units = "px", pointsize = 12)
ggplot(get(paste0("group",i)),aes(x=Discharge,y=Suspended.Load/10000))+
  geom_point(aes(shape=group,color=group,size=group))+
  geom_line(aes(color=group))+
  scale_x_continuous(guide = "prism_minor", #x軸副刻度
                     limits = c(min(get(paste0("group",i))$Discharge), max(get(paste0("group",i))$Discharge)),
                     minor_breaks = seq(0, max(get(paste0("group",i))$Discharge), 2.5))+
  scale_y_continuous(guide = "prism_minor", #y軸副刻度
                     limits = c(min(get(paste0("group",i))$Suspended.Load)/10000, max(get(paste0("group",i))$Suspended.Load)/10000),
                     minor_breaks = seq(min(get(paste0("group",i))$Suspended.Load)/10000, max(get(paste0("group",i))$Suspended.Load)/10000, 1))+
  scale_color_manual(values=c("black","red","orange","darkgreen","blue","purple"))+
  scale_shape_manual(values=c(19,0,1,2,4,6))+
  scale_size_manual(values=c(4,3,3,3,3,3))+
  labs(x=TeX("$日流量Q \\,(m^3/s)$"),y=TeX("$日懸浮載輸砂量L \\,(10^4 \\,ton/d)$")) + # 座標軸名稱
  theme_bw()+ #白底
  theme(legend.position = c(.15,.75))+ #圖例位置座標
  theme(legend.title=element_blank())+ #隱藏圖例標題
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(prism.ticks.length=unit(.7,"lines"))+
  theme(axis.ticks.length=unit(1.2,"lines"))+
  theme(text=element_text(size=30,color="black"))+  # 字體大小
  theme(axis.text = element_text(colour = "black"))
dev.off()




i <- 0 # 0 or 3
setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站討論(5種推估方法)第",i,"組.png"),width = 800, height = 600, units = "px", pointsize = 12)
ggplot(get(paste0("group",i)),aes(x=Discharge,y=Suspended.Load/10000))+
  geom_point(aes(shape=group,color=group,size=group))+
  geom_line(aes(color=group))+
  scale_x_continuous(guide = "prism_minor", #x軸副刻度
                     limits = c(0, max(get(paste0("group",i))$Discharge)),
                     minor_breaks = seq(0, max(get(paste0("group",i))$Discharge), 100))+
  scale_y_continuous(guide = "prism_minor", #y軸副刻度
                     limits = c(0, max(get(paste0("group",i))$Suspended.Load)/10000),
                     minor_breaks = seq(0, max(get(paste0("group",i))$Suspended.Load)/10000, 20))+
  scale_color_manual(values=c("black","red","orange","darkgreen","blue","purple"))+
  scale_shape_manual(values=c(19,0,1,2,4,6))+
  scale_size_manual(values=c(4,3,3,3,3,3))+
  labs(x=TeX("$日流量Q \\,(m^3/s)$"),y=TeX("$日懸浮載輸砂量L \\,(10^4 \\,ton/d)$")) + # 座標軸名稱
  theme_bw()+ #白底
  theme(legend.position = c(.15,.75))+ #圖例位置座標
  theme(legend.title=element_blank())+ #隱藏圖例標題
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(prism.ticks.length=unit(.7,"lines"))+
  theme(axis.ticks.length=unit(1.2,"lines"))+
  theme(text=element_text(size=30,color="black"))+  # 字體大小
  theme(axis.text = element_text(colour = "black"))
dev.off()

for (i in 0:3){
  l1 <- summary(subset(get(paste0("group",i)),group=="觀測資料")$Suspended.Load)
  l2 <- summary(subset(get(paste0("group",i)),group=="率定曲線")$Suspended.Load)
  l3 <- summary(subset(get(paste0("group",i)),group=="方法 1")$Suspended.Load)
  l4 <- summary(subset(get(paste0("group",i)),group=="方法 2")$Suspended.Load)
  l5 <- summary(subset(get(paste0("group",i)),group=="方法 3")$Suspended.Load)
  l6 <- summary(subset(get(paste0("group",i)),group=="方法 4")$Suspended.Load)
  stat.table <- cbind(l1,l2,l3,l4,l5,l6)
  colnames(stat.table) <- c("觀測資料","率定曲線","方法1","方法2","方法3","方法4")
  file <- paste("F:/R_output/",station,"/vinecopula/",
                year[1],"到", year[y],station_ch,"率定+驗證_敘述統計(第",i,"組).csv", sep="") #存檔路徑
  write.csv(stat.table,file)
}

# 2.比較相同或相近流量產生之輸砂量
dis <- subset(all.table,Discharge<9.4&Discharge>9)
dis$Discharge <- paste0(dis$Discharge,"cms")
setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站討論(5種推估方法)9cms比較.png"),width = 1000, height = 600, units = "px", pointsize = 12)
ggplot(dis,aes(x=group,y=Suspended.Load,fill=Discharge,label=Suspended.Load))+
  geom_bar(position = "dodge",stat="identity")+
  geom_text(aes(label=round(Suspended.Load)), position=position_dodge(width=1), size=6,vjust=-0.25)+
  #geom_label_repel(min.segment.length = 0, box.padding = 0.5)
  labs(x="",y=TeX("$日懸浮載輸砂量(ton/d)$")) + # 座標軸名稱
  theme_bw()+ #白底
  theme(legend.position = c(0.7,0.75))+ #圖例位置座標
  #theme(legend.text=element_blank())+ #隱藏圖例標題
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(prism.ticks.length=unit(.7,"lines"))+
  theme(axis.ticks.length=unit(.3,"lines"))+
  theme(text=element_text(size=30,color="black"))+  # 字體大小
  theme(axis.text = element_text(colour = "black"))
dev.off()

dis <- subset(all.table,Discharge>14&Discharge<15)
dis$Discharge <- paste0(dis$Discharge,"cms")
#dis$Suspended.Load <- dis$Suspended.Load/10000
setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站討論(5種推估方法)15cms比較.png"),width = 1000, height = 600, units = "px", pointsize = 12)
ggplot(dis,aes(x=group,y=Suspended.Load,fill=Discharge,label=Suspended.Load))+
  geom_bar(position = "dodge",stat="identity")+
  geom_text(aes(label=round(Suspended.Load)), position=position_dodge(width=1), size=5.5,vjust=-0.25)+
  #geom_label_repel(min.segment.length = 0, box.padding = 0.5)+
  labs(x="",y=TeX("$日懸浮載輸砂量( ton/d)$")) + # 座標軸名稱
  theme_bw()+ #白底
  theme(legend.position = c(0.7,0.75))+ #圖例位置座標
  #theme(legend.text=element_blank())+ #隱藏圖例標題
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(prism.ticks.length=unit(.7,"lines"))+
  theme(axis.ticks.length=unit(.3,"lines"))+
  theme(text=element_text(size=30,color="black"))+  # 字體大小
  theme(axis.text = element_text(colour = "black"))
dev.off()

dis <- subset(all.table,Discharge>35&Discharge<36)
dis$Discharge <- paste0(dis$Discharge,"cms")
#dis$Suspended.Load <- dis$Suspended.Load/10000
setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站討論(5種推估方法)35cms比較.png"),width = 1000, height = 600, units = "px", pointsize = 12)
ggplot(dis,aes(x=group,y=Suspended.Load,fill=Discharge,label=Suspended.Load))+
  geom_bar(position = "dodge",stat="identity")+
  geom_text(aes(label=round(Suspended.Load)), position=position_dodge(width=1), size=5.5,vjust=-0.25)+
  #geom_label_repel(min.segment.length = 0, box.padding = 0.5)+
  labs(x="",y=TeX("$日懸浮載輸砂量( ton/d)$")) + # 座標軸名稱
  theme_bw()+ #白底
  theme(legend.position = c(0.7,0.75))+ #圖例位置座標
  #theme(legend.text=element_blank())+ #隱藏圖例標題
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(prism.ticks.length=unit(.7,"lines"))+
  theme(axis.ticks.length=unit(.3,"lines"))+
  theme(text=element_text(size=30,color="black"))+  # 字體大小
  theme(axis.text = element_text(colour = "black"))
dev.off()

dis <- subset(all.table,Discharge>67&Discharge<68)
dis$Discharge <- paste0(dis$Discharge,"cms")
#dis$Suspended.Load <- dis$Suspended.Load/10000
setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站討論(5種推估方法)70cms比較.png"),width = 1000, height = 600, units = "px", pointsize = 12)
ggplot(dis,aes(x=group,y=Suspended.Load,fill=Discharge,label=Suspended.Load))+
  geom_bar(position = "dodge",stat="identity")+
  geom_text(aes(label=round(Suspended.Load)), position=position_dodge(width=1), size=5.5,vjust=-0.25)+
  #geom_label_repel(min.segment.length = 0, box.padding = 0.5)+
  labs(x="",y=TeX("$日懸浮載輸砂量( ton/d)$")) + # 座標軸名稱
  theme_bw()+ #白底
  theme(legend.position = c(0.7,0.75))+ #圖例位置座標
  #theme(legend.text=element_blank())+ #隱藏圖例標題
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(prism.ticks.length=unit(.7,"lines"))+
  theme(axis.ticks.length=unit(.3,"lines"))+
  theme(text=element_text(size=30,color="black"))+  # 字體大小
  theme(axis.text = element_text(colour = "black"))
dev.off()

dis <- subset(all.table,Discharge<98&Discharge>95)
dis$Discharge <- paste0(dis$Discharge,"cms")
setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站討論(5種推估方法)100cms比較.png"),width = 1000, height = 600, units = "px", pointsize = 12)
ggplot(dis,aes(x=group,y=Suspended.Load,fill=Discharge,label=Suspended.Load))+
  geom_bar(position = "dodge",stat="identity")+
  geom_text(aes(label=round(Suspended.Load)), position=position_dodge(width=1), size=5.5,vjust=-0.25)+
  #geom_label_repel(min.segment.length = 0, box.padding = 0.5)
  labs(x="",y=TeX("$日懸浮載輸砂量(ton/d)$")) + # 座標軸名稱
  theme_bw()+ #白底
  theme(legend.position = c(0.7,0.75))+ #圖例位置座標
  #theme(legend.text=element_blank())+ #隱藏圖例標題
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(prism.ticks.length=unit(.7,"lines"))+
  theme(axis.ticks.length=unit(.3,"lines"))+
  theme(text=element_text(size=30,color="black"))+  # 字體大小
  theme(axis.text = element_text(colour = "black"))
dev.off()
#
#
# 3. 觀測輸砂量as histogram；推估輸砂量as density
setwd(paste0("F:/R_output/",station,"/vinecopula")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站輸砂量觀測分布與推估分布.png"),width = 1000, height = 700, units = "px", pointsize = 12)
ggplot(data=all.table,aes(x=Suspended.Load))+
  geom_histogram(data=subset(all.table,group=="觀測資料"),aes(y=..density..,fill=group),binwidth = 100)+
  stat_density(data=subset(all.table,group!="觀測資料"),aes(color=group),geom="line",position="identity",size=1.2)+
  scale_x_continuous(guide = "prism_minor", #y軸副刻度
                     limits = c(0,10000),
                     minor_breaks = seq(0, 10000, 500))+
  scale_y_continuous(guide = "prism_minor", #y軸副刻度
                     limits = c(0,0.001),
                     minor_breaks = seq(0, 0.001,0.00005))+
  labs(y="Frequency",x=TeX("日懸浮載輸砂量 ($ton/d$)")) + # 座標軸名稱
  theme_bw()+ #白底
  scale_fill_manual(values=c("gray"))+
  scale_color_discrete(breaks =c("觀測資料","率定曲線", "方法 1","方法 2","方法 3","方法 4"))+
  theme(legend.position = c(0.8,0.7))+ #圖例位置座標
  theme(legend.title=element_blank())+ #隱藏圖例標題
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(prism.ticks.length=unit(.7,"lines"))+
  theme(axis.ticks.length=unit(1.2,"lines"))+
  theme(text=element_text(size=40,color="black"))+  # 字體大小
  theme(axis.text = element_text(colour = "black"))
dev.off()
