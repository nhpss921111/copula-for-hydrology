# 邊際分布：parametric method
# copula函數：IFM method
# 開始撰寫日期：2020/02/16
# 完成撰寫日期：2021/03/03
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
# ===========
# 家源橋CHIA-YUANG year <- c(1974:2009,2012:2019)
# 彰雲橋CHUNYUN BRIDGE year <- c(1987:2019)
# 內茅埔(NEI-MAO-PU)：year <- c(1972:2001,2003:2019)
# ===========
station <- c("CHUNYUN BRIDGE") # 測站名稱
station_ch <-c("彰雲橋")
group.number <- c(9) # 分組的組數
log.group.number <- c(9) # 分組的組數
set.seed(101)
perc.mis <- 0.3 # 多少%的資料當成NA
year <- c(1987:2019) # 年分
MD.input <- c(paste0(year,"QandQs.csv"))
output <- c(paste0(year,"imp.csv"))
ob.data <- c()
# ---- 誤差指標公式 ----
e <- function(actual, predicted) { # 相對誤差
  (actual - predicted)/actual
}
mse <- function(actual, predicted) { # 均方誤差
  mean((actual - predicted) ^ 2)
}
rmse <- function(actual, predicted) { # 均方根誤差
  sqrt(mean((actual - predicted) ^ 2))
}
nmse <- function(actual, predicted) { # 正歸化均方根誤差
  mean((actual - predicted) ^ 2)/mean((actual - mean(actual)) ^ 2)
}
mape <- function(actual, predicted) { # 平均絕對百分誤差
  100*mean(abs((actual - predicted)/actual))
}
# ----------------
# 建立margin參數估計表格(參數估計方法：mle)(候選分布：norm, lnorm, gumbel, weibull, gamma)
fitmargin.par <- matrix(nrow=12,ncol=20)
colnames(fitmargin.par) <- c("norm.Q.par1","norm.Q.par2","norm.Qs.par1","norm.Qs.par2",
                             "lnorm.Q.par1","lnorm.Q.par2","lnorm.Qs.par1","lnorm.Qs.par2",
                             "gumbel.Q.par1","gumbel.Q.par2","gumbel.Qs.par1","gumbel.Qs.par2",
                             "weibull.Q.par1","weibull.Q.par2","weibull.Qs.par1","weibull.Qs.par2",
                             "gamma.Q.par1","gamma.Q.par2","gamma.Qs.par1","gamma.Qs.par2")

# 建立copula參數估計表格(參數估計方法：itau, irho, mpl, ml)
fitcopula.par <- matrix(nrow=12,ncol=12)
colnames(fitcopula.par) <- c("gumbel.itau","gumbel.irho","gumbel.mpl","gumbel.ml",
                             "frank.itau","frank.irho","frank.mpl","frank.ml",
                             "clayton.itau","clayton.irho","clayton.mpl","clayton.ml")
# 新增最後一行放選擇的參數值(預設值為0)
fitcopula.par <- data.frame(fitcopula.par, copula.parameter=0)

# 建立p-value表格(參數估計方法：itau, irho, mpl, ml)
gof.pvalue <- matrix(nrow=12,ncol=12)
colnames(gof.pvalue) <- c("gumbel.itau","gumbel.irho","gumbel.mpl","gumbel.ml",
                          "frank.itau","frank.irho","frank.mpl","frank.ml",
                          "clayton.itau","clayton.irho","clayton.mpl","clayton.ml")
# 新增最後一行放選擇的聯結函數(預設值為0)
gof.pvalue <- data.frame(gof.pvalue, resultcopula=0)
# --------------
# 主要迴圈(以年份為底)
for( y in 1:length(year)){
  # 讀取有缺失的資料
  setwd(paste0("F:/R_reading/",station,"/missingdata"))
  data <- read.csv(file.path(getwd(),MD.input[y]))
  data <- data[,-1]
  ob.data <- rbind(ob.data,data)
}
ob.data <- subset(ob.data,ob.data$Discharge>0) # 把流量為0(無觀測資料)刪除
ob.data[ob.data==0] <- NA #將輸砂量0的資料當成NA
log.data <- cbind(ob.data[,1:3],log10(ob.data[,4:5]))
rm.ob.data <- ob.data[complete.cases(ob.data), ] # 移除原始觀測資料中全部NA
file <- paste("F:/R_output/",station,"/parametric&coimp/",
              year[1],"到", year[y],station_ch,"流量&輸砂量(100%觀測資料).csv", sep="") #存檔路徑
write.csv(rm.ob.data,file)
rm.log.data <- log.data[complete.cases(log.data), ] # 移除觀測資料取對數後全部NA
#ob.data <- subset(ob.data,ob.data[,4]>20 & ob.data[,5]>1000)

# ----------- 將同時有Q與Qs的資料分兩組 (80%資料總數建模，剩下20%當成驗證) -------------
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
file <- paste("F:/R_output/",station,"/parametric&coimp/",
              year[1],"到", year[y],station_ch,"率定模式(70%觀測資料).csv", sep="") #存檔路徑
write.csv(MD.rmNA,file)

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
  # group 8：98% ~  99%
  # group 9：99% ~ 100%
  
  rank.data <- cbind(ob.data,rank(ob.data$Discharge))
  persent <- (rank.data$`rank(ob.data$Discharge)`) / length(rank.data$Discharge)
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
  
  group.BC <- c(0,max(data.1$Discharge),max(data.2$Discharge),max(data.3$Discharge),
                max(data.4$Discharge),max(data.5$Discharge),max(data.6$Discharge),
                max(data.7$Discharge),max(data.8$Discharge),max(data.9$Discharge))
  persent.BC <- c(0,20,40,60,80,90,95,98,99,100)
}

# 分組的目的：讓CoImp計算用
if (log.group.number==9){
  #  -------------- 決定流量分組範圍 (來源：log.data) ---------------
  # group 1： 0% ~  20%
  # group 2：20% ~  40%
  # group 3：40% ~  60%
  # group 4：60% ~  80%
  # group 5：80% ~  90%
  # group 6：90% ~  95%
  # group 8：95% ~  98%
  # group 8：98% ~  99%
  # group 9：99% ~ 100%
  
  log.rank.data <- cbind(log.data,rank(log.data$Discharge))
  log.persent <- (log.rank.data$`rank(log.data$Discharge)`) / length(log.rank.data$Discharge)
  log.per.data <- cbind(log.rank.data,persent)
  
  log.data.1 <- data.frame(subset(log.per.data, persent<=0.2),group="group1")
  log.data.2 <- data.frame(subset(log.per.data, persent>0.2 & persent<=0.4),group="group2")
  log.data.3 <- data.frame(subset(log.per.data, persent>0.4 & persent<=0.6),group="group3")
  log.data.4 <- data.frame(subset(log.per.data, persent>0.6 & persent<=0.8),group="group4")
  log.data.5 <- data.frame(subset(log.per.data, persent>0.8 & persent<=0.9),group="group5")
  log.data.6 <- data.frame(subset(log.per.data, persent>0.9 & persent<=0.96),group="group6")
  log.data.7 <- data.frame(subset(log.per.data, persent>0.96 & persent<=0.98),group="group7")
  log.data.8 <- data.frame(subset(log.per.data, persent>0.98 & persent<=0.99),group="group8")
  log.data.9 <- data.frame(subset(log.per.data, persent>0.99),group="group9")
  
  log.data.group <- rbind(log.data.1, log.data.2, log.data.3,
                      log.data.4, log.data.5, log.data.6,
                      log.data.7, log.data.8, log.data.9)
  
  log.group.BC <- c(0,max(log.data.1$Discharge),max(log.data.2$Discharge),max(log.data.3$Discharge),
                max(log.data.4$Discharge),max(log.data.5$Discharge),max(log.data.6$Discharge),
                max(log.data.7$Discharge),max(log.data.8$Discharge),max(log.data.9$Discharge))
  log.persent.BC <- c(0,20,40,60,80,90,95,98,99,100)
}

ggplot(data=MD.rmNA)+
  geom_point(aes(x=Discharge,y=Suspended.Load))+
  scale_color_discrete(name="年",labels=c(""))+
  labs(x="流量(cms)",y="輸砂量Qs (公噸)") + # 座標軸名稱
  #geom_mark_ellipse(expand = 0,aes(fill=group))+ # 橢圓形圈
  #geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=group))+ # 多邊形
  #geom_mark_hull(expand=0.01,aes(fill=group))+ # 凹凸多邊形
  ggtitle(paste0(year[1],"到",year[y],"年",station_ch,"測站建立模型")) +
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=30))  # 字體大小

# =========== 不分組 (建立 rating curve)===============
# 從"率定曲線初始值.csv"找到a,b的初始值再帶入
file <- paste("F:/R_output/",station,"/parametric&coimp/",
              year[1],"到", year[y],"年不分組率定曲線初始值查找.csv", sep="") #存檔路徑
write.csv(MD.rmNA,file)
# 從線性模型計算log10a,b 初始值
rating <- lm(Suspended.Load ~ Discharge,data=MD.log.rmNA)
summary(rating) # 觀察

# 將線性計算的係數當成初始值，再帶入非線性模型中
log10.a <- rating$coefficients[1] # 迴歸係數log10(a)
a.start <- 10^log10.a # 線性迴歸係數log10(a) -> 非線性迴歸係數a
b.start <- rating$coefficients[2] # 線性迴歸係數b -> 非線性迴歸係數b

rating <- nls(Suspended.Load ~ a*Discharge^b, algorithm="port",
              control = list(maxiter = 200,minFactor = 1/2^300,warnOnly = TRUE),
              start=list(a=a.start,b=b.start), data=MD.rmNA,trace=T) # 初始值要給好!!!
a <- environment(rating[["m"]][["resid"]])[["env"]][["a"]]
b <- environment(rating[["m"]][["resid"]])[["env"]][["b"]]
rating.par <- cbind(a,b)
file <- paste("F:/R_output/",station,"/parametric&coimp/",
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
setwd(paste0("F:/R_output/",station,"/parametric&coimp")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站率定曲線(不分組).png"),width = 1250, height = 700, units = "px", pointsize = 12)
Rating.all <- ggplot(data=rating.all)+
  geom_point(aes(x=Discharge,y=asNA,color="觀測資料")) +
  geom_line(aes(x=Discharge,y=ratingSSL,color="率定曲線")) +
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年",station_ch,"測站率定曲線(不分組)")) +
  theme(text=element_text(size=20))  # 字體大小
plot(Rating.all)
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
# 建立copula參數估計表格(參數估計方法：itau, irho, mpl, ml)
fitcopula.par <- matrix(nrow=1,ncol=12)
colnames(fitcopula.par) <- c("gumbel.itau","gumbel.irho","gumbel.mpl","gumbel.ml",
                             "frank.itau","frank.irho","frank.mpl","frank.ml",
                             "clayton.itau","clayton.irho","clayton.mpl","clayton.ml")
# 新增最後一行放選擇的參數值(預設值為0)
fitcopula.par <- data.frame(fitcopula.par, copula.parameter=0)

# 建立p-value表格(參數估計方法：itau, irho, mpl, ml)
gof.pvalue <- matrix(nrow=1,ncol=12)
colnames(gof.pvalue) <- c("gumbel.itau","gumbel.irho","gumbel.mpl","gumbel.ml",
                          "frank.itau","frank.irho","frank.mpl","frank.ml",
                          "clayton.itau","clayton.irho","clayton.mpl","clayton.ml")
# 新增最後一行放選擇的聯結函數(預設值為0)
gof.pvalue <- data.frame(gof.pvalue, resultcopula=0)
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

# 儲存 margins 相關資料
file <- paste("F:/R_output/",station,"/parametric&coimp/",
              year[1],"到", year[y],station_ch," margins parameter.csv", sep="") #存檔路徑
write.csv(margin.par,file)
file <- paste("F:/R_output/",station,"/parametric&coimp/",
              year[1],"到", year[y],station_ch," margins ks_test.csv", sep="") #存檔路徑
write.csv(margin.ks,file)
file <- paste("F:/R_output/",station,"/parametric&coimp/",
              year[1],"到", year[y],station_ch," margins AIC.csv", sep="") #存檔路徑
write.csv(margin.aic,file)
#
# ================== Determine copula function ==========================
#
# ----- fitCopula (Inference For Margin estimation) ------ 
candidate.copula <- c("gumbel","frank","clayton")
#Q與Qs的機率
var_a <- pobs(MD.anal$Discharge)
var_b <- pobs(MD.anal$Suspended.Load)
data.probs <- cbind(var_a, var_b)
dist <- 1
for(dist in c(1:length(candidate.copula))){
  print(paste0(candidate.copula[dist],"copula"))
  copula.char <- c(paste0(candidate.copula[dist],"Copula"))
  copula.func <- get(copula.char)(2) # initial theta = 2
  # mvd2 <- mvdc(copula.func, margin.dist,
  #             paramMargins = list(list(margin.par[margin.num[1],1], margin.par[margin.num[1],2]),
  #                                 list(margin.par[margin.num[2],3], margin.par[margin.num[2],4])))
  # copula參數估計
  # loglikMvdc(c(1,1,1,1), MD.anal, mvd2)
  # fitMvdc(MD.anal,mvd2,start=)
  fit.tau <- fitCopula(copula.func, data.probs, method="itau")
  fit.rho <- fitCopula(copula.func, data.probs, method="irho")
  fit.mpl <- fitCopula(copula.func, data.probs, method="mpl",optim.method="BFGS")
  fit.ml <- fitCopula(copula.func, data.probs, method="ml",optim.method="BFGS")
  #儲存估計參數到表格中
  if(dist==1){ #gumbelcopula參數估計
    fitcopula.par[1] <- fit.tau@estimate
    fitcopula.par[2] <- fit.rho@estimate
    fitcopula.par[3] <- fit.mpl@estimate
    fitcopula.par[4] <- fit.ml@estimate
  }
  if(dist==2){ # frankcopula參數估計
    fitcopula.par[5] <- fit.tau@estimate
    fitcopula.par[6] <- fit.rho@estimate
    fitcopula.par[7] <- fit.mpl@estimate
    fitcopula.par[8] <- fit.ml@estimate
  }
  if(dist==3){ # claytoncopula參數估計
    fitcopula.par[9] <- fit.tau@estimate
    fitcopula.par[10] <- fit.rho@estimate
    fitcopula.par[11] <- fit.mpl@estimate
    fitcopula.par[12] <- fit.ml@estimate
  }
  #
  # Goodness of fit test 
  #
  #if (gof=="y"){
  #print(paste0(m,"月的適合度檢定"))
  gof.tau <- gofCopula(get(copula.char)(fit.tau@estimate, dim=2), pobs(MD.anal),N = 2000,
                       method = "Sn", estim.method = "itau", simulation = "pb",ties=TRUE,
                       optim.method = "BFGS") # BFGS：牛頓法
  gof.rho <- gofCopula(get(copula.char)(fit.rho@estimate, dim=2), pobs(MD.anal),N = 2000,
                       method = "Sn", estim.method = "irho", simulation = "pb",ties=TRUE,
                       optim.method = "BFGS")
  gof.mpl <- gofCopula(get(copula.char)(fit.mpl@estimate, dim=2), pobs(MD.anal),N = 2000,
                       method = "Sn",start=gof.tau$parameter, estim.method = "mpl", 
                       simulation = "pb",ties=TRUE,optim.method = "BFGS")
  gof.ml <- gofCopula(get(copula.char)(fit.ml@estimate, dim=2), pobs(MD.anal),N = 2000,
                      method = "Sn", estim.method = "ml", simulation = "pb",ties=TRUE,
                      optim.method = "BFGS")
  if(copula.char=="gumbelCopula"){
    gof.pvalue[1] <- gof.tau$p.value
    gof.pvalue[2] <- gof.rho$p.value
    gof.pvalue[3] <- gof.mpl$p.value
    gof.pvalue[4] <- gof.ml$p.value
  }
  if(copula.char=="frankCopula"){
    gof.pvalue[5] <- gof.tau$p.value
    gof.pvalue[6] <- gof.rho$p.value
    gof.pvalue[7] <- gof.mpl$p.value
    gof.pvalue[8] <- gof.ml$p.value
  }
  if(copula.char=="claytonCopula"){
    gof.pvalue[9] <- gof.tau$p.value
    gof.pvalue[10] <- gof.rho$p.value
    gof.pvalue[11] <- gof.mpl$p.value
    gof.pvalue[12] <- gof.ml$p.value
  }
  #}  
}
# 選擇pvalue最大的值，並找出是哪個聯結函數(每個聯結函數使用四種參數估計法)
gof.pvalue[13] <- floor(which.max(gof.pvalue)%/%(4+.1)+1) #每四個為一組
# 最佳copula function的參數
fitcopula.par[13] <- fitcopula.par[which.max(gof.pvalue[1:12])]

# 儲存copula function相關資料
file <- paste("F:/R_output/",station,"/parametric&coimp/",
              year[1],"到", year[y],station_ch," copula parameter.csv", sep="") #存檔路徑
write.csv(fitcopula.par,file)
file <- paste("F:/R_output/",station,"/parametric&coimp/",
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
setwd(paste0("F:/R_output/",station,"/parametric&coimp")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站流量之機率密度函數.png"),width = 1250, height = 700, units = "px", pointsize = 12)
ggplot()+
  geom_line(aes(x=MD.anal$Discharge,y=f_Q))+
  labs(x="流量Q (cms)",y="PDF函數值") + # 座標軸名稱
  ggtitle(paste0(station_ch,"測站流量之機率密度函數")) +
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=30))  # 字體大小
dev.off()
# f_L
setwd(paste0("F:/R_output/",station,"/parametric&coimp")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站輸砂量之機率密度函數.png"),width = 1250, height = 700, units = "px", pointsize = 12)
ggplot()+
  geom_line(aes(x=MD.anal$Suspended.Load,y=f_L))+
  labs(x="輸砂量Qs (公噸)",y="PDF函數值") + # 座標軸名稱
  xlim(0,2000) +
  ggtitle(paste0(station_ch,"測站輸砂量之機率密度函數")) +
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=30))  # 字體大小
dev.off()
# F_Q
setwd(paste0("F:/R_output/",station,"/parametric&coimp")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站流量之累積分布函數.png"),width = 1250, height = 700, units = "px", pointsize = 12)
ggplot()+
  geom_line(aes(x=MD.anal$Discharge,y=F_Q))+
  labs(x="流量Q (cms)",y="CDF機率值") + # 座標軸名稱
  ggtitle(paste0(station_ch,"測站流量之累積分布函數")) +
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=30))  # 字體大小
dev.off()
# F_L
setwd(paste0("F:/R_output/",station,"/parametric&coimp")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站輸砂量之累積分布函數.png"),width = 1250, height = 700, units = "px", pointsize = 12)
ggplot()+
  geom_line(aes(x=MD.anal$Suspended.Load,y=F_L))+
  labs(x="輸砂量Qs (公噸)",y="CDF機率值") + # 座標軸名稱
  #xlim(0,2000) +
  ggtitle(paste0(station_ch,"測站輸砂量之累積分布函數")) +
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=30))  # 字體大小
dev.off()
# ---- 3. joint probability distribution ----
# joint CDF of L and Q：
# F_{L,Q}(l,q) = C(F_L(l),F_Q(q)), C：copula
# joint PDF of L and Q：
# f_{L,Q}(l,q) = c(F_L(l),F_Q(q))*f_L(l)*f_Q(q) , c：copula density
#
copula.func <- get(paste0(candidate.copula[as.integer(gof.pvalue[13])],"Copula"))(fitcopula.par$copula.parameter,dim=2)
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

setwd(paste0("F:/R_output/",station,"/parametric&coimp")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站 mvd 散佈圖.png"),width = 1250, height = 700, units = "px", pointsize = 12)
par(mfrow = c(1, 2))
scatterplot3d(v[,1],v[,2], mvd.pdf, color="red", main= paste0("joint PDF 散佈圖"), xlab = "Q", ylab="Qs", zlab="pMvdc",pch=".")
scatterplot3d(v[,1],v[,2], mvd.cdf, color="red", main=paste0("joint CDF 散佈圖"), xlab = "Q", ylab="Qs", zlab="pMvdc",pch=".")
dev.off()

setwd(paste0("F:/R_output/",station,"/parametric&coimp")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站joint PDF.png"),width = 1250, height = 700, units = "px", pointsize = 12)
par(mfrow = c(1, 2))
persp(mymvdc, dMvdc, xlim = c(0, max(MD.anal[,1])), ylim=c(0, max(MD.anal[,2])), main = paste0("joint PDF 透視圖"), xlab = "Q", ylab="Qs")
contour(mymvdc, dMvdc, xlim = c(0, max(MD.anal[,1])), ylim=c(0, max(MD.anal[,2])), main =  paste0("joint PDF 等高線圖"), xlab = "Q", ylab="Qs")
dev.off()

setwd(paste0("F:/R_output/",station,"/parametric&coimp")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站joint CDF.png"),width = 1250, height = 700, units = "px", pointsize = 12)
par(mfrow = c(1, 2))
persp(mymvdc, pMvdc, xlim = c(0, max(MD.anal[,1])), ylim=c(0, max(MD.anal[,2])), main =  paste0("joint CDF 透視圖"),xlab = "Q", ylab="Qs")
contour(mymvdc, pMvdc, xlim = c(0, max(MD.anal[,1])), ylim=c(0, max(MD.anal[,2])), main =  paste0("joint CDF 等高線圖"), xlab = "Q", ylab="Qs")
dev.off()
# ---- 4. conditional probability distribution ----

## conditional pdf of L given q0
## conditional PDF of L given the observation discharge q_0：
## f_{L|q_0}(l) = f_L(l)*c(F_L(l),F_Q(q))



copula.pdf <- function(u,v,theta){
  exp(-((-log(u))^theta*(-log(v))^theta)^(1/theta)) * 
  ((-log(u))*(-log(v)))^(theta-1)/(u*v) *
  ((-log(u))^(theta)+(-log(v))^(theta))^(2/theta-2) *
  ((theta-1)*((-log(u))^(theta)+(-log(v))^(theta))^(-1/theta)+1)
}

givenQ <- c(600) # 設定條件流量大小：

small.Qs <- seq(from=0.01,to=10,by=0.01)
all.Qs <- append(small.Qs,c(11:1000000))
## caculate conditional PDF of L given the q0
F_Q <- get(paste0("p",margin.dist[1]))(givenQ, margin.par[margin.num[1],1], margin.par[margin.num[1],2])
F_L <- get(paste0("p",margin.dist[2]))(all.Qs, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
f_L <- get(paste0("d",margin.dist[2]))(all.Qs, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
con.pdf <-copula.pdf(F_Q,F_L,fitcopula.par$copula.parameter)*f_L

# plot conditional PDF of L given the q0
setwd(paste0("F:/R_output/",station,"/parametric&coimp")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站Q=",givenQ,"cms之conditional PDF.png"),width = 1250, height = 700, units = "px", pointsize = 12)
ggplot()+
  geom_line(aes(x=all.Qs,y=con.pdf)) +
  labs(x="輸砂量Qs (公噸)",y="PDF函數值") + # 座標軸名稱
  xlim(0,1000000) +
  ggtitle(paste0(station_ch,"測站Q=",givenQ,"cms之條件機率密度函數")) +
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=30))  # 字體大小
dev.off()

## caculate conditional CDF of L given the q0
F_Q <- get(paste0("p",margin.dist[1]))(givenQ, margin.par[margin.num[1],1], margin.par[margin.num[1],2])
F_L <- get(paste0("p",margin.dist[2]))(all.Qs, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
con.cdf <- cCopula(cbind(F_Q,F_L), copula = copula.func,indices = 1:dim(copula.func), inverse = F)

# plot conditional CDF of L given the q0
setwd(paste0("F:/R_output/",station,"/parametric&coimp")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站Q=",givenQ,"cms之conditional CDF.png"),width = 1250, height = 700, units = "px", pointsize = 12)
ggplot()+
  geom_line(aes(x=all.Qs,y=con.cdf[,2])) +
  labs(x="輸砂量Qs (公噸)",y="CDF累積機率值") + # 座標軸名稱
  ggtitle(paste0(station_ch,"測站Q=",givenQ,"cms之條件累積分布函數")) +
  xlim(0,1000000)+
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=30))  # 字體大小
dev.off()



#  =========== 輸砂量推估 ============
# 1. mode(中位數)
# 2. (Peng et al. 2020)
# 3. (Bezak et al. 2017)
# 4. Suspended load rating curve
# 5. (Di Lascio et al. 2015)原始尺度
# 6. (Di Lascio et al. 2015)對數尺度
# ====================================
MD.onlyNA <- MD[complete.cases(MD)==F, ]
MD.onlyNA.rmNASSL <- MD.onlyNA[,-5]
ori.data.vad <- left_join(MD.onlyNA.rmNASSL,rm.ob.data)
givenQ.table <- MD.onlyNA$Discharge
estimate.SSL <- c()
for (q in 1:length(givenQ.table)){
  
  # 1. Suspended load rating curve 率定曲線
  # by rating curve coefficient a, b
  est.L1 <- a * (givenQ.table[q])^b
  
  # 2.(Peng et al. 2020) 從CDF隨機取一個點
  #set.seed(100)
  r2 <- runif(1) #隨機產生1個均勻分布
  F_Q <- get(paste0("p",margin.dist[1]))(givenQ.table[q], margin.par[margin.num[1],1], margin.par[margin.num[1],2])
  con.cdf2 <- cCopula(cbind(F_Q,r2), copula = copula.func,indices = 2, inverse = T)
  est.L2 <- get(paste0("q",margin.dist[2]))(con.cdf2,margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  
  # 3.(Bezak et al. 2017) 從CDF隨機取1000個點，再取中位數
  #set.seed(100)
  r3 <- runif(100) #隨機產生1000個均勻分布
  F_Q <- get(paste0("p",margin.dist[1]))(givenQ.table[q], margin.par[margin.num[1],1], margin.par[margin.num[1],2])
  con.cdf3 <- cCopula(cbind(F_Q,r3), copula = copula.func,indices = 2, inverse = T)
  est.L3 <- get(paste0("q",margin.dist[2]))(con.cdf3,margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  med.est.L3 <- median(est.L3) # 取中位數

  # 4. mode(眾數)
  F_Q <- get(paste0("p",margin.dist[1]))(givenQ.table[q], margin.par[margin.num[1],1], margin.par[margin.num[1],2])
  F_L <- get(paste0("p",margin.dist[2]))(MD.anal$Suspended.Load, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  f_L <- get(paste0("d",margin.dist[2]))(MD.anal$Suspended.Load, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  con.pdf <-copula.pdf(F_Q,F_L,fitcopula.par$copula.parameter)*f_L
  max(con.pdf) # 最大PDF值
  est.L4 <- MD.anal$Suspended.Load[which.max(con.pdf)] # PDF最大值所推估的輸砂量

  
  each.est <- cbind(est.L1,est.L2,med.est.L3,est.L4)
  estimate.SSL <- rbind(estimate.SSL,each.est)
  print(paste0("完成",q,"次，還剩",length(givenQ.table)-q,"次"))
}

# 合併前4個方法的推估值
colnames(estimate.SSL) <- c("est.SSL1","est.SSL2","est.SSL3","est.SSL4")
#subset(ori.data.vad,Discharge==givenQ.table)

SSL.result <- cbind(ori.data.vad,estimate.SSL)

# 5.(Di Lascio et al. 2015) 正常尺度

# ---- CoImp補遺 分九組 (source：MD) ----

imp.table <- c() # 輸砂量總表
coimp.copula.parameter.table <- c() # 各組聯結函數及其參數
loss.table <- c() # 各組損失函數比較

# ---- 分組迴圈  ----
# 將一開始移除的 % 加回來

for (g in 1:(length(group.BC)-1)){
  #maxBC1 <- max(data.1$Discharge)
  print(paste0("第",g,"組開始補遺"))
  MD.bygroup <- as.matrix(subset(MD[,4:5],Discharge>group.BC[g] & Discharge<=group.BC[g+1])) #提取Q與Qs出來，並限制流量範圍
  colnames(MD.bygroup) <- c("Discharge","Suspended.Load")
  rm.MD.bygroup <- MD.bygroup[complete.cases(MD.bygroup), ] # 移除全部NA 
  #upper <- c(max(rm.MD.bygroup[,1]),max(rm.MD.bygroup[,2]))
  #lower <- c(min(rm.MD.bygroup[,1]),min(rm.MD.bygroup[,2]))
  n.marg <- 2 # 兩個變數(Q、Qs)
  imp <- CoImp(MD.bygroup, n.marg=n.marg, smoothing = c(0.7,0.7),
               plot=T, q.lo=c(0.01,0.01), q.up=c(0.99,0.99),
               model=list(gumbelCopula(),frankCopula(),claytonCopula())) # 補遺計算
  # 邊際分布
  #setwd(paste0("F:/R_output/",station,"/parametric&coimp")) # 請修改儲存路徑：
  #png(paste0(year[1],"到",year[y],"年",station_ch,"測站group",g,"邊際分布.png"),width = 1250, height = 700, units = "px", pointsize = 12)
  plot(imp)
  #dev.off()
  
  coimp.copula.func <- imp@Estimated.Model.Imp[["model"]]
  coimp.copula.para <- imp@Estimated.Model.Imp[["parameter"]]
  coimp.copula.list <- cbind(coimp.copula.func, coimp.copula.para)
  rownames(coimp.copula.list) <- paste0("group",g)
  coimp.copula.parameter.table <- rbind(coimp.copula.parameter.table,
                                        coimp.copula.list) # 傳承
  
  imp.group <- cbind(subset(rm.ob.data,Discharge>group.BC[g] & Discharge<=group.BC[g+1]),
                     subset(MD[,4:5],Discharge>group.BC[g] & Discharge<=group.BC[g+1])[,2],
                     imp@Imputed.data.matrix[,2]) # 提取補遺值，並合併成結果表格
  
  colnames(imp.group) <- c("Year","Month","Day","Discharge",
                           "Suspended.Load","asNA","Imp.SSL") # 為框架的行命名
  
  imp.table <- rbind(imp.table,imp.group) # 傳承
  
  # # 計算損失函數
  # imp_vs_rating.group <- imp.group[(is.na(imp.group$asNA)),] #保留 na
  # 
  # # rmse (root mean square error)
  # rmse.imp <- rmse(imp_vs_rating.group$Suspended.Load, imp_vs_rating.group$Imp.SSL) #原本移除的觀測值與補遺值
  # rmse.rating <- rmse(imp_vs_rating.group$Suspended.Load, imp_vs_rating.group$ratingSSL_group) #原本移除的觀測值與率定值
  # 
  # # mse (mean absolute error)
  # mse.imp <- mse(imp_vs_rating.group$Suspended.Load, imp_vs_rating.group$Imp.SSL) #原本移除的觀測值與補遺值
  # mse.rating <- mse(imp_vs_rating.group$Suspended.Load, imp_vs_rating.group$ratingSSL_group) #原本移除的觀測值與率定值
  # 
  # # mape (mean absolute persentage error)
  # mape.imp <- mape(imp_vs_rating.group$Suspended.Load, imp_vs_rating.group$Imp.SSL) #原本移除的觀測值與補遺值
  # mape.rating <- mape(imp_vs_rating.group$Suspended.Load, imp_vs_rating.group$ratingSSL_group) #原本移除的觀測值與率定值
  # 
  # # loss function
  # loss <- cbind(mse.imp, mse.rating,mape.imp, mape.rating,rmse.imp, rmse.rating)
  # rownames(loss) <- paste0("group",g)
  # loss.table <- rbind(loss.table,loss)
  # 補遺資料出圖
  setwd(paste0("F:/R_output/",station,"/parametric&coimp")) # 請修改儲存路徑：
  png(paste0(year[1],"到",year[y],"年",station_ch,"測站補遺驗證group",g,".png"),width = 1250, height = 700, units = "px", pointsize = 12)
  Imp.group <- ggplot(data=imp.group)+
    geom_point(aes(x=Discharge,y=Imp.SSL,color="補遺值")) +
    geom_point(aes(x=Discharge,y=Suspended.Load,color="觀測值")) +
    scale_color_discrete(name="圖例") + #圖例名稱
    ggtitle(paste0(year[1],"到",year[y],"年",station_ch,"測站補遺驗證",persent.BC[g],"% ~ ",persent.BC[g+1],"%")) +
    theme(text=element_text(size=20))  # 字體大小
  plot(Imp.group)
  dev.off()
}
print("CoImp 分組計算完成")



# 6.(Di Lascio et al. 2015) 對數尺度

# ---- CoImp補遺 分九組 (source：MD) ----

imp.log.table <- c() # 輸砂量總表
coimp.log.copula.parameter.table <- c() # 各組聯結函數及其參數
loss.table <- c() # 各組損失函數比較

# ---- 分組迴圈  ----
# 將一開始移除的 % 加回來

for (g in 1:(length(log.group.BC)-1)){
  #maxBC1 <- max(data.1$Discharge)
  print(paste0("第",g,"組開始補遺"))
  MD.log.bygroup <- as.matrix(subset(MD.log[,4:5],Discharge>log.group.BC[g] & Discharge<=log.group.BC[g+1])) #提取Q與Qs出來，並限制流量範圍
  colnames(MD.log.bygroup) <- c("Discharge","Suspended.Load")
  rm.MD.log.bygroup <- MD.log.bygroup[complete.cases(MD.log.bygroup), ] # 移除全部NA 
  #upper <- c(max(rm.MD.bygroup[,1]),max(rm.MD.bygroup[,2]))
  #lower <- c(min(rm.MD.bygroup[,1]),min(rm.MD.bygroup[,2]))
  n.marg <- 2 # 兩個變數(Q、Qs)
  imp.log <- CoImp(MD.log.bygroup, n.marg=n.marg, smoothing = c(0.7,0.7),
               plot=T, q.lo=c(0.01,0.01), q.up=c(0.99,0.99),
               model=list(gumbelCopula(),frankCopula(),claytonCopula())) # 補遺計算
  # 邊際分布
  #setwd(paste0("F:/R_output/",station,"/parametric&coimp")) # 請修改儲存路徑：
  #png(paste0(year[1],"到",year[y],"年",station_ch,"測站group",g,"邊際分布.png"),width = 1250, height = 700, units = "px", pointsize = 12)
  plot(imp)
  #dev.off()
  
  coimp.copula.func <- imp.log@Estimated.Model.Imp[["model"]]
  coimp.copula.para <- imp.log@Estimated.Model.Imp[["parameter"]]
  coimp.copula.list <- cbind(coimp.copula.func, coimp.copula.para)
  rownames(coimp.copula.list) <- paste0("group",g)
  coimp.log.copula.parameter.table <- rbind(coimp.log.copula.parameter.table,
                                        coimp.copula.list) # 傳承
  
  imp.log.group <- cbind(subset(rm.log.data,Discharge>log.group.BC[g] & Discharge<=log.group.BC[g+1]),
                     subset(MD.log[,4:5],Discharge>log.group.BC[g] & Discharge<=log.group.BC[g+1])[,2],
                     imp.log@Imputed.data.matrix[,2]) # 提取補遺值，並合併成結果表格
  
  colnames(imp.log.group) <- c("Year","Month","Day","Discharge",
                           "Suspended.Load","asNA","Imp.log.SSL") # 為框架的行命名
  
  imp.log.table <- rbind(imp.log.table,imp.log.group) # 傳承
  
  # 補遺資料出圖
  setwd(paste0("F:/R_output/",station,"/parametric&coimp")) # 請修改儲存路徑：
  png(paste0(year[1],"到",year[y],"年",station_ch,"測站補遺驗證(對數)group",g,".png"),width = 1250, height = 700, units = "px", pointsize = 12)
  Imp.log.group <- ggplot(data=imp.log.group)+
    geom_point(aes(x=Discharge,y=Imp.log.SSL,color="補遺值")) +
    geom_point(aes(x=Discharge,y=Suspended.Load,color="觀測值")) +
    scale_color_discrete(name="圖例") + #圖例名稱
    ggtitle(paste0(year[1],"到",year[y],"年",station_ch,"測站補遺驗證",persent.BC[g],"% ~ ",persent.BC[g+1],"%")) +
    theme(text=element_text(size=20))  # 字體大小
  plot(Imp.log.group)
  dev.off()
}
print("CoImp 分組計算完成")
# 將資料返還成原始尺度(10^)
imp.log.final.table <- cbind(imp.log.table[,1:3],10^imp.log.table[,4:7])


# 合併5種推估輸砂量方法的值
imp.log.final.onlyNA <- imp.log.final.table[complete.cases(imp.log.final.table)==F,]
imp.log.final.onlyNA <- imp.log.final.onlyNA[,-4:-6] # 只留日期和補遺的輸砂量
validation.table1 <- left_join(SSL.result,imp.table)
validation.table <- left_join(validation.table1,imp.log.final.onlyNA)
validation.table <- validation.table[,-10]
validation.table[is.na(validation.table)] <-0 
colnames(validation.table) <- c("Year","Month","Day","Discharge",
                            "Suspended.Load","est.SSL1","est.SSL2","est.SSL3","est.SSL4","est.SSL5","est.SSL6")

file <- paste("F:/R_output/",station,"/parametric&coimp/",
              year[1],"到", year[y],station_ch,"驗證總表(30%觀測資料).csv", sep="") #存檔路徑
write.csv(validation.table,file)


# 1. mse
mse1 <- mse(validation.table$Suspended.Load, validation.table$est.SSL1)
mse2 <- mse(validation.table$Suspended.Load, validation.table$est.SSL2)
mse3 <- mse(validation.table$Suspended.Load, validation.table$est.SSL3)
mse4 <- mse(validation.table$Suspended.Load, validation.table$est.SSL4)
mse5 <- mse(validation.table$Suspended.Load, validation.table$est.SSL5)
mse6 <- mse(validation.table$Suspended.Load, validation.table$est.SSL6)
mse.table <- cbind(mse1,mse2,mse3,mse4,mse5,mse6)
# 2. rmse
rmse1 <- rmse(validation.table$Suspended.Load, validation.table$est.SSL1)
rmse2 <- rmse(validation.table$Suspended.Load, validation.table$est.SSL2)
rmse3 <- rmse(validation.table$Suspended.Load, validation.table$est.SSL3)
rmse4 <- rmse(validation.table$Suspended.Load, validation.table$est.SSL4)
rmse5 <- rmse(validation.table$Suspended.Load, validation.table$est.SSL5)
rmse6 <- rmse(validation.table$Suspended.Load, validation.table$est.SSL6)
rmse.table <- cbind(rmse1,rmse2,rmse3,rmse4,rmse5,rmse6)
# 3. nmse
nmse1 <- nmse(validation.table$Suspended.Load, validation.table$est.SSL1)
nmse2 <- nmse(validation.table$Suspended.Load, validation.table$est.SSL2)
nmse3 <- nmse(validation.table$Suspended.Load, validation.table$est.SSL3)
nmse4 <- nmse(validation.table$Suspended.Load, validation.table$est.SSL4)
nmse5 <- nmse(validation.table$Suspended.Load, validation.table$est.SSL5)
nmse6 <- nmse(validation.table$Suspended.Load, validation.table$est.SSL6)
nmse.table <- cbind(nmse1,nmse2,nmse3,nmse4,nmse5,nmse6)
# 4. mape
mape1 <- mape(validation.table$Suspended.Load, validation.table$est.SSL1)
mape2 <- mape(validation.table$Suspended.Load, validation.table$est.SSL2)
mape3 <- mape(validation.table$Suspended.Load, validation.table$est.SSL3)
mape4 <- mape(validation.table$Suspended.Load, validation.table$est.SSL4)
mape5 <- mape(validation.table$Suspended.Load, validation.table$est.SSL5)
mape6 <- mape(validation.table$Suspended.Load, validation.table$est.SSL6)
mape.table <- cbind(mape1,mape2,mape3,mape4,mape5,mape6)
error.table <- rbind(mse.table,rmse.table,nmse.table,mape.table)
colnames(error.table) <- c("est.SSL1","est.SSL2","est.SSL3","est.SSL4","est.SSL5","est.SSL6")
rownames(error.table) <- c("mse","rmse","nmse","mape")
file <- paste("F:/R_output/",station,"/parametric&coimp/",
              year[1],"到", year[y],station_ch,"驗證_誤差指標(30%觀測資料).csv", sep="") #存檔路徑
write.csv(error.table,file)

setwd(paste0("F:/R_output/",station,"/parametric&coimp")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站驗證(6種推估方法).png"),width = 1250, height = 700, units = "px", pointsize = 12)
ggplot(validation.table)+
  geom_point(aes(x=Discharge,y=Suspended.Load,color="真實值"),size=5)+
  geom_point(aes(x=Discharge,y=est.SSL1,color="率定曲線"),size=3)+
  geom_point(aes(x=Discharge,y=est.SSL2,color="CDF取1點"),size=3)+
  geom_point(aes(x=Discharge,y=est.SSL3,color="CDF取100點之中位數"),size=3)+
  geom_point(aes(x=Discharge,y=est.SSL4,color="PDF眾數"),size=3)+
  geom_point(aes(x=Discharge,y=est.SSL5,color="HitorMiss"),size=3)+
  geom_point(aes(x=Discharge,y=est.SSL6,color="log(HitorMiss)"),size=3)+
  # xlim(0,200)+
  # ylim(0,10000)+
  labs(x="流量Q(cms)",y="輸砂量Qs (公噸)") + # 座標軸名稱
  scale_color_discrete(name="圖例")+  #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年",station_ch,"測站驗證(30%觀測資料)"))+
  theme(text=element_text(size=30))  # 字體大小
dev.off()

a1 <- subset(validation.table,Discharge==5)
a2 <- subset(validation.table,Discharge==10)
a3 <- subset(validation.table,Discharge==15)
a4 <- subset(validation.table,Discharge==20)
a5 <- subset(validation.table,Discharge==30)
check.table <- rbind(a1,a2,a3,a4,a5)
file <- paste("F:/R_output/",station,"/parametric&coimp/",
              year[1],"到", year[y],station_ch,"驗證特定流量(30%觀測資料).csv", sep="") #存檔路徑
write.csv(check.table,file)

