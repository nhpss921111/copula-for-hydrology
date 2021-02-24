# 邊際分布：parametric method
# copula函數：IFM method
# 開始撰寫日期：2020/02/16
# 完成撰寫日期：2021/02/24
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
#library(aomisc)
#library(Metrics) # 評估指標
# ===========
# 家源橋CHIA-YUANG year <- c(1974:2009,2012:2019)
# 彰雲橋CHUNYUN BRIDGE year <- c(1987:2019)
# ===========
station <- c("CHIA-YUANG") # 測站名稱
station_ch <-c("家源橋")
group.number <- c(9) # 分組的組數
set.seed(100)
perc.mis    <- 0.3 # 多少%的資料當成NA
year <- c(1974:2009,2012:2019) # 年分
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
ob.data <- subset(ob.data,ob.data[,4]>0) # 把流量為0(無觀測資料)刪除
log.data <- cbind(ob.data[,1:3],log10(ob.data[,4:5]))
rm.ob.data <- ob.data[complete.cases(ob.data), ] # 移除原始觀測資料中全部NA
file <- paste("F:/R_output/",station,"/parametric method/",
              year[1],"到", year[y],station_ch,"流量&輸砂量(100%觀測資料).csv", sep="") #存檔路徑
write.csv(rm.ob.data,file)
rm.log.data <- log.data[complete.cases(log.data), ] # 移除觀測資料取對數後全部NA
#ob.data <- subset(ob.data,ob.data[,4]>20 & ob.data[,5]>1000)

x.samp <- as.matrix(rm.ob.data[,4:5])
miss.row    <- sample(1:length(rm.ob.data$Discharge), perc.mis*length(rm.ob.data$Discharge), replace=FALSE)
miss.col    <- rep(2,perc.mis*length(rm.ob.data$Discharge))
miss        <- cbind(miss.row,miss.col) # NA的欄位座標
samp.miss <- replace(x.samp,miss,NA) # 將 ?% 觀測資料轉換成NA
MD <- cbind(rm.ob.data[,1:3],samp.miss) 
MD.withNA <- MD[,4:5]
MD.rmNA <- MD[complete.cases(MD), ] # 移除全部NA (剩餘資料)
file <- paste("F:/R_output/",station,"/parametric method/",
              year[1],"到", year[y],station_ch,"率定模式(70%觀測資料).csv", sep="") #存檔路徑
write.csv(MD.rmNA,file)
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
file <- paste("F:/R_output/",station,"/parametric method/",
              year[1],"到", year[y],station_ch," margins parameter.csv", sep="") #存檔路徑
write.csv(margin.par,file)
file <- paste("F:/R_output/",station,"/parametric method/",
              year[1],"到", year[y],station_ch," margins ks_test.csv", sep="") #存檔路徑
write.csv(margin.ks,file)
file <- paste("F:/R_output/",station,"/parametric method/",
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
  fit.mpl <- fitCopula(copula.func, data.probs, method="mpl")
  fit.ml <- fitCopula(copula.func, data.probs, method="ml")
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
file <- paste("F:/R_output/",station,"/parametric method/",
              year[1],"到", year[y],station_ch," copula parameter.csv", sep="") #存檔路徑
write.csv(fitcopula.par,file)
file <- paste("F:/R_output/",station,"/parametric method/",
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
setwd(paste0("F:/R_output/",station,"/parametric method")) # 請修改儲存路徑：
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
setwd(paste0("F:/R_output/",station,"/parametric method")) # 請修改儲存路徑：
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
setwd(paste0("F:/R_output/",station,"/parametric method")) # 請修改儲存路徑：
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
setwd(paste0("F:/R_output/",station,"/parametric method")) # 請修改儲存路徑：
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

setwd(paste0("F:/R_output/",station,"/parametric method")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站 mvd 散佈圖.png"),width = 1250, height = 700, units = "px", pointsize = 12)
par(mfrow = c(1, 2))
scatterplot3d(v[,1],v[,2], mvd.pdf, color="red", main= paste0("joint PDF 散佈圖"), xlab = "Q", ylab="Qs", zlab="pMvdc",pch=".")
scatterplot3d(v[,1],v[,2], mvd.cdf, color="red", main=paste0("joint CDF 散佈圖"), xlab = "Q", ylab="Qs", zlab="pMvdc",pch=".")
dev.off()

setwd(paste0("F:/R_output/",station,"/parametric method")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站joint PDF.png"),width = 1250, height = 700, units = "px", pointsize = 12)
par(mfrow = c(1, 2))
persp(mymvdc, dMvdc, xlim = c(0, max(MD.anal[,1])), ylim=c(0, max(MD.anal[,2])), main = paste0("joint PDF 透視圖"), xlab = "Q", ylab="Qs")
contour(mymvdc, dMvdc, xlim = c(0, max(MD.anal[,1])), ylim=c(0, max(MD.anal[,2])), main =  paste0("joint PDF 等高線圖"), xlab = "Q", ylab="Qs")
dev.off()

setwd(paste0("F:/R_output/",station,"/parametric method")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站joint CDF.png"),width = 1250, height = 700, units = "px", pointsize = 12)
par(mfrow = c(1, 2))
persp(mymvdc, pMvdc, xlim = c(0, max(MD.anal[,1])), ylim=c(0, max(MD.anal[,2])), main =  paste0("joint CDF 透視圖"),xlab = "Q", ylab="Qs")
contour(mymvdc, pMvdc, xlim = c(0, max(MD.anal[,1])), ylim=c(0, max(MD.anal[,2])), main =  paste0("joint CDF 等高線圖"), xlab = "Q", ylab="Qs")
dev.off()
# ---- 4. conditional probability distribution ----

## conditional pdf of L given q0
## conditional PDF of L given the observation discharge q_0：
## f_{L|q_0}(l) = f_L(l)*c(F_L(l),F_Q(q))

givenQ <- 20 # 設定條件流量大小：
copula.pdf <- function(u,v,theta){
  exp(-((-log(u))^theta*(-log(v))^theta)^(1/theta)) * 
  ((-log(u))*(-log(v)))^(theta-1)/(u*v) *
  ((-log(u))^(theta)+(-log(v))^(theta))^(2/theta-2) *
  ((theta-1)*((-log(u))^(theta)+(-log(v))^(theta))^(-1/theta)+1)
}

## caculate conditional PDF of L given the q0
F_Q <- get(paste0("p",margin.dist[1]))(givenQ, margin.par[margin.num[1],1], margin.par[margin.num[1],2])
F_L <- get(paste0("p",margin.dist[2]))(MD.anal$Suspended.Load, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
f_L <- get(paste0("d",margin.dist[2]))(MD.anal$Suspended.Load, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
con.pdf <-copula.pdf(F_Q,F_L,fitcopula.par$copula.parameter)*f_L

# plot conditional PDF of L given the q0
setwd(paste0("F:/R_output/",station,"/parametric method")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站conditional PDF.png"),width = 1250, height = 700, units = "px", pointsize = 12)
ggplot()+
  geom_line(aes(x=MD.anal$Suspended.Load,y=con.pdf)) +
  labs(x="輸砂量Qs (公噸)",y="PDF函數值") + # 座標軸名稱
  #xlim(0,2000) +
  ggtitle(paste0(station_ch,"測站Q=",givenQ,"cms之條件機率密度函數")) +
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=30))  # 字體大小
dev.off()

## caculate conditional CDF of L given the q0
F_Q <- get(paste0("p",margin.dist[1]))(givenQ, margin.par[margin.num[1],1], margin.par[margin.num[1],2])
F_L <- get(paste0("p",margin.dist[2]))(MD.anal$Suspended.Load, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
con.cdf <- cCopula(cbind(F_Q,F_L), copula = copula.func,indices = 1:dim(copula.func), inverse = F)

# plot conditional CDF of L given the q0
setwd(paste0("F:/R_output/",station,"/parametric method")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站conditional CDF.png"),width = 1250, height = 700, units = "px", pointsize = 12)
ggplot()+
  geom_line(aes(x=MD.anal$Suspended.Load,y=con.cdf[,2])) +
  labs(x="輸砂量Qs (公噸)",y="CDF累積機率值") + # 座標軸名稱
  ggtitle(paste0(station_ch,"測站Q=",givenQ,"cms之條件累積分布函數")) +
  #xlim(0,1000)+
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=30))  # 字體大小
dev.off()

#  =========== 輸砂量推估 ============
# 1. mode(中位數)
# 2. (Peng et al. 2020)
# 3. (Bezak et al. 2017)
# ====================================
MD.onlyNA <- MD[complete.cases(MD)==F, ]
givenQ.table <- MD.onlyNA$Discharge
estimate.SSL <- c()
for (q in 1:length(givenQ.table)){
  # 1. mode(眾數)
  F_Q <- get(paste0("p",margin.dist[1]))(givenQ.table[q], margin.par[margin.num[1],1], margin.par[margin.num[1],2])
  F_L <- get(paste0("p",margin.dist[2]))(MD.anal$Suspended.Load, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  f_L <- get(paste0("d",margin.dist[2]))(MD.anal$Suspended.Load, margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  con.pdf <-copula.pdf(F_Q,F_L,fitcopula.par$copula.parameter)*f_L
  max(con.pdf) # 最大PDF值
  est.L1 <- MD.anal$Suspended.Load[which.max(con.pdf)] # PDF最大值所推估的輸砂量
  
  # 2.(Peng et al. 2020)
  set.seed(100)
  r2 <- runif(1) #隨機產生1個均勻分布
  F_Q <- get(paste0("p",margin.dist[1]))(givenQ.table[q], margin.par[margin.num[1],1], margin.par[margin.num[1],2])
  con.cdf2 <- cCopula(cbind(F_Q,r2), copula = copula.func,indices = 2, inverse = T)
  est.L2 <- get(paste0("q",margin.dist[2]))(con.cdf2,margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  
  # 3.(Bezak et al. 2017)
  set.seed(100)
  r3 <- runif(1000) #隨機產生1000個均勻分布
  F_Q <- get(paste0("p",margin.dist[1]))(givenQ.table[q], margin.par[margin.num[1],1], margin.par[margin.num[1],2])
  con.cdf3 <- cCopula(cbind(F_Q,r3), copula = copula.func,indices = 2, inverse = T)
  est.L3 <- get(paste0("q",margin.dist[2]))(con.cdf3,margin.par[margin.num[2],3], margin.par[margin.num[2],4])
  med.est.L3 <- median(est.L3) # 取中位數
  each.est <- cbind(est.L1,est.L2,med.est.L3)
  estimate.SSL <- rbind(estimate.SSL,each.est)
  
  print(paste0("第",q,"個完成"))
}
colnames(estimate.SSL) <- c("est.SSL1","est.SSL2","est.SSL3")
SSL.result <- cbind(MD.onlyNA[,1:4],estimate.SSL)
est.table <- left_join(rm.ob.data,SSL.result)
file <- paste("F:/R_output/",station,"/parametric method/",
              year[1],"到", year[y],station_ch,"驗證(30%觀測資料).csv", sep="") #存檔路徑
write.csv(est.table,file)

rm.est.table <- est.table[complete.cases(est.table),]
rmse1 <- rmse(rm.est.table$Suspended.Load, rm.est.table$est.SSL1)
rmse2 <- rmse(rm.est.table$Suspended.Load, rm.est.table$est.SSL2)
rmse3 <- rmse(rm.est.table$Suspended.Load, rm.est.table$est.SSL3)
rmse.table <- cbind(rmse1,rmse2,rmse3)
file <- paste("F:/R_output/",station,"/parametric method/",
              year[1],"到", year[y],station_ch,"驗證_誤差指標(30%觀測資料).csv", sep="") #存檔路徑
write.csv(rmse.table,file)

setwd(paste0("F:/R_output/",station,"/parametric method")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站驗證(3種推估方法).png"),width = 1250, height = 700, units = "px", pointsize = 12)
ggplot(rm.est.table)+
  geom_line(aes(x=Discharge,y=Suspended.Load,color="真實值"))+
  geom_line(aes(x=Discharge,y=est.SSL1,color="推估值1"))+
  geom_line(aes(x=Discharge,y=est.SSL2,color="推估值2"))+
  geom_line(aes(x=Discharge,y=est.SSL3,color="推估值3"))+
  scale_color_discrete(name="圖例")+  #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年",station_ch,"測站驗證(30%觀測資料)"))+
  theme(text=element_text(size=30))  # 字體大小
dev.off()

