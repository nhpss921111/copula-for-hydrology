# Missing data impution
# 開始撰寫日期：2020/11/12
# 完成撰寫日期：2020/11/12
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



month <- c(1:12)
MD.input <- c(paste0(month,"monthQandQs.csv"))
D.input <-  c(paste0(month,"month.csv"))
output <- c(paste0(month,"month_imp.csv"))
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

# 主要迴圈(以月分為底)
for( m in 1:length(month)){
  # 讀取有缺失的資料
  setwd("F:/R_reading/CHIA-YUANG/missingdata")
  MD <- read.csv(file.path(getwd(),MD.input[m]))
  MD <- MD[,-1]
  # 讀取有完整Q與Qs的資料
  setwd("F:/R_reading/CHIA-YUANG")
  D <- read.csv(file.path(getwd(),D.input[m])) 
  D <- D[,-1]
  # -------------------------------- Remove missing value --------------------------------
  #
  #summary(data)
  #describe(data)
  # 計算missing value個數
  sum(is.na(D$Discharge)) #data have missing value
  sum(is.na(D$Suspended.Load)) #data have missing value
  # 移除missing value
  complete.cases(D)
  rm.data <- D[complete.cases(D), ]
  # 重新定義變數
  Q <- rm.data$Discharge
  S <- rm.data$Suspended.Load
  #
  # -------------------------------------------------------------------------------------
  #
  print("邊際分布之參數估計")
  candidate <- c("norm","lnorm","gumbel","weibull","gamma")
  dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b)) # a:location, b:scale
  pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
  qgumbel <- function(p, a, b) a-b*log(-log(p))
  #
  variable <- cbind(Q,S)
  # 邊際分布參數表格
  margin.par <- matrix(nrow=length(candidate),ncol=2*2)
  rownames(margin.par) <- c(candidate)
  colnames(margin.par) <- c("par1","par2","par1","par2") 
  # margin.ks 放置p-value
  margin.ks <- matrix(nrow=length(candidate)+1,ncol=dim(variable)[2]) # +1 是為了要放最佳分布的欄位
  rownames(margin.ks) <- c(candidate,"good dist")
  colnames(margin.ks) <- c(colnames(variable))
  # margin.aic  放置AIC value
  margin.aic <- matrix(nrow=length(candidate)+1,ncol=dim(variable)[2]) # +1 是為了要放最小AIC的欄位
  rownames(margin.aic) <- c(candidate,"good dist")
  colnames(margin.aic) <- c(colnames(variable))
  # 邊際分布編號
  margin.num <- matrix(ncol=2)
  #
  for(i in 1:dim(variable)[2]){ # Q ,QS
    var <- variable[,i]
    print(paste0("第",i,"個變數：",colnames(variable)[i])) #顯示第幾個及變數名稱
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
  
  #
  # ---------------------- estimate copula parameter ---------------------
  candidate.copula <- c("gumbel","frank","clayton")

  for(dist in c(1:length(candidate.copula))){
    print(paste0(candidate.copula[dist],"copula"))
    copula.char <- c(paste0(candidate.copula[dist],"Copula"))
    copula.func <- get(copula.char)(2)
    #Q與Qs的機率
    var_a <- pobs(Q)
    var_b <- pobs(S)
    data.probs <- cbind(var_a, var_b)
    #建立Mvdc
    Mvd2 <- mvdc(copula.func,margin.dist,
                 param =list(list(margin.par[margin.num[1],1], margin.par[margin.num[1],2]),
                             list(margin.par[margin.num[2],3], margin.par[margin.num[2],4])))
    # copula參數估計
    fit.tau <- fitCopula(Mvd2@copula, data.probs, method="itau")
    fit.rho <- fitCopula(Mvd2@copula, data.probs, method="irho")
    fit.mpl <- fitCopula(Mvd2@copula, data.probs, method="mpl")
    fit.ml <- fitCopula(Mvd2@copula, data.probs, method="ml")
    #儲存估計參數到表格中
    if(dist==1){ #gumbelcopula參數估計
      fitcopula.par[m,1] <- fit.tau@estimate
      fitcopula.par[m,2] <- fit.rho@estimate
      fitcopula.par[m,3] <- fit.mpl@estimate
      fitcopula.par[m,4] <- fit.ml@estimate
    }
    if(dist==2){ # frankcopula參數估計
      fitcopula.par[m,5] <- fit.tau@estimate
      fitcopula.par[m,6] <- fit.rho@estimate
      fitcopula.par[m,7] <- fit.mpl@estimate
      fitcopula.par[m,8] <- fit.ml@estimate
    }
    if(dist==3){ # claytoncopula參數估計
      fitcopula.par[m,9] <- fit.tau@estimate
      fitcopula.par[m,10] <- fit.rho@estimate
      fitcopula.par[m,11] <- fit.mpl@estimate
      fitcopula.par[m,12] <- fit.ml@estimate
    }
    #
    # Goodness of fit test 
    #
    #if (gof=="y"){
      print(paste0(month[m],"月的適合度檢定"))
      gof.tau <- gofCopula(get(copula.char)(fit.tau@estimate, dim=2), pobs(variable),N = 2000,
                           method = "Sn", estim.method = "itau", simulation = "pb",ties=TRUE,
                           optim.method = "BFGS") # BFGS：牛頓法
      gof.rho <- gofCopula(get(copula.char)(fit.rho@estimate, dim=2), pobs(variable),N = 2000,
                           method = "Sn", estim.method = "irho", simulation = "pb",ties=TRUE,
                           optim.method = "BFGS")
      gof.mpl <- gofCopula(get(copula.char)(fit.mpl@estimate, dim=2), pobs(variable),N = 2000,
                           method = "Sn",start=gof.tau$parameter, estim.method = "mpl", 
                           simulation = "pb",ties=TRUE,optim.method = "BFGS")
      gof.ml <- gofCopula(get(copula.char)(fit.ml@estimate, dim=2), pobs(variable),N = 2000,
                          method = "Sn", estim.method = "ml", simulation = "pb",ties=TRUE,
                          optim.method = "BFGS")
      if(copula.char=="gumbelCopula"){
        gof.pvalue[m,1] <- gof.tau$p.value
        gof.pvalue[m,2] <- gof.rho$p.value
        gof.pvalue[m,3] <- gof.mpl$p.value
        gof.pvalue[m,4] <- gof.ml$p.value
      }
      if(copula.char=="frankCopula"){
        gof.pvalue[m,5] <- gof.tau$p.value
        gof.pvalue[m,6] <- gof.rho$p.value
        gof.pvalue[m,7] <- gof.mpl$p.value
        gof.pvalue[m,8] <- gof.ml$p.value
      }
      if(copula.char=="claytonCopula"){
        gof.pvalue[m,9] <- gof.tau$p.value
        gof.pvalue[m,10] <- gof.rho$p.value
        gof.pvalue[m,11] <- gof.mpl$p.value
        gof.pvalue[m,12] <- gof.ml$p.value
      }
    #}  
  }
  # 選擇pvalue最大的值，並找出是哪個聯結函數(每個聯結函數使用四種參數估計法)
  gof.pvalue[m,13] <- floor(which.max(gof.pvalue[m,])%/%(4+.1)+1) #每四個為一組
  # 最佳連結函數的參數
  fitcopula.par[m,13] <- fitcopula.par[m,which.max(gof.pvalue[m,1:12])]
  
  # choice copula function and its parameter
  mycopula <- get(paste0(candidate.copula[gof.pvalue[m,13]],"Copula"))(param = fitcopula.par[m,13])
  
  # imputation 
  n.marg <- 2
  samp.miss <- as.matrix(MD[,4:5])
  imp <- CoImp(samp.miss, n.marg=n.marg, smoothing = rep(0.6,n.marg), 
               plot=TRUE,type.data="continuous", 
               model=list(get(paste0(candidate.copula[gof.pvalue[m,13]],"Copula"))(param = fitcopula.par[m,13], dim=n.marg)))
  show(imp)
  plot(imp)
  imp@Imputed.data.matrix[,2]
  imp.data <- cbind(MD,imp@Imputed.data.matrix[,2])
  imp.data <- imp.data[,-5]
  colnames(imp.data) <- c("Year","Month","Day","Discharge","Imp.SSL")
  
  # 存檔
  file <- paste("F:/R_output/CHIA-YUANG/imputation/", output[m], sep="") #存檔路徑  
  write.csv(imp.data,file)
}
