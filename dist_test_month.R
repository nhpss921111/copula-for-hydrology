# 開始撰寫日期：2020/06/13
# 完成撰寫日期：2020/07/06
# 統計檢定
# 以"蘭陽溪-家源橋"為例
# 候選分布：norm,lnorm,gumbel,weibull,gamma
# 用最大概似法
# K-S適合度檢定
# 等機率卡方適合度檢定
# AIC準則選最佳分布
# gamma的初始值要調整!!!
# 結果請參考：par.table & ks.table & chi.table & aic.table
# By連育成 

rm(list = ls())
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
library(gumbel)
library(ggplot2) #繪圖用
#
# Read data from excel flie
month<- c(1:12) # 請輸入月分：
input <- c(paste0(month,"month.csv"))
for (m in month){
  setwd("E:/R_reading/CHIA-YUANG")
  data <- read.csv(file.path(getwd(),input[m]),header = T) # 請輸入月分：
  data <- data[,-1]
  attach(data)
  #
  # ----------------------------- Remove missing value ------------------------------
  #
  #summary(data)
  #describe(data)
  # 計算missing value個數
  sum(is.na(data$Discharge)) #data have missing value
  sum(is.na(data$Suspended.Load)) #data have missing value
  # 移除missing value
  complete.cases(data)
  rm.data <- data[complete.cases(data), ]
  # 重新定義變數
  Q <- rm.data$Discharge
  S <- rm.data$Suspended.Load
  #
  # ------------------------------------------------------------------------------------
  #
  print("參數估計")
  candidate <- c("norm","lnorm","gumbel","weibull","gamma")
  dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b)) # a:location, b:scale
  pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
  qgumbel <- function(p, a, b) a-b*log(-log(p))
  #
  variable <- cbind(Q,S)
  #
  # par.table
  par.table <- matrix(nrow=length(candidate),ncol=2*2)
  rownames(par.table) <- c(candidate)
  colnames(par.table) <- c("par1","par2","par1","par2")
  # ks.table 放置p-value
  ks.table <- matrix(nrow=length(candidate)+1,ncol=dim(variable)[2]) # +1 是為了要放最佳分布的欄位
  rownames(ks.table) <- c(candidate,"good dist")
  colnames(ks.table) <- c(colnames(variable))
  # chi.table 放置p-value
  chi.table <- matrix(nrow=length(candidate)+1,ncol=dim(variable)[2]) # +1 是為了要放最佳分布的欄位  
  rownames(chi.table) <- c(candidate,"good dist")
  colnames(chi.table) <- c(colnames(variable))
  # aic.table  放置AIC value
  aic.table <- matrix(nrow=length(candidate)+1,ncol=dim(variable)[2]) # +1 是為了要放最小AIC的欄位
  rownames(aic.table) <- c(candidate,"good dist")
  colnames(aic.table) <- c(colnames(variable))
  #
  for(i in 1:dim(variable)[2]){
    var <- variable[,i]
    print(paste0("第",i,"個變數：",colnames(variable)[i])) #顯示第幾個及變數名稱
    # By Maximun Likelihood Estimate Method
    for(dist in c(1:length(candidate))){
      # -------------------------- parameter estimate -----------------------
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
  
      #表示要儲存成png的格式
      setwd("E:/R_output/CHIA-YUANG/result")
      png(paste(month[m],"month",colnames(variable)[i],candidate[dist],".png",sep=""), width=1200, height=800)
      cdfcomp(md,)
      dev.off() #最後要關掉輸出圖片裝置 
      # parameter 擺放設定
      if(i==1){
        par.table[dist,1] <- as.numeric(par1) 
        par.table[dist,2] <- as.numeric(par2)}
      if(i==2){
        par.table[dist,3] <- as.numeric(par1) 
        par.table[dist,4] <- as.numeric(par2)}
      #
      # ------------------------------- K-S test ----------------------
      #
      print("KS test")
      if(candidate[dist] == "norm"){xi.cdf <- get(dist.char[3])(var, par1,par2)}
      if(candidate[dist] == "lnorm"){xi.cdf <- get(dist.char[3])(var, meanlog=par1, sdlog=par2)}
      if(candidate[dist] == "gumbel"){xi.cdf <- get(dist.char[3])(var, par1, par2)}
      if(candidate[dist] == "weibull"){xi.cdf <- get(dist.char[3])(var, par1, par2)}
      if(candidate[dist] == "gamma"){xi.cdf <- get(dist.char[3])(var, par1, par2)}
      #
      result <- ks.test(var, dist.char[3], par1, par2)
      print(paste0(candidate[dist], "的KS檢定P-value: ", result$p.value))
      #
      # 將P-value整理成表格
      ks.table[dist,i] <- result$p.value
      #
      # ------------------------- chi square test ------------------------
      #
      print("equal prob. chi square test")
      k <- ceiling(1+3.3*log10(length(var))) # K：組數，無條件進位
      fv <- rep(1/k, time=k) # 分成K組，每組的機率密度 fx(x:Q,S)
      Fv <- cumsum(fv)
      oi <- c()
      for(j in c(1:k)){
        o <- get(dist.char[4])((Fv[j]),par1,par2)
        oi <- append(oi,o)
      }
      count <- cut(var,oi)
      Oi <- table(count)
      result <- chisq.test(as.numeric(Oi), p = fv[1:(k-1)], rescale.p = TRUE, simulate.p.value = TRUE)
      print(paste0(candidate[dist], "的卡方檢定P-value: ", result$p.value))
      chi.table[dist,i] <- result$p.value
      # --------------------------------- AIC -----------------------------
      print("AIC")
      if(candidate[dist] != "gumbel"){
        aic <- gofstat(md)
        aic.table[dist,i] <- aic$aic # 將AIC-value整理成表格
        print(paste0(candidate[dist], "AIC值: ", aic$aic))}
      if(candidate[dist] == "gumbel"){
        aic <- abic.gumbel(var,par1,par2)
        aic.table[dist,i] <- aic$AIC # 將AIC-value整理成表格
        print(paste0(candidate[dist], "AIC值: ", aic$AIC))}
    }
    # 每個延時P-value排序，變成數值再排序
    ks.choice <- rank(as.numeric(ks.table[(1:dist),i])) 
    ks.table[length(candidate)+1,i] <- candidate[which.max(ks.choice)] #最大的P-value對應的機率分布
    # 每個延時P-value排序，變成數值再排序
    chi.choice <- rank(as.numeric(chi.table[(1:dist),i])) 
    chi.table[length(candidate)+1,i] <- candidate[which.max(chi.choice)] #最大的P-value對應的機率分布  
    # 每個延時AIC-value排序，變成數值再排序
    aic.choice <- rank(as.numeric(aic.table[(1:dist),i])) 
    aic.table[length(candidate)+1,i] <- candidate[which.min(aic.choice)] #最小的AIC-value對應的機率分布
  }
  #
  # ------------------- export table ------------------------------------
  #
  file <- paste("E:/R_output/CHIA-YUANG/result/", month[m], "month_ks.csv", sep="")
  write.csv(ks.table,file)
  file <- paste("E:/R_output/CHIA-YUANG/result/", month[m], "month_chi.csv", sep="")
  write.csv(chi.table,file)
  file <- paste("E:/R_output/CHIA-YUANG/result/", month[m], "month_aic.csv", sep="")
  write.csv(aic.table,file)
}