# copula
# 開始日期：2020/07/08
# 完成日期：2020/07/13
# By 連育成
#
rm(list=ls())
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
library(VineCopula)
library(copula)
# 輸入邊際分布
margin.dist <-c("lnorm","lnorm") 
# Read data from excel flie
month<- c(1) # 請輸入月分：
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
  # 參數表格
  par.table <- matrix(nrow=length(candidate),ncol=2*2)
  rownames(par.table) <- c(candidate)
  colnames(par.table) <- c("par1","par2","par1","par2")
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
      # parameter 擺放設定
      if(i==1){
        par.table[dist,1] <- as.numeric(par1) 
        par.table[dist,2] <- as.numeric(par2)}
      if(i==2){
        par.table[dist,3] <- as.numeric(par1) 
        par.table[dist,4] <- as.numeric(par2)}
    }
  }
  var_a <- pobs(Q)
  var_b <- pobs(S)
  data.probs <- cbind(var_a, var_b)
  plot(var_a,var_b,col="red")
  #selectedCopula <- BiCopSelect(var_a, var_b, familyset = NA)
  # 建立copula
  g3 <- gumbelCopula(1,dim=2)
  f3 <- frankCopula(1,dim=2)
  c3 <- claytonCopula(1,dim=2)
  # 建立mvdc
  gMvd2 <- mvdc(g3,c("lnorm","lnorm"),
                param =list(list(meanlog=par.table[2,1], sdlog=par.table[2,2]), 
                            list(meanlog=par.table[2,3], sdlog=par.table[2,4])))
  fMvd2 <- mvdc(f3,c("lnorm","lnorm"),
                param =list(list(meanlog=par.table[2,1], sdlog=par.table[2,2]), 
                            list(meanlog=par.table[2,3], sdlog=par.table[2,4])))
  cMvd2 <- mvdc(c3,c("lnorm","lnorm"),
                param =list(list(meanlog=par.table[2,1], sdlog=par.table[2,2]), 
                            list(meanlog=par.table[2,3], sdlog=par.table[2,4])))
  # fit mvdc
  #loglikMvdc(c(2, 2, 2, 2, 3), variable, gMvd2)
  mm <- apply(variable, 2, mean)
  vv <- apply(variable, 2, var)
  b1.0 <- c(mm[1]^2/vv[1], vv[1]/mm[1])
  b2.0 <- c(mm[2]^2/vv[2], vv[2]/mm[2])
  a.0 <- sin(cor(variable[, 1], variable[, 2], method = "kendall") * pi/2)
  start <- c(b1.0, b2.0, a.0)
  fit <- fitMvdc(variable, gMvd2)# c(par.table[2,1],par.table[2,2],
                                            #par.table[2,3],par.table[2,4], a.0))
                 #optim.control = list(trace = TRUE, maxit = 1000))
  # mvdc(archmCopula(family="clayton",param=2),
  #     margins=c("lnorm","lnorm"),
  #    paramMargins=list(list(meanlog=par.table[2,1], sdlog=par.table[2,2]), 
  #                     list(meanlog=par.table[2,3], sdlog=par.table[2,4])))
}
