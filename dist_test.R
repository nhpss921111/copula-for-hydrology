# 撰寫日期：2020/06/13
# 統計檢定
# 候選分布：norm,lnorm,gumbel,weibull,gamma,lgamma
# 用最大概似法
# K-S檢定
# AIC準則選最佳分布
# gamma和lgamma的初始值要調整!!!
# 結果請參考：par.table & dist.table & aic.table
# By連育成

rm(list = ls())
library(xlsx) #讀取excel
library(stats) #機率分布
library(actuar) #機率分布
library(dplyr) #資料整理
library(FAdist) # 水文上常用機率分布
library(vcd) # 估計gumbel參數
library(fitdistrplus) # 估計參數
library(EnvStats) # Environmental Statistics, Including US EPA Guidance
library(reliaR) # AIC of gumbel
library(lestat) #inverse CDF
library(goftest) #適合度檢定
library(goft) #適合度檢定
library(gumbel)
#
# Read data from excel flie
data <- read.xlsx(file="E:\\R_reading\\test_2.xlsx",
                  sheetIndex="DATA",startRow = 1,endRow = 960,
                  header = T,colIndex =2:3,
                  encoding = "UTF-8")
attach(data)
#
# ----------------------------- Remove missing value ------------------------------
#
summary(data)
describe(data)
# 計算missing value個數
sum(is.na(data$Q)) #data have missing value
sum(is.na(data$S)) #data have missing value
# 移除missing value
complete.cases(data)
rm.data <- data[complete.cases(data), ]
# 重新定義變數
Q <- rm.data$Q
S <- rm.data$S
#
# ------------------------------------------------------------------------------------
#
print("參數估計")
candidate <- c("norm","lnorm","gumbel","weibull","gamma","lgamma")
dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))
#
variable <- cbind(Q,S)
#
# par.table
par.table <- matrix(nrow=length(candidate),ncol=2*2)
rownames(par.table) <- c(candidate)
colnames(par.table) <- c("par1","par2","par1","par2")
# dist.table 放置p-value
dist.table <- matrix(nrow=length(candidate)+1,ncol=dim(variable)[2]) # +1 是為了要放最佳分布的欄位
rownames(dist.table) <- c(candidate,"good dist")
colnames(dist.table) <- c(colnames(variable))
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
      md <- fitdist(var, dist = dist.char[1])
      par1 <- md$estimate[1] #fitting參數1
      par2 <- md$estimate[2] #參數2
      print(c(par1, par2))}
    if(candidate[dist] == "lnorm"){
      md <- fitdist(var, dist = dist.char[1], start = list(meanlog=1, sdlog=1))
      par1 <- md$estimate[1] #fitting參數1
      par2 <- md$estimate[2] #參數2
      print(c(par1, par2))}
    if(candidate[dist] == "gumbel"){
      md <- eevd(var,method = "mle") 
      par1 <- md$parameters[1] #fitting參數1
      par2 <- md$parameters[2] #參數2
      print(c(par1, par2))}
    if(candidate[dist] == "weibull"){
      md <- fitdist(var, dist = dist.char[1])
      par1 <- md$estimate[1] #fitting參數1
      par2 <- md$estimate[2] #參數2
      print(c(par1, par2))}
    if(candidate[dist] == "gamma"){
      md <- fitdist(var, dist = dist.char[1],lower = 0, upper=Inf) 
      # start=list(shape=1,scale=1) # control=list(trace=1, REPORT=1)
      par1 <- md$estimate[1] #fitting參數1
      par2 <- md$estimate[2] #參數2
      print(c(par1, par2))}
    if(candidate[dist] == "lgamma"){
      if(i==1){
        md <- fitdist(var, dist = dist.char[1],start=list(shapelog=1,ratelog=1),
                      lower=c(0,0),upper=c(Inf,Inf)) # control=list(trace=1, REPORT=1)
        par1 <- md$estimate[1] #fitting參數1
        par2 <- md$estimate[2] #參數2
        print(c(par1, par2))}
      if(i==2){
        next # 若無法估計就跳過
        prefit(var,dist = dist.char[1],"mle",list(shapelog=0.1,ratelog=0.5),lower=0,upper=Inf)
        md <- fitdist(var, dist = dist.char[1],start=list(shapelog=0.1,ratelog=0.5),
                      lower=c(0,0),upper=c(Inf,Inf)) # control=list(trace=1, REPORT=1)
        par1 <- md$estimate[1] #fitting參數1
        par2 <- md$estimate[2] #參數2
        print(c(par1, par2))}
    }
    if(i==1){
      par.table[dist,1] <- as.numeric(par1) 
      par.table[dist,2] <- as.numeric(par2)}
    if(i==2){
      par.table[dist,3] <- as.numeric(par1) 
      par.table[dist,4] <- as.numeric(par2)}
    
    # ------------------------------- K-S test ----------------------
    print("KS test")
    if(candidate[dist] == "norm"){xi.cdf <- get(dist.char[3])(var, par1,par2)}
    if(candidate[dist] == "lnorm"){xi.cdf <- get(dist.char[3])(var, meanlog=par1, sdlog=par2)}
    if(candidate[dist] == "gumbel"){xi.cdf <- get(dist.char[3])(var, par1, par2)}
    if(candidate[dist] == "weibull"){xi.cdf <- get(dist.char[3])(var, par1, par2)}
    if(candidate[dist] == "gamma"){xi.cdf <- get(dist.char[3])(var, par1, par2)}
    if(candidate[dist] == "lgamma"){xi.cdf <- get(dist.char[3])(var, par1, par2)}
    
    result <- ks.test(var, dist.char[3], par1,par2)
    print(paste0(candidate[dist], "的KS檢定P-value: ", result$p.value))
    
    # 將P-value整理成表格
    dist.table[dist,i] <- result$p.value
    
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
  if(i==2 | candidate[dist] == "lgamma"){dist.table[6,2] <- -Inf; aic.table[6,2] <- Inf}
  # 每個延時P-value排序，變成數值再排序
  dist.choice <- rank(as.numeric(dist.table[(1:dist),i])) 
  dist.table[length(candidate)+1,i] <- candidate[which.max(dist.choice)] #最大的P-value對應的機率分布
  # 每個延時AIC-value排序，變成數值再排序
  aic.choice <- rank(as.numeric(aic.table[(1:dist),i])) 
  aic.table[length(candidate)+1,i] <- candidate[which.min(aic.choice)] #最小的AIC-value對應的機率分布
}