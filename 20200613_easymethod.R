# 撰寫日期：2020/06/13
# 統計檢定
# 候選分布：norm,lnorm,gumbel,gamma3,lgamma3
# 用最大概似法
# K-S檢定
# AIC準則選最佳分布
# By連育成

rm(list = ls())
library("xlsx") #讀取excel
require("xlsx")
library("stats") #機率分布
library("grDevices")
library("actuar") #機率分布
library("psych")
library("dplyr") #資料整理
library("FAdist") # 水文上常用機率分布
require("grid")
library("vcd") # 估計gumbel參數
library("fitdistrplus")
require("survival")
require("npsurv")
require("lsei")
library("EnvStats") 
library("ismev")
library("mgcv")
library("nlme")
library("reliaR")
library("lestat") #inverse CDF
library("goftest") #適合度檢定
library("goft") #適合度檢定
library("gumbel")
library(PearsonDS)


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

# ------------------------------------------------------------------------------------

print("參數估計")
candidate <- c("norm","lnorm","gumbel","weibull","gamma3")

dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))

variable <- cbind(Q,S)
dist.table <- matrix(nrow=length(candidate)+1,ncol=dim(variable)[2])
rownames(dist.table) <- c("norm","lnorm","gumbel","weibull","gamma3","good dist")



for(i in 1:dim(variable)[2]){
  var <- variable[,i]
  if (i==1){print("Q")}
  if (i==2){print("S")}
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
      print(c(par1, par2))
    }else if(candidate[dist] == "lnorm"){
      md <- fitdist(var, dist = dist.char[1], start = list(meanlog=1, sdlog=1))
      par1 <- md$estimate[1] #fitting參數1
      par2 <- md$estimate[2] #參數2
      print(c(par1, par2))
    }else if(candidate[dist] == "gumbel"){
      md <- eevd(var,method = "mle") 
      par1 <- md$parameters[1] #fitting參數1
      par2 <- md$parameters[2] #參數2
      print(c(par1, par2))
    }else if(candidate[dist] == "weibull"){
      md <- fitdist(var, dist = dist.char[1])
      par1 <- md$estimate[1] #fitting參數1
      par2 <- md$estimate[2] #參數2
      print(c(par1, par2))
    }else if(candidate[dist] == "gamma3"){
      #prefit(var,"gamma3","mle",list(shape=0.1,scale=200,thres=1),lower=0, upper=Inf)
      md <- fitdist(var, dist = dist.char[1], start=list(shape=1,scale=50,thres=1),
                    lower = c(0,0,0),upper=c(Inf,Inf,min(var))) #control=list(trace=1, REPORT=1)
      par1 <- md$estimate[1] #fitting參數1
      par2 <- md$estimate[2] #參數2
      par3 <- md$estimate[3] #參數3
      print(c(par1, par2, par3))
      #plot(md)
    # }else if(candidate[dist] == "lgamma3"){
    #   md <- fitdist(var, dist = dist.char[1],start=list(shape=10,scale=10,thres=1)
    #                 ,lower=c(0.1,0.1,0.1),upper=c(Inf,Inf,min(var)))
    #   par1 <- md$estimate[1] #fitting參數1
    #   par2 <- md$estimate[2] #參數2
    #   par3 <- md$estimate[3] #參數3
    #   print(c(par1, par2, par3))
      #plot(md)
    }else{break}
    
    # ------------------------------- K-S test ----------------------
    print("KS test")
    if(candidate[dist] == "norm"){xi.cdf <- get(dist.char[3])(var, par1,par2)}
    if(candidate[dist] == "lnorm"){xi.cdf <- get(dist.char[3])(var, meanlog = par1, sdlog = par2)}
    if(candidate[dist] == "gumbel"){xi.cdf <- get(dist.char[3])(var, par1, par2)}
    if(candidate[dist] == "weibull"){xi.cdf <- get(dist.char[3])(var, par1, par2)}
    if(candidate[dist] == "gamma3"){xi.cdf <- get(dist.char[3])(var, shape=par1, scale=par2, thres=par3)}
    if(candidate[dist] == "lgamma3"){xi.cdf <- get(dist.char[3])(var, shape=par1, scale=par2, thres=par3)}
    
    result <- ks.test(var, dist.char[3], par1,par2)
    print(paste0(candidate[dist], "的KS檢定P-value: ", result$p.value))
    
    # ------------- 將P-value整理成表格-------------------
    
    dist.table[dist,i] <- result$p.value
    
  }
  # 每個延時P-value排序，變成數值再排序
  dist.choice <- rank(as.numeric(dist.table[(1:dist),i])) 
  dist.table[length(candidate)+1,i] <- candidate[which.max(dist.choice)] #最大的P-value對應的機率分布
  
}

