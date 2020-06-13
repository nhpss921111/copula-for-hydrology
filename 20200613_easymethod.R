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
candidate <- c("norm","lnorm","gumbel","gamma3","lgamma3")

dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))

variable <- cbind(Q,S)
dist.table <- matrix(nrow=length(candidate)+1,ncol=dim(variable))
rownames(dist.table) <- c("norm","lnorm","gumbel","gamma3","lgamma3","good dist")



for(i in 1:dim(variable)[2]){
  var <- variable[,i]
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
    }else if(candidate[dist] == "gamma3"){
      prefit(var,"gamma3","mle",list(shape=1,scale=1,thres=1),lower=-Inf, upper=Inf)
      md <- fitdist(var, dist = dist.char[1], start=list(shape=1,scale=1,thres=0.1),
                    lower = c(0,0,0),upper=c(Inf,Inf,min(var)),control=list(trace=1, REPORT=1))
      par1 <- md$estimate[1] #fitting參數1
      par2 <- md$estimate[2] #參數2
      par3 <- md$estimate[3] #參數3
      print(c(par1, par2, par3))
      #plot(md)
    }else if(candidate[dist] == "lgamma3"){
      md <- fitdist(var, dist = dist.char[1],start=list(shape=0.1,scale=0.1,thres=0.1),lower=c(0.1,0.1,0.1),
                    control=list(trace=1, REPORT=1))
      par1 <- md$estimate[1] #fitting參數1
      par2 <- md$estimate[2] #參數2
      par3 <- md$estimate[3] #參數3
      print(c(par1, par2, par3))
      #plot(md)
    }else{break}
    
    
    
    
    # ------------------------------- K-S test ----------------------
    print("KS test")
    if(candidate[dist] == "norm"){xi.cdf <- get(dist.char[3])(var, par1,par2)}
    if(candidate[dist] == "lnorm"){xi.cdf <- get(dist.char[3])(var, meanlog = par1, sdlog = par2)}
    if(candidate[dist] == "gumbel"){xi.cdf <- get(dist.char[3])(var, par1, par2)}
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



# weibull distribution
Q.par.weibull   <- fitdist(Q, distr = "weibull", method = "mle") # Q parameter
Q.shape.weibull <- Q.par.weibull[["estimate"]][["shape"]]
Q.scale.weibull <- Q.par.weibull[["estimate"]][["scale"]]
S.par.weibull   <- fitdist(S, distr = "weibull", method = "mle",lower=c(0, 0))# s parameter
S.shape.weibull <- S.par.weibull[["estimate"]][["shape"]]
S.scale.weibull <- S.par.weibull[["estimate"]][["scale"]]
## gamma distribution
Q.par.gamma   <- fitdist(Q, distr = "gamma", method = "mle",lower=c(0, 0)) # Q parameter
Q.shape.gamma <- Q.par.gamma[["estimate"]][["shape"]]
Q.rate.gamma  <- Q.par.gamma[["estimate"]][["rate"]]
S.par.gamma   <- fitdist(S, distr = "gamma", method = "mle",lower=c(0, 0))# s parameter
S.shape.gamma <- S.par.gamma[["estimate"]][["shape"]]
S.rate.gamma  <- S.par.gamma[["estimate"]][["rate"]]
## log-normal distribution
Q.par.lnorm       <- fitdist(Q,"lnorm",method = "mle") # Q parameter
Q.meanlog.lnorm   <- Q.par.lnorm[["estimate"]][["meanlog"]]
Q.sdlog.lnorm     <- Q.par.lnorm[["estimate"]][["sdlog"]]
S.par.lnorm       <- fitdist(S,"lnorm",method = "mle") # S parameter
S.meanlog.lnorm   <- S.par.lnorm[["estimate"]][["meanlog"]]
S.sdlog.lnorm     <- S.par.lnorm[["estimate"]][["sdlog"]]
# pearson tpye three distribution
# Q.par.lgamma      <- fitdist(Q,"lgamma",method = "mle") # Q parameter
# Q.shapelog.lgamma <- Q.par.lgamma[["estimate"]][["shapelog"]]
# Q.ratelog.lgamma  <- Q.par.lgamma[["estimate"]][["ratelog"]]
# S.par.lgamma      <- fitdist(S,"lgamma",method = "mle") # S parameter
# S.shapelog.lgamma <- S.par.lgamma[["estimate"]][["shapelog"]]
# S.ratelog.lgamma  <- S.par.lgamma[["estimate"]][["ratelog"]]
# gumbel distribution
Q.par.gumbel      <- eevd(Q, method = "mle") #EnvStats # Q parameter
Q.location.gumbel <- Q.par.gumbel[["parameters"]][["location"]]
Q.scale.gumbel    <- Q.par.gumbel[["parameters"]][["scale"]]
S.par.gumbel      <- eevd(S, method = "mle") #EnvStats # S parameter
S.location.gumbel <- S.par.gumbel[["parameters"]][["location"]]
S.scale.gumbel    <- S.par.gumbel[["parameters"]][["scale"]]

plot(Q.par.weibull)
plot(Q.par.gamma)
plot(Q.par.lnorm)
#plot(Q.par.lgamma)


plot(S.par.weibull)
plot(S.par.gamma)
plot(S.par.lnorm)
#plot(S.par.lgamma)