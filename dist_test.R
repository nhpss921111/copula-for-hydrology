# 撰寫日期：2020/06/13
# 統計檢定
# 以"家源橋"為例
# 候選分布：norm,lnorm,gumbel,weibull,gamma
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
setwd("E:/R_reading/CHIA-YUANG")
data <- read.xlsx(file.path(getwd(),"final_data.xlsx"),sheetIndex="Sheet1",
                  startRow = 1,endRow = 1324,
                  header = T,colIndex =2:3,encoding = "UTF-8")
attach(data)
#
# ----------------------------- Remove missing value ------------------------------
#
#summary(data)
#describe(data)
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
candidate <- c("norm","lnorm","gumbel","weibull","gamma")
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
# pvalue.table 放置p-value
pvalue.table <- matrix(nrow=length(candidate)+1,ncol=dim(variable)[2]) # +1 是為了要放最佳分布的欄位
rownames(pvalue.table) <- c(candidate,"good dist")
colnames(pvalue.table) <- c(colnames(variable))
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
      next
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
    pvalue.table[dist,i] <- result$p.value
    #
    
    # -------------------------- Chi-square equal probability method ---------------------
    print("Equal Probability method")
    ## Q的組數:k
    k <- round(1+3.3*log10(length(var)))
    print(paste0("理論組數",k))
    
    adjust <- 0
    k.after <- k
    k.ini <- k
    Oi <- c(0)
    
    for (i in c(1:k)){
      #if(any(Oi<5)){
      
      if(k.ini==k.after+adjust){
        
        ## Q分成K組，每組的累積機率 Fx (x:Q,S)
        Fx <- c()
        i <- 1
        while(i <= k){
          b <- 1/k*i
          Fx <- append(Fx,b)
          i = i + 1
        }
        #print(paste("累積機率: ", Fx))
        
        
        ## 分成K組，每組的機率密度 fx (x:Q,S)
        
        fx <- rep(1/k, time=k)
        
        ## 累積機率對應到的水文量 xi( Qi 和 Si )
        ### 5種機率分布套配 (用累積機率回推 xi)
        
        if (candidate[dist]=="norm"){
          xi <- get(dist.char[4])(Fx, par1, par2)
          xi[k] <- max(var)
          
        }else if (candidate[dist]=="lnorm"){
          xi <- get(dist.char[4])(Fx, meanlog = par1, sdlog = par2)
          xi[k] <- max(var)
          
        }else if (candidate[dist]=="gumbel"){
          xi <- get(dist.char[4])(Fx, par1, par2)
          xi[k] <- max(var)
          
        }else if (candidate[dist]=="weibull"){
          xi <- get(dist.char[4])(Fx,  par1, par2)
          xi[k] <- max(var)
          
        }else if (candidate[dist]=="gamma"){
          xi <- get(dist.char[4])(Fx, par1,  par2)
          xi[k] <- max(var)
          
        }else{break}
        
        ## 每個等機率區間的個數 Oi
        O <- c()
        i <- 1
        for(i in c(1:k)){
          a <- xi[i]>var
          table(a)["TRUE"]
          n <- length(a[a == TRUE])
          O <- append(O,n)
        }
        #print(O)
        
        Oi <- diff(O) #累積次數換成區間次數
        if (length(Oi)!= k){
          Oi <- append(Oi, O[1], after=0) # 把第一個值補回去
        }else{break}
        #print(cbind(xi,Oi))
        
        # 檢查每組至少要5個 Oi[]>=5
        
        adjust <- adjust +1
        if(any(Oi<5)){
          k.after <- k.ini - adjust #實際組數
        }else{break}
        
        k <- k.after
      
      }else{break}
      
    }
    
    print(paste0("實際組數", k))
    # 合併Qi與Oi
    #print(cbind(xi,Oi))
    
    # 卡方檢定
    # function explain:
    #   chisq.test(x,p) x:每組裡面的個數, p:每個區間的機率密度
    if( k > 2 ){
      chisq.result <- chisq.test(Oi, p = fx, rescale.p = TRUE, simulate.p.value = TRUE)
      print(paste0(candidate[dist], "的等機率卡方檢定P-value: ", chisq.result$p.value))
    }else{
      print("skip equal probability method")
    }
    
    
    # ----------------------------- chi square test ------------------------
    #
    print("Equal Width method")
    #
    # 分組
    #
    k <- round(1+3.3*log10(length(var)))
    print(paste0("理論分組數:", k))
    adjust <- 0
    k.after <- k
    k.ini <- k
    Oi <- c(0)
    
    for(i in c(1:k)){
      
      #if(any(Oi<5)){
      
      if(k.ini==k.after+adjust){
        
        #等間距的 xi
        interval <- (ceiling(max(var))-floor(min(var)))/k
        
        xi <- seq(floor(min(var+interval)), ceiling(max(var+1)), by = interval ) #floor(比給定值小的整數) # "+1"為了組數可以維持k
        
        ## 每個等間距區間的個數 Oi
        O <- c()
        i <- 1
        for(i in c(1:k)){
          a <- xi[i]>var
          table(a)["TRUE"]
          n <- length(a[a == TRUE])
          O <- append(O,n)
        }
        Oi <- diff(O) #累積次數換成區間次數
        if (length(Oi)!= k){
          Oi <- append(Oi, O[1], after=0) # 把第一個值補回去
        }else{break}
        # 檢查每組至少要5個 Oi[]>=5
        adjust <- adjust +1
        if(any(Oi<5)){
           k.after <- k.ini - adjust #實際組數
        }else{break}
        k <- k.after
      }else{break}
    }
    print(paste0("實際分組數:", k))
    #print(cbind(xi,Oi))
    
    #每個間距的累積機率
    library(zoo)
    #brks.cdf <- get(dist.char[3])(xi, par1, par2)
    if(candidate[dist] == "norm"){xi.cdf <- get(dist.char[3])(xi, par1, par2)}
    if(candidate[dist] == "lnorm"){xi.cdf <- get(dist.char[3])(xi, meanlog = par1, sdlog = par2)}
    if(candidate[dist] == "gumbel"){xi.cdf <- get(dist.char[3])(xi, par1, par2)}
    if(candidate[dist] == "weibull"){xi.cdf <- get(dist.char[3])(xi, par1, par2)}
    if(candidate[dist] == "gamma"){xi.cdf <- get(dist.char[3])(xi, par1, par2)}
    if(candidate[dist] == "lgamma"){xi.cdf <- get(dist.char[3])(xi, par1, par2)}
    
    null.probs <- rollapply(xi.cdf, 2, function(x) x[2]-x[1])
    
    null.probs <- append(null.probs, xi.cdf[1], after = 0)
    
    if( k > 2 ){
      result <- chisq.test(Oi, p = null.probs, rescale.p = TRUE, simulate.p.value = TRUE)
      print(paste0(candidate[dist], "的等間距卡方檢定P-value: ", result$p.value))
    }else{
      print("skip equal width method")
      next
    }
    
    #
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
  if(i==2 | candidate[dist] == "lgamma"){pvalue.table[6,2] <- -Inf; aic.table[6,2] <- Inf}
  # 每個延時P-value排序，變成數值再排序
  pvalue.choice <- rank(as.numeric(pvalue.table[(1:dist),i])) 
  pvalue.table[length(candidate)+1,i] <- candidate[which.max(pvalue.choice)] #最大的P-value對應的機率分布  
  # 每個延時AIC-value排序，變成數值再排序
  aic.choice <- rank(as.numeric(aic.table[(1:dist),i])) 
  aic.table[length(candidate)+1,i] <- candidate[which.min(aic.choice)] #最小的AIC-value對應的機率分布
}

