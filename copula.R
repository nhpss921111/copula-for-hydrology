# copula
# 開始日期：2020/07/08
# 完成日期：2020/08/19
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
library(scatterplot3d)
# 輸出表格
export <- c("n")
# 輸入邊際分布
margin.dist <-c("lnorm","lnorm") 
# Read data from excel flie
month<- c(1:12) # 請輸入月分：
input <- c(paste0(month,"month.csv"))
#
pvalue.table <- matrix(nrow=12,ncol=3)
colnames(pvalue.table) <- c("gumbel","frank","clayton")
for (m in month){
  setwd("E:/R_reading/CHIA-YUANG")
  data <- read.csv(file.path(getwd(),input[m]),header = T) # 請輸入月分：
  data <- data[,-1]
  attach(data)
  print(paste0(m,"月"))
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
  #
  # ------------------- ML method -----------------
  #
  var_a <- pobs(Q)
  var_b <- pobs(S)
  data.probs <- cbind(var_a, var_b)
  plot(var_a,var_b,col="red")
  # 建立copula
  g3 <- gumbelCopula(1,use.indepC="FALSE")
  f3 <- frankCopula(1)
  c3 <- claytonCopula(1)
  a3 <- amhCopula(1)
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
  aMvd2 <- mvdc(a3,c("lnorm","lnorm"),
                param =list(list(meanlog=par.table[2,1], sdlog=par.table[2,2]), 
                            list(meanlog=par.table[2,3], sdlog=par.table[2,4])))
  # 估計 mvdc參數：margins的參數 和 copula的參數
  mm <- apply(variable, 2, mean)
  vv <- apply(variable, 2, var)
  b1.0 <- c(mm[1]^2/vv[1], vv[1]/mm[1])
  b2.0 <- c(mm[2]^2/vv[2], vv[2]/mm[2])
  a.0 <- sin(cor(variable[, 1], variable[, 2], method = "kendall") * pi/2)
  # start <- c(b1.0, b2.0, a.0)
  #fit.ml.g <- fitMvdc(variable, gMvd2, start=as.numeric(start))
  #fit.ml.f <- fitMvdc(variable, fMvd2, start=as.numeric(start))
  #fit.ml.c <- fitMvdc(variable, cMvd2, start=as.numeric(start))
  #fit.ml.a <- fitMvdc(variable, aMvd2, start=as.numeric(start))
  #optim.control = list(trace = TRUE, maxit = 1000))
  # 
  # ------------------------ IFM method -------------------------
  # 
  loglik.marg <- function(b, x) sum(dlnorm(x, meanlog = b[1], sdlog = b[2], log = TRUE))
  ctrl <- list(fnscale = -1)
  b1hat <- optim(b1.0, fn = loglik.marg, x = variable[, 1], control = ctrl)$par
  b2hat <- optim(b2.0, fn = loglik.marg, x = variable[, 2], control = ctrl)$par
  udat <- cbind(plnorm(variable[, 1], meanlog = b1hat[1], sdlog = b1hat[2]),
                plnorm(variable[, 2], meanlog = b2hat[1], sdlog = b2hat[2]))
  fit.ifm.g <- fitCopula(gMvd2@copula, udat, start = a.0)
  summary(fit.ifm.g)
  fit.ifm.f <- fitCopula(fMvd2@copula, udat, start = a.0)
  fit.ifm.c <- fitCopula(cMvd2@copula, udat, start = a.0)
  fit.ifm.a <- fitCopula(aMvd2@copula, udat, start = a.0)
  #
  # ----------------------- 適合度檢定 -----------------------------
  #
  print(paste0("第",m,"個月的適合度檢定"))
  gfg <- gofCopula(gumbelCopula(fit.ifm.g@estimate, dim=2), pobs(variable),N = 2000
                   ,method = "Sn", estim.method = "mpl", simulation = "pb")
  gff <- gofCopula(frankCopula(fit.ifm.f@estimate, dim=2), pobs(variable),N = 2000
                   ,method = "Sn", estim.method = "mpl", simulation = "pb")
  gfc <- gofCopula(claytonCopula(fit.ifm.c@estimate, dim=2), pobs(variable),N = 1000
                   ,method = "Sn", estim.method = "mpl", simulation = "pb", ties=TRUE
                   ,optim.method = "BFGS")
  #gfa <- gofCopula(amhCopula(dim=2), pobs(variable),N = 2000)
  pvalue.table[m,1] <- gfg$p.value
  pvalue.table[m,2] <- gff$p.value
  pvalue.table[m,3] <- gfc$p.value
  
  #aic.choice <- rank(as.numeric(aic.table[(1:dist),i])) 
  #aic.table[length(candidate)+1,i] <- candidate[which.max(aic.choice)] #最大的P-value對應的機率分布
  #
  # # ------------------------- plotting --------------------------------
  # #
  # #
  # # Build the bivariate distribution
  # my_dist <- mvdc(gumbelCopula(param = fit.ifm.g@estimate, dim = 2), margins = c("lnorm","lnorm"), 
  #                               paramMargins=list(list(meanlog=par.table[2,1], sdlog=par.table[2,2]), 
  #                                            list(meanlog=par.table[2,3], sdlog=par.table[2,4])))
  # 
  # # Generate random sample observations from the multivariate distribution
  # #v <- rMvdc(5000, my_dist)
  # # Compute the density
  # pdf_mvd <- dMvdc(data.probs, my_dist)
  # # Compute the CDF
  # cdf_mvd <- pMvdc(data.probs, my_dist)
  # 
  # # 3D plain scatterplot of the generated bivariate distribution
  # par(mfrow = c(1, 2))
  # scatterplot3d(var_a, var_b, pdf_mvd, color="red", main="Density", xlab = "Q", ylab="Qs", zlab="pMvdc",pch=".")
  # scatterplot3d(var_a, var_b, cdf_mvd, color="red", main="CDF", xlab = "Q", ylab="Qs", zlab="pMvdc",pch=".")
  # persp(my_dist, dMvdc, xlim = c(0, 1), ylim=c(0, 1), main = "Density",xlab = "Q", ylab="Qs")
  # contour(my_dist, dMvdc, xlim = c(0, 1), ylim=c(0, 1), main = "Contour plot",xlab = "Q", ylab="Qs")
  # persp(my_dist, pMvdc, xlim = c(0, 1), ylim=c(0, 1), main = "CDF",xlab = "Q", ylab="Qs")
  # contour(my_dist, pMvdc, xlim = c(0, 1), ylim=c(0, 1), main = "Contour plot",xlab = "Q", ylab="Qs")
  
  # # ------------------------------------------------------------------------------
  # # Compute the density
  # pdf_ <- dCopula(Q, mycopula)
  # # Compute the CDF
  # cdf <- pCopula(Q, mycopula)
  # # Compute the density
  # pdf_mvd <- dMvdc(v, multivariate_dist)
  # # Compute the CDF
  # cdf_mvd <- pMvdc(v, multivariate_dist)
  # #
  # # 3D plain scatterplot of the density, plot of the density and contour plot
  # par(mfrow = c(1, 3))
  # scatterplot3d(variable[, 1], variable[, 2], pdf_, color="red", main="Density", xlab ="Q", ylab="Qs", zlab="dCopula", pch=".")
  # persp(mycopula, dCopula, main ="Density")
  # contour(mycopula, dCopula, xlim = c(0, 1), ylim=c(0, 1), main = "Contour plot")
  # #
  # # 3D plain scatterplot of the CDF, plot of the CDF and contour plot
  # par(mfrow = c(1, 3))
  # scatterplot3d(u[,1], u[,2], cdf, color="red", main="CDF", xlab = "u1", ylab="u2", zlab="pCopula",pch=".")
  # persp(mycopula, pCopula, main = "CDF")
  # contour(mycopula, pCopula, xlim = c(0, 1), ylim=c(0, 1), main = "Contour plot")
  # #
  # # 3D plain scatterplot of the multivariate distribution
  # par(mfrow = c(1, 2))
  # scatterplot3d(v[,1],v[,2], pdf_mvd, color="red", main="Density", xlab = "u1", ylab="u2", zlab="pMvdc",pch=".")
  # scatterplot3d(v[,1],v[,2], cdf_mvd, color="red", main="CDF", xlab = "u1", ylab="u2", zlab="pMvdc",pch=".")
  # persp(multivariate_dist, dMvdc, xlim = c(-4, 4), ylim=c(0, 2), main = "Density")
  # contour(multivariate_dist, dMvdc, xlim = c(-4, 4), ylim=c(0, 2), main = "Contour plot")
  # persp(multivariate_dist, pMvdc, xlim = c(-4, 4), ylim=c(0, 2), main = "CDF")
  # contour(multivariate_dist, pMvdc, xlim = c(-4, 4), ylim=c(0, 2), main = "Contour plot")
}
# export table
if (export=="y"){
  colnames(pvalue.table)<-c("gumbelcopula","frankcopula","claytoncopula")
  file <- paste("E:/R_output/CHIA-YUANG/result/copula_pvalue.csv", sep="")
  write.csv(pvalue.table,file)
}

