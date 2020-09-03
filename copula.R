# copula
# 開始日期：2020/07/08
# 完成日期：2020/09/03
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
library(gumbel) # for gumbel copula
library(ggplot2) # 繪圖用
library(VineCopula)
library(copula) # 聯結函數
library(scatterplot3d) #畫3D圖
library(ggplot2) #繪圖用

# ======================== 執行前請先設定以下參數 =========================
#
# Read data from csv flie
month<- c(5:8) # 請輸入月分： (連續輸入、單獨輸入、跳著輸入都可以)
input <- c(paste0(month,"month.csv"))
#
# 1.相同流量不同月份之情況(只能擇一)
sameq <- c("y")  # "n" or "y"
if (sameq == "y"){
  qn <- c(50) # 輸入流量：
}
# 2.相同月份不同流量之情況(只能擇一)
samemonth <- c("n") # "n" or "y"
if (samemonth=="y"){
  qn <- seq(from=25, to=100, length.out=4) # 輸入流量:
}
# 輸入邊際分布
margin.dist <-c("lnorm","lnorm") # 請輸入邊際分布：
#
# 執行適合度檢定
gof <- c("n") # "n" or "y"
#
# 輸出表格
export <- c("n") # "n" or "y"
#
# copula function plotting
cfp <- c("n") # "n" or "y"

# ============================ 主程式 ==================================
#
# 建立copula參數表格
ifm.table <- matrix(nrow=12,ncol=3)
colnames(ifm.table) <- c("gumbel","frank","clayton")

# 建立p-value表格
pvalue.table <- matrix(nrow=12,ncol=3)
colnames(pvalue.table) <- c("gumbel","frank","clayton")

pdf.new <- c() # 不同月份下，相同流量之PDF
j <- 0 # 計算月份次數

# 主要執行迴圈
for (m in month){
  j <- j + 1 # 計算月份次數
  setwd("E:/R_reading/CHIA-YUANG")
  data <- read.csv(file.path(getwd(),input[j])) 
  data <- data[,-1]
  attach(data)
  print(paste0(m,"月"))
  #
  # ------------------------ Remove missing value ------------------------------
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
  # 邊際分布參數表格
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
      if(i==1){ #Q
        par.table[dist,1] <- as.numeric(par1) 
        par.table[dist,2] <- as.numeric(par2)}
      if(i==2){ #S
        par.table[dist,3] <- as.numeric(par1) 
        par.table[dist,4] <- as.numeric(par2)}
    }
  
  }
  #file <- paste0("E:/R_output/CHIA-YUANG/result/",m,"margin_par.csv", sep="")
  #write.csv(par.table,file)
  #
  # ------------------------ IFM method -------------------------
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
  # 對數概似函數
  loglik.marg <- function(b, x) sum(dlnorm(x, meanlog = b[1], sdlog = b[2], log = TRUE))
  ctrl <- list(fnscale = -1)
  b1hat <- optim(b1.0, fn = loglik.marg, x = variable[, 1], control = ctrl)$par
  b2hat <- optim(b2.0, fn = loglik.marg, x = variable[, 2], control = ctrl)$par
  udat <- cbind(plnorm(variable[, 1], meanlog = b1hat[1], sdlog = b1hat[2]),
                plnorm(variable[, 2], meanlog = b2hat[1], sdlog = b2hat[2]))
  if(m!=12){
    fit.ifm.g <- fitCopula(gMvd2@copula, udat, start = a.0)
    fit.ifm.f <- fitCopula(fMvd2@copula, udat, start = a.0)
    fit.ifm.c <- fitCopula(cMvd2@copula, udat, start = a.0)}
  if (m==12){
    fit.ifm.g <- fitCopula(gMvd2@copula, udat, start = a.0)
    fit.ifm.f <- fitCopula(fMvd2@copula, udat, start = a.0)}
  #fit.ifm.a <- fitCopula(aMvd2@copula, udat, start = a.0)
  #儲存參數到ifm.table
  if(m!=12){
    ifm.table[m,1] <- fit.ifm.g@estimate
    ifm.table[m,2] <- fit.ifm.f@estimate
    ifm.table[m,3] <- fit.ifm.c@estimate}
  if (m==12){
    ifm.table[m,1] <- fit.ifm.g@estimate
    ifm.table[m,2] <- fit.ifm.f@estimate} #12月無法使用
  #
  # ------------------------ Goodness of fit test -----------------------------
  #
  if (gof=="y"){
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
    #
  }
  #aic.choice <- rank(as.numeric(aic.table[(1:dist),i]))
  #aic.table[length(candidate)+1,i] <- candidate[which.max(aic.choice)] #最大的P-value對應的機率分布
  #
  # # ------------------------ copula function plotting --------------------------------
  # #
  # 尚未寫選擇出圖之copula的迴圈!!!
  ## Generate the gunbel copula and sample some observations
  #
  if (cfp == "y"){
    mycopula <- gumbelCopula(param = fit.ifm.g@estimate, dim = 2)
    u <- rCopula(2000, mycopula)
    # Compute the density
    pdf_ <- dCopula(u, mycopula)
    # Compute the CDF
    cdf <- pCopula(u, mycopula)
  
    ## 3D plain scatterplot of the density, plot of the density and contour plot
    par(mfrow = c(1, 3))
    scatterplot3d(u[,1], u[,2], pdf_, color="red", main="Density",xlab ="u1", ylab="u2", zlab="dCopula", pch=".")
    persp(mycopula, dCopula, main ="Density")
    contour(mycopula, dCopula, xlim = c(0, 1), ylim=c(0, 1), main = "Contour plot")
  
    ## 3D plain scatterplot of the CDF, plot of the CDF and contour plot
    par(mfrow = c(1, 3))
    scatterplot3d(u[,1], u[,2], cdf, color="red", main="CDF", xlab = "u1", ylab="u2", zlab="pCopula",pch=".")
    persp(mycopula, pCopula, main = "CDF")
    contour(mycopula, pCopula, xlim = c(0, 1), ylim=c(0, 1), main = "Contour plot")
  
    # Build the bivariate distribution
  
    my_dist <- mvdc(mycopula, margin.dist,
                   paramMargins=list(list(meanlog=par.table[2,1], sdlog=par.table[2,2]),
                                     list(meanlog=par.table[2,3], sdlog=par.table[2,4])))
    # Generate random sample observations from the multivariate distribution
    v <- rMvdc(5000, my_dist)
    # Compute the density
    pdf_mvd <- dMvdc(v, my_dist)
    # Compute the CDF
    cdf_mvd <- pMvdc(v, my_dist)
    ## 3D plain scatterplot of the generated bivariate distribution
    par(mfrow = c(1, 2))
    scatterplot3d(v[,1],v[,2], pdf_mvd, color="red", main= paste0("第",m,"個月的Density"), xlab = "Q", ylab="Qs", zlab="pMvdc",pch=".")
    scatterplot3d(v[,1],v[,2], cdf_mvd, color="red", main=paste0("第",m,"個月的CDF"), xlab = "Q", ylab="Qs", zlab="pMvdc",pch=".")
    persp(my_dist, dMvdc, xlim = c(0, 50), ylim=c(0, 1500), main = paste0("第",m,"個月的Density"), xlab = "Q", ylab="Qs")
    contour(my_dist, dMvdc, xlim = c(0, 50), ylim=c(0, 1500), main =  paste0("第",m,"個月的Contour plot"), xlab = "Q", ylab="Qs")
    persp(my_dist, pMvdc, xlim = c(0, 50), ylim=c(0, 1500), main =  paste0("第",m,"個月的CDF"),xlab = "Q", ylab="Qs")
    contour(my_dist, pMvdc, xlim = c(0, 50), ylim=c(0, 1500), main =  paste0("第",m,"個月的Contour plot"), xlab = "Q", ylab="Qs")
  }
  #
  # ----------------- 建立 conditional copula distribution function ------------------
  #
  #gumbel copula parameter
  print("開始建立conditional copula distribution function")
  gum.cop <- gumbelCopula(param=fit.ifm.g@estimate)
  # 建立雙變數分布函數
  Mvdc <- mvdc(gum.cop, margin.dist,
               param =list(list(meanlog=par.table[2,1], sdlog=par.table[2,2]),
                           list(meanlog=par.table[2,3], sdlog=par.table[2,4])))
  # 建立Q與Qs表格: x.samp
  test.samp <- matrix(ncol=2)
  qs <- seq(from=1, to=250, by=1) # 調整Qs大小
  test.samp <- matrix(nrow=250,ncol=1) # 調整Qs大小
  n <- 1
  # 全部範圍出圖
  # while (q<max(Q)){ # 調整Q的大小
  #   test.samp[,1] <- q
  #   x.samp <- cbind(test.samp,qs)
  #   copula_density <- dMvdc(x.samp, Mvdc, log=FALSE)
  #   fqs <- dlnorm(qs,meanlog=par.table[2,3], sdlog=par.table[2,4])
  #   fx <- copula_density*fqs
  #   par(mfrow = c(1, 1))
  #   # 輸出PDF的圖
  #   setwd("C:/Users/user/Desktop/PDF")
  #   png(paste0(m,"月流量",n,"cms.png"),width = 1250, height = 700, units = "px", pointsize = 12)
  #   plot(qs,fx,type="line",xlab="Qs",ylab="probs",main=paste0(m,"月時流量為",n,"cms之輸砂量PDF"))
  #   dev.off()
  #   n <- n+1
  # }
  #
  # 同月份，不同流量之PDF疊在一起
  #
  if (samemonth == "y"){
    n <- 1
    par(mfrow = c(1, 1))
    
    pdf.new <- c()
    for (q in qn){
      test.samp[,1] <- q
      x.samp <- cbind(test.samp,qs)
      copula_density <- dMvdc(x.samp, Mvdc, log=FALSE)
      fqs <- dlnorm(qs,meanlog=par.table[2,3], sdlog=par.table[2,4])
      fx <- copula_density*fqs
      q <- paste0(q,"cms")
      pdf<- data.frame(q,qs,fx)
      pdf.new <- rbind(pdf.new,pdf)
    }
    setwd("C:/Users/user/Desktop/PDF")
    png(paste0(m,"月流量",n,"cms.png"),width = 1250, height = 700, units = "px", pointsize = 12)
    sameMonth <- ggplot(pdf.new) +
      geom_line(aes(x = qs, y = fx, color = q))+
      labs(title=paste0(m,"月之PDF"),
           x="輸砂量(公噸/日)",
           y="PDF") +
      theme_bw() + # 白底
      theme(panel.grid.major = element_blank()) + # 隱藏主要格線
      theme(panel.grid.minor = element_blank())  # 隱藏次要格線
    dev.off()
    print(sameMonth)
  }
  #
  # 同流量，不同月份之PDF疊在一起
  #
  if (sameq=="y"){
    n <- 1
    par(mfrow = c(1, 1))
    for (q in qn){
      test.samp[,1] <- q
      x.samp <- cbind(test.samp,qs)
      copula_density <- dMvdc(x.samp, Mvdc, log=FALSE)
      fqs <- dlnorm(qs,meanlog=par.table[2,3], sdlog=par.table[2,4])
      fx <- copula_density*fqs
      
      q <- paste0(q,"cms")
      mon <- paste0(m,"month")
      pdf<- data.frame(mon, q, qs, fx)
      pdf.new <- rbind(pdf.new, pdf)
    }
  }

}

if (sameq == "y"){
sameQ <- ggplot(pdf.new) +
  geom_line(aes(x = qs, y = fx, color = mon))+
  labs(title=paste0("不同月下",q,"之PDF"),
       x="輸砂量(公噸/日)",
       y="PDF") +
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank())  # 隱藏次要格線
print(sameQ)
}
# export table
if (export=="y"){
  colnames(pvalue.table)<-c("gumbelcopula","frankcopula","claytoncopula")
  file <- paste("E:/R_output/CHIA-YUANG/result/copula_pvalue.csv", sep="")
  write.csv(pvalue.table,file)
}

