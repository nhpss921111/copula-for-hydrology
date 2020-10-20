# 以聯結函數建立流量及輸砂量之雙變數機率分布模式
# 開始日期：2020/07/08
# 完成日期：2020/10/19
# By 連育成
# ----------------- Setting --------------------- 
# 邊際分布：norm, lnorm, gumbel, weibull, gamma
# 邊際分布參數估計：mle
# 邊際分布適合度檢定：K-S test
# 聯結函數：gumbel, frank, clayton
# 聯結函數估計方法：itau, irho, mpl, ml
# 聯結函數適合度檢定： method=Sn, optim=BFGS
# -----------------------------------------------
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
library(VineCopula) # for high dimension copula function
library(copula) # 聯結函數
library(scatterplot3d) #畫3D圖
#
# ============================= 執行前請先設定以下參數 =================================
#
# 1. Read data from csv flie
month<- c(5) # 請輸入月分： (連續輸入、單獨輸入、跳著輸入都可以)
input <- c(paste0(month,"month.csv"))
#
# Case(1)相同流量不同月份之情況
sameq <- c("n")  # "n" or "y" (case(1),(2) 只能擇一執行)
if (sameq == "y"){
  qn <- c(20) # 輸入流量：(只能輸入一個值)
}
# Case(2)相同月份不同流量之情況
samemonth <- c("y") # "n" or "y" (case(1),(2) 只能擇一執行)
if (samemonth=="y"){
  qn <-c(5,10,15,20)  # 輸入流量：seq(from=10.30, to=10.30, length.out=1);(起始值,最大值,總分組數)
}
observation.qs <- c("n")
if (observation.qs =="y"){
  ob.qs <- c(12245.21,71363.72,36069.72,28545.90) # 輸入觀測輸砂量：
  rc.qs <- c(8901.20,8158.45,9328.86,8727.26) # 輸入率定曲線推估輸砂量：
  qs.table <- data.frame(ob.qs,0)
}
#
# 2. 輸入邊際分布
margin.dist <-c("lnorm","lnorm") # 請輸入邊際分布：
margin.num <- c(2,2) # 邊際分布代號(norm=1, lnorm=2, gumbel=3, weibull=4, gamma=5)
# 3. 輸入聯結函數
 copula <-c("gumbelCopula") # 請輸入聯結函數：
# 4. 執行適合度檢定
gof <- c("y") # "n" or "y"
# 5. 輸出表格
export <- c("n") # "n" or "y"
# 6. copula function plotting
cfp <- c("n") # "n" or "y"
#
# 7. PDF儲存路徑請至line 313、line 348 附近修改
#
# 8. PDF之點位資料存在pdf.new裡面
#
# 9.# 建立Q表格: q.samp
qs <- seq(from=1, to=1000, by=1) # 調整Qs範圍
q.samp <- matrix(nrow=1000,ncol=1) # qs的範圍決定q的數量
#
# ===================================== 主程式 ===========================================
#
# 建立margin參數估計表格(參數估計方法：mle)(候選分布：norm, lnorm, gumbel, weibull, gamma)
fitmargin.par <- matrix(nrow=12,ncol=20)
colnames(fitmargin.par) <- c("norm.Q.par1","norm.Q.par2","norm.Qs.par1","norm.Qs.par2",
                             "lnorm.Q.par1","lnorm.Q.par2","lnorm.Qs.par1","lnorm.Qs.par2",
                             "gumbel.Q.par1","gumbel.Q.par2","gumbel.Qs.par1","gumbel.Qs.par2",
                             "weibull.Q.par1","weibull.Q.par2","weibull.Qs.par1","weibull.Qs.par2",
                             "gamma.Q.par1","gamma.Q.par2","gamma.Qs.par1","gamma.Qs.par2")
# 新增最後一行放選擇分布的參數值(預設值為0)
# fitmargin.par <- data.frame(fitmargin.par, Q.par1=0, Q.par2=0, Qs.par1=0, Qs.par2=0)
# # 建立margin的K-S檢定
# margin.ks <- matrix(nrow=12,ncol=10)
# colnames(margin.ks) <- c("norm.Q.ks","norm.Qs.ks",
#                          "lnorm.Q.ks","lnorm.Qs.ks",
#                          "gumbel.Q.ks","gumbel.Qs.ks",
#                          "weibull.Q.ks","weibull.Qs.ks",
#                          "gamma.Q.ks","gamma.Qs.ks")
# 新增最後一行放選擇分布的參數值(預設值為0)
# margin.ks <- data.frame(margin.ks, Q.par1=0, Q.par2=0, Qs.par1=0, Qs.par2=0)
# # 建立margin的AIC
# margin.aic <- matrix(nrow=12,ncol=10)
# colnames(margin.aic) <- c("norm.Q.aic","norm.Qs.aic",
#                          "lnorm.Q.aic","lnorm.Qs.aic",
#                          "gumbel.Q.aic","gumbel.Qs.aic",
#                          "weibull.Q.aic","weibull.Qs.aic",
#                          "gamma.Q.aic","gamma.Qs.aic")
# # 新增最後一行放選擇分布的參數值(預設值為0)
# margin.aic <- data.frame(margin.aic, choice.Q.dist=0, choice.Qs.dist=0)

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

pdf.new <- c() # 不同月份下，相同流量之PDF
j <- 0 # 計算月份次數

# 資料放在 pdf.new

# 主要執行迴圈
for (m in month){
  j <- j + 1 # 計算月份次數
  setwd("E:/R_reading/NEI-MAO-PU")
  data <- read.csv(file.path(getwd(),input[j])) 
  data <- data[,-1]
  attach(data)
  print(paste0(m,"月"))
  #
  # -------------------------------- Remove missing value --------------------------------
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
    }
  }
  #
  # ---------------- copula parameter estimate & Goodness of fit test ------------------
  #
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
    if (gof=="y"){
      print(paste0(m,"月的適合度檢定"))
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
    }
  }
  # 選擇pvalue最大的值，並找出是哪個聯結函數(每個聯結函數使用四種參數估計法)
  gof.pvalue[m,13] <- floor(which.max(gof.pvalue[m,])%/%(4+.1)+1) #每四個為一組
  # 最佳連結函數的參數
  fitcopula.par[m,13] <- fitcopula.par[m,which.max(gof.pvalue[m,1:12])]
  #
  # ------------------------ copula function plotting --------------------------------
  # 
  # Generate the gunbel copula and sample some observations
  #
  if (cfp == "y"){
    mycopula <- get(paste0(candidate.copula[gof.pvalue[m,13]],"Copula"))(param = fitcopula.par[m,13], dim = 2)
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
                   paramMargins=list(list(margin.par[margin.num[1],1], margin.par[margin.num[1],2]),
                                     list(margin.par[margin.num[2],3], margin.par[margin.num[2],4])))
    # Generate random sample observations from the multivariate distribution
    v <- rMvdc(5000, my_dist) # 隨機生成5000點的模型
    # Compute the density
    pdf_mvd <- dMvdc(v, my_dist)
    # Compute the CDF
    cdf_mvd <- pMvdc(v, my_dist)
    
    ## 3D plain scatterplot of the generated bivariate distribution
    par(mfrow = c(1, 2))
    scatterplot3d(v[,1],v[,2], pdf_mvd, color="red", main= paste0(m,"月的Density"), xlab = "Q", ylab="Qs", zlab="pMvdc",pch=".")
    scatterplot3d(v[,1],v[,2], cdf_mvd, color="red", main=paste0(m,"月的CDF"), xlab = "Q", ylab="Qs", zlab="pMvdc",pch=".")
    persp(my_dist, dMvdc, xlim = c(0, 50), ylim=c(0, 1500), main = paste0(m,"月的Density"), xlab = "Q", ylab="Qs")
    contour(my_dist, dMvdc, xlim = c(0, 50), ylim=c(0, 1500), main =  paste0(m,"月的Contour plot"), xlab = "Q", ylab="Qs")
    persp(my_dist, pMvdc, xlim = c(0, 50), ylim=c(0, 1500), main =  paste0(m,"月的CDF"),xlab = "Q", ylab="Qs")
    contour(my_dist, pMvdc, xlim = c(0, 50), ylim=c(0, 1500), main =  paste0(m,"月的Contour plot"), xlab = "Q", ylab="Qs")
  }
  #
  # ----------------- 建立 conditional copula distribution function ------------------
  #
  #gumbel copula parameter
  print("開始建立conditional copula distribution function")
  mycopula <- get(paste0(candidate.copula[gof.pvalue[m,13]],"Copula"))(param = fitcopula.par[m,13])
  # 建立雙變數分布函數
  Mvdc <- mvdc(mycopula, margin.dist,
               param =list(list(margin.par[margin.num[1],1], margin.par[margin.num[1],2]),
                           list(margin.par[margin.num[2],3], margin.par[margin.num[2],4])))
  v <- rMvdc(20000,Mvdc)
  plot(sort(v[,1]),sort(v[,2]))
  pdf_mvd <- dMvdc(v, Mvdc)
  # Compute the CDF
  cdf_mvd <- pMvdc(v, Mvdc)
  contour(Mvdc,pMvdc,xlim = c(0, 80), ylim=c(0, 2000),main=paste0(m,"月雙變數機率分布(CDF)"),
          xlab="流量Q (cms)",ylab="輸砂量Qs (公噸)",labcex = 1.2,lwd=1.5,drawlabels = TRUE)
  scatterplot3d(x=v[,1], y=v[,2],z=pdf_mvd,color="red", main=paste0(m,"月的PDF"), xlab = "Q", ylab="Qs", zlab="pMvdc",pch=".")
  scatterplot3d(x=v[,1], y=v[,2],z=cdf_mvd,color="red", main=paste0(m,"月雙變數機率分布(CDF)"), xlab = "流量Q (cms)", ylab="輸砂量Qs (公噸)", zlab="累積機率",pch=".")
  #scatterplot3d(x=v[,1], y=v[,2],z=checkcdf,color="red", main=paste0(m,"月的checkCDF"), xlab = "Q", ylab="Qs", zlab="pMvdc",pch=".")
  #
  # 同流量，不同月份之PDF疊在一起
  #
  if (sameq=="y"){
    n <- 1
    par(mfrow = c(1, 1))
    for (q in qn){
      q.samp[,1] <- q
      data.samp <- cbind(q.samp,qs)
      mvdc.density <- dMvdc(data.samp, Mvdc, log=FALSE)
      q.pdf <- get(paste0("d",margin.dist[1]))(qn,margin.par[margin.num[1],1], margin.par[margin.num[1],2])
      condition.pdf <- mvdc.density/q.pdf
      q <- paste0(q,"cms")
      mon <- paste0(m,"月")
      pdf<- data.frame(mon, q, qs, condition.pdf)
      pdf.new <- rbind(pdf.new, pdf)
    }
    # PDF出圖
    setwd("C:/Users/user/Desktop/PDF") # 請修改儲存路徑：
    png(paste0("流量",q,"之PDF.png"),width = 1250, height = 700, units = "px", pointsize = 12)
    sameQ <- ggplot(pdf.new) +
      geom_line(aes(x = qs, y = condition.pdf, color = mon),size=1.3) + # 畫線圖
      labs(x="輸砂量Qs (公噸)", y="PDF") + # 坐標軸單位
      scale_color_discrete(name="月份") + # 圖例名稱
      theme_bw() + # 白底
      theme(panel.grid.major = element_blank()) + # 隱藏主要格線
      theme(panel.grid.minor = element_blank()) +  # 隱藏次要格線
      theme(text=element_text(size=50)) + # 調整字型大小
      theme(legend.position = c(0.9,0.7)) # 調整圖例位置
    print(sameQ)
    dev.off()
    # CDF驗證
  }
  #
  # 同月份，不同流量之PDF疊在一起
  #
  if (samemonth == "y"){
    n <- 1
    par(mfrow = c(1, 1))
    pdf.new <- c()
    q.group <- 1
    for (q in qn){
      q.samp[,1] <- q
      data.samp <- cbind(q.samp,qs)
      mvdc.density <- dMvdc(data.samp, Mvdc, log=FALSE)
      q.pdf <- get(paste0("d",margin.dist[1]))(qn,margin.par[margin.num[1],1], margin.par[margin.num[1],2])
      condition.pdf <- mvdc.density/q.pdf[q.group]
      q <- paste0(q,"cms")
      if (observation.qs=="y"){
        colnames(qs.table) <- c("Qs","PDF")
        qs.obser <- qs.table
        for(i in c(1:length(ob.qs))){
          plus <- c(0:3)
          qs.obser <-  add_row(qs.obser,Qs=ob.qs[i],PDF=condition.pdf[ob.qs][i],.after=i+plus[i])
        }
      }
      pdf<- data.frame(q,qs,condition.pdf)
      pdf.new <- rbind(pdf.new,pdf)
      q.group <- q.group + 1
    }
    # PDF出圖
    setwd("C:/Users/user/Desktop/PDF") # 請修改儲存路徑：
    png(paste0(m,"月流量之PDF.png"),width = 1250, height = 700, units = "px", pointsize = 12)
    sameMonth <- ggplot() +
      geom_line(data=pdf.new,aes(x = qs, y = condition.pdf, color = q),size=1.3)+  # 畫線圖
      labs(x="輸砂量Qs (公噸)",y="PDF") + #座標軸名稱
      scale_color_discrete(name="流量") + #圖例名稱
      theme_bw() + # 白底
      theme(panel.grid.major = element_blank()) + # 隱藏主要格線
      theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
      theme(text=element_text(size=50)) + # 字體大小
      theme(legend.position = c(0.85,0.7)) # 調整圖例位置
    print(sameMonth)
    dev.off()
    if(observation.qs=="y"){
      setwd("C:/Users/user/Desktop/PDF") # 請修改儲存路徑：
      png(paste0(m,"月流量之PDF.png"),width = 1250, height = 700, units = "px", pointsize = 12)
      sameMonth <- ggplot() +
        geom_line(data=pdf.new,aes(x = qs, y = condition.pdf, color = q),size=1.3)+  # 畫線圖
        geom_vline(data=qs.table,aes(xintercept=Qs),color="blue",size=1.3)+ # 觀測輸砂量
        geom_vline(aes(xintercept=rc.qs),color="red",size=1.3)+ # 率定曲線的推估輸砂量
        labs(x="輸砂量Qs (公噸)",y="PDF") + #座標軸名稱
        scale_color_discrete(name="流量") + #圖例名稱
        theme_bw() + # 白底
        theme(panel.grid.major = element_blank()) + # 隱藏主要格線
        theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
        theme(text=element_text(size=50)) + # 字體大小
        theme(legend.position = c(0.85,0.7)) # 調整圖例位置
      print(sameMonth)
      dev.off()
    }
    # CDF驗證(單一流量)
    #plot(qs,cumsum(condition.pdf),type="l",ylim=c(0,1),xlab="輸砂量Qs(公噸)",ylab="probability")
    #cumsum(condition.pdf)[2500]
  }
}
# export table
if (export=="y"){
  colnames(pvalue.table)<-c("gumbelcopula","frankcopula","claytoncopula")
  file <- paste("E:/R_output/CHIA-YUANG/result/copula_gof_pvalue.csv", sep="")
  write.csv(pvalue.table,file)
}