# Missing data impution 驗證
# 開始撰寫日期：2020/12/24
# 完成撰寫日期：2021/01/14
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
library(locfit) # local polynomial estimation
#library(Metrics) # for RMSE
library(nls2)
library(devtools)
#library(aomisc)
#library(Metrics) # 評估指標
# ===========
# 家源橋CHIA-YUANG year <- c(1974:2009,2012:2019)
# 彰雲橋CHUNYUN BRIDGE year <- c(1987:2019)
# ===========
station <- c("CHIA-YUANG") # 測站名稱
station_ch <-c("家源橋") 
group.number <- c(9) # 分組的組數
perc.mis    <- 0.3 # 多少%的資料當成NA 
year <- c(1974:2009,2012:2019) # 年分
MD.input <- c(paste0(year,"QandQs.csv"))
output <- c(paste0(year,"imp.csv"))
ob.data <- c()
seednum <- 100
e <- function(actual, predicted) { # 相對誤差
  (actual - predicted)/actual
}
mse <- function(actual, predicted) { # 均方誤差
  mean((actual - predicted) ^ 2)
}
rmse <- function(actual, predicted) { # 均方根誤差
  sqrt(mean((actual - predicted) ^ 2))
}
nmse <- function(actual, predicted) { # 正歸化均方根誤差
  mean((actual - predicted) ^ 2)/mean((actual - mean(actual)) ^ 2)
}
mape <- function(actual, predicted) { # 平均絕對百分誤差
  100*mean(abs((actual - predicted)/actual))
}

#while (seednum<101){
  set.seed(100)
  
  #
  # 主要迴圈(以年份為底)
  for( y in 1:length(year)){
    # 讀取有缺失的資料
    setwd(paste0("F:/R_reading/",station,"/missingdata"))
    data <- read.csv(file.path(getwd(),MD.input[y]))
    data <- data[,-1]
    ob.data <- rbind(ob.data,data)
  }
  ob.data <- subset(ob.data,ob.data[,4]>0) # 把流量為0(無觀測資料)刪除
  log.data <- cbind(ob.data[,1:3],log10(ob.data[,4:5]))
  rm.ob.data <- ob.data[complete.cases(ob.data), ] # 移除原始觀測資料中全部NA 
  rm.log.data <- log.data[complete.cases(log.data), ] # 移除觀測資料取對數後全部NA 
  #ob.data <- subset(ob.data,ob.data[,4]>20 & ob.data[,5]>1000)
  #
  # ---- 所有觀測資料(原始) ----
  #
  ggplot(data=ob.data,aes(x=Discharge,y=Suspended.Load))+
    geom_point()+
    scale_color_discrete(name="年",labels=c(""))+
    labs(x="流量(cms)",y="輸砂量Qs (公噸)") + # 座標軸名稱
    #geom_mark_ellipse(expand = 0,aes(fill=group))+ # 橢圓形圈
    #geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=group))+ # 多邊形
    #geom_mark_hull(expand=0.01,aes(fill=group))+ # 凹凸多邊形
    ggtitle(paste0(year[1],"到",year[y],"年",station_ch,"測站觀測資料")) +
    theme_bw() + # 白底
    theme(panel.grid.major = element_blank()) + # 隱藏主要格線
    theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
    theme(text=element_text(size=30))  # 字體大小
  #
  # ---- 所有觀測資料(對數) ----
  #
  ggplot(data=log.data,aes(x=Discharge,y=Suspended.Load))+
    geom_point()+
    scale_color_discrete(name="年",labels=c(""))+
    labs(x="log10 流量(cms)",y="log10 輸砂量Qs (公噸)") + # 座標軸名稱
    #geom_mark_ellipse(expand = 0,aes(fill=group))+ # 橢圓形圈
    #geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=group))+ # 多邊形
    #geom_mark_hull(expand=0.01,aes(fill=group))+ # 凹凸多邊形
    ggtitle(paste0(year[1],"到",year[y],"年",station_ch,"測站觀測資料(對數)")) +
    theme_bw() + # 白底
    theme(panel.grid.major = element_blank()) + # 隱藏主要格線
    theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
    theme(text=element_text(size=30))  # 字體大小
  #
  # ------- 分九組
  if (group.number==9){
    #  -------------- 決定流量分組範圍 (來源：ob.data) ---------------
    # group 1： 0% ~  20%
    # group 2：20% ~  40%
    # group 3：40% ~  60%
    # group 4：60% ~  80%
    # group 5：80% ~  90%
    # group 6：90% ~  95%
    # group 8：95% ~  98%
    # group 8：98% ~  99%
    # group 9：99% ~ 100%
    
    rank.data <- cbind(log.data,rank(log.data$Discharge))
    persent <- (rank.data$`rank(log.data$Discharge)`) / length(rank.data$Discharge)
    per.data <- cbind(rank.data,persent)
    
    data.1 <- data.frame(subset(per.data, persent<=0.2),group="group1")
    data.2 <- data.frame(subset(per.data, persent>0.2 & persent<=0.4),group="group2")
    data.3 <- data.frame(subset(per.data, persent>0.4 & persent<=0.6),group="group3")
    data.4 <- data.frame(subset(per.data, persent>0.6 & persent<=0.8),group="group4")
    data.5 <- data.frame(subset(per.data, persent>0.8 & persent<=0.9),group="group5")
    data.6 <- data.frame(subset(per.data, persent>0.9 & persent<=0.96),group="group6")
    data.7 <- data.frame(subset(per.data, persent>0.96 & persent<=0.98),group="group7")
    data.8 <- data.frame(subset(per.data, persent>0.98 & persent<=0.99),group="group8")
    data.9 <- data.frame(subset(per.data, persent>0.99),group="group9")
    
    data.group <- rbind(data.1, data.2, data.3,
                        data.4, data.5, data.6,
                        data.7, data.8, data.9)
    
    group.BC <- c(0,max(data.1$Discharge),max(data.2$Discharge),max(data.3$Discharge),
                  max(data.4$Discharge),max(data.5$Discharge),max(data.6$Discharge),
                  max(data.7$Discharge),max(data.8$Discharge),max(data.9$Discharge))
    persent.BC <- c(0,20,40,60,80,90,95,98,99,100)
  }
  
  if (group.number==5){
    #  -------------- 決定流量分組範圍 (來源：log.data) ---------------
    # group 1： 0% ~  20%
    # group 2：20% ~  40%
    # group 3：40% ~  60%
    # group 4：60% ~  80%
    # group 5：80% ~  100%
    
    rank.data <- cbind(log.data,rank(log.data$Discharge))
    persent <- (rank.data$`rank(log.data$Discharge)`) / length(rank.data$Discharge)
    per.data <- cbind(rank.data,persent)
    
    data.1 <- data.frame(subset(per.data, persent<=0.2),group="group1")
    data.2 <- data.frame(subset(per.data, persent>0.2 & persent<=0.4),group="group2")
    data.3 <- data.frame(subset(per.data, persent>0.4 & persent<=0.6),group="group3")
    data.4 <- data.frame(subset(per.data, persent>0.6 & persent<=0.8),group="group4")
    data.5 <- data.frame(subset(per.data, persent>0.8),group="group5")
    
    data.group <- rbind(data.1, data.2, data.3,
                        data.4, data.5)
    
    group.BC <- c(0,max(data.1$Discharge),max(data.2$Discharge),max(data.3$Discharge),
                  max(data.4$Discharge),max(data.5$Discharge))
    persent.BC <- c(0,20,40,60,80,100)
  }
  
  
  #
  # ----------- 將同時有Q與Qs的資料分兩組 (80%資料總數建模，剩下20%當成驗證) -------------
  #
  
  x.samp <- as.matrix(rm.log.data[,4:5])
  miss.row    <- sample(1:length(rm.log.data$Discharge), perc.mis*length(rm.log.data$Discharge), replace=FALSE)
  miss.col    <- rep(2,perc.mis*length(rm.log.data$Discharge))
  miss        <- cbind(miss.row,miss.col)
  samp.miss <- replace(x.samp,miss,NA) # NA的欄位座標
  MD <- cbind(rm.log.data[,1:3],samp.miss) # 將 ?% 觀測資料轉換成NA
  
  MD.rmNA <- MD[complete.cases(MD), ] # 移除全部NA (剩餘資料)
  
  # ----- 建模(同時有Q與Qs)的80%觀測資料 -------
  ggplot(data=MD.rmNA,aes(x=Discharge,y=Suspended.Load))+
    geom_point()+
    scale_color_discrete(name="年",labels=c(""))+
    labs(x="流量(cms)",y="輸砂量Qs (公噸)") + # 座標軸名稱
    #geom_mark_ellipse(expand = 0,aes(fill=group))+ # 橢圓形圈
    #geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=group))+ # 多邊形
    #geom_mark_hull(expand=0.01,aes(fill=group))+ # 凹凸多邊形
    ggtitle(paste0(year[1],"到",year[y],"年",station_ch,"測站建立模型")) +
    theme_bw() + # 白底
    theme(panel.grid.major = element_blank()) + # 隱藏主要格線
    theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
    theme(text=element_text(size=30))  # 字體大小
  #
  # ------------------------------- imputation_validation ----------------------------------------
  # =========== 不分組 (建立 rating curve)===============
  file <- paste("F:/R_output/",station,"/imputation_validation(",perc.mis,"NA)(",group.number,"groups)/",
                year[1],"到", year[y],"年不分組率定曲線初始值查找.csv", sep="") #存檔路徑
  write.csv(MD.rmNA,file)
  # 從"率定曲線初始值.csv"找到a,b的初始值再帶入
  # SSpower <- selfStart(~ A*x^B,
  #                      function(mCall, data, LHS)
  #                      {
  #                        xy <- sortedXyData(mCall[["x"]], LHS, data)
  #                        if(nrow(xy) < 3) {
  #                          stop("Too few distinct x values to fit a power function")
  #                        }
  #                        z <- xy[["y"]]
  #                        xy[["logx"]] <- log(xy[["x"]])     
  #                        xy[["logy"]] <- log(xy[["y"]])  
  #                        aux <- coef(lm(logy ~ logx, xy))
  #                        pars <- c(exp(aux[[1]]), aux[[2]])
  #                        setNames(pars,
  #                                 mCall[c("A", "B")])
  #                      }, c("A", "B"))
  # rating <- nls(Suspended.Load ~SSpower(Discharge,a,b),data=MD.rmNA)
  # rating <- nls(Suspended.Load ~ a*Discharge^b, algorithm="port",
  #               control = list(maxiter = 200,minFactor = 1/2^300,warnOnly = TRUE),
  #               start=list(a=1,b=1), data=MD.rmNA,trace=T) # 初始值要給好!!!
  rating <- lm(Suspended.Load ~ Discharge,data=MD.rmNA)
  summary(rating) # 觀察
  log10.a <- rating$coefficients[1] # 迴歸係數log10(a)
  b <- rating$coefficients[2] # 迴歸係數b
  rating.par <- cbind(log10.a,b)
  file <- paste("F:/R_output/",station,"/imputation_validation(",perc.mis,"NA)(",group.number,"groups)/",
                year[1],"到", year[y],"年不分組率定曲線係數.csv", sep="") #存檔路徑
  write.csv(rating.par,file)
  ratingSSL <- log10.a + b* (MD$Discharge) # 計算率定曲線推估的輸砂量
  rating.all <- cbind(rm.log.data,MD[,5],ratingSSL)
  colnames(rating.all) <- c("Year","Month","Day","Discharge",
                           "Suspended.Load","asNA","ratingSSL_all")
  rmse <- function(actual, predicted) { # 均方根誤差公式
    sqrt(mean((actual - predicted) ^ 2))
  }
  # 計算損失函數
  rating.all_vs_observation <- rating.all[(is.na(rating.all$asNA)),] #保留 na
  # rmse (root mean square error)
  
  rmse.rating <- rmse(rating.all_vs_observation$Suspended.Load, rating.all_vs_observation$ratingSSL_all) #原本移除的觀測值與率定值
  
  # nmse (normalize mean squared error)
  
  nmse.rating <- nmse(rating.all_vs_observation$Suspended.Load, rating.all_vs_observation$ratingSSL_all) #原本移除的觀測值與率定值
  # mse (mean absolute error)
  
  mse.rating <- mse(rating.all_vs_observation$Suspended.Load, rating.all_vs_observation$ratingSSL_all) #原本移除的觀測值與率定值
  
  # mape (mean absolute persentage error)
  
  mape.rating <- mape(rating.all_vs_observation$Suspended.Load, rating.all_vs_observation$ratingSSL_all) #原本移除的觀測值與率定值
  
  rating.table <- cbind(rmse.rating,
                        nmse.rating,
                        mse.rating,
                        mape.rating) #沒分組的率定曲線vs觀測資料
  # 率定曲線(不分組)資料出圖
  setwd(paste0("F:/R_output/",station,"/imputation_validation(",perc.mis,"NA)(",group.number,"groups)")) # 請修改儲存路徑：
  png(paste0(year[1],"到",year[y],"年",station_ch,"測站率定曲線(不分組).png"),width = 1250, height = 700, units = "px", pointsize = 12)
  Rating.all <- ggplot(data=rating.all)+
    geom_point(aes(x=Discharge,y=asNA,color="觀測資料")) +
    geom_line(aes(x=Discharge,y=ratingSSL,color="率定曲線")) +
    scale_color_discrete(name="圖例") + #圖例名稱
    ggtitle(paste0(year[1],"到",year[y],"年",station_ch,"測站率定曲線(不分組)")) +
    theme(text=element_text(size=20))  # 字體大小
  plot(Rating.all)
  dev.off()

  
  # =========== 分九組 or 分五組 (source：MD)===============
  # 分組數
  
  imp.table <- c() # 輸砂量總表
  copula.parameter.table <- c() # 各組聯結函數及其參數
  loss.table <- c() # 各組損失函數比較
  rating.par.table <- c() #各組率定曲線係數
  # ---- 分組迴圈  ----
  # 將一開始移除的 % 加回來
  for (g in 1:(length(group.BC)-1)){
    #maxBC1 <- max(data.1$Discharge)
    print(paste0("第",g,"組開始補遺"))
    MD.bygroup <- as.matrix(subset(MD[,4:5],Discharge>group.BC[g] & Discharge<=group.BC[g+1])) #提取Q與Qs出來，並限制流量範圍
    colnames(MD.bygroup) <- c("Discharge","Suspended.Load")
    rm.MD.bygroup <- MD.bygroup[complete.cases(MD.bygroup), ] # 移除全部NA 
    #upper <- c(max(rm.MD.bygroup[,1]),max(rm.MD.bygroup[,2]))
    #lower <- c(min(rm.MD.bygroup[,1]),min(rm.MD.bygroup[,2]))
    n.marg <- 2 # 兩個變數(Q、Qs)
    f <- rep(0.9,n.marg)
    imp <- CoImp(MD.bygroup, n.marg=n.marg,smoothing = c(1,1),
                  plot=F, q.lo=c(0.01,0.01), q.up=c(0.99,0.99),
                 model=list(gumbelCopula(),frankCopula(),claytonCopula())) # 補遺計算
    # 邊際分布
    setwd(paste0("F:/R_output/",station,"/imputation_validation(",perc.mis,"NA)(",group.number,"groups)")) # 請修改儲存路徑：
    png(paste0(year[1],"到",year[y],"年",station_ch,"測站group",g,"邊際分布.png"),width = 1250, height = 700, units = "px", pointsize = 12)
    plot(imp)
    dev.off()
    
    copula.func <- imp@Estimated.Model.Imp[["model"]]
    copula.para <- imp@Estimated.Model.Imp[["parameter"]]
    copula.list <- cbind(copula.func, copula.para)
    rownames(copula.list) <- paste0("group",g)
    copula.parameter.table <- rbind(copula.parameter.table,copula.list) # 傳承
    
    imp.group <- cbind(subset(rm.log.data,Discharge>group.BC[g] & Discharge<=group.BC[g+1]),
                   subset(MD[,4:5],Discharge>group.BC[g] & Discharge<=group.BC[g+1])[,2],
                   imp@Imputed.data.matrix[,2]) # 提取補遺值，並合併成結果表格
    
    colnames(imp.group) <- c("Year","Month","Day","Discharge",
                             "Suspended.Load","asNA","Imp.SSL") # 為框架的行命名
    require(minpack.lm)
    rating <- nlsLM(Suspended.Load ~ a*Discharge^b, 
                    data=as.data.frame(MD.bygroup[complete.cases(MD.bygroup), ]),
                    start = list(a=1,b=1))
    # rating <- nls(Suspended.Load ~SSpower(Discharge,a,b),
    #               data=as.data.frame(MD.bygroup[complete.cases(MD.bygroup), ]),
    #               control = list(maxiter = 5000, minFactor = 1/2^20, warnOnly = T),
    #               trace=T)
    # rating <- nls(Suspended.Load ~ a*Discharge^b, algorithm="port",
    #               control = list(maxiter = 200, minFactor = 1/2^300, warnOnly = F),weights,
    #               start=list(a=0.1,b=2), data=as.data.frame(MD.bygroup[complete.cases(MD.bygroup), ]),
    #               trace=T) # 乘冪迴歸式
    rating <- lm(Suspended.Load ~ Discharge,data=as.data.frame(MD.bygroup[complete.cases(MD.bygroup), ]))
    summary(rating) # 觀察
    log10.a <- rating$coefficients[1] # 迴歸係數log10(a)
    b <- rating$coefficients[2] # 迴歸係數b
    rating.par.group <- cbind(log10.a,b)
    rownames(rating.par.group) <- paste0("group",g)
    rating.par.table <- rbind(rating.par.table,rating.par.group)
    ratingSSL <- log10.a + b *(imp.group$Discharge) # 計算率定曲線推估的輸砂量
    imp.group <- cbind(imp.group,ratingSSL)
    colnames(imp.group) <- c("Year","Month","Day","Discharge","Suspended.Load",
                             "asNA","Imp.SSL","ratingSSL_group")
    imp.table <- rbind(imp.table,imp.group) # 傳承
   
    # 計算損失函數
    imp_vs_rating.group <- imp.group[(is.na(imp.group$asNA)),] #保留 na
    
    # rmse (root mean square error)
    rmse.imp <- rmse(imp_vs_rating.group$Suspended.Load, imp_vs_rating.group$Imp.SSL) #原本移除的觀測值與補遺值
    rmse.rating <- rmse(imp_vs_rating.group$Suspended.Load, imp_vs_rating.group$ratingSSL_group) #原本移除的觀測值與率定值
    
    # mse (mean absolute error)
    mse.imp <- mse(imp_vs_rating.group$Suspended.Load, imp_vs_rating.group$Imp.SSL) #原本移除的觀測值與補遺值
    mse.rating <- mse(imp_vs_rating.group$Suspended.Load, imp_vs_rating.group$ratingSSL_group) #原本移除的觀測值與率定值
    
    # mape (mean absolute persentage error)
    mape.imp <- mape(imp_vs_rating.group$Suspended.Load, imp_vs_rating.group$Imp.SSL) #原本移除的觀測值與補遺值
    mape.rating <- mape(imp_vs_rating.group$Suspended.Load, imp_vs_rating.group$ratingSSL_group) #原本移除的觀測值與率定值
    
    # loss function
    loss <- cbind(mse.imp, mse.rating,mape.imp, mape.rating,rmse.imp, rmse.rating)
    rownames(loss) <- paste0("group",g)
    loss.table <- rbind(loss.table,loss)
    # 補遺資料出圖
    setwd(paste0("F:/R_output/",station,"/imputation_validation(",perc.mis,"NA)(",group.number,"groups)")) # 請修改儲存路徑：
    png(paste0(year[1],"到",year[y],"年",station_ch,"測站補遺驗證group",g,".png"),width = 1250, height = 700, units = "px", pointsize = 12)
    Imp.group <- ggplot(data=imp.group)+
      geom_point(aes(x=Discharge,y=Imp.SSL,color="補遺值")) +
      geom_point(aes(x=Discharge,y=Suspended.Load,color="觀測值")) +
      geom_line(aes(x=Discharge,y=ratingSSL_group,color="率定曲線(各組)")) +
      scale_color_discrete(name="圖例") + #圖例名稱
      ggtitle(paste0(year[1],"到",year[y],"年",station_ch,"測站補遺驗證",persent.BC[g],"% ~ ",persent.BC[g+1],"%")) +
      theme(text=element_text(size=20))  # 字體大小
    plot(Imp.group)
    dev.off()
  }
  print("分組完成")
  # ------------- 把分組合併 -------------------
  
  # 組數合併資料出圖
  setwd(paste0("F:/R_output/",station,"/imputation_validation(",perc.mis,"NA)(",group.number,"groups)")) # 請修改儲存路徑：
  png(paste0(year[1],"到",year[y],"年",station_ch,"測站補遺驗證(all group).png"),width = 1250, height = 700, units = "px", pointsize = 12)
  Imp.sum <- ggplot(data=imp.table)+
    geom_point(aes(x=Discharge,y=Imp.SSL,color="補遺值"))+
    geom_point(aes(x=Discharge,y=Suspended.Load,color="觀測值"))+
    geom_line(aes(x=Discharge,y=ratingSSL_group,color="率定曲線(合併各組)")) + #
    #geom_vline(xintercept=group.BC) +
    scale_color_discrete(name="圖例") + #圖例名稱
    ggtitle(paste0(year[1],"到",year[y],"年",station_ch,"測站補遺驗證(all group)"))+
    theme(text=element_text(size=20))  # 字體大小
  plot(Imp.sum)
  dev.off()
  
  # 率定曲線不分組 vs 率定曲線有分組
  result.table <- full_join(imp.table,rating.all)
  colnames(result.table) <- c("Year","Month","Day","Discharge",
                          "Suspended.Load","asNA","Imp.SSL","ratingSSL_group","ratingSSL_all")
  setwd(paste0("F:/R_output/",station,"/imputation_validation(",perc.mis,"NA)(",group.number,"groups)")) # 請修改儲存路徑：
  png(paste0(year[1],"到",year[y],"年",station_ch,"測站率定曲線(有分組 vs 沒分組).png"),width = 1250, height = 700, units = "px", pointsize = 12)
  Loss.cf <- ggplot(data=result.table)+
    geom_point(aes(x=Discharge,y=Imp.SSL,color="補遺值"))+
    geom_point(aes(x=Discharge,y=Suspended.Load,color="觀測值"))+
    geom_line(aes(x=Discharge,y=ratingSSL_group,color="率定曲線(有分組)")) + #
    geom_line(aes(x=Discharge,y=ratingSSL_all,color="率定曲線(不分組)")) + #
    scale_color_discrete(name="圖例") + #圖例名稱
    ggtitle(paste0(year[1],"到",year[y],"年",station_ch,"測站率定曲線(有分組 vs 沒分組)"))+
    theme(text=element_text(size=20))  # 字體大小
  plot(Loss.cf)
  dev.off()
  #
  # ------ 計算損失函數 (補遺 and 率定分組)--------------
  imp_vs_rating.allgroup <- result.table[(is.na(result.table$asNA)),] #保留 na
  
  # rmse (root mean square error)
  
  rmse.rating <- rmse(imp_vs_rating.allgroup$Suspended.Load, imp_vs_rating.allgroup$ratingSSL_group) #原本移除的觀測值與率定值
  rmse.imp <- rmse(imp_vs_rating.allgroup$Suspended.Load, imp_vs_rating.allgroup$Imp.SSL) #原本移除的觀測值與補遺值
  
  # nmse (normalize mean square error)
  
  nmse.rating <- nmse(imp_vs_rating.allgroup$Suspended.Load, imp_vs_rating.allgroup$ratingSSL_group) #原本移除的觀測值與率定值
  nmse.imp <- nmse(imp_vs_rating.allgroup$Suspended.Load, imp_vs_rating.allgroup$Imp.SSL) #原本移除的觀測值與補遺值
  
  # mse (mean absolute error)
  
  mse.rating <- mse(imp_vs_rating.allgroup$Suspended.Load, imp_vs_rating.allgroup$ratingSSL_group) #原本移除的觀測值與率定值
  mse.imp <- mse(imp_vs_rating.allgroup$Suspended.Load, imp_vs_rating.allgroup$Imp.SSL) #原本移除的觀測值與補遺值
  
  # mape (mean absolute persentage error)
  
  mape.rating <- mape(imp_vs_rating.allgroup$Suspended.Load, imp_vs_rating.allgroup$ratingSSL_group) #原本移除的觀測值與率定值
  mape.imp <- mape(imp_vs_rating.allgroup$Suspended.Load, imp_vs_rating.allgroup$Imp.SSL) #原本移除的觀測值與補遺值
  
  rating.bygroup <- cbind(rmse.rating,nmse.rating,mse.rating,mape.rating)
  imp.total <- cbind(rmse.imp,nmse.imp,mse.imp,mape.imp)
  loss.cf <- rbind(rating.table,rating.bygroup,imp.total)
  rownames(loss.cf) <- c("rating不分組","rating有分組","copula補遺")
  # ------- 結果儲存 -------
  file <- paste("F:/R_output/",station,"/imputation_validation(",perc.mis,"NA)(",group.number,"groups)/",
                year[1],"到", year[y],"年各組率定曲線參數.csv", sep="") #存檔路徑
  write.csv(rating.par.table,file)
  
  file <- paste("F:/R_output/",station,"/imputation_validation(",perc.mis,"NA)(",group.number,"groups)/",
                year[1],"到", year[y],"年率定曲線vs補遺法(取對數).csv", sep="") #存檔路徑
  write.csv(loss.cf,file)
  
  file <- paste("F:/R_output/",station,"/imputation_validation(",perc.mis,"NA)(",group.number,"groups)/",
                year[1],"到", year[y],"年",station_ch,"流量輸砂量總表.csv", sep="") #存檔路徑
  write.csv(result.table,file)
  
  file <- paste("F:/R_output/",station,"/imputation_validation(",perc.mis,"NA)(",group.number,"groups)/",
                year[1],"到", year[y],"年各組聯結函數&參數.csv", sep="") #存檔路徑
  write.csv(copula.parameter.table,file)
  
  file <- paste("F:/R_output/",station,"/imputation_validation(",perc.mis,"NA)(",group.number,"groups)/",
                year[1],"到", year[y],"年各組損失函數比較.csv", sep="") #存檔路徑
  write.csv(loss.table,file)

#}


setwd(paste0("F:/R_output/",station,"/imputation_validation(",perc.mis,"NA)(",group.number,"groups)")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站補遺驗證(all-對數值).png"),width = 1250, height = 700, units = "px", pointsize = 12)
All.cf <- ggplot(data=result.table)+
  geom_point(aes(x=Discharge,y=Imp.SSL,color="補遺值"))+
  geom_point(aes(x=Discharge,y=Suspended.Load,color="觀測值"))+
  geom_line(aes(x=Discharge,y=ratingSSL_group,color="率定曲線(有分組)")) + #
  geom_line(aes(x=Discharge,y=ratingSSL_all,color="率定曲線(不分組)")) + #
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年",station_ch,"測站補遺驗證(all-對數值)"))+
  theme(text=element_text(size=20))  # 字體大小
plot(All.cf)
dev.off()
#
# ------- 返回成原本尺度 ------
#
final.table <- cbind(result.table[,1:3],10^result.table[,4:9])
final.cf.table <- final.table[(is.na(final.table$asNA)),] #保留 na
# rmse (root mean square error)
rmse.rating_a <- rmse(final.cf.table$Suspended.Load, final.cf.table$ratingSSL_all) #原本移除的觀測值與率定值
rmse.rating_g <- rmse(final.cf.table$Suspended.Load, final.cf.table$ratingSSL_group) #原本移除的觀測值與率定值
rmse.imp <- rmse(final.cf.table$Suspended.Load, final.cf.table$Imp.SSL) #原本移除的觀測值與補遺值

# nmse (normalize mean square error)
nmse.rating_a <- nmse(final.cf.table$Suspended.Load, final.cf.table$ratingSSL_all) #原本移除的觀測值與率定值
nmse.rating_g <- nmse(final.cf.table$Suspended.Load, final.cf.table$ratingSSL_group) #原本移除的觀測值與率定值
nmse.imp <- nmse(final.cf.table$Suspended.Load, final.cf.table$Imp.SSL) #原本移除的觀測值與補遺值

# mse (mean absolute error)
mse.rating_a <- mse(final.cf.table$Suspended.Load, final.cf.table$ratingSSL_all) #原本移除的觀測值與率定值
mse.rating_g <- mse(final.cf.table$Suspended.Load, final.cf.table$ratingSSL_group) #原本移除的觀測值與率定值
mse.imp <- mse(final.cf.table$Suspended.Load, final.cf.table$Imp.SSL) #原本移除的觀測值與補遺值

# mape (mean absolute persentage error)
mape.rating_a <- mape(final.cf.table$Suspended.Load, final.cf.table$ratingSSL_all) #原本移除的觀測值與率定值
mape.rating_g <- mape(final.cf.table$Suspended.Load, final.cf.table$ratingSSL_group) #原本移除的觀測值與率定值
mape.imp <- mape(final.cf.table$Suspended.Load, final.cf.table$Imp.SSL) #原本移除的觀測值與補遺值

rating.allgroup <- cbind(rmse.rating_a,nmse.rating_a,mse.rating_a,mape.rating_a)
rating.bygroup <- cbind(rmse.rating_g,nmse.rating_g,mse.rating_g,mape.rating_g)
imp.bygroup <- cbind(rmse.imp,nmse.imp,mse.imp,mape.imp)
loss.cf <- rbind(rating.allgroup,rating.bygroup ,imp.bygroup)
rownames(loss.cf) <- c("rating不分組","rating有分組","copula補遺")
colnames(loss.cf) <- c("RMSE","NMES","MES","MAPE")
file <- paste("F:/R_output/",station,"/imputation_validation(",perc.mis,"NA)(",group.number,"groups)/",
              year[1],"到", year[y],"年率定曲線vs補遺法(實際值).csv", sep="") #存檔路徑
write.csv(loss.cf,file)

setwd(paste0("F:/R_output/",station,"/imputation_validation(",perc.mis,"NA)(",group.number,"groups)")) # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年",station_ch,"測站補遺驗證(all-實際值).png"),width = 1250, height = 700, units = "px", pointsize = 12)
All.cf <- ggplot(data=final.table)+
  geom_point(aes(x=Discharge,y=Imp.SSL,color="補遺值"))+
  geom_point(aes(x=Discharge,y=Suspended.Load,color="觀測值"))+
  geom_line(aes(x=Discharge,y=ratingSSL_group,color="率定曲線(有分組)")) + #
  geom_line(aes(x=Discharge,y=ratingSSL_all,color="率定曲線(不分組)")) + #
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年",station_ch,"測站補遺驗證(all-實際值)"))+
  theme(text=element_text(size=20))  # 字體大小
plot(All.cf)
dev.off()
