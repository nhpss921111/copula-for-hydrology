# Missing data impution
# 開始撰寫日期：2020/11/12
# 完成撰寫日期：2020/11/30
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

year <- c(2019) 
MD.input <- c(paste0(year,"QandQs.csv"))
output <- c(paste0(year,"imp.csv"))
MD <- c()

# 主要迴圈(以年份為底)
for( y in 1:length(year)){
  # 讀取有缺失的資料
  setwd("F:/R_reading/CHIA-YUANG/missingdata")
  data <- read.csv(file.path(getwd(),MD.input[y]))
  data <- data[,-1]
  #data$Discharge <- data$Discharge %>% log() # 取自然對數
  #data$Suspended.Load <- data$Suspended.Load %>% log() # 取自然對數
  MD <- rbind(MD,data)
}
MD <- subset(MD,MD[,4]>0) # 把流量為0(無觀測資料)刪除

#
# ------------------------------- Imputation ----------------------------------------
n.marg <- 2
samp.miss <- as.matrix(MD[,4:5])
imp <- CoImp(samp.miss, n.marg=n.marg, q.up=rep(0.99,n.marg), q.lo=rep(0.01,n.marg), 
             smoothing = c(0.2,0.2),type.data="continuous", 
             model=list(gumbelCopula(),frankCopula(),claytonCopula()))
show(imp)

setwd("F:/R_output/CHIA-YUANG/imputation") # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"年分布-0.2.png"),width = 1250, height = 700, units = "px", pointsize = 12)
plot(imp)
dev.off()

imp@Imputed.data.matrix[,2]
imp.data <- cbind(MD,imp@Imputed.data.matrix[,2])
imp.data <- imp.data[,-5]
colnames(imp.data) <- c("Year","Month","Day","Discharge","Imp.SSL")
# 資料還原
#imp.data$Discharge<- imp.data$Discharge %>% exp() # 取自然對數
#imp.data$Imp.SSL <- imp.data$Imp.SSL %>% exp() # 取自然對數
#MD$Discharge<- MD$Discharge %>% exp() # 取自然對數
#MD$Suspended.Load <- MD$Suspended.Load %>% exp() # 取自然對數

# 補遺資料出圖
setwd("F:/R_output/CHIA-YUANG/imputation") # 請修改儲存路徑：
png(paste0(year[1],"到",year[y],"-0.2.png"),width = 1250, height = 700, units = "px", pointsize = 12)
Imp.data <- ggplot()+
  geom_point(data=imp.data,aes(x=Discharge,y=Imp.SSL,color="補遺值"))+
  geom_point(data=MD,aes(x=Discharge,y=Suspended.Load,color="觀測值"))+
  scale_color_discrete(name="圖例") + #圖例名稱
  ggtitle(paste0(year[1],"到",year[y],"年家園橋測站補遺成果"))+
  theme(text=element_text(size=20))  # 字體大小
plot(Imp.data)
dev.off()
# 存檔
file <- paste("F:/R_output/CHIA-YUANG/imputation/",year[1],"到", year[y],"imp-0.2.csv", sep="") #存檔路徑  
write.csv(imp.data,file)


