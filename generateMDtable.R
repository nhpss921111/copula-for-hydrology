# Generate Missing Data Table 
# 開始撰寫時間：2020/11/10
# 完成撰寫時間：2021/03/08
# 將"流量資料"與"流量及輸砂量資料"合併
# ============================================
# 家源橋(CHIA-YUANG)：year <- c(1974:2009,2012:2019)
# 彰雲橋(CHUNYUN BRIDGE)：year <- c(1987:2019)
# 內茅埔(NEI-MAO-PU)：year <- c(1972:2001,2003:2019)
# 仁壽橋(JEN-SHOU BRIDGE)：year <- c(1960:2019)
# ===========================================
library(dplyr)
rm(list=ls())
station <- c("JEN-SHOU BRIDGE")
year <- c(1960:2019)
month <- c(2:12) # loop從2月開始跑
input <- c(paste0(year,".csv"))

for (y in 1:length(input)){ 
  setwd(paste0("F:/R_reading/",station,"/discharge"))
  data <- read.csv(file.path(getwd(),input[y])) # 第y年所有流量資料
  data <- data[,-1] # 整理多餘的行數
  data <- as.data.frame(data)
  data[is.na(data)] <- 0
  data[data == "*" ] <- 0
  # 先做完第一月
  Q <- data[,(2+1)] # 選取第1月的流量
  month.data <- as.data.frame(cbind(data[,1],1,data[,2],Q,0)) # 整理表格
  colnames(month.data) <- c("Year","Month","Day","Discharge","Suspended.Load")
  setwd(paste0("F:/R_reading/",station,"/SSL_arranged"))
  SSL <- read.csv(file.path(getwd(),paste0(1,"month.csv"))) # 讀取第1月流量及輸砂量資料
  SSL <- SSL[,-1] # 整理表格
  SSL <- SSL[,-5] # 整理表格
  subSSL <- subset(SSL,Year==year[y]) # 選取y年第1月的Q與Qs資料
  month.data$Discharge <- replace(month.data$Discharge,c(subSSL$Day),c(subSSL$Discharge))
  cb.data <- full_join(month.data, subSSL, by = c("Year","Month","Day","Discharge"))
  cb.data <- cb.data[,-5]
  colnames(cb.data) <- c("Year","Month","Day","Discharge","Suspended.Load")
  
  # 從第二月開始做
  for(m in month){
    Q <- data[,(2+m)] # 選取第m月的流量
    month.data <- as.data.frame(cbind(data[,1],m,data[,2],Q,0)) # 整理表格
    colnames(month.data) <- c("Year","Month","Day","Discharge","Suspended.Load")
    setwd(paste0("F:/R_reading/",station,"/SSL_arranged"))
    SSL <- read.csv(file.path(getwd(),paste0(m,"month.csv"))) # 讀取第m月流量及輸砂量資料
    SSL <- SSL[,-1] # 整理表格
    SSL <- SSL[,-5] # 整理表格
    subSSL <- subset(SSL,Year==year[y]) # 選取y年m月的Q與Qs資料
    month.data$Discharge <- replace(month.data$Discharge,c(subSSL$Day),c(subSSL$Discharge))
    cb.data.new <- full_join(month.data, subSSL, by = c("Year","Month","Day","Discharge"))
    cb.data.new <- cb.data.new[,-5]
    colnames(cb.data.new) <- c("Year","Month","Day","Discharge","Suspended.Load")
    cb.data <- rbind(cb.data,cb.data.new)
  }
  cb.data$Discharge[is.na(cb.data$Discharge)] <- 0
  file <- paste0("F:/R_output/",station,"/missingdata/", paste0(year[y],"QandQs.csv"), sep="") #存檔路徑  
  write.csv(cb.data,file) #結果寫到csv裡面
}