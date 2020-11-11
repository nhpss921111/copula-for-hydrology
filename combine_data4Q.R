# Reseach
# 把流量資料中"*"轉換成"0"
# 開始撰寫日期：2020/11/10
# 完成撰寫日期：2020/11/10
# By連育成

rm(list=ls())
year <- c(2012:2019) # 請輸入年分： 第2年資料到第n年資料
setwd("F:/R_reading/CHIA-YUANG/discharge")
input <- c(paste0(year,".csv")) 
for(i in c(1:length(input))){
  data <- read.csv(file.path(getwd(),input[i]),stringsAsFactors = FALSE)
  data <- data[,-1]
  data[data== "*" ] <- 0
  file <- paste("F:/R_output/CHIA-YUANG/discharge/", input[i], sep="") #存檔路徑  
  write.csv(data,file) #結果寫到csv裡面
}

