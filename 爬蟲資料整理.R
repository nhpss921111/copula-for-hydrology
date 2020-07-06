# 爬蟲資料整理
# 開始撰寫日期：2020/07/05
# 完成撰寫日期：2020/07/05
#
rm(list=ls())
library(dplyr)
setwd("E:/R_reading/CHIA-YUANG")
data <- read.csv(file.path(getwd(),"2560H017.csv"), header = T)
data <- data[,-1]
data[1,2:7]
year <- c(1974:2019) # 請輸入年分：

new.data <- c()

for (n in c(1:length(year))){
  ini.data <- data[n,2:7]
  colnames(ini.data) <- c("No.","Month","Day","Q","Sc","Sl")
  i <- 1
  while (i <145 ){
    slice.data <- data[n,(i+1):(i+6)]
    colnames(slice.data) <- c("No.","Month","Day","Q","Sc","Sl")
    new.data <- rbind(new.data,slice.data)
    i <- i+6
    new.data <- arrange(new.data,No.) #依照No.欄位，由小到大排序
    
  }
}
