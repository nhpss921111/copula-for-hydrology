# 爬蟲資料整理
# 開始撰寫日期：2020/07/05
# 完成撰寫日期：2020/07/06
#
rm(list=ls())
library(dplyr)
# 參數給定----
year <- c(1974:2000) # 請輸入年分：
input.file <- c(paste0(year,".csv"))
# -----------
setwd("E:/R_reading/CHIA-YUANG/") # 讀取檔案之路徑

for (n in c(1:length(year))){
  data <- read.csv(file.path(getwd(),input.file[n]), header = T) #按年份讀取資料
  data <- data[,-1]
  ini.data <- data[2:7]
  arrange.data <- c()
  i <- 1
  while (i <300){
    slice.data <- data[(i+1):(i+6)] #每一天的資料
    arrange.data <- rbind(arrange.data,slice.data)
    colnames(arrange.data) <- c("No.","Month","Day","Discharge",
                                "Sediment Content","Suspended Load")
    i <- i+6
  }
  rm.data <- arrange.data[complete.cases(arrange.data), ]
  new.data <- arrange(as.data.frame(rm.data),No.) #依照No.欄位，由小到大排序
  new.data[,1] <- year[n]
  colnames(new.data) <- c("Year","Month","Day","Discharge",
                          "Sediment Content","Suspended Load")
  write.csv(new.data, file = paste0("E:/R_output/CHIA-YUANG/",year[n], ".csv"))
}
