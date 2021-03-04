# 爬蟲資料整理
# 整理同時有流量和輸砂量的資料(使用crawling.R爬下來的資料)
# 開始撰寫日期：2020/07/05
# 完成撰寫日期：2021/03/04
# =============================================
# 更改路徑就能做不同的測站
#station <- c("CHIA-YUANG")
#station.num <- c("2560H017") #蘭陽溪-家源橋
#year <- c(1974:2000) # 1974才開始有SSL的紀錄
#
#station <- c("CHUNYUN BRIDGE")
#station.num <- c("1510H057") #濁水溪-彰雲橋
#year <- c(1987:2000) # 1987才開始有SSL的紀錄
#
#station <- c("NEI-MAO-PU")
#station.num <- c("1510H049") #濁水溪-內茅埔
#year <- c(1972:2000) # 1972才開始有SSL的紀錄
# =============================================
rm(list=ls())
library(dplyr)
# 參數給定----
station <- c("NEI-MAO-PU")
year <- c(1972:2000) # 請輸入年分：
input.file <- c(paste0(year,".csv"))
# -----------
setwd(paste0("F:/R_output/crawling/",station,"/discharge+SSL")) # 讀取檔案之路徑


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
  write.csv(new.data, file = paste0("F:/R_output/",station,"/SSL_arranged/",year[n], ".csv"))
}
