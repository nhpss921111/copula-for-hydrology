# Reseach
# 以"家源橋"測站為例，將所有流量資料合併
# 先輸入第一年資料
# 第二年開始用迴圈寫進來，直到最後一年
# 輸出csv格式
# 開始撰寫日期：2020/11/10
# 完成撰寫日期：2020/11/10
# By連育成

rm(list=ls())
year <- c(2013:2019) # 請輸入年分： 第2年資料到第n年資料
setwd("F:/R_reading/CHIA-YUANG/discharge")
ori.data <- read.csv(file.path(getwd(),"2012.csv"),stringsAsFactors = FALSE) # 第一年資料
new.data <- ori.data
input <- c(paste0(year,".csv")) 
for(i in c(1:length(input))){
  data <- read.csv(file.path(getwd(),input[i]),stringsAsFactors = FALSE)

  new.data <- rbind(new.data,data) # 把資料新增進來
  
}
new.data <- new.data[,-1]
new.data[new.data== "*" ] <- 0
file <- paste("F:/R_output/CHIA-YUANG/discharge/", "2012-2019.csv", sep="") #存檔路徑
write.csv(new.data,file) #結果寫到csv裡面
