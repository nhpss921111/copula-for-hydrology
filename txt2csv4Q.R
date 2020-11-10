# Research
# 將流量Q之txt轉換成csv
# 撰寫日期：2020/11/10
# 完成日期：2020/11/10
# 事前準備：每年的流量
# Input："年分".txt
# Ouput："年分".csv

rm(list = ls()) # 清除資料
setwd("F:/R_reading/CHIA-YUANG/discharge/")
year <- c(2012:2019) # 請輸入年分：
input <- c(paste0(year,".txt")) 
output <- c(paste0(year,".csv")) 
for(i in c(1:length(input))){ 
  data <- read.table(file.path(getwd(),input[i]),fill=TRUE) #按年份讀取資料
  # 資料解析
  #data[data == 0] <- NA # 把"0"當成NA
  data <- cbind(year[i],data)
  colnames(data) <- c("年","日","一月","二月","三月","四月","五月","六月",
                      "七月","八月","九月","十月","十一","十二月")
  file <- paste("F:/R_output/CHIA-YUANG/discharge/", output[i], sep="") #存檔路徑
  write.csv(data,file) #結果寫到xlsx裡面
}
