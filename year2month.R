# 年資料分成月資料
# 以濁水溪內茅埔測站為例
# 2010年~2019年
# Input：1974.csv, 1975.csv,...,2019.csv
# Output：1month.csv, 2month.csv,...,12month.csv
# 開始撰寫日期：2020/07/06
# 完成撰寫日期：2020/10/20
# By連育成

rm(list=ls())
setwd("E:/R_reading/NEI-MAO-PU/")
data <- as.data.frame(read.csv(file.path(getwd(),"2010-2019.csv"),header = T)) 
month <- c(1:12) #分組月份
output <- c(paste0(month,"month.csv"))
for(m in month){
  month.data <- subset(data,Month==m)
  file <- paste("E:/R_output/NEI-MAO-PU/", output[m], sep="") #存檔路徑
  write.csv(month.data,file) #結果寫到csv裡面
}
