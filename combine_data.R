# Reseach
# 把憲宗老師提供資料(~2008年)更新到最新年分
# 以"家源橋"測站為例
# 先輸入第一年資料
# 第二年開始用迴圈寫進來，直到最後一年
# 輸出csv格式
# 開始撰寫日期：2020/06/15
# 完成撰寫日期：2020/07/06
# By連育成

rm(list=ls())
year <- c(2011:2019) # 請輸入年分： 第2年資料到第n年資料
setwd("E:/R_reading/NEI-MAO-PU/")
ori.data <- read.csv(file.path(getwd(),"2010.csv"),header = T) # 第一年資料
new.data <- ori.data 
input <- c(paste0(year,".csv")) 
for(i in c(1:length(input))){
  data <- read.csv(file.path(getwd(),input[i]),header = T)
  new.data <- rbind(new.data,data) # 把資料新增進來
}
new.data <- new.data[,-1]
file <- paste("E:/R_output/NEI-MAO-PU/", "2010-2019.csv", sep="") #存檔路徑
write.csv(new.data,file) #結果寫到csv裡面
