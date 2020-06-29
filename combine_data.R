# Reseach
# 把憲宗老師提供資料(~2008年)更新到最新年分
# 以"家源橋"測站為例
# 先輸入第一年資料
# 第二年開始用迴圈寫進來，直到最後一年
# 開始撰寫日期：2020/06/15
# 完成撰寫日期：2020/06/29
# By連育成

rm(list=ls())

library(xlsx)

year <- c(99:108) # 請輸入年分： 第2年資料到第n年資料
setwd("E:/R_reading/CHIA-YUANG/")
ori.data <- read.xlsx(file.path(getwd(),"98.xlsx"),sheetIndex="Sheet1", # 第一年資料                      startRow = 1,endRow = 1005,
                      header = T,colIndex =c(5,7),encoding = "UTF-8")
new.data <- ori.data 
colnames(new.data) <- c("Q","S")
input <- c(paste0(year,".xlsx")) 
for(i in c(1:length(input))){
  data <- read.xlsx(file.path(getwd(),input[i]),sheetIndex="Sheet1",
                    startRow = 1,endRow = 50,
                    header = T,colIndex =c(5,7),encoding = "UTF-8")
  colnames(data) <- c("Q","S") # 統一colnames
  new.data <- rbind(new.data,data) # 把資料新增進來
}
file <- paste("E:/R_output/CHIA-YUANG/", "2009-2019.xlsx", sep="") #存檔路徑
write.xlsx(new.data,file) #結果寫到xlsx裡面