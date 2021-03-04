# (SSL)年資料分成月資料
# 以濁水溪內茅埔測站為例
# 2010年~2019年
# Input：1974.csv, 1975.csv,...,2019.csv
# Output：1month.csv, 2month.csv,...,12month.csv
# 開始撰寫日期：2020/07/06
# 完成撰寫日期：2021/03/04
# By連育成
# ============================================
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
# ===========================================
rm(list=ls())
station <- c("NEI-MAO-PU")
year <- c(1972:2001,2003:2019)
setwd(paste0("F:/R_reading/",station,"/SSL_arranged/"))
input <- c(paste0(year,".csv"))
data <- c() #放所有資料的表格
for(y in 1:length(year)){
  year.data <- read.csv(file.path(getwd(),input[y]),header = T)
  data <- rbind(data,year.data)
}
month <- c(1:12) #分組月份
output <- c(paste0(month,"month.csv"))
for(m in month){
  month.data <- subset(data,Month==m)
  month.data <- month.data[,-1]
  file <- paste0("F:/R_output/",station,"/SSL_arranged/",output[m]) #存檔路徑
  write.csv(month.data,file) #結果寫到csv裡面
}