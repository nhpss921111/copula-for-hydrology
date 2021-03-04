# Reseach
# 把流量資料中"*"轉換成"0"
# 開始撰寫日期：2020/11/10
# 完成撰寫日期：2021/03/04
# By連育成
# ============================================
# 家源橋(CHIA-YUANG)：year <- c(1974:2009,2012:2019)
# 彰雲橋(CHUNYUN BRIDGE)：year <- c(1987:2019)
# 內茅埔(NEI-MAO-PU)：year <- c(1972:2001,2003:2019)
# ===========================================
rm(list=ls())
station <- c("NEI-MAO-PU")
year <- c(1972:2001,2003:2019) # 請輸入年分： 第2年資料到第n年資料
setwd(paste0("F:/R_reading/",station,"/discharge"))
input <- c(paste0(year,".csv")) 
for(i in c(1:length(input))){
  data <- read.csv(file.path(getwd(),input[i]),stringsAsFactors = FALSE)
  data <- data[,-1]
  data[data== "*" ] <- 0
  file <- paste0("F:/R_output/",station,"/discharge/", input[i]) #存檔路徑  
  write.csv(data,file) #結果寫到csv裡面
}

