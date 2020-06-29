# Research
# txt轉換成xlsx
# 撰寫日期：2020/06/15
# 事前準備：每年的輸砂量與流量
# Input："年分".txt
# Ouput："年分".xlsx

rm(list = ls()) # 清除資料
library(xlsx)
setwd("E:/R_reading/CHIA-YUANG/suspend")
year <- c(98:108) # 請輸入年分：
input <- c(paste0(year,".txt")) 
output <- c(paste0(year,".xlsx")) 
for(i in c(1:length(input))){ 
  data <- read.table(file.path(getwd(),input[i])) #按年份讀取資料
  # 資料解析
  sus1 <- data[,1:6] #把資料分成sus1和sus2
  sus2 <- data[,7:12]
  colnames(sus2) <- c(paste0("V",1:6))
  sus.table <- rbind(sus1,sus2) #合併sus1和sus2
  complete.cases(sus.table) # 回傳是否有na值
  sus.table <- sus.table[complete.cases(sus.table), ] #移除na值
  colnames(sus.table) <- c("No.","Month","Day","Discharge",
                         "Sediment Content","Suspended Load") #新增欄位名稱
  file <- paste("E:/R_output/CHIA-YUANG/", output[i], sep="") #存檔路徑
  write.xlsx(sus.table,file) #結果寫到xlsx裡面
}
