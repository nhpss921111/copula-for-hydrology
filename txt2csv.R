# Research
# from：txt2xlsx
# txt轉換成csv
# 撰寫日期：2020/07/06
# 完成日期：2020/12/21
# 事前準備：每年的輸砂量與流量
# 更改路徑就能做不同的測站
# Input："年分".txt
# Ouput："年分".csv
# CHIA-YUANG：家源橋
# CHUNYUN BRIDGE：彰雲橋

rm(list = ls()) # 清除資料
setwd("F:/R_reading/CHUNYUN BRIDGE/form_pdf/SSL/") #濁水溪-彰雲橋
year <- c(2001:2019) # 請輸入年分：
input <- c(paste0(year,".txt")) 
output <- c(paste0(year,".csv")) 
for(i in c(1:length(input))){ 
  data <- read.table(file.path(getwd(),input[i])) #按年份讀取資料
  # 資料解析
  sus1 <- data[,1:6] #把資料分成sus1和sus2
  sus2 <- data[,7:12]
  colnames(sus2) <- c(paste0("V",1:6))
  sus.table <- rbind(sus1,sus2) #合併sus1和sus2
  complete.cases(sus.table) # 回傳是否有na值
  sus.table <- sus.table[complete.cases(sus.table), ] #移除na值
  sus.table[,1] <- year[i]
  colnames(sus.table) <- c("Year","Month","Day","Discharge",
                         "Sediment Content","Suspended Load") #新增欄位名稱
  file <- paste("F:/R_output/CHUNYUN BRIDGE/SSL_arranged/", output[i], sep="") #存檔路徑
  write.csv(sus.table,file) #結果寫到csv裡面
}
