# Research
# 將流量Q之txt轉換成csv
# 撰寫日期：2020/11/10
# 完成日期：2020/11/10
# 事前準備：每年的流量
# Input："年分".txt
# Ouput："年分".csv
# ========================================
# 更改路徑就能做不同的測站
#station <- c("CHIA-YUANG") #家源橋
#year <- c(2001:2019) 
#
#station <- c("CHUNYUN BRIDGE") #彰雲橋
#year <- c(2001:2019) 
#
#station <- c("NEI-MAO-PU") #內茅埔
#year <- c(2001,2003:2019) 
# ========================================
rm(list = ls()) # 清除資料
station <- c("NEI-MAO-PU")
year <- c(2001,2003:2019) # 請輸入年分：
setwd(paste0("F:/R_reading/",station,"/from_pdf/discharge"))
input <- c(paste0(year,".txt")) 
output <- c(paste0(year,".csv")) 
for(i in c(1:length(input))){ 
  data <- read.table(file.path(getwd(),input[i]),fill=TRUE) #按年份讀取資料
  # 資料解析
  #data[data == 0] <- NA # 把"0"當成NA
  data <- cbind(year[i],data)
  colnames(data) <- c("年","日","一月","二月","三月","四月","五月","六月",
                      "七月","八月","九月","十月","十一","十二月")
  file <- paste0("F:/R_output/",station,"/discharge/", output[i], sep="") #存檔路徑
  write.csv(data,file) #結果寫到csv裡面
}
