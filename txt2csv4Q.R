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
#
#station <- c("JEN-SHOU BRIDGE") #花蓮溪-仁壽橋
#year <- c(2001:2019) # 
#
#station <- c("LIU-KWEI")#高屏溪-六龜
#year <- c(2001:2009,2011:2019) 
#
#station <- c("NEI-WAN")# 頭前溪-內灣
#year <- c(2001:2019) 
#
#station <- c("JEIN-KUO BRIDGE")頭前溪-經國橋
#year <- c(2001:2019) 
#
#station <- c("LAN-YANG BRIDGE")蘭陽溪-蘭陽大橋
#year <- c(2002:2019)  
#
#station <- c("HENG CHI") #三峽河-橫溪
#year <- c(2001:2005,2007:2019) 
#
#station <- c("HSIU-LUNG") #新店溪-秀朗
#year <- c(2001,2004,2005,2007:2019) 
# ========================================
rm(list = ls()) # 清除資料

# =============== 參數設定區 ======================================
station <- c("JEN-SHOU BRIDGE") #花蓮溪-仁壽橋
year <- c(2001:2019) 
input_file_path <- paste0("F:/copula/",station,"/from_pdf/discharge")
output_file_path <- paste0("F:/copula/",station,"/discharge/")
# ==================================================================
input <- c(paste0(year,".txt")) 
output <- c(paste0(year,".csv")) 
for(i in c(1:length(input))){ 
  setwd(input_file_path)
  data <- read.table(file.path(getwd(),input[i]),fill=TRUE) #按年份讀取資料
  # 資料解析
  #data[data == 0] <- NA # 把"0"當成NA
  data <- cbind(year[i],data)
  colnames(data) <- c("年","日","一月","二月","三月","四月","五月","六月",
                      "七月","八月","九月","十月","十一","十二月")
  file <- paste0(output_file_path, output[i], sep="") #存檔路徑
  write.csv(data,file) #結果寫到csv裡面
}
