# Research
# from：txt2xlsx
# txt轉換成csv
# 撰寫日期：2020/07/06
# 完成日期：2021/03/04
# 事前準備：每年的輸砂量與流量
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
#year <- c(2001:2019) # 1960才開始有SSL的紀錄
#
#station <- c("LIU-KWEI")#高屏溪-六龜
#year <- c(2001:2009,2011:2019) # 
#
#station <- c("NEI-WAN")# 頭前溪-內灣
#year <- c(2001:2005,2009:2019) 
#
#station <- c("JEIN-KUO BRIDGE")頭前溪-經國橋
#year <- c(2001:2005,2009:2019) # 1990才開始有SSL的紀錄
#
#station <- c("LAN-YANG BRIDGE")蘭陽溪-蘭陽大橋
#year <- c(2001:2017,2019)  
#
#station <- c("HENG CHI") #三峽河-橫溪
#year <- c(2001:2004,2006:2019)
#
#station <- c("HSIU-LUNG") #新店溪-秀朗
#year <- c(2001:2004,2006:2019) 
#
#station <- c("I-LI") #大安溪-義里
#year <- c(2001:2003,2006:2019) 
#
#station <- c("LAO-NUNG")#荖濃(新發大橋)
#year <- c(2001:2009) 
#
#station <- c("LI-LIN BRIDGE")#里嶺大橋
#year <- c(2001:2004,2007:2019) 
# ========================================
rm(list = ls()) # 清除資料
station <- c("LI-LIN BRIDGE")#里嶺大橋
year <- c(2001:2004,2007:2019) 
setwd(paste0("F:/R_reading/",station,"/from_pdf/discharge+SSL/")) 
input <- c(paste0(year,".txt")) 
output <- c(paste0(year,".csv")) 
for(i in c(1:length(input))){ 
  data <- read.table(file.path(getwd(),input[i]),fill=TRUE) #按年份讀取資料
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
  file <- paste(paste0("F:/R_output/",station,"/SSL_arranged/", output[i], sep="")) #存檔路徑
  write.csv(sus.table,file) #結果寫到csv裡面
}
