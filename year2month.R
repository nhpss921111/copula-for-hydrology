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
#year <- c(1974:2019) # 1974才開始有SSL的紀錄
#
#station <- c("CHUNYUN BRIDGE")
#station.num <- c("1510H057") #濁水溪-彰雲橋
#year <- c(1987:2019) # 1987才開始有SSL的紀錄
#
#station <- c("NEI-MAO-PU")
#station.num <- c("1510H049") #濁水溪-內茅埔
#year <- c(1972:2001,2003:2019) # 1972才開始有SSL的紀錄
#
#station <- c("JEN-SHOU BRIDGE")
#station.num <- c("2420H019") #花蓮溪-仁壽橋
#year <- c(1960:2019) # 1960才開始有SSL的紀錄
#
#station <- c("LIU-KWEI")
#station.num <- c("1730H039") #高屏溪-六龜
#year <- c(1982:2009,2011:2019) # 1982才開始有SSL的紀錄
#
#station <- c("NEI-WAN")
#station.num <- c("1300H013")# 頭前溪-內灣
#year <- c(1971:2005,2009:2019) # 1971才開始有SSL的紀錄
#
#station <- c("JEIN-KUO BRIDGE")
#station.num <- c("1300H017")# 頭前溪-經國橋
#year <- c(1990:2005,2009:2019) # 1990才開始有SSL的紀錄
#
#station.num <- c("2560H006")
#station <- c("LAN-YANG BRIDGE") # 蘭陽溪-蘭陽大橋
#year <- c(1949:2017,2019) # 1949才開始有SSL的紀錄
#
#station.num <- c("1140H049")
#station <- c("HENG CHI") #三峽河-橫溪
#year <- c(1974:2004,2006:2019) 
#
#station.num <- c("1140H066")
#station <- c("HSIU-LUNG") #新店溪-秀朗
#year <- c(1970:2004,2006:2019) 
#
#station.num <- c("1400H009")
#station <- c("I-LI") #大安溪-義里
#year <- c(1966:2003,2006:2010,2013:2019) 
#
#station.num <- c("1730H031")
#station <- c("LAO-NUNG")#荖濃(新發大橋)
#year <- c(1956:2009) 
#
#station <- c("LI-LIN BRIDGE")#里嶺大橋
#year <- c(1991:2004,2007:2019) 

rm(list=ls())
# ============== 參數設定區 ====================
station <- c("JEN-SHOU BRIDGE")
station.num <- c("2420H019") #花蓮溪-仁壽橋
year <- c(1960:2019) # 1960才開始有SSL的紀錄
input_file_path <- paste0("F:/copula/",station,"/discharge+SSL/")
output_file_path <- paste0("F:/copula/",station,"/discharge+SSL/")
# ==================================================
setwd(input_file_path)
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
  file <- paste0(output_file_path, output[m]) #存檔路徑
  write.csv(month.data,file) #結果寫到csv裡面
}