# 來源：爬蟲，水文年報(2000年以前)  2019/11/02 BY奕廷學長
# 爬取水文年報之測站的日平均流量(cms)
# 全部抓下來後，再清洗資料
# 開始撰寫日期：2020/07/02
# 完成撰寫日期：2021/03/04
rm(list = ls())
library(rvest)
library(magrittr)
library(jsonlite)
output <- "y"
# ========================================
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
#
#station <- c("JEN-SHOU BRIDGE")
#station.num <- c("2420H019") #花蓮溪-仁壽橋
#year <- c(1960:2000) # 1960才開始有SSL的紀錄
#
#station <- c("LIU-KWEI")
#station.num <- c("1730H039") #高屏溪-六龜
#year <- c(1982:2000) # 1960才開始有SSL的紀錄
#
#station <- c("NEI-WAN")
#station.num <- c("1300H013")# 頭前溪-內灣
#year <- c(1971:2000) # 1971才開始有SSL的紀錄
#
#station <- c("JEIN-KUO BRIDGE")
#station.num <- c("1300H017")# 頭前溪-經國橋
#year <- c(1990:2000) # 1990才開始有SSL的紀錄
#
#station.num <- c("2560H006")
#station <- c("LAN-YANG BRIDGE") # 蘭陽溪-蘭陽大橋
#year <- c(1949:2000) # 1949才開始有SSL的紀錄
#
#station.num <- c("1140H049")
#station <- c("HENG CHI") #三峽河-橫溪
#year <- c(1974:2000) 
#
#station.num <- c("1140H066")
#station <- c("HSIU-LUNG") #新店溪-秀朗
#year <- c(1970:2000) 
#
#station.num <- c("1340H009")
#station <- c("YUN HSIN BRIDGE") #中港溪-永興橋
#year <- c(1986:2000) 
#
#station.num <- c("1400H009")
#station <- c("I-LI") #大安溪-義里
#year <- c(1966:2000) 
#
#station.num <- c("1730H031")
#station <- c("LAO-NUNG")#荖濃(新發大橋)
#year <- c(1956:2000) 
# ========================================
#----起始資料
#測站編號
station <- c("JEN-SHOU BRIDGE")
station.num <- c("2420H019") #花蓮溪-仁壽橋
year <- c(1960:2000) # 1960才開始有SSL的紀錄
#----
# 先建立錯誤函數：
# 1.空白值將錯誤，建立tryCatch將空白或缺測值轉為0
get.rec.err.func<- function(x){
  tryCatch({x},
           warning = function(msg){
             return(matrix(0, nrow = 1, ncol = 600))# 回傳0值
           },
           error = function(msg){
             return(matrix(0, nrow = 1, ncol = 600))# 回傳0值
           }
  )
}
# 2.測站缺測時值都是0
daily.rec.err.func <- function(x){
  tryCatch({x},
           warning = function(msg){
             return(0)
           },
           error = function(msg){
             return(0)
           }
  )
}
#----
sum.data <- c()
# 抓取每年資料
for(s in c(1:length(station.num))){
  
  for(yr in year){ #資料年分：
    url <- paste0("https://gweb.wra.gov.tw/wrhygis/ebooks/ebook/ebook/hyb", yr, "/", station.num[s], ".HTM")
    get.rec.err.func(
      data <- read_html(url) %>% html_nodes(xpath = "//td") %>% html_text()
    )
    t = 1
    while(t < length (data)){
      data[t] <- daily.rec.err.func(as.numeric(data[t]))
      t = t + 1
    }
    sum.data <- cbind(sum.data, c(yr, data))
    sum.data <- as.data.frame(sum.data)
  } 
}
rm.data <- sum.data[28:430,] # 移除多餘的欄位
colnames(rm.data) <- year # 列的名稱

# 清洗資料
for (i in 1:length(year)){
  year.line <- t(data.frame(rm.data[paste0(year[i])])) #選出一年的資料
  
  year.table <- t(data.frame(matrix(year.line,ncol=31,nrow=13))) # 將數列資料轉成(月*日)資料
  
  year.data<- as.data.frame(year.table) # (月*日)資料轉成 data.frame
  year.data <- cbind(year[i],year.data)
  colnames(year.data) <- c("年","日","一月","二月","三月","四月","五月","六月",
                           "七月","八月","九月","十月","十一","十二月")
  #輸出csv
  if(output == "y"){
    write.csv(year.data, file = paste0("F:/R_output/",station,"/discharge/",year[i], ".csv"))}
}

