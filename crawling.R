# 來源：爬蟲，水文年報(2000年以前)  2019/11/02 BY奕廷學長
# 因為每年觀測次數不同，所以無法一次全部爬完再存
# 所以需要每爬完一年就儲存
# 存好的檔案再丟到"爬蟲資料整理.R"清洗資料
# 開始撰寫日期：2020/07/02
# 完成撰寫日期：2021/03/08
rm(list = ls())
library(rvest)
library(magrittr)
library(jsonlite)
library(dplyr)
output <- "y" 
# ========================================
#station <- c("CHIA-YUANG")
#station.num <- c("2560H017") #蘭陽溪-家源橋
#year <- c(1974:2000)
#
#station <- c("CHUNYUN BRIDGE")
#station.num <- c("1510H057") #濁水溪-彰雲橋
#year <- c(1987:2000) 
#
#station <- c("NEI-MAO-PU")
#station.num <- c("1510H049") #濁水溪-內茅埔
#year <- c(1972:2000)
#
#station <- c("JEN-SHOU BRIDGE")
#station.num <- c("2420H019") #花蓮溪-仁壽橋 
#year <- c(1960:2000) 
#
#station <- c("LIU-KWEI")
#station.num <- c("1730H039") #高屏溪-六龜
#year <- c(1982:2000) 
#
#station <- c("NEI-WAN")
#station.num <- c("1300H013")# 頭前溪-內灣
#year <- c(1971:2000) 

#station <- c("JEIN-KUO BRIDGE")
#station.num <- c("1300H017")# 頭前溪-經國橋
#year <- c(1990:2000) 
#
#station.num <- c("2560H006")
#station <- c("LAN-YANG BRIDGE") # 蘭陽溪-蘭陽大橋
#year <- c(1949:2000) 
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
#
#station.num <- c("1730H043")
#station <- c("LI-LIN BRIDGE")#里嶺大橋
#year <- c(1991:2000) 
# ========================================
#----起始資料
#測站編號
station <- c("JEN-SHOU BRIDGE")
station.num <- c("2420H019") #花蓮溪-仁壽橋 
year <- c(1960:2000) 
#----
# 先建立錯誤函數：
# 1.空白值將錯誤，建立tryCatch將空白或缺測值轉為0
get.rec.err.func<- function(x){
  tryCatch({x}, 
           warning = function(msg){
             return(matrix(0, nrow = 1, ncol = 600))# 回傳-50值
           }, 
           error = function(msg){
             return(matrix(0, nrow = 1, ncol = 600))# 回傳-50值
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
year.data <- c()

#抓取每年資料
for(s in c(1:length(station.num))){
  
  for(yr in year){ #資料年分：
    url <- paste0("https://gweb.wra.gov.tw/wrhygis/ebooks/ebook/ebook/hyb", yr, "/", station.num[s], ".HTML")
    get.rec.err.func(
      data <- read_html(url) %>% html_nodes(xpath = "//td") %>% html_text()
    )
    t = 1
    while(t < length (data)){
      data[t] <- daily.rec.err.func(as.numeric(data[t]))
      t = t + 1
    }
    year.data <- c(yr, data)
    year.data <- year.data[-2:-24] # 刪除多餘欄位
    #year.data <- as.data.frame(year.data)
    
    # #輸出csv
    # if(output == "y"){
    #   file <- paste0("F:/R_output/crawling/",station,"/discharge+SSL/",yr,".csv")
    #   write.csv(year.data, file)
    # }
  }
}
for (n in c(1:length(year))){
  ini.data <- year.data[2:7]
  arrange.data <- c()
  i <- 1
  while (i <300){
    slice.data <- year.data[(i+1):(i+6)] #每一天的資料
    arrange.data <- rbind(arrange.data,slice.data)
    colnames(arrange.data) <- c("No.","Month","Day","Discharge",
                                "Sediment Content","Suspended Load")
    i <- i+6
  }
  rm.data <- arrange.data[complete.cases(arrange.data), ]
  new.data <- arrange(as.data.frame(rm.data),No.) #依照No.欄位，由小到大排序
  new.data[,1] <- year[n]
  colnames(new.data) <- c("Year","Month","Day","Discharge",
                          "Sediment Content","Suspended Load")
  write.csv(new.data, file = paste0("F:/R_output/",station,"/SSL_arranged/",year[n], ".csv"))
}