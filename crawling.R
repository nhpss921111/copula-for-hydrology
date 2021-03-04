# 來源：爬蟲，水文年報(2000年以前)  2019/11/02 BY奕廷學長
# 因為每年觀測次數不同，所以無法一次全部爬完再存
# 所以需要每爬完一年就儲存
# 存好的檔案再丟到"爬蟲資料整理.R"清洗資料
# 開始撰寫日期：2020/07/02
# 完成撰寫日期：2021/03/04
rm(list = ls())
library(rvest)
library(magrittr)
library(jsonlite)
library(dplyr)
output <- "y"
# ---- 測站資料放置區 ----
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

#----起始資料
#測站編號
station <- c("NEI-MAO-PU")
station.num <- c("1510H049") #濁水溪-內茅埔
year <- c(1972:2000) # 1972才開始有SSL的紀錄
#
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
    year.data <- as.data.frame(year.data)
    #輸出csv
    if(output == "y"){
      file <- paste0("F:/R_output/crawling/",station,"/discharge+SSL/",yr,".csv")
      write.csv(year.data, file)
    }
  }
}