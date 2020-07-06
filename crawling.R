# 來源：爬蟲，水文年報(2000年以前)  2019/11/02 BY奕廷學長
# 
# 開始撰寫日期：2020/07/02
# 完成撰寫日期：2020/07/02
rm(list = ls())
library(rvest)
output <- "n"
#
#----起始資料
#測站編號
station <- c("2560H017")
#
#----
# 先建立錯誤函數：
# 1.空白值將錯誤，建立tryCatch將空白或缺測值轉為-50
get.rec.err.func<- function(x){
  tryCatch({x}, 
           warning = function(msg){
             return(matrix(-50, nrow = 1, ncol = 600))# 回傳-50值
           }, 
           error = function(msg){
             return(matrix(-50, nrow = 1, ncol = 600))# 回傳-50值
           }
  )
}
# 2.測站缺測時值都是-50
daily.rec.err.func <- function(x){
  tryCatch({x}, 
           warning = function(msg){
             return(-50)
           }, 
           error = function(msg){
             return(-50)
           }
  )
}
#----

#抓取每年資料
for(s in c(1:length(station))){
  all.data <- c(NULL)
  
  for(yr in c(1974:2000)){ #資料年分
    url <- paste0("https://gweb.wra.gov.tw/wrhygis/ebooks/ebook/ebook/hyb", yr, "/", station[s], ".HTML")
     #get.rec.err.func(
    data <- read_html(url) %>% html_nodes(".Section1 td font") %>% html_text
    #)
    t = 1
    while(t < length (data)){
      data[t] <- daily.rec.err.func(as.numeric(data[t]))
      t = t + 1
    }
    
    all.data <- rbind(all.data, c(yr, data))
    break
    
  }
  #輸出csv
  #all.data <- all.data[,-2:-18]
  if(output == "y"){
    write.csv(all.data, file = paste0("E:/R_output/crawling/",station[s], ".csv")) 
  }
}

#接著匯入"爬蟲資料整理.R"整理資料