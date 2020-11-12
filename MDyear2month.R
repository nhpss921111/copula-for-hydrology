# Missing Data year convert to month
# 
rm(list=ls())
year <- c(2012:2019)
month <- c(1:12)
setwd("F:/R_reading/CHIA-YUANG/missingdata")
input <- paste0(year,"QandQs.csv")
output <- paste0(month,"monthQandQs.csv")
cb.data <- data.frame() #設定初始空的表格
for(y in 1:length(year)){
  data <- read.csv(file.path(getwd(),input[y]))
  cb.data <- rbind(cb.data,data)
}

for(m in month){
  month.data <- subset(cb.data,Month==month[m])
  month.data <- month.data[,-1]
  file <- paste("F:/R_output/CHIA-YUANG/missingdata/", output[m], sep="") #存檔路徑  
  write.csv(month.data,file) #結果寫到csv裡面
}
