# 數據呈現
# 開始撰寫日期：2020/06/24
# 完成撰寫日期：2020/06/29
#
#
rm(list = ls())
library(xlsx) # 讀取excel
library(ggplot2) #繪圖用
#
# Read data from excel flie
#
setwd("E:/R_reading/CHIA-YUANG")
data <- read.xlsx(file.path(getwd(),"1974-2019.xlsx"),sheetIndex="Sheet1",
                  startRow = 1,endRow = 2000,
                  header = T,colIndex =2:3,encoding = "UTF-8")
data1 <- read.xlsx(file.path(getwd(),"1974-1986.xlsx"),sheetIndex="Sheet1",
                   startRow = 1,endRow = 2000,
                   header = T,colIndex =2:3,encoding = "UTF-8")
data2 <- read.xlsx(file.path(getwd(),"1987-1997.xlsx"),sheetIndex="Sheet1",
                   startRow = 1,endRow = 2000,
                   header = T,colIndex =2:3,encoding = "UTF-8")
data3 <- read.xlsx(file.path(getwd(),"1998-2008.xlsx"),sheetIndex="Sheet1",
                   startRow = 1,endRow = 2000,
                   header = T,colIndex =2:3,encoding = "UTF-8")
data4 <- read.xlsx(file.path(getwd(),"2009-2019.xlsx"),sheetIndex="Sheet1",
                   startRow = 1,endRow = 2000,
                   header = T,colIndex =2:3,encoding = "UTF-8")
Q1 <- data1$Q
S1 <- data1$S
Q2 <- data2$Q
S2 <- data2$S
Q3 <- data3$Q
S3 <- data3$S
Q4 <- data4$Q
S4 <- data4$S
#
# ----------------------------------- plot --------------------------
#
ggplot(data) + 
  geom_point(aes(x = Q, y = S),color="black") +
  geom_smooth(aes(x = Q, y =S),color="black")+
  xlim(0,350)+
  ylim(0,800000)+
  labs(title="逕流和輸砂量的分布1974-2019",x="Q(cms)",y="Qs(公噸/日)") 
ggplot(data1) + 
  geom_point(aes(x = Q, y = S),color="red") +
  geom_smooth(aes(x = Q, y =S),color="red")+
  xlim(0,350)+
  ylim(0,800000)+
  labs(title="逕流和輸砂量的分布1974-1986",x="Q(cms)",y="Qs(公噸/日)") 
ggplot(data2) + 
  geom_point(aes(x = Q, y = S),color="purple") +
  geom_smooth(aes(x = Q, y =S),color="purple")+
  xlim(0,350)+
  ylim(0,800000)+
  labs(title="逕流和輸砂量的分布1987-1997",x="Q(cms)",y="Qs(公噸/日)") 
ggplot(data3) + 
  geom_point(aes(x = Q, y = S),color="green")+
  geom_smooth(aes(x = Q, y =S),color="green")+
  xlim(0,350)+
  ylim(0,800000)+
  labs(title="逕流和輸砂量的分布1998-2008",x="Q(cms)",y="Qs(公噸/日)") 
ggplot(data4) + 
  geom_point(aes(x = Q, y = S),color="blue") +
  geom_smooth(aes(x = Q, y =S),color="blue")+
  xlim(0,350)+
  ylim(0,800000)+
  labs(title="逕流和輸砂量的分布2009-2019",x="Q(cms)",y="Qs(公噸/日)") 


ggplot()+
  geom_point(aes(x=Q1, y=S1),color="red") +
  geom_smooth(aes(x=Q1, y=S1),color="red") +
  geom_point(aes(x=Q2, y=S2),color="purple") +
  geom_smooth(aes(x=Q2, y=S2),color="purple") +
  geom_point(aes(x=Q3, y=S3),color="green") +
  geom_smooth(aes(x=Q3, y=S3),color="green") +
  geom_point(aes(x=Q4, y=S4),color="blue") +
  geom_smooth(aes(x=Q4, y=S4),color="blue") +
  labs(title="逕流和輸砂量的分布1974-2019分4組",x="Q(cms)",y="Qs(公噸/日)") 
  #theme(legend.position = "right")

