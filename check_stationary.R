# 檢查流量與輸砂量資料是否為定常性
# 單位根檢定中的ADF檢定(擴張DF檢定)
# 開始撰寫時間：2020/12/02
# 完成撰寫日期：2020/12/07
# By連育成
library(ggplot2)
library(magrittr)
library(ggforce)
library(ggalt)
library(concaveman)
rm(list=ls())
setwd("F:/R_reading/CHIA-YUANG/")
data <- read.csv(file.path(getwd(),"1974-2019.csv"))

data$Discharge<- data$Discharge %>% log() # 取自然對數
data$Suspended.Load<- data$Suspended.Load %>% log() # 取自然對數
# ---- 分成兩組資料畫成 scatter plot ----
data.2.1 <- data.frame(subset(data, Year<=1997),group="group1")
data.2.2 <- data.frame(subset(data, Year>1997),group="group2")
data.2 <- rbind(data.2.1, data.2.2)

ggplot(data=data.2,aes(x=Discharge,y=Suspended.Load,colour=group))+
  geom_point()+
  scale_color_discrete(name="分組(23年)",labels=c("1974-1997","1998-2019"))+
  labs(x="ln流量(cms)",y="ln輸砂量Qs (公噸)") + # 座標軸名稱
  geom_mark_ellipse(expand = 0,aes(fill=group))+ # 橢圓形圈
  #geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=group))+ # 多邊形
  #geom_mark_hull(expand=0.01,aes(fill=group))+ # 凹凸多邊形
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=30))  # 字體大小

# ---- 分成三組資料畫成 scatter plot ----
data.3.1 <- data.frame(subset(data, Year<=1989),group="group1")
data.3.2 <- data.frame(subset(data, Year>1989 & Year<=2004),group="group2")
data.3.3 <- data.frame(subset(data, Year>2004),group="group3")
data.3 <- rbind(data.3.1, data.3.2, data.3.3)

ggplot(data=data.3,aes(x=Discharge,y=Suspended.Load,colour=group))+
  geom_point()+
  scale_color_discrete(name="分組(15年)",labels=c("1974-1989","1990-2004","2005-2019"))+
  labs(x="ln流量(cms)",y="ln輸砂量Qs (公噸)") + # 座標軸名稱
  geom_mark_ellipse(expand = 0,aes(fill=group))+ # 橢圓形圈
  #geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=group))+ # 多邊形
  #geom_mark_hull(expand=0.01,aes(fill=group))+ # 凹凸多邊形
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=30))  # 字體大小
  
# ---- 分成四組資料畫成 scatter plot ----
data.4.1 <- data.frame(subset(data, Year<=1985),group="group1")
data.4.2 <- data.frame(subset(data, Year>1985 & Year<=1996),group="group2")
data.4.3 <- data.frame(subset(data, Year>1996 & Year<=2007),group="group3")
data.4.4 <- data.frame(subset(data, Year>2007),group="group4")
data.4 <- rbind(data.4.1, data.4.2, data.4.3, data.4.4)

ggplot(data=data.4,aes(x=Discharge,y=Suspended.Load,colour=group))+
  geom_point()+
  scale_color_discrete(name="分組(11年)",labels=c("1974-1985","1986-1996",
                                                  "1997-2007","2008-2019"))+
  labs(x="ln流量(cms)",y="ln輸砂量Qs (公噸)") + # 座標軸名稱
  geom_mark_ellipse(expand = 0,aes(fill=group))+ # 橢圓形圈
  #geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=group))+ # 多邊形
  #geom_mark_hull(expand=0.01,aes(fill=group))+ # 凹凸多邊形
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=30))  # 字體大小

# ---- 分成九組資料畫成 scatter plot ----
data.9.1 <- data.frame(subset(data, Year<=1979),group="group1")
data.9.2 <- data.frame(subset(data, Year>1979 & Year<=1984),group="group2")
data.9.3 <- data.frame(subset(data, Year>1984 & Year<=1989),group="group3")
data.9.4 <- data.frame(subset(data, Year>1989 & Year<=1994),group="group4")
data.9.5 <- data.frame(subset(data, Year>1994 & Year<=1999),group="group5")
data.9.6 <- data.frame(subset(data, Year>1999 & Year<=2004),group="group6")
data.9.7 <- data.frame(subset(data, Year>2004 & Year<=2009),group="group7")
data.9.8 <- data.frame(subset(data, Year>2009 & Year<=2014),group="group8")
data.9.9 <- data.frame(subset(data, Year>2014),group="group9")
data.9 <- rbind(data.9.1, data.9.2, data.9.3, data.9.4, data.9.5,
                data.9.6, data.9.7, data.9.8, data.9.9)

ggplot(data=data.9,aes(x=Discharge,y=Suspended.Load,colour=group))+
  geom_point()+
  scale_color_discrete(name="分組(5年)",
                       labels=c("1974-1979","1980-1984","1985-1989",
                                "1990-1994","1995-1999","2000-2004",
                                "2005-2009","2010-2014","2015-2019"))+
  labs(x="ln流量(cms)",y="ln輸砂量Qs (公噸)") + # 座標軸名稱
  geom_mark_ellipse(expand = 0,aes(fill=group))+ # 橢圓形圈
  #geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=group))+ # 多邊形
  #geom_mark_hull(expand=0.01,aes(fill=group))+ # 凹凸多邊形
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=30))  # 字體大小

#---- 單變數定常性檢定(增強DF檢定) ----

# ggplot(data, aes(x=Discharge, fill=Month)) +
#   geom_density() +
#   facet_grid(Month~.)+
#   xlim(0,50)
# ggplot(data, aes(x=Suspended.Load, fill=Month)) +
#   geom_density() +
#   facet_grid(Month~.)

library(tseries)
library(forecast)
# 擴張DF檢定(檢定流量是否為定常性)
adf.test(data$Discharge) 
# 流量的時間序列
ggplot(data,aes(x=X,y=Discharge))+
  geom_line(aes(x=X,y=Discharge,color=Month))+
  geom_smooth(method="gam",formula = y ~ poly(x,2),color="red")+
  labs(x="時間序列(day)",y="流量(cms)") + #座標軸名稱
  ggtitle("1974~2019年家園橋觀測流量")+
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=20))  # 字體大小

# 擴張DF檢定(檢定輸砂量是否為定常性)
adf.test(data$Suspended.Load)
# 輸砂量的時間序列
ggplot(data,aes(x=X,y=Suspended.Load))+
  geom_line(aes(x=X,y=Suspended.Load,color=Month))+
  geom_smooth(method="gam",formula = y ~ poly(x,2),color="red")+
  labs(x="時間序列(day)",y="輸砂量(t)") + #座標軸名稱
  ggtitle("1974~2019年家園橋觀測輸砂量")+
  theme_bw() + # 白底
  theme(panel.grid.major = element_blank()) + # 隱藏主要格線
  theme(panel.grid.minor = element_blank()) + # 隱藏次要格線
  theme(text=element_text(size=20))  # 字體大小

