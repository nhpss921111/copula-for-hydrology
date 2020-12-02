# 檢查流量與輸砂量資料是否為定常性
# 單位根檢定中的ADF檢定(擴張DF檢定)
# 開始撰寫時間：2020/12/02
# 完成撰寫日期：2020/12/02
# By連育成
library(ggplot2)
rm(list=ls())
setwd("F:/R_reading/CHIA-YUANG/")
data <- read.csv(file.path(getwd(),"1974-2019.csv"))
ggplot(data=data)+
  geom_point(aes(x=Discharge,y=Suspended.Load,color=Month,shape=))

ggplot(data=data)+
  geom_point(aes(x=Discharge,y=Suspended.Load,color=Year))

ggplot(data, aes(x=Discharge, fill=Month)) +
  geom_density() +
  facet_grid(Month~.)+
  xlim(0,50)
ggplot(data, aes(x=Suspended.Load, fill=Month)) +
  geom_density() +
  facet_grid(Month~.)

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

