# 以聯結函數補遺日懸浮載輸砂量

Published in water of MDPI : [Copula-Based Infilling Methods for Daily Suspended Sediment Loads](https://www.mdpi.com/2073-4441/13/12/1701/htm)

Master's thesis：[以聯結函數補遺日懸浮載輸砂量](https://thesis.lib.ncku.edu.tw/thesis/detail/bee9602ce9debe703eaa908b5075e30b/?seq=1)

## 研究資料來源：[水利署_水文年報_流量與輸砂量之觀測資料](https://gweb.wra.gov.tw/wrhygis/)

## 觀察原始資料格式

* 註：2000年以前資料為網頁形式，2001年以後資料為PDF形式

### 流量資料(日觀測資料)
1. [民國89年 日平均流量](https://gweb.wra.gov.tw/wrhygis/ebooks/ebook/ebook/hyb2000/2420H019.HTM)
2. 民國108年 日平均流量
![](https://i.imgur.com/eUn5Oy3.png)
### 同時有流量及輸砂量資料(一個月約只有2~4筆觀測資料)
1. [民國89年 懸移質實測紀錄](https://gweb.wra.gov.tw/wrhygis/ebooks/ebook/ebook/hyb2000/2420H019.HTML)
2. 民國108年 懸移質實測紀錄
![](https://i.imgur.com/WpKNNwN.png)

# 程式碼使用說明
## 註：以下程式碼之讀檔路徑與寫檔路徑需自行更改
---
## 一、蒐集&清洗資料
### 蒐集流量-執行步驟
1. 2000年之前的資料為網頁格式(使用爬蟲擷取資料)使用 [crawling4discharge.R](https://github.com/nhpss921111/copula-for-hydrology/blob/master/crawling4discharge.R)
* 寫檔路徑  F:/copula/JEN-SHOU BRIDGE/discharge/2000.csv
2. 2001年之後的資料為PDF檔(一年一年複製資料到記事本，缺測的觀測要資料要補0)使用 [txt2csv4Q.R](https://github.com/nhpss921111/copula-for-hydrology/blob/master/txt2csv4Q.R)
* 讀檔路徑 F:/copula/JEN-SHOU BRIDGE/from_pdf/discharge/2019.txt
* 寫檔路徑 F:/copula/JEN-SHOU BRIDGE/discharge/2019.csv
* PDF檔複製範圍：
![](https://i.imgur.com/aTcKLfC.png)

### 蒐集流量-結果
* 寫檔路徑 F:/copula/JEN-SHOU BRIDGE/discharge/

* excel結果範例圖(2000.csv)
![](https://i.imgur.com/GkBwRGp.png)

* excel結果範例圖(2019.csv)
![](https://i.imgur.com/2Q5gduE.png)
### 蒐集流量及輸砂量-執行步驟
1. 2000年之前的資料為網頁格式(使用爬蟲擷取資料)使用 [crawling.R](https://github.com/nhpss921111/copula-for-hydrology/blob/master/crawling.R)
* 寫檔路徑 F:/copula/JEN-SHOU BRIDGE/discharge+SSL/2000.csv
2. 2001年之後的資料為PDF檔(一年一年複製資料到記事本，缺測的觀測要資料要補0)使用 [txt2csv.R](https://github.com/nhpss921111/copula-for-hydrology/blob/master/txt2csv.R)
* 讀檔路徑 F:/copula/JEN-SHOU BRIDGE/from_pdf/discharge+SSL/2019.txt
* 寫檔路徑 F:/copula/JEN-SHOU BRIDGE/discharge+SSL/2019.csv
* PDF檔複製範圍：

![](https://i.imgur.com/Jg8GdWV.png)

### 蒐集流量及輸砂量-結果
* 寫檔路徑 F:/copula/JEN-SHOU BRIDGE/discharge+SSL/

* excel結果範例圖(2000.csv)
![](https://i.imgur.com/wrLpAUG.png)

* excel結果範例圖(2019.csv)
![](https://i.imgur.com/9P2elKC.png)

### 將同時有流量及輸砂量的年資料換成月資料
* 此步驟是為了將"**流量**"與"**流量及輸砂量**"資料合併，使用[year2month.R](https://github.com/nhpss921111/copula-for-hydrology/blob/master/year2month.R)將"**流量及輸砂量**"的年份資料轉換成月份資料
* 讀檔路徑 F:/copula/JEN-SHOU BRIDGE/discharge+SSL/2019.csv
* 寫檔路徑 F:/copula/JEN-SHOU BRIDGE/discharge+SSL/1month.csv
* [結果範例，以一月份為例](https://drive.google.com/file/d/16H0rNzcdCvAIkyQ_62rAWEPZfniaVfj5/view?usp=sharing)

---

## 二、合併流量及輸砂量資料


## 三、Copula分析
