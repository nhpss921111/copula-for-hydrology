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



# 註：以下程式碼之讀檔路徑與寫檔路徑需自行更改

## 蒐集資料


### 流量的資料
1. 2000年之前的資料為網頁格式(使用爬蟲擷取資料)使用 [crawling4discharge.R](https://github.com/nhpss921111/copula-for-hydrology/blob/master/crawling4discharge.R)
存檔路徑 F:/R_output/crawling/JEN-SHOU BRIDGE/discharge/2000.csv
2. 2001年之後的資料為PDF檔(一年一年複製資料到記事本)使用 [txt2csv4Q.R](https://github.com/nhpss921111/copula-for-hydrology/blob/master/txt2csv4Q.R)
存檔路徑 F:/R_output/JEN-SHOU BRIDGE/discharge/2019.csv

* PDF檔複製範圍：
![](https://i.imgur.com/aTcKLfC.png)

### 流量及輸砂量的資料

1. 2000年之前的資料為網頁格式(使用爬蟲擷取資料)使用 [crawling.R](https://github.com/nhpss921111/copula-for-hydrology/blob/master/crawling.R)
存檔路徑 F:/R_output/crawling/JEN-SHOU BRIDGE/discharge+SSL/2000.csv
2. 2001年之後的資料為PDF檔(一年一年複製資料到記事本)使用 [txt2csv.R](https://github.com/nhpss921111/copula-for-hydrology/blob/master/txt2csv.R)
* PDF檔複製範圍：
![](https://i.imgur.com/Jg8GdWV.png)


## 清洗資料

1. 2000年之後的需要整理資料[爬蟲資料整理.R](https://github.com/nhpss921111/copula-for-hydrology/blob/master/%E7%88%AC%E8%9F%B2%E8%B3%87%E6%96%99%E6%95%B4%E7%90%86.R)，1999年以前的不用 
    
2. 原始資料為的檔名為【年分.csv】

3. 如果要把觀測年分的資料全部合成一個檔案使用 [combine_data.R](https://github.com/nhpss921111/copula-for-hydrology/blob/master/combine_data.R)
   
4. 觀測資料以"月分"分組，使用excel操作即可檔名為【1month.csv】【2month.csv】...以此類推

### 統計檢定 
1. 邊際分布之參數估計、適合度檢定、最佳機率分布選擇 
   [dist_test_month.R](https://github.com/nhpss921111/copula-for-hydrology/blob/master/dist_test_month.R)
   
2. 使用copula建立聯合分布函數
   把 [dist_test_month.R](https://github.com/nhpss921111/copula-for-hydrology/blob/master/dist_test_month.R)算出來的邊際分布名稱代到 [copula.R](https://github.com/nhpss921111/copula-for-hydrology/blob/master/copula.R)
