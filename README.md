# research
論文程式碼

以聯結函數建立流量和輸砂量雙變量之機率模型

## code使用方法
1. 蒐集水文年報裡流量和輸砂量的資料
    1999年之前的資料為網頁格式(使用爬蟲擷取資料)
    使用 [crawling.R](https://github.com/nhpss921111/research/blob/master/crawling.R)
    2000年之後的資料為PDF檔(一年一年複製資料到記事本)
    使用 [txt2csv.R](https://github.com/nhpss921111/research/blob/master/txt2csv.R)
   
2. 2000年之後的需要整理資料，1999年以前的不用 
    [爬蟲資料整理.R](https://github.com/nhpss921111/research/blob/master/%E7%88%AC%E8%9F%B2%E8%B3%87%E6%96%99%E6%95%B4%E7%90%86.R)

3. 原始資料為的檔名為【年分.csv】

4. 如果要把觀測年分的資料全部合成一個檔案使用 [combine_data.R](https://github.com/nhpss921111/research/blob/master/combine_data.R)
   
5. 觀測資料以 月分 分組，使用excel操作即可檔名為【1month.csv】【2month.csv】...以此類推
   
5. 邊際分布之參數估計、適合度檢定、最佳機率分布選擇 
   [dist_test_month.R](https://github.com/nhpss921111/research/blob/master/dist_test_month.R)
   
6. 使用copula建立聯合分布函數
   把 [dist_test_month.R](https://github.com/nhpss921111/research/blob/master/dist_test_month.R)算出來的邊際分布名稱代到 [copula.R](https://github.com/nhpss921111/research/blob/master/copula.R)


