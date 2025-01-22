# BMP 檔案壓縮專案-程式使用方法
專案程式連結：https://github.com/maxjiang323/BMP_Compression_Report
本專案使用Makefile以結構化的方式管理程式的編譯與執行。下面以 Makefile 中的 make 指令進行分點說明

## 1. make build
內容如下：
```
    gcc encoder.c $(SOURCES_ENCODER) -o encoder.exe -lm -Iheaders
    gcc decoder.c $(SOURCES_DECODER) -o decoder.exe -lm -Iheaders
```
- 專案中包含多個副程式，導致指令過於冗長，故以變數的方式進行簡化。
- 在此使用GCC編譯器分別編譯encoder.c和decoder.c兩個程式，並將生成的執行檔命名為encoder.exe和decoder.exe。
- 指令尾端的 -lm用來告訴GCC編譯器鏈接數學函式庫 <math.h>，而 -Iheaders則指定編譯器在在處理 #include指令時，若路徑中未找到相應的標頭檔（header file），會額外在headers這個目錄下搜尋。

## 2. make ascii
內容如下：
```  
./encoder.exe x.bmp ascii codebook.txt huffman_code.txt
./decoder.exe QResa_x.bmp ascii codebook.txt huffman_code.txt
```
###
- x.bmp：輸入的測試 BMP 檔案
- codebook.txt：以ASCII表示的 codebook，儲存以下內容：
    - 符號（Symbol）。
    - 符號對應的數量（count）。
    - 符號對應的編碼（codeword）。
- huffman_code.txt: 以ASCII方式儲存的 Encoding，包含：
    - 圖片的header資訊。
    - Symbol串（string）對應的bitstream，就是Y, Cb, Cr三個頻道（channel）中數值對應的Huffman編碼。
- QResa_x.bmp：使用ASCII方式儲存的huffman_code.txt及   codebook.txt還原的 BMP 檔案。

## 3. make binary
內容如下：
```
./encoder.exe x.bmp binary codebook.txt huffman_code.bin
./decoder.exe QResb_x.bmp binary codebook.txt huffman_code.bin
```
- x.bmp：輸入的測試 BMP 檔案。
- codebook.txt：以ASCII表示的 codebook，儲存以下內容：
    - 符號（Symbol）。
    - 符號對應的數量（count）。
    - 符號對應的編碼（codeword）。
- huffman_code.bin: 以Binary方式儲存的Encoding，包含：
    - 圖片的header資訊。
    - Symbol串對應的bitstream。也就是Y, Cb, Cr三個頻道中數值對應的Huffman編碼。
- QResb_x.bmp：使用Binary方式儲存的huffman_code.bin及 codebook.txt還原的 BMP 檔案。


## 4. make clean
內容如下：
```
rm -f *.txt *.exe QRes*.bmp *.bin
```
使用rm -f 指令強制清除執行後產生的各種檔案。

