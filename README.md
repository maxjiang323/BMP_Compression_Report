# BMP 檔案壓縮專案-程式使用方法 (BMP Compression Report-Usage)
- This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.
- 本專案使用 Makefile 以結構化的方式管理程式的編譯與執行。下面以 Makefile 中的 make 指令進行分點說明

## 1. make build
內容如下：
```
    gcc encoder.c $(SOURCES_ENCODER) -o encoder.exe -lm -Iheaders
    gcc decoder.c $(SOURCES_DECODER) -o decoder.exe -lm -Iheaders
```
- 專案中包含多個副程式，導致指令過於冗長，故以變數的方式進行簡化。
- 在此使用 GCC 編譯器分別編譯 encoder.c 和 decoder.c 兩個程式，並將生成的執行檔命名為 encoder.exe 和 decoder.exe。
- 指令尾端的 -lm 用來告訴 GCC 編譯器鏈接數學函式庫 <math.h>，而 -Iheaders 則指定編譯器在在處理 #include指令時，若路徑中未找到相應的標頭檔（header file），會額外在headers這個目錄下搜尋。

## 2. make ascii
內容如下：
```  
./encoder.exe x.bmp ascii codebook.txt huffman_code.txt
./decoder.exe QResa_x.bmp ascii codebook.txt huffman_code.txt
```
###
- x.bmp：輸入的測試 BMP 檔案
- codebook.txt：以 ASCII 表示的 codebook，儲存以下內容：
    - 符號（Symbol）。
    - 符號對應的數量（count）。
    - 符號對應的編碼（codeword）。
- huffman_code.txt: 以 ASCII 方式儲存的 Encoding，包含：
    - 圖片的 header 資訊。
    - Symbol串（string）對應的 bitstream，就是Y, Cb, Cr三個頻道（channel）中數值對應的 Huffman 編碼。
- QResa_x.bmp：使用 ASCII 方式儲存的 huffman_code.txt 及    codebook.txt 還原的 BMP 檔案。

## 3. make binary
內容如下：
```
./encoder.exe x.bmp binary codebook.txt huffman_code.bin
./decoder.exe QResb_x.bmp binary codebook.txt huffman_code.bin
```
- x.bmp：輸入的測試 BMP 檔案。
- codebook.txt：以 ASCII 表示的 codebook，儲存以下內容：
    - 符號（Symbol）。
    - 符號對應的數量（count）。
    - 符號對應的編碼（codeword）。
- huffman_code.bin: 以 Binary 方式儲存的 Encoding，包含：
    - 圖片的 header 資訊。
    - Symbol 串對應的 bitstream，也就是Y, Cb, Cr三個頻道中數值對應的 Huffman 編碼。
- QResb_x.bmp：使用 Binary 方式儲存的 huffman_code.bin 及  codebook.txt 還原的 BMP 檔案。


## 4. make clean
內容如下：
```
rm -f *.txt *.exe QRes*.bmp *.bin
```
使用 rm -f 指令強制清除執行後產生的各種檔案。

