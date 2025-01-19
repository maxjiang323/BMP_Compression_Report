# build
build: encoder.c decoder.c
	gcc encoder.c headerInfo.c pixeldata.c testouput.c mathCalculation.c sort.c huffmanEncoding.c -o encoder.exe -lm
	gcc decoder.c headerInfo.c pixeldata.c testouput.c mathCalculation.c sort.c huffmanDecoding.c -o decoder.exe -lm

# ascii
ascii: encoder.exe decoder.exe x.bmp
	./encoder.exe x.bmp ascii codebook.txt huffman_code.txt
	./decoder.exe QResa_x.bmp ascii codebook.txt huffman_code.txt 

# binary
binary: encoder.exe decoder.exe x.bmp
	./encoder.exe x.bmp binary codebook.txt huffman_code.bin 
	./decoder.exe QResb_x.bmp binary codebook.txt huffman_code.bin 

# clean
clean:
	rm -f *.txt *.exe QRes*.bmp *.bin
