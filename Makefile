# build
build: encoder.c decoder.c
	gcc encoder.c -o encoder.exe -lm
	gcc decoder.c -o decoder.exe -lm

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
