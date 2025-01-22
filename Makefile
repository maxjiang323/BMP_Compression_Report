# build

CC = gcc 
CFLAGS = -lm -Iheaders

SOURCES_ENCODER = sources/headerInfo.c sources/pixeldata.c \
                  sources/testouput.c sources/mathCalculation.c \
                  sources/sort.c sources/runLengthCoding.c \
                  sources/huffmanEncoding.c sources/encoderMemory.c \
				  sources/SQNR.c

SOURCES_DECODER = sources/headerInfo.c sources/pixeldata.c \
                  sources/testouput.c sources/mathCalculation.c \
                  sources/sort.c sources/runLengthCoding.c \
                  sources/huffmanDecoding.c sources/decoderMemory.c

build: encoder.c decoder.c
	@$(CC) encoder.c $(SOURCES_ENCODER) -o encoder.exe $(CFLAGS)
	@echo "gcc encoder.c (SOURCES_ENCODER) -o encoder.exe -lm -Iheaders"

	@$(CC) decoder.c $(SOURCES_DECODER) -o decoder.exe $(CFLAGS)
	@echo "gcc decoder.c (SOURCES_DECODER) -o decoder.exe -lm -Iheaders"


# ascii
ascii: encoder.exe decoder.exe x.bmp
	./encoder.exe x.bmp ascii codebook.txt huffman_code.txt SQNR_a.txt
	./decoder.exe QResa_x.bmp ascii codebook.txt huffman_code.txt 

# binary
binary: encoder.exe decoder.exe x.bmp
	./encoder.exe x.bmp binary codebook.txt huffman_code.bin SQNR_b.txt
	./decoder.exe QResb_x.bmp binary codebook.txt huffman_code.bin 

# clean
clean:
	rm -f *.txt *.exe QRes*.bmp *.bin

.PHONY: build ascii binary clean
