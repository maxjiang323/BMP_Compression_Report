# build
SOURCES_ENCODER = sources/headerInfo.c sources/pixeldata.c \
                  sources/testouput.c sources/mathCalculation.c \
                  sources/sort.c sources/runLengthCoding.c \
                  sources/huffmanEncoding.c sources/encoderMemory.c

SOURCES_DECODER = sources/headerInfo.c sources/pixeldata.c \
                  sources/testouput.c sources/mathCalculation.c \
                  sources/sort.c sources/runLengthCoding.c \
                  sources/huffmanDecoding.c sources/decoderMemory.c

build: encoder.c decoder.c
	gcc encoder.c $(SOURCES_ENCODER) -o encoder.exe -lm -Iheaders
	gcc decoder.c $(SOURCES_DECODER) -o decoder.exe -lm -Iheaders


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

.PHONY: build ascii binary clean
