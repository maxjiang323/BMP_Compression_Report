#include <stdio.h>
#include "headers/headerInfo.h"

void read_dim(FILE *dim_txt, Bmpheader *header){
    fscanf(dim_txt, 
        "identifier: %c%c, filesize: %u, reserved: %hu, reserved2: %hu, "
        "bitmap_dataoffset: %u, bitmap_headersize: %u, width: %u, height: %u, "
        "planes: %hu, bits_perpixel: %hu, compression: %u, bitmap_datasize: %u, "
        "hresolution: %u, vresolution: %u, usedcolors: %u, importantcolors: %u",
        &header->identifier[0], &header->identifier[1],
        &header->filesize,
        &header->reserved, &header->reserved2,
        &header->bitmap_dataoffset,
        &header->bitmap_headersize,
        &header->width, &header->height,
        &header->planes, &header->bits_perpixel,
        &header->compression, &header->bitmap_datasize,
        &header->hresolution, &header->vresolution,
        &header->usedcolors, &header->importantcolors); 
}

void write_dim(FILE *fp, Bmpheader header){
    fprintf(fp, 
        "identifier: %c%c, filesize: %u, reserved: %hu, reserved2: %hu, "
        "bitmap_dataoffset: %u, bitmap_headersize: %u, width: %u, height: %u, "
        "planes: %hu, bits_perpixel: %hu, compression: %u, bitmap_datasize: %u, "
        "hresolution: %u, vresolution: %u, usedcolors: %u, importantcolors: %u",
        header.identifier[0], header.identifier[1],
        header.filesize,
        header.reserved, header.reserved2,
        header.bitmap_dataoffset,
        header.bitmap_headersize,
        header.width, header.height,
        header.planes, header.bits_perpixel,
        header.compression, header.bitmap_datasize,
        header.hresolution, header.vresolution,
        header.usedcolors, header.importantcolors);
}

void readheader(FILE* fp, Bmpheader *header) {
	fread(&header->identifier, sizeof(header->identifier), 1, fp);
	fread(&header->filesize, sizeof(header->filesize), 1, fp);
	fread(&header->reserved, sizeof(header->reserved), 1, fp);
	fread(&header->reserved2, sizeof(header->reserved2), 1, fp);
	fread(&header->bitmap_dataoffset, sizeof(header->bitmap_dataoffset), 1, fp);
	fread(&header->bitmap_headersize, sizeof(header->bitmap_headersize), 1, fp);
	fread(&header->width, sizeof(header->width), 1, fp);
	fread(&header->height, sizeof(header->height), 1, fp);
	fread(&header->planes, sizeof(header->planes), 1, fp);
	fread(&header->bits_perpixel, sizeof(header->bits_perpixel), 1, fp);
	fread(&header->compression, sizeof(header->compression), 1, fp);
	fread(&header->bitmap_datasize, sizeof(header->bitmap_datasize), 1, fp);
	fread(&header->hresolution, sizeof(header->hresolution), 1, fp);
	fread(&header->vresolution, sizeof(header->vresolution), 1, fp);
	fread(&header->usedcolors, sizeof(header->usedcolors), 1, fp);
	fread(&header->importantcolors, sizeof(header->importantcolors), 1, fp);
}

void write_header(FILE* outf, Bmpheader header1){
	fwrite(&header1.identifier, sizeof(char), 2, outf);
	fwrite(&header1.filesize, sizeof(unsigned int), 1, outf);
	fwrite(&header1.reserved, sizeof(unsigned short), 1, outf);
	fwrite(&header1.reserved2, sizeof(unsigned short), 1, outf);
	fwrite(&header1.bitmap_dataoffset, sizeof(unsigned int), 1, outf);
	fwrite(&header1.bitmap_headersize, sizeof(unsigned int), 1, outf);
	fwrite(&header1.width, sizeof(unsigned int), 1, outf);
	fwrite(&header1.height, sizeof(unsigned int), 1, outf);
	fwrite(&header1.planes, sizeof(unsigned short), 1, outf);
	fwrite(&header1.bits_perpixel, sizeof(unsigned short), 1, outf);
	fwrite(&header1.compression, sizeof(unsigned int), 1, outf);
	fwrite(&header1.bitmap_datasize, sizeof(unsigned int), 1, outf);
	fwrite(&header1.hresolution, sizeof(unsigned int), 1, outf);
	fwrite(&header1.vresolution, sizeof(unsigned int), 1, outf);
	fwrite(&header1.usedcolors, sizeof(unsigned int), 1, outf);
	fwrite(&header1.importantcolors, sizeof(unsigned int), 1, outf);
}
