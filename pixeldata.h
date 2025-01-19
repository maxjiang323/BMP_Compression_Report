#ifndef PIXELDATA_H
#define PIXELDATA_H

#include <stdio.h>
#include "definitions.h"

// encoder.c
void read_rgb_data(FILE *fp, Bmpheader header, int *width, int *height, int *row_length, unsigned char **pix_rgb);
YCbCr RGB2YCbCr(unsigned char R, unsigned char G, unsigned char B);
void put_RGB2YCbCr_into_block(int height, int width, int block_size, int row_length, unsigned char *pix_rgb, float ***f_Y, float ***f_Cb, float ***f_Cr);

// decoder.c
RGB YCbCr2RGB(float Y, float Cb, float Cr);
void set_rgb_data(int width, int height, int *row_length, unsigned char **pix_rgb);

#endif // PIXELDATA_H
