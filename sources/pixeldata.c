#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "headers/pixeldata.h"

#define W 8 // dimension of basis vector (width)
#define H 8 // dimension of basis vector (height)


// encoder.c
void read_rgb_data(FILE *fp, Bmpheader header, int *width, int *height, int *row_length, unsigned char **pix_rgb){
    *width = header.width;  
    *height = header.height;
    int width_rgb = *width * 3;  // rgb -> *3
    int padding = 4 - (width_rgb % 4);
    if ((width_rgb % 4) == 0){
        padding = 0;
    }
    *row_length = width_rgb + padding; 
    // add the padding to make (the modified width, the length of a row) % 4 == 0
    *pix_rgb = (unsigned char *)calloc((*row_length) * (*height), sizeof(unsigned char));

    // read data
    fseek(fp, header.bitmap_dataoffset, SEEK_SET); // jump to the start of pixel data  
    fread(*pix_rgb, sizeof(unsigned char), (*row_length) * (*height), fp);
}

YCbCr RGB2YCbCr(unsigned char R, unsigned char G, unsigned char B) { // Rec. 709
    YCbCr result;
    result.Y = 0.2126 * R + 0.7152 * G + 0.0722 * B;
    result.Cb = -0.1146 * R - 0.3854 * G + 0.5 * B;
    result.Cr = 0.5 * R - 0.4542 * G - 0.0458 * B;
    return result;
}

void put_RGB2YCbCr_into_block(int height, int width, int block_size, int row_length, unsigned char *pix_rgb, float ***f_Y, float ***f_Cb, float ***f_Cr){
    // add FILE *en_f_Y to parameter if want to see the output
    int rM = H - (height % H), rN = W - (width % W); // the remainder of the height and width
    // may generate non-complete H*W block    
    if (rM == H) rM = 0; // if the height is a multiple of H, no extra row is needed 
    if (rN == W) rN = 0; // if the width is a multiple of W, no extra column is needed 

    int block_x, block_y, x, y;
    for (block_x = 0; block_x < height+rM; block_x += block_size) {
        for (block_y = 0; block_y < width+rN; block_y += block_size) {
            for (x = 0; x < block_size; x++) {
                for (y = 0; y < block_size; y++) {
                    int index = (block_x + x) * row_length + (block_y + y) * 3;
                    // Check if we are within the image bounds
                    if ((block_x + x) < height && (block_y + y) < width) {
                        unsigned char b = pix_rgb[index];
                        unsigned char g = pix_rgb[index + 1];
                        unsigned char r = pix_rgb[index + 2];

                        YCbCr result = RGB2YCbCr(r, g, b);

                        (*f_Y)[block_x + x][block_y + y] = result.Y;
                        (*f_Cb)[block_x + x][block_y + y] = result.Cb;
                        (*f_Cr)[block_x + x][block_y + y] = result.Cr;
                    } else {
                            (*f_Y)[block_x + x][block_y + y] = 0.0;
                            (*f_Cb)[block_x + x][block_y + y] = 0.0;
                            (*f_Cr)[block_x + x][block_y + y] = 0.0;
                    }

                    // fprintf(en_f_Y, "%f ", (*f_Y)[block_x + x][block_y + y]);
                }
                // fprintf(en_f_Y, "\n");
            }
            // fprintf(en_f_Y, "\n\n");
        }
    }
}


// decoder.c
RGB YCbCr2RGB(float Y, float Cb, float Cr) { // Rec. 709
    RGB result;
    float R = Y + 1.5748*Cr;
    float G = Y - 0.1873*Cb - 0.4681*Cr;
    float B = Y + 1.8556*Cb;

    // Clamp values to 0-255 range
    if (R < 0){
        result.R = 0;
    } else if (R > 255){
        result.R = 255;
    } else{
        result.R = (unsigned char)(round(R));
    }

    if (G < 0){
        result.G = 0;
    } else if (G > 255){
        result.G = 255;
    } else{
        result.G = (unsigned char)(round(G));
    }

    if (B < 0){
        result.B = 0;
    } else if (B > 255){
        result.B = 255;
    } else{
        result.B = (unsigned char)(round(B));
    }

    return result;
}

void set_rgb_data(int width, int height, int *row_length, unsigned char **pix_rgb){
    int width_rgb = width * 3;  // rgb -> *3
    int padding = 4 - (width_rgb % 4);
    if ((width_rgb % 4) == 0){
        padding = 0;
    }
    *row_length = width_rgb + padding; 
    // add the padding to make (the modified width, the length of a row) % 4 == 0
    *pix_rgb = (unsigned char *)calloc((*row_length) * (height), sizeof(unsigned char));
}
