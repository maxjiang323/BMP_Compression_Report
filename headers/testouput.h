#ifndef TESTOUT_H
#define TESTOUT_H

#include <stdio.h>
#include "headers/definitions.h"

// encoder.c
void test_output_of_dpcm(FILE *fp, int height, int width, short ****dpcm);
void test_output_of_zz(FILE *fp, int height, int width, short ***output);

// decoder.c
void test_output_of_each_block(FILE *fp, int width, int height, int block_size, float **output);
void test_output_of_rle_txt(FILE *rle_txt, FILE* rle_out, long *read_pos, long *pos, short **rle_code);
void test_output_of_un_rle(FILE *fp, int height, int width, short ***zigzag);
void test_output_of_un_zigzag_dpcm(FILE *fp, int height, int width, short ****output);
void test_output_of_un_quantizaton(FILE *fp, int height, int width, float ****F);

#endif // TESTOUT_H
