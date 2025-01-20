// runLength.h
#ifndef RUNLENGTHCODING_H
#define RUNLENGTHCODING_H

#include "definitions.h"
void run_length_encoding(int height, int width, short ***zz, short **rle_code, long *size);
void un_rle(short *rle_code, short *zz_mn, long *pos);

#endif // RUNLENGTHCODING_H
