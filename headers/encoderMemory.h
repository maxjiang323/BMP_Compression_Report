#ifndef ENCODERMEMORY_H
#define ENCODERMEMORY_H

#include "definitions.h"

// allocate memory
void calloc_f(int h, int w, float ***f);
void malloc_F(int M, int N, float *****F);
void malloc_qF(int M, int N, short *****qF);
void malloc_zz(int M, int N, short ****zz);
void malloc_dpcm(int M, int N, short *****dpcm);

// free memory
void free_F(int M, int N, float ****F);
void free_qF(int M, int N, short ****qF);
void free_f(int k, float **f);
void free_zz(int M, short ***zz);

#endif // ENCODERMEMORY_H
