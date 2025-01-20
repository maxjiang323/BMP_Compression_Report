#ifndef DECODERMEMORY_H
#define DECODERMEMORY_H

#include "headers/definitions.h"

// allocate memory
void malloc_zz(int M, int N, short ****zz);
void malloc_dpcm(int M, int N, short *****dpcm);
void malloc_F(int M, int N, float *****F) ;
void calloc_f(int h, int w, float ***f);;

// free memory
void free_F(int M, int N, float ****F);
void free_f(int k, float **f);;
void free_dev(int k, float **dev);;
void free_zz(int M, int N, short ***zz);;
void free_dpcm(int M, int N, short ****dpcm);

#endif // DECODERMEMORY_H
