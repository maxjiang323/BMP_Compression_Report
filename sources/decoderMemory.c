#include <stdio.h>
#include <stdlib.h>
#include "headers/decoderMemory.h"

void malloc_zz(int M, int N, short ****zz){
    *zz = (short ***)malloc(M * sizeof(short **));
    if (!(*zz)) {
        fprintf(stderr, "Memory allocation failed for outer array\n");
        exit(EXIT_FAILURE);
    }
        
    int i, j;
    for (i = 0; i < M; i++) {
        (*zz)[i] = (short **)malloc(N * sizeof(short *));
        if (!(*zz)[i]) {
            fprintf(stderr, "Memory allocation failed for row %d\n", i);
            exit(EXIT_FAILURE);
        }
        for (j = 0; j < N; j++) {
            (*zz)[i][j] = (short *)malloc(64 * sizeof(short));
            if (!(*zz)[i][j]) {
                fprintf(stderr, "Memory allocation failed for element (%d, %d)\n", i, j);
                exit(EXIT_FAILURE);
            }
        }
    }
}

void malloc_dpcm(int M, int N, short *****dpcm){
    *dpcm = malloc(H * sizeof(short ***));

    int u, v, m;
    for (u = 0; u < H; u++) {
        (*dpcm)[u] = malloc(W * sizeof(short **));
        for (v = 0; v < W; v++) {
            (*dpcm)[u][v] = malloc(M * sizeof(short *));
            for (m = 0; m < M; m++) {
                (*dpcm)[u][v][m] = malloc(N * sizeof(short));
            }
        }
    }
}

void malloc_F(int M, int N, float *****F) {
    int u, v, m;
    *F = (float ****)malloc(H * sizeof(float ***));
    for (u = 0; u < H; u++) {
        (*F)[u] = (float ***)malloc(W * sizeof(float **));
        for (v = 0; v < W; v++) {
            (*F)[u][v] = (float **)malloc(M * sizeof(float *));
            for (m = 0; m < M; m++) {
                (*F)[u][v][m] = (float *)malloc(N * sizeof(float));
            }
        }
    }
}

void calloc_f(int h, int w, float ***f){
    *f = (float **)calloc(h, sizeof(float *));
    int k;
    for (k = 0; k < h; k++) {
        (*f)[k] = (float *)calloc(w, sizeof(float));
    } 
}

void free_F(int M, int N, float ****F){
    int u, v, m;
    for (u = 0; u < H; u++) {
        for (v = 0; v < W; v++) {
            for (m = 0; m < M; m++) {
                free(F[u][v][m]);
            }
            free(F[u][v]);
        }
        free(F[u]);
    }
    free(F);
}

void free_f(int k, float **f){
    int fr;
    for (fr = 0; fr < k; fr++) {
        free(f[fr]);
    }
    free(f);  
}

void free_dev(int k, float **dev){
    int fr;
    for (fr = 0; fr < k; fr++) {
        free(dev[fr]);
    }
    free(dev);
}

void free_zz(int M, int N, short ***zz){
    int i, j;
    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            free(zz[i][j]);
        }
        free(zz[i]);  
    }
    free(zz); 
}

void free_dpcm(int M, int N, short ****dpcm){
    int u, v, m;
    for (u = 0; u < H; u++) {
        for (v = 0; v < W; v++) {
            for (m = 0; m < M; m++) {
                free(dpcm[u][v][m]);
            }
            free(dpcm[u][v]);
        }
        free(dpcm[u]);
    }
    free(dpcm);
}
