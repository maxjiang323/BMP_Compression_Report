#include <stdio.h>
#include <stdlib.h>
#include "../headers/encoderMemory.h"

// allocate memory
void calloc_f(int h, int w, float ***f){
    *f = (float **)calloc(h, sizeof(float *));
    int k;
    for (k = 0; k < h; k++) {
        (*f)[k] = (float *)calloc(w, sizeof(float));
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

void malloc_qF(int M, int N, short *****qF) {
    int u, v, m;
    *qF = (short ****)malloc(H * sizeof(short ***));
    for (u = 0; u < H; u++) {
        (*qF)[u] = (short ***)malloc(W * sizeof(short **));
        for (v = 0; v < W; v++) {
            (*qF)[u][v] = (short **)malloc(M * sizeof(short *));
            for (m = 0; m < M; m++) {
                (*qF)[u][v][m] = (short *)malloc(N * sizeof(short));
            }
        }
    }
}

void malloc_zz(int M, int N, short ****zz){
    *zz = (short ***)malloc(64 * sizeof(short **));
    if (!(*zz)) {
        fprintf(stderr, "Memory allocation failed for outer array\n");
        exit(EXIT_FAILURE);
    }
        
    int i, j;
    for (i = 0; i < 64; i++) {
        (*zz)[i] = (short **)malloc(M * sizeof(short *));
        if (!(*zz)[i]) {
            fprintf(stderr, "Memory allocation failed for row %d\n", i);
            exit(EXIT_FAILURE);
        }
        for (j = 0; j < M; j++) {
            (*zz)[i][j] = (short *)malloc(N * sizeof(short));
            if (!(*zz)[i][j]) {
                fprintf(stderr, "Memory allocation failed for element (%d, %d)\n", i, j);
                exit(EXIT_FAILURE);
            }
        }
    }
}


// free memory
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

void free_qF(int M, int N, short ****qF){
    int u, v, m;
    for (u = 0; u < H; u++) {
        for (v = 0; v < W; v++) {
            for (m = 0; m < M; m++) {
                free(qF[u][v][m]);
            }
            free(qF[u][v]);
        }
        free(qF[u]);
    }
    free(qF);
}

void free_f(int k, float **f){
    int fr;
    for (fr = 0; fr < k; fr++) {
        free(f[fr]);
    }
    free(f);  
}

void free_zz(int M, short ***zz){
    int i, j;
    for (i = 0; i < 64; i++) {
        for (j = 0; j < M; j++) {
            free(zz[i][j]);
        }
        free(zz[i]);  
    }
    free(zz); 
}
