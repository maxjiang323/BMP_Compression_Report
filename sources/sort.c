#include <stdio.h>
#include <stdlib.h>
#include "../headers/sort.h"

#define W 8 // dimension of basis vector (width)
#define H 8 // dimension of basis vector (height)

const int zz_matrix[8][8] = {
    { 1,  2,  6,  7, 15, 16, 28, 29},
    { 3,  5,  8, 14, 17, 27, 30, 43},
    { 4,  9, 13, 18, 26, 31, 42, 44},
    {10, 12, 19, 25, 32, 41, 45, 54},
    {11, 20, 24, 33, 40, 46, 53, 55},
    {21, 23, 34, 39, 47, 52, 56, 61},
    {22, 35, 38, 48, 51, 57, 60, 62},
    {36, 37, 49, 50, 58, 59, 63, 64}
};


// encoder.c
void dpcm(int height, int width, short *****qF){
    int ind_M, ind_N, ind_H, ind_W;
    int rM = H - (height % H), rN = W - (width % W); // the remainder of the height and width
    // may generate non-complete H*W block    
    if (rM == H) rM = 0; // if the height is a multiple of H, no extra row is needed 
    if (rN == W) rN = 0; // if the width is a multiple of W, no extra column is needed 
    int M = (height+rM)/H;  // Number of blocks vertically
    int N = (width+rN)/W;   // Number of blocks horizontally

    for (ind_M = 0; ind_M < M; ind_M++) {
        for (ind_N = N-1; ind_N > 0; ind_N--) {
            // perform horizontal DPCM
            (*qF)[0][0][ind_M][ind_N] -= (*qF)[0][0][ind_M][ind_N-1];
        }
    }

    // perform vertical DPCM
    for (ind_M = M-1; ind_M > 0; ind_M--) {
        (*qF)[0][0][ind_M][0] -= (*qF)[0][0][ind_M-1][0];
    }
}

void apply_zigzag(int height, int width, short ****dpcm, short ***zz){
    int rM = H - (height % H), rN = W - (width % W);
    if (rM == H) rM = 0; 
    if (rN == W) rN = 0; 
    int M = (height + rM) / H; 
    int N = (width + rN) / W; 

    // printf("M: %d, N: %d\n", M, N); // Debugging output

    int m, n, r, c;
    for (m = 0; m < M; m++) {
        for (n = 0; n < N; n++) {
            for (r = 0; r < H; r++) {
                for (c = 0; c < W; c++) {
                    int index = zz_matrix[r][c] - 1;
                    if (index < 0 || index >= 64) {
                        fprintf(stderr, "Index out of bounds: %d\n", index);
                        exit(EXIT_FAILURE); // Prevent further access
                    }
                    zz[index][m][n] = dpcm[r][c][m][n];
                    // printf("%hd ", dpcm[r][c][m][n]);
                }
            }
        }
    }
}


// decoder.c
void un_zigzag(int height, int width, short ***zz, short *****dpcm) {
    int rM = H - (height % H), rN = W - (width % W);
    if (rM == H) rM = 0;
    if (rN == W) rN = 0;
    int M = (height + rM) / H;
    int N = (width + rN) / W;

    int m, n, r, c;
    for (m = 0; m < M; m++) {
        for (n = 0; n < N; n++) {
            for (r = 0; r < H; r++) {
                for (c = 0; c < W; c++) {
                    int index = zz_matrix[r][c] - 1; // Ensure zz_matrix is within bounds
                    if (index < 0 || index >= 64) {
                        fprintf(stderr, "Index out of bounds: %d\n", index);
                        exit(EXIT_FAILURE);
                    }
                    (*dpcm)[r][c][m][n] = zz[m][n][index];
                }
            }
        }
    }
}

void un_dpcm(int height, int width, short *****dpcm){
    // dpcm dimension: (H, W, M, N)
    int ind_M, ind_N, ind_H, ind_W;
    int rM = H - (height % H), rN = W - (width % W); // the remainder of the height and width
    // may generate non-complete H*W block    
    if (rM == H) rM = 0; // if the height is a multiple of H, no extra row is needed 
    if (rN == W) rN = 0; // if the width is a multiple of W, no extra column is needed 
    int M = (height+rM)/H;  // Number of blocks vertically
    int N = (width+rN)/W;   // Number of blocks horizontally

    // perform vertical un-DPCM
    for (ind_M = 1; ind_M < M; ind_M++) {
        (*dpcm)[0][0][ind_M][0] += (*dpcm)[0][0][ind_M-1][0];
    }       


    for (ind_M = 0; ind_M < M; ind_M++) {
        for (ind_N = 1; ind_N < N; ind_N++) {
            // perform horizontal un-DPCM
            (*dpcm)[0][0][ind_M][ind_N] += (*dpcm)[0][0][ind_M][ind_N-1];
        }
    }
}
