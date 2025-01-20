#include <stdio.h>
#include <stdlib.h>
#include "headers/runLengthCoding.h"

// encoder.c
void run_length_encoding(int height, int width, short ***zz, short **rle_code, long *size) {
    int rM = H - (height % H), rN = W - (width % W); // the remainder of the height and width
    // may generate non-complete H*W block    
    if (rM == H) rM = 0; // if the height is a multiple of H, no extra row is needed 
    if (rN == W) rN = 0; // if the width is a multiple of W, no extra column is needed 
    int M = (height+rM)/H;  // Number of blocks vertically
    int N = (width+rN)/W;   // Number of blocks horizontally
    int max_rle_length = H*W*2+2; // Maximum size for each block (worst-case scenario)
    int m, n, i;
    // zz dimension: (8x8, M, N); rle dimension: (M, N, H*W+2 --> idx)
    short ***temp_rle = (short ***)malloc(M * sizeof(short **)); // RLE blocks for each (m, n)

    for (m = 0; m < M; m++) {
        temp_rle[m] = (short **)malloc(N * sizeof(short *));
        for (n = 0; n < N; n++) {
            temp_rle[m][n] = (short *)malloc(max_rle_length * sizeof(short));

            int idx = 0; // Index for RLE array
            int nzero = 0; // Count of consecutive zeros

            for (i = 0; i < H * W; i++) {
                if (zz[i][m][n] != 0) {
                    temp_rle[m][n][idx] = nzero; // Count of consecutive zeros
                    idx += 1;

                    temp_rle[m][n][idx] = zz[i][m][n];  // non-zero value
                    idx += 1;

                    nzero = 0; // reset
                } else {
                    nzero += 1;
                }
            }
            temp_rle[m][n][idx] = 0; 
            idx += 1;
            temp_rle[m][n][idx] = 0; // End of block marker --> 0 0
            idx += 1;
            temp_rle[m][n] = realloc(temp_rle[m][n], idx * sizeof(short)); // Adjust to actual size
        }
    }

    // Flatten RLE blocks into a single array
    int total_size = 0;
    for (m = 0; m < M; ++m) {
        for (n = 0; n < N; ++n) {
            int idx = 0;
            while (temp_rle[m][n][idx] != 0 || temp_rle[m][n][idx + 1] != 0) { // Find block length
                idx += 2;
            }
            idx += 2; // Include end-of-block marker
            total_size += idx;
        }
    }

    *size = total_size;
    (*rle_code) = malloc(total_size * sizeof(short)); 
    long pos = 0;
    for (m = 0; m < M; ++m) {
        for (n = 0; n < N; ++n) {
            int idx = 0;
            while (temp_rle[m][n][idx] != 0 || temp_rle[m][n][idx + 1] != 0) { // not "0 0"
                (*rle_code)[pos] = temp_rle[m][n][idx]; 
                // (*(&rle_code))[pos], &rle_code --> input
                pos += 1;
                idx += 1;
                (*rle_code)[pos] = temp_rle[m][n][idx];
                pos += 1;
                idx += 1;
            }
            (*rle_code)[pos] = temp_rle[m][n][idx]; // Copy end-of-block markers
            pos += 1;
            idx += 1;
            (*rle_code)[pos] = temp_rle[m][n][idx];
            pos += 1;
            idx += 1;
        }
    }

    // Free "temp_rle" dynamically 
    for (m = 0; m < M; ++m) {
        for (n = 0; n < N; ++n) {
            free(temp_rle[m][n]);
        }
        free(temp_rle[m]);
    }
    free(temp_rle);
}


// decoder.c
void un_rle(short *rle_code, short *zz_mn, long *pos) {
    // zz dimension: (M, N, H*W)
    int nzero = 64; // the remaining count of zeros
    int index = 0; // the index in each (m, n) block
    short skip, value;
    while (rle_code[(*pos)] != 0 || rle_code[(*pos) + 1] != 0) {
        skip = rle_code[(*pos)];
        value =  rle_code[(*pos) + 1];
        nzero -= (skip + 1);

        int i, j;
        for (i = 0; i < skip; i++){
            zz_mn[index] = 0;
            index += 1;
        }
        zz_mn[index] = value;
        index += 1;
        (*pos) += 2;
    }
    (*pos) += 2;  // end-of-block marker 0 0

    int k;
    for (k = 0; k < nzero; k++){
        zz_mn[index] = 0;
        index += 1;
    }
}
