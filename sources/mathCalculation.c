#include <math.h>
#include "headers/mathCalculation.h"

#define W 8 // dimension of basis vector (width)
#define H 8 // dimension of basis vector (height)
#define Pi 3.14159265359

// Const quantization table 
// Y channel
const int Qt_Y[8][8] = {
    {16, 11, 10, 16, 24, 40, 51, 61},
    {12, 12, 14, 19, 26, 58, 60, 55},
    {14, 13, 16, 24, 40, 57, 69, 56},
    {14, 17, 22, 29, 51, 87, 80, 62},
    {18, 22, 37, 56, 68, 109, 103, 77},
    {24, 35, 55, 64, 81, 104, 113, 92},
    {49, 64, 78, 87, 103, 121, 120, 101},
    {72, 92, 95, 98, 112, 100, 103, 99}
};

// Cb, Cr channel
const int Qt_C[8][8] = {
    {17, 18, 24, 47, 99, 99, 99, 99},
    {18, 21, 26, 66, 99, 99, 99, 99},
    {24, 26, 56, 99, 99, 99, 99, 99},
    {47, 66, 99, 99, 99, 99, 99, 99},
    {99, 99, 99, 99, 99, 99, 99, 99},
    {99, 99, 99, 99, 99, 99, 99, 99},
    {99, 99, 99, 99, 99, 99, 99, 99},
    {99, 99, 99, 99, 99, 99, 99, 99},
};


// encoder.c
void generate_basis_vector(float basis_vector[H][W][H][W]){
    int u, v, r, c;
    for (u = 0; u < H; u++) {
        for (v = 0; v < W; v++) {
            for (r = 0; r < H; r++) {
                for (c = 0; c < W; c++) {
                    basis_vector[r][c][u][v] = cos(Pi*u*(2*r+1)/2.0/H)*cos(Pi*v*(2*c+1)/2.0/W);
                }
            }
        }
    }
}

void compute_2D_DCT(int height, int width, float **f, float ****F, float basis_vector[H][W][H][W]){
    int rM = H - (height % H), rN = W - (width % W); // the remainder of the height and width
    // may generate non-complete H*W block    
    if (rM == H) rM = 0; // if the height is a multiple of H, no extra row is needed 
    if (rN == W) rN = 0; // if the width is a multiple of W, no extra column is needed 
    int M = (height+rM)/H;  // Number of blocks vertically
    int N = (width+rN)/W;   // Number of blocks horizontally
    
    float beta[H];
    int i;
    for (i = 0; i < H; i++) {
        beta[i] = 1.0;
    }
    beta[0] = 1.0 / sqrt(2);

    int m, n, j, u, v;
    for (m = 0; m < M; m++) {  // loop for # of blocks vertically
        for (n = 0; n < N; n++) {  // loop for # of blocks horizontally
            float target_8x8[H][W];
            // Extract target 8x8 block
            for (i = 0; i < H; i++) {
                for (j = 0; j < W; j++) {
                    target_8x8[i][j] = f[(m*H)+i][(n*W)+j] - 128;  // subtract 128
                }
            }

            // Loop through the frequency components
            for (u = 0; u < H; u++) {  // frequencies in the vertical direction
                for (v = 0; v < W; v++) {  // frequencies in the horizontal direction
                    F[u][v][m][n] = 2.0 / sqrt(H * W) * beta[u] * beta[v];

                    // Compute the DCT projection on the basis vector
                    float sum = 0.0;
                    int r, c;
                    for (r = 0; r < H; r++) {
                        for (c = 0; c < W; c++) {
                            sum += target_8x8[r][c] * basis_vector[r][c][u][v];
                        }
                    }
                    F[u][v][m][n] *= sum;  // Store the result
                }
            }
        }
    }
}

void quantization(int height, int width, short *****qF, float ****F, char channel){
    // add FILE *qF_en to parameter if want to see the output
    int ind_M, ind_N, ind_H, ind_W;
    int rM = H - (height % H), rN = W - (width % W); // the remainder of the height and width
    // may generate non-complete H*W block    
    if (rM == H) rM = 0; // if the height is a multiple of H, no extra row is needed 
    if (rN == W) rN = 0; // if the width is a multiple of W, no extra column is needed 
    int M = (height+rM)/H;  // Number of blocks vertically
    int N = (width+rN)/W;   // Number of blocks horizontally
    
    for (ind_M = 0; ind_M < M; ind_M++) {
        for (ind_N = 0; ind_N < N; ind_N++) {
            for (ind_H = 0; ind_H < H; ind_H ++) {
                for (ind_W = 0; ind_W < W; ind_W++) {
                    if (channel == 'Y'){
                        (*qF)[ind_H][ind_W][ind_M][ind_N] = (short)round(F[ind_H][ind_W][ind_M][ind_N]/Qt_Y[ind_H][ind_W]);
                    }
                    else if (channel == 'C'){
                        (*qF)[ind_H][ind_W][ind_M][ind_N] = (short)round(F[ind_H][ind_W][ind_M][ind_N]/Qt_C[ind_H][ind_W]);
                    }                   
                    // test output
                    // fprintf(qF_en, "%hd ", qF[ind_H][ind_W][ind_M][ind_N]);
                }
                // fprintf(qF_en, "\n");
            }
            // fprintf(qF_en, "\n\n");
        }
    }
}


// decoder.c
void un_quantization(int height, int width, short ****qF, float *****F, char channel){
    int ind_M, ind_N, ind_H, ind_W;
    int rM = H - (height % H), rN = W - (width % W); // the remainder of the height and width
    // may generate non-complete H*W block    
    if (rM == H) rM = 0; // if the height is a multiple of H, no extra row is needed 
    if (rN == W) rN = 0; // if the width is a multiple of W, no extra column is needed 
    int M = (height+rM)/H;  // Number of blocks vertically
    int N = (width+rN)/W;   // Number of blocks horizontally
    
    for (ind_M = 0; ind_M < M; ind_M++) {
        for (ind_N = 0; ind_N < N; ind_N++) {
            for (ind_H = 0; ind_H < H; ind_H ++) {
                for (ind_W = 0; ind_W < W; ind_W++) {
                    if (channel == 'Y'){
                        (*F)[ind_H][ind_W][ind_M][ind_N] = qF[ind_H][ind_W][ind_M][ind_N]*Qt_Y[ind_H][ind_W];
                    }
                    else if (channel == 'C'){
                        (*F)[ind_H][ind_W][ind_M][ind_N] = qF[ind_H][ind_W][ind_M][ind_N]*Qt_C[ind_H][ind_W];
                    }                   
                }
            }
        }
    }
}

void compute_2D_IDCT(int height, int width, float **f, float ****F, float basis_vector[H][W][H][W]) {
    int rM = H - (height % H), rN = W - (width % W); // the remainder of the height and width
    // may generate non-complete H*W block    
    if (rM == H) rM = 0; // if the height is a multiple of H, no extra row is needed 
    if (rN == W) rN = 0; // if the width is a multiple of W, no extra column is needed 
    int M = (height+rM)/H;  // Number of blocks vertically
    int N = (width+rN)/W;   // Number of blocks horizontally
    
    float beta[H];
    int i;
    for (i = 0; i < H; i++) {
        beta[i] = 1.0;
    }
    beta[0] = 1.0 / sqrt(2);

    int m, n, j, u, v;
    for (m = 0; m < M; m++) {  // loop for # of blocks vertically
        for (n = 0; n < N; n++) {  // loop for # of blocks horizontally
            float reconstructed_8x8[H][W] = {0};

            // Reconstruct the spatial domain block from frequency domain coefficients
            int r, c;
            for (r = 0; r < H; r++) {  // Rows in the spatial domain
                for (c = 0; c < W; c++) {  // Columns in the spatial domain
                    float sum = 0.0;

                    // Perform the summation over all frequency components
                    for (u = 0; u < H; u++) {
                        for (v = 0; v < W; v++) {
                            sum += beta[u] * beta[v] * F[u][v][m][n] * basis_vector[r][c][u][v];
                        }
                    }

                    // Normalize and store the result
                    reconstructed_8x8[r][c] = 2.0 / sqrt(H * W) * sum;
                }
            }

            // Write back the reconstructed block to the original image array
            for (i = 0; i < H; i++) {
                for (j = 0; j < W; j++) {
                    f[(m * H) + i][(n * W) + j] = reconstructed_8x8[i][j] + 128;  // Add back 128
                }
            }
        }
    }
}
