#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SQNR.h"     
#include "definitions.h"


void calculate_Ps_Pe(int height, int width, float Ps[H][W], float Pe[H][W], float ****F, char channel) {      
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
                    // qF_ and eF_ value
                    short qF;
                    float eF;
                    if (channel == 'Y'){
                        qF = (short)round(F[ind_H][ind_W][ind_M][ind_N]/Qt_Y[ind_H][ind_W]);
                        eF = F[ind_H][ind_W][ind_M][ind_N] - qF * Qt_Y[ind_H][ind_W];
                    }
                    else if (channel == 'C'){
                        qF = (short)round(F[ind_H][ind_W][ind_M][ind_N]/Qt_C[ind_H][ind_W]);
                        eF = F[ind_H][ind_W][ind_M][ind_N] - qF * Qt_C[ind_H][ind_W];
                    }
                    else{
                        printf("Error: channel must be 'Y' or 'C'\n");
                        exit(1);
                    }
                        
                    // calculate Ps_ and Pe_
                    Ps[ind_H][ind_W] += F[ind_H][ind_W][ind_M][ind_N] * F[ind_H][ind_W][ind_M][ind_N];
                    Pe[ind_H][ind_W] += eF * eF;                     
                }
            }
        }
    }
}

void write_SQNR(FILE *fp, float sqnr[H][W], float Ps[H][W], float Pe[H][W]){
    int ind_H, ind_W;
    for (ind_H = 0; ind_H < H; ind_H++) {
        for (ind_W = 0; ind_W < W; ind_W++) {
            sqnr[ind_H][ind_W] = 10 * log10(Ps[ind_H][ind_W] / Pe[ind_H][ind_W]);
            fprintf(fp, "%f ", sqnr[ind_H][ind_W]);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n\n");
}
