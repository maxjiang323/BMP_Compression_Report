#include <stdio.h>
#include "../headers/testouput.h"

// encoder.c
void test_output_of_dpcm(FILE *fp, int height, int width, short ****dpcm){
    int rM = H - (height % H), rN = W - (width % W); // the remainder of the height and width
    // may generate non-complete H*W block    
    if (rM == H) rM = 0; // if the height is a multiple of H, no extra row is needed 
    if (rN == W) rN = 0; // if the width is a multiple of W, no extra column is needed 
    int M = (height+rM)/H;  // Number of blocks vertically
    int N = (width+rN)/W;   // Number of blocks horizontally

    int ind_H, ind_W, ind_M, ind_N; // index of H/W/M/N
    for (ind_M = 0; ind_M < M; ind_M++) {
        for (ind_N = 0; ind_N < N; ind_N++) {
            for (ind_H = 0; ind_H < H; ind_H++) {
                for (ind_W = 0; ind_W < W; ind_W++) {
                    fprintf(fp, "%hd ", dpcm[ind_H][ind_W][ind_M][ind_N]);
                }
                fprintf(fp, "\n");
            }
            fprintf(fp, "\n\n");
        }
    }
}

void test_output_of_zz(FILE *fp, int height, int width, short ***output){
    int rM = H - (height % H), rN = W - (width % W); // the remainder of the height and width
    // may generate non-complete H*W block    
    if (rM == H) rM = 0; // if the height is a multiple of H, no extra row is needed 
    if (rN == W) rN = 0; // if the width is a multiple of W, no extra column is needed 
    int M = (height+rM)/H;  // Number of blocks vertically
    int N = (width+rN)/W;   // Number of blocks horizontally

    int ind_H_W, ind_M, ind_N;
    for (ind_M = 0; ind_M < M; ind_M++){
        for (ind_N = 0; ind_N < N; ind_N++){
            for (ind_H_W = 0; ind_H_W < 64; ind_H_W++){
                fprintf(fp, "%hd ", output[ind_H_W][ind_M][ind_N]);
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n\n");
    }
}


// decoder.c
void test_output_of_each_block(FILE *fp, int width, int height, int block_size, float **output){
    int rM = H - (height % H), rN = W - (width % W); // the remainder of the height and width
    // may generate non-complete H*W block    
    if (rM == H) rM = 0; // if the height is a multiple of H, no extra row is needed 
    if (rN == W) rN = 0; // if the width is a multiple of W, no extra column is needed 

    int block_x, block_y, x, y;
    for (block_x = 0; block_x < height+rM; block_x += block_size) {
        for (block_y = 0; block_y < width+rN; block_y += block_size) {
            for (x = 0; x < block_size; x++) {
                for (y = 0; y < block_size; y++) {
                    // Check if we are within the image bounds
                    if ((block_x + x) < height && (block_y + y) < width) {

                    } else {
                        output[block_x + x][block_y + y] = 0.0;
                    }
                    fprintf(fp, "%f ", output[block_x + x][block_y + y]);
                }
                fprintf(fp, "\n");
            }
            fprintf(fp, "\n\n");
        }
    }
}

void test_output_of_rle_txt(FILE *rle_txt, FILE* rle_out, long *read_pos, long *pos, short **rle_code){ 
    *read_pos = ftell(rle_txt); // the realistic position of the skips and values

    char buffer[100];
    fgets(buffer, 100, rle_txt);
    // printf("%s", buffer);
    if (buffer[0] == 'e'){
        // empty
    }else{
        fseek(rle_txt, *read_pos, SEEK_SET);  // back to the position before the buffer read
        while (fscanf(rle_txt, "skip %hd value %hd ", &(*rle_code)[(*pos)], &(*rle_code)[(*pos)+1]) == 2) {
            fprintf(rle_out, "%hd %hd ", (*rle_code)[(*pos)], (*rle_code)[(*pos)+1]);
            (*pos) += 2;
        }
    }
    fscanf(rle_txt, "\n");
    (*rle_code)[(*pos)] = 0;
    (*rle_code)[(*pos) + 1] = 0;
    fprintf(rle_out, "%hd %hd ", (*rle_code)[(*pos)], (*rle_code)[(*pos)+1]);
    fprintf(rle_out, "\n");
    (*pos) += 2;
}

void test_output_of_rle_bin(FILE *rle_bin, FILE *rle_bin_out, long *pos, short **rle_code){
    int flag = 0;  // 0: continue; 1: end of a block
    while (flag == 0){
        // check the position of the file (debug)
        // long current_pos = ftell(rle_bin);
        // fprintf(stderr, "Current file position: %ld\n", current_pos);
                        
        fread(&(*rle_code)[(*pos)], sizeof(short), 1, rle_bin);
        fread(&(*rle_code)[(*pos)+1], sizeof(short), 1, rle_bin);
        // test output to debug
        // fprintf(stderr, "Y data (%d): %hd %hd\n", pos, (*rle_code)[(*pos)], (*rle_code)[(*pos)+1]);
        // fprintf(stderr, "Processing block (%d, %d)\n", m, n);

        fprintf(rle_bin_out, "%hd %hd ", (*rle_code)[(*pos)], (*rle_code)[(*pos)+1]);

        if ((*rle_code)[(*pos)] == 0 && (*rle_code)[(*pos)+1] == 0){
            flag = 1;
            fprintf(rle_bin_out, "\n");
        }
        (*pos) += 2;
    }
}

void test_output_of_un_rle(FILE *fp, int height, int width, short ***zigzag){
    int rM = H - (height % H), rN = W - (width % W); // the remainder of the height and width
    // may generate non-complete H*W block    
    if (rM == H) rM = 0; // if the height is a multiple of H, no extra row is needed 
    if (rN == W) rN = 0; // if the width is a multiple of W, no extra column is needed 
    int M = (height+rM)/H;  // Number of blocks vertically
    int N = (width+rN)/W;   // Number of blocks horizontally

    int ind_H_W, ind_M, ind_N;
    for (ind_M = 0; ind_M < M; ind_M++){
        for (ind_N = 0; ind_N < N; ind_N++){
            for (ind_H_W = 0; ind_H_W < 64; ind_H_W++){
                if (zigzag[ind_M][ind_N][ind_H_W] < -255 || zigzag[ind_M][ind_N][ind_H_W] > 255) { 
                    fprintf(stderr, "Error: Value out of range: %hd\n", zigzag[ind_M][ind_N][ind_H_W]); 
                }
                fprintf(fp, "%hd ", zigzag[ind_M][ind_N][ind_H_W]);
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n\n");
    }
}

void test_output_of_un_zigzag_dpcm(FILE *fp, int height, int width, short ****output){
    int rM = H - (height % H), rN = W - (width % W); // the remainder of the height and width
    // may generate non-complete H*W block    
    if (rM == H) rM = 0; // if the height is a multiple of H, no extra row is needed 
    if (rN == W) rN = 0; // if the width is a multiple of W, no extra column is needed 
    int M = (height+rM)/H;  // Number of blocks vertically
    int N = (width+rN)/W;   // Number of blocks horizontally

    int ind_H, ind_W, ind_M, ind_N; // index of H/W/M/N
    for (ind_M = 0; ind_M < M; ind_M++) {
        for (ind_N = 0; ind_N < N; ind_N++) {
            for (ind_H = 0; ind_H < H; ind_H++) {
                for (ind_W = 0; ind_W < W; ind_W++) {
                    fprintf(fp, "%hd ", output[ind_H][ind_W][ind_M][ind_N]);
                }
                fprintf(fp, "\n");
            }
            fprintf(fp, "\n\n");
        }
    }
}

void test_output_of_un_quantizaton(FILE *fp, int height, int width, float ****F){
    int rM = H - (height % H), rN = W - (width % W); // the remainder of the height and width
    // may generate non-complete H*W block    
    if (rM == H) rM = 0; // if the height is a multiple of H, no extra row is needed 
    if (rN == W) rN = 0; // if the width is a multiple of W, no extra column is needed 
    int M = (height+rM)/H;  // Number of blocks vertically
    int N = (width+rN)/W;   // Number of blocks horizontally

    int ind_H, ind_W, ind_M, ind_N; // index of H/W/M/N
    for (ind_M = 0; ind_M < M; ind_M++) {
        for (ind_N = 0; ind_N < N; ind_N++) {
            for (ind_H = 0; ind_H < H; ind_H++) {
                for (ind_W = 0; ind_W < W; ind_W++) {
                    fprintf(fp, "%f ", F[ind_H][ind_W][ind_M][ind_N]);
                }
                fprintf(fp, "\n");
            }
            fprintf(fp, "\n\n");
        }
    }
}
