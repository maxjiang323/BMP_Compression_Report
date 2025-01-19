#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "headerInfo.h"
#include "pixeldata.h"
#include "testouput.h"
#include "mathCalculation.h"
#include "sort.h"
#include "huffmanDecoding.h"

// #define W 8 // dimension of basis vector (width)
// #define H 8 // dimension of basis vector (height)


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

// Function to free memory of the Huffman tree
void freeHuffmanTree(Node *root) {
    if (root == NULL) return;

    // Recursively free left and right subtrees
    freeHuffmanTree(root->left);
    freeHuffmanTree(root->right);

    // Free the current node
    free(root);
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

int main(int argc, char *argv[]){
    Bmpheader header;
    int height, width, M, N, max_size, rM, rN;
    long pos_Y = 0, pos_Cb = 0, pos_Cr = 0;  // pos_ : index of Y, Cb, Cr
    int m, n; 
    short *rle_code_Y = NULL, *rle_code_Cb = NULL, *rle_code_Cr = NULL;

    FILE *codebook = fopen(argv[3], "r");
    Node *root = buildHuffmanTreeFromCodebook(codebook);

    // ascii: read huffman data from hf_txt
    if (strcmp(argv[2], "ascii") == 0){        
        FILE *hf_txt = fopen(argv[4], "r");
        read_dim(hf_txt, &header);
        // printf("%u\n", header.height);

        width = header.width;  
        height = header.height;
        rM = H - (height % H);
        rN = W - (width % W); // the remainder of the height and width
        // may generate non-complete H*W block    
        if (rM == H) rM = 0; // if the height is a multiple of H, no extra row is needed 
        if (rN == W) rN = 0; // if the width is a multiple of W, no extra column is needed 
        M = (height+rM)/H;  // Number of blocks vertically
        N = (width+rN)/W;   // Number of blocks horizontally

        fscanf(hf_txt, "\n\n"); 
        fscanf(hf_txt, "Bitstream: \n");  

        // // read all bitstream
        // char *bitstream = (char *)malloc((datasize + 1)*sizeof(char));
        // fread(bitstream, 1, datasize + 1, hf_txt);
        // bitstream[datasize] = '\0';
        // // printf("%s", bitstream);

        max_size = (((H*W)*2)+2)*M*N; // max size of "rle_code" arrays
        // (skip, value)*H*W + 0 0 (end-of-block marker) in (m, n) block
        rle_code_Y = (short *)malloc(max_size * sizeof(short));
        rle_code_Cb = (short *)malloc(max_size * sizeof(short));
        rle_code_Cr = (short *)malloc(max_size * sizeof(short));
        long posY = 0, posCb = 0, posCr = 0;
                
        decodeHuffmanTxt(root, hf_txt, &rle_code_Y, &rle_code_Cb, &rle_code_Cr, &pos_Y, &pos_Cb, &pos_Cr);

        // free(bitstream);
        fclose(hf_txt);
    }
    // binary: read huffman data from hf_bin
    else if (strcmp(argv[2], "binary") == 0){     
        FILE *hf_bin = fopen(argv[4], "rb");
        readheader(hf_bin, &header);

        width = header.width;  
        height = header.height;
        rM = H - (height % H);
        rN = W - (width % W); // the remainder of the height and width
        // may generate non-complete H*W block    
        if (rM == H) rM = 0; // if the height is a multiple of H, no extra row is needed 
        if (rN == W) rN = 0; // if the width is a multiple of W, no extra column is needed 
        M = (height+rM)/H;  // Number of blocks vertically
        N = (width+rN)/W;   // Number of blocks horizontally

        max_size = (((H*W)*2)+2)*M*N; // max size of "rle_code" arrays
        // (skip, value)*H*W + 0 0 (end-of-block marker) in (m, n) block
        rle_code_Y = (short *)malloc(max_size * sizeof(short));
        rle_code_Cb = (short *)malloc(max_size * sizeof(short));
        rle_code_Cr = (short *)malloc(max_size * sizeof(short));
        long posY = 0, posCb = 0, posCr = 0;

        if (rle_code_Y == NULL || rle_code_Cb == NULL || rle_code_Cr == NULL) {
                fprintf(stderr, "Memory allocation failed.\n");
                exit(EXIT_FAILURE);
        }

        decodeHuffmanBin(root, hf_bin, &rle_code_Y, &rle_code_Cb, &rle_code_Cr, &pos_Y, &pos_Cb, &pos_Cr);              
        fclose(hf_bin);
    }
    freeHuffmanTree(root);
    fclose(codebook);

    // un-rle
    // zz dimension: (M, N, H*W)
    short ***zz_Y = NULL, ***zz_Cb = NULL, ***zz_Cr = NULL;
    malloc_zz(M, N, &zz_Y);
    malloc_zz(M, N, &zz_Cb);
    malloc_zz(M, N, &zz_Cr);

    //perform un-rle
    pos_Y = 0;
    pos_Cb = 0;
    pos_Cr = 0;
    // printf("Debug: M = %d, N = %d\n", M, N);
    for (m = 0; m < M; m++){
        for (n = 0; n < N; n++){
            un_rle(rle_code_Y, zz_Y[m][n], &pos_Y);
            un_rle(rle_code_Cb, zz_Cb[m][n], &pos_Cb);
            un_rle(rle_code_Cr, zz_Cr[m][n], &pos_Cr);
        }
    }

    // test output of un_rle_Y
    // FILE *un_rle_Y_de_txt = fopen("un_rle_Y_de.txt", "w");
    // test_output_of_un_rle(un_rle_Y_de_txt, height, width, zz_Y);
    // fclose(un_rle_Y_de_txt);


    // un-zigzag
    // dpcm dimension: (H, W, M, N)
    short ****dpcm_Y = NULL, ****dpcm_Cb = NULL, ****dpcm_Cr = NULL;
    malloc_dpcm(M, N, &dpcm_Y);
    malloc_dpcm(M, N, &dpcm_Cb);
    malloc_dpcm(M, N, &dpcm_Cr);

    // perform un-zigzag
    un_zigzag(height, width, zz_Y, &dpcm_Y);
    un_zigzag(height, width, zz_Cb, &dpcm_Cb);
    un_zigzag(height, width, zz_Cr, &dpcm_Cr);

    // test output of un_zz_Y
    // FILE *un_zz_Y_de_txt = fopen("un_zz_Y_de.txt", "w");

    // test_output_of_un_zigzag_dpcm(un_zz_Y_de_txt, height, width, dpcm_Y);
    // fclose(un_zz_Y_de_txt);


    // un-DPCM
    // dpcm --> qF
    // qF dimension: (H, W, M, N)
    un_dpcm(height, width, &dpcm_Y);
    un_dpcm(height, width, &dpcm_Cb);
    un_dpcm(height, width, &dpcm_Cr);

    // test output of un_dpcm_Y
    // FILE *un_dpcm_Y_de_txt = fopen("un_dpcm_Y_de.txt", "w");
    // test_output_of_un_zigzag_dpcm(un_dpcm_Y_de_txt, height, width, dpcm_Y);
    // fclose(un_dpcm_Y_de_txt);


    // un-quantization
    // Allocate memory for F_Y, F_Cb, F_Cr dynamically
    // F dimension: (H, W, M, N)
    float ****F_Y = NULL, ****F_Cb = NULL, ****F_Cr = NULL;
    malloc_F(M, N, &F_Y);
    malloc_F(M, N, &F_Cb);
    malloc_F(M, N, &F_Cr);

    // perform un-quantization
    un_quantization(height, width, dpcm_Y, &F_Y, 'Y');
    un_quantization(height, width, dpcm_Cb, &F_Cb, 'C');
    un_quantization(height, width, dpcm_Cr, &F_Cr, 'C');

    // test output of un_qF_Y
    // FILE *un_qF_Y_de = fopen("un_qF_Y_de.txt", "w");
    // test_output_of_un_quantizaton(un_qF_Y_de2, height, width, F_Y);
    // fclose(un_qF_Y_de);

        
    // [IDCT] F_Y -> f_Y; ... 
    // "f" arrays; dimension:(height+rM, width+rN)
    float basis_vector[H][W][H][W];
    float **f_Y = NULL, **f_Cb = NULL, **f_Cr = NULL;
    calloc_f(height+rM, width+rN, &f_Y);
    calloc_f(height+rM, width+rN, &f_Cb);
    calloc_f(height+rM, width+rN, &f_Cr);


    // perform IDCT 
    generate_basis_vector(basis_vector);
    compute_2D_IDCT(height, width, f_Y, F_Y, basis_vector);
    compute_2D_IDCT(height, width, f_Cb, F_Cb, basis_vector);
    compute_2D_IDCT(height, width, f_Cr, F_Cr, basis_vector);

    // test output of IDCT_Y, ...
    // FILE *de_IDCT_Y = fopen("de_IDCT_Y.txt", "w");
    // int block_size = 8;
    // test_output_of_each_block(de_IDCT_Y, width, height, block_size, f_Y);
    // fclose(de_IDCT_Y);


    // [YCbCr2RGB] : R, G, B 
    int row_length;
    unsigned char *pix_rgb = NULL;
    set_rgb_data(width, height, &row_length, &pix_rgb);    

    // test output of G, ...
    // FILE *de_G = fopen("de_G.txt", "w");
    int ind_H, ind_W;
    for (ind_H = 0; ind_H < height; ind_H ++) {
        for (ind_W = 0; ind_W < width; ind_W++) {        
            int index = ind_H * row_length + ind_W * 3;              
            // "f" arrays; dimension:(height+rM, width+rN)
            RGB result = YCbCr2RGB(f_Y[ind_H][ind_W], f_Cb[ind_H][ind_W], f_Cr[ind_H][ind_W]);

            pix_rgb[index] = result.B;
            pix_rgb[index+1] = result.G;
            pix_rgb[index+2] = result.R; 
                
            // test output
            // fprintf(de_G, "%d ", result.G);
        }
        // fprintf(de_G, "\n");
    }
    // fclose(de_G);
        
    // write header information into the bmp file
    FILE *bmp = fopen(argv[1], "wb");
    write_header(bmp, header);

    // write rgb pixel data into bmp file
    fseek(bmp, header.bitmap_dataoffset, SEEK_SET); // jump to the start of pixel data  
    fwrite(pix_rgb, sizeof(unsigned char), row_length * height, bmp);

    free(pix_rgb);
    fclose(bmp);


    // Free allocated memory for f_Y, f_Cb, f_Cr
    free_f(height+rM, f_Y);
    free_f(height+rM, f_Cb);
    free_f(height+rM, f_Cr);

    // Free allocated memory for F_Y, F_Cb, F_Cr
    free_F(M, N, F_Y);
    free_F(M, N, F_Cb);
    free_F(M, N, F_Cr);

    // Free allocated memory for dpcm_Y, dpcm_Cb, dpcm_Cr
    free_dpcm(M, N, dpcm_Y);
    free_dpcm(M, N, dpcm_Cb);
    free_dpcm(M, N, dpcm_Cr);

    // Free allocated memory for zz_Y, zz_Cb, zz_Cr
    free_zz(M, N, zz_Y);
    free_zz(M, N, zz_Cb);
    free_zz(M, N, zz_Cr);

    free(rle_code_Y);
    free(rle_code_Cb); 
    free(rle_code_Cr);

    return 0;
}
