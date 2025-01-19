#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "headerInfo.h"
#include "pixeldata.h"
#include "testouput.h"
#include "mathCalculation.h"

#define W 8 // dimension of basis vector (width)
#define H 8 // dimension of basis vector (height)
#define Pi 3.14159265359


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

void write_rle_txt(FILE *rle_txt, long *pos, short **rle_code){
    if ((*rle_code)[(*pos)] == 0 && (*rle_code)[(*pos) + 1] == 0) {
        fprintf(rle_txt, "empty\n");
    } else {
        while ((*rle_code)[(*pos)] != 0 || (*rle_code)[(*pos) + 1] != 0) {
            fprintf(rle_txt, "skip %hd value %hd ", (*rle_code)[(*pos)], (*rle_code)[(*pos) + 1]);
            (*pos) += 2;
        }
        fprintf(rle_txt, "\n");
    }
    (*pos) += 2; // skip end-of-block marker (0 0)
}

void write_rle_bin(FILE *rle_bin, long *pos, short **rle_code, long *size_bin, long *size_bin_channel){
    while ((*rle_code)[(*pos)] != 0 || (*rle_code)[(*pos) + 1] != 0) {
        fwrite(&(*rle_code)[(*pos)], sizeof(short), 1, rle_bin); // skip
        fwrite(&(*rle_code)[(*pos) + 1], sizeof(short), 1, rle_bin); // value
        (*size_bin) += sizeof(short)*2;
        (*size_bin_channel) += sizeof(short)*2;
        (*pos) += 2;
    }
    short end_marker = 0;
    fwrite(&end_marker, sizeof(short), 1, rle_bin); // end-of-block marker
    fwrite(&end_marker, sizeof(short), 1, rle_bin); // end-of-block marker
    (*pos) += 2;
    (*size_bin) += sizeof(short)*2;
    (*size_bin_channel) += sizeof(short)*2;
}


#define MAX_RANGE 512 // -255 to 255 + one extra symbol for "0 0"

// Define the node structure for the Huffman tree
typedef struct Node {
    int value;         // Frequency value of the node
    int symbol;        // Corresponding symbol (-255 ~ 255 or 511 for "0 0")
    struct Node *left; // Left child node
    struct Node *right; // Right child node
} Node;

// Create a new node
Node *createNode(int value, int symbol, Node *left, Node *right) {
    Node *newNode = (Node *)malloc(sizeof(Node));
    newNode->value = value;
    newNode->symbol = symbol;
    newNode->left = left;
    newNode->right = right;
    return newNode;
}

void countFrequency(long size, int (*frequencies)[512], short *data){
    int symbol, prevSymbol = 256;
    int i;
    for (i = 0; i < size; i++) {
        symbol = data[i];
        if (prevSymbol == 0 && symbol == 0) {
            (*frequencies)[511]++;
            prevSymbol = 256; // Reset
        } else {
            if (prevSymbol != 256) {
                (*frequencies)[prevSymbol + 255]++;
            }
            prevSymbol = symbol;
        }
    }
    if (prevSymbol != 256) {
        (*frequencies)[prevSymbol + 255]++;
    }
}

// Find the two nodes with the smallest frequencies
void findTwoSmallest(Node *nodes[], int size, int *first, int *second) {
    *first = -1;
    *second = -1;
    int i;
    for (i = 0; i < size; i++) {
        if (nodes[i] == NULL) continue;
        if (*first == -1 || nodes[i]->value < nodes[*first]->value) {
            *second = *first;
            *first = i;
        } else if (*second == -1 || nodes[i]->value < nodes[*second]->value) {
            *second = i;
        }
    }
}

// Build the Huffman tree
Node *buildHuffmanTree(int frequencies[], int size) {
    Node *nodes[MAX_RANGE] = {NULL};

    // Initialize nodes
    int i;
    for (i = 0; i < size; i++) {
        if (frequencies[i] > 0) {
            nodes[i] = createNode(frequencies[i], i - 255, NULL, NULL);
        }
    }

    // Merge nodes until only one tree remains
    int first;
    while (1) {
        int second;
        findTwoSmallest(nodes, size, &first, &second);

        if (second == -1) break; // Only one tree left

        // Merge the two smallest nodes
        Node *newNode = createNode(
            nodes[first]->value + nodes[second]->value, 
            -1, 
            nodes[first], 
            nodes[second]
        );
        nodes[first] = newNode;
        nodes[second] = NULL;
    }

    return nodes[first];
}

// Populate the Huffman code table
void generateCodes(Node *tree, char *codeTable[], char *currentCode, int depth) {
    if (!tree) return;

    // If it's a leaf node, store the code
    if (tree->left == NULL && tree->right == NULL) {
        currentCode[depth] = '\0';
        codeTable[tree->symbol + 255] = strdup(currentCode);
        // codeTable[tree->symbol + 255] --> pointer
        return;
    }

    // Recur for the left subtree
    currentCode[depth] = '0';
    generateCodes(tree->left, codeTable, currentCode, depth + 1);

    // Recur for the right subtree
    currentCode[depth] = '1';
    generateCodes(tree->right, codeTable, currentCode, depth + 1);
}

// Recursively free the nodes of the Huffman tree
void freeHuffmanTree(Node *tree) {
    if (!tree) return;

    // First, free the child nodes
    freeHuffmanTree(tree->left);
    freeHuffmanTree(tree->right);

    // Finally, free the current node
    free(tree);
}

// Free the Huffman code table
void freeCodeTable(char *codeTable[]) {
    int i;
    for (i = 0; i < MAX_RANGE; i++) {
        if (codeTable[i] != NULL) {
            free(codeTable[i]); // Free the dynamically allocated string
            codeTable[i] = NULL; // Avoid dangling pointers
        }
    }
}


void writeCodebook(FILE *fp, char *codeTable[], int frequencies[]) {
    if (!fp) {
        printf("Unable to open codebook file\n");
        return;
    }
    int i;
    for (i = 0; i < MAX_RANGE; i++) {
        if (codeTable[i] != NULL) {
            if (i == 511) {
                fprintf(fp, "Symbol: \"0 0\", Frequency: %d, Code: %s\n", frequencies[i], codeTable[i]);
            } else {
                fprintf(fp, "Symbol: %d, Frequency: %d, Code: %s\n", i - 255, frequencies[i], codeTable[i]);
            }
        }
    }
}

void write_hf_txt(short **rle_code, FILE *fp, char *codeTable[], long *pos) {
    if ((*rle_code)[(*pos)] == 0 && (*rle_code)[(*pos)+1] == 0) {
        // empty
    } else {
        while ((*rle_code)[(*pos)] != 0 || (*rle_code)[(*pos)+1] != 0) {
            fprintf(fp, "%s%s", codeTable[(*rle_code)[(*pos)] + 255], codeTable[(*rle_code)[(*pos)+1] + 255]);            
            (*pos) += 2;
        }
    }
    fprintf(fp, "%s", codeTable[511]); // Handle "0 0"
    (*pos) += 2; // skip end-of-block marker (0 0)
}

void bit_alignment(FILE *fp, char *code, unsigned char *byte, int *temp, long *bitCount, long *byteCount){
    int j;
    for (j = 0; code[j] != '\0'; j++){
        (*byte) |= (code[j] - '0') << (7 - (*temp));
        (*temp) += 1;
        (*bitCount) += 1;

        if ((*temp) == 8) {
            fwrite(&(*byte), sizeof(unsigned char), 1, fp);
            (*temp) = 0;
            (*byte) = 0;
            (*byteCount) += 1;
        }    
    }
}

void write_hf_bin(short **rle_code, FILE *fp, char *codeTable[], long *pos, long *bitCount, unsigned char *byte, long *byteCount, int *temp) {
    char *code1, *code2;

    if ((*rle_code)[(*pos)] == 0 && (*rle_code)[(*pos)+1] == 0) { 
        // empty           
    }
    else {
        while ((*rle_code)[(*pos)] != 0 || (*rle_code)[(*pos)+1] != 0) {
            code1 = codeTable[(*rle_code)[(*pos)] + 255];
            code2 = codeTable[(*rle_code)[(*pos)+1] + 255]; 

            bit_alignment(fp, code1, &(*byte), &(*temp), &(*bitCount), &(*byteCount));     
            bit_alignment(fp, code2, &(*byte), &(*temp), &(*bitCount), &(*byteCount));    

            (*pos) += 2;
        }
    }
    bit_alignment(fp, codeTable[511], &(*byte), &(*temp), &(*bitCount), &(*byteCount));  // 0 0     
    (*pos) += 2; // skip end-of-block marker (0 0)

}


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


int main(int argc, char *argv[]){
    FILE *bmp = fopen(argv[1], "rb");                   
    Bmpheader header;
    readheader(bmp, &header); // read header from bmp files


    // RGB2YCbCr
    int width = header.width;  
    int height = header.height;
    unsigned char *pix_rgb = NULL;
    int row_length;
    read_rgb_data(bmp, header, &width, &height, &row_length, &pix_rgb);      

    // write YCbCr data to the txt files
    // and record the spatial domain data to "f" arrays
    int block_x, block_y, x, y;

    int rM = H - (height % H), rN = W - (width % W); // the remainder of the height and width
    // may generate non-complete H*W block    
    if (rM == H) rM = 0; // if the height is a multiple of H, no extra row is needed 
    if (rN == W) rN = 0; // if the width is a multiple of W, no extra column is needed 

    // "f" arrays; dimension:(height+rM, width+rN)
    float **f_Y = NULL, **f_Cb = NULL, **f_Cr = NULL;
    calloc_f(height+rM, width+rN, &f_Y);
    calloc_f(height+rM, width+rN, &f_Cb);
    calloc_f(height+rM, width+rN, &f_Cr);

    // test output of f_Y, ...
    // FILE *en_f_Y = fopen("en_f_Y.txt", "w");
    int block_size = 8;
    // may generate non-complete H*W block 
    put_RGB2YCbCr_into_block(height, width, block_size, row_length, pix_rgb, &f_Y, &f_Cb, &f_Cr);
    // fclose(en_f_Y);


    // generate basis vector for 2D-DCT
    float basis_vector[H][W][H][W];
    generate_basis_vector(basis_vector);

    // 2D-DCT
    // "F" arrays dimension: (H, W, M, N)
    // Allocate memory for F_Y, F_Cb, F_Cr, qF_Y, qF_Cb, qF_Cr dynamically
    int M = (height+rM)/H;  // Number of blocks vertically
    int N = (width+rN)/W;   // Number of blocks horizontally

    float ****F_Y = NULL, ****F_Cb = NULL, ****F_Cr = NULL;
    malloc_F(M, N, &F_Y);
    malloc_F(M, N, &F_Cb);
    malloc_F(M, N, &F_Cr);

    compute_2D_DCT(height, width, f_Y, F_Y, basis_vector);
    compute_2D_DCT(height, width, f_Cb, F_Cb, basis_vector);
    compute_2D_DCT(height, width, f_Cr, F_Cr, basis_vector);


    // Quantization
    short ****qF_Y = NULL, ****qF_Cb = NULL, ****qF_Cr = NULL;
    malloc_qF(M, N, &qF_Y);
    malloc_qF(M, N, &qF_Cb);
    malloc_qF(M, N, &qF_Cr);


    // perform quantization
    quantization(height, width, &qF_Y, F_Y, 'Y');
    quantization(height, width, &qF_Cb, F_Cb, 'C');
    quantization(height, width, &qF_Cr, F_Cr, 'C');
    // test output of qF_Y
    // FILE *qF_Y_en = fopen("qF_Y_en.txt", "w");
    // test_output_of_dpcm(qF_Y_en, height, width, qF_Y);
    // fclose(qF_Y_en);


    // DPCM
    dpcm(height, width, &qF_Y);
    dpcm(height, width, &qF_Cb);
    dpcm(height, width, &qF_Cr);

    // test output of dpcm_Y
    // FILE *dpcm_Y_en = fopen("dpcm_Y_en.txt", "w");
    // test_output_of_dpcm(dpcm_Y_en, height, width, qF_Y);
    // fclose(dpcm_Y_en);

    // zigzag
    // zz dimension: (H*W, M, N)
    short ***zz_Y = NULL, ***zz_Cb = NULL, ***zz_Cr = NULL;
    malloc_zz(M, N, &zz_Y);
    malloc_zz(M, N, &zz_Cb);
    malloc_zz(M, N, &zz_Cr);

    // perform zigzag
    apply_zigzag(height, width, qF_Y, zz_Y);
    apply_zigzag(height, width, qF_Cb, zz_Cb);
    apply_zigzag(height, width, qF_Cr, zz_Cr);

        
    // test output of zz_Y
    // FILE *zz_Y_en_txt = fopen("zz_Y_en.txt", "w");
    // test_output_of_zz(zz_Y_en_txt, height, width, zz_Y);
    // fclose(zz_Y_en_txt);        

    // Run Length Encoding (RLE)
    // declare "rle_code" array
    // rle dimension: ((((H*W)*2)+2)*M*N)
    int max_size = (((H*W)*2)+2)*M*N; // max size of "rle_code" arrays
    short *rle_code_Y = NULL, *rle_code_Cb = NULL, *rle_code_Cr = NULL;
    long size_Y, size_Cb, size_Cr;
    // perform RLE
    run_length_encoding(height, width, zz_Y, &rle_code_Y, &size_Y);
    run_length_encoding(height, width, zz_Cb, &rle_code_Cb, &size_Cb);
    run_length_encoding(height, width, zz_Cr, &rle_code_Cr, &size_Cr);

   
    // Part 3: write the codebook.txt
    int frequencies[MAX_RANGE] = {0}; // Store the frequency of each number (-255 to 255 + "0 0")
    char *codeTable[MAX_RANGE] = {NULL}; // string array, ex: {"a, b", "c", "de"}
    char currentCode[256];
    Node *huffmanTree;

    // count
    countFrequency(size_Y, &frequencies, rle_code_Y);
    countFrequency(size_Cb, &frequencies, rle_code_Cb);
    countFrequency(size_Cr, &frequencies, rle_code_Cr);

    // Build the Huffman tree
    huffmanTree = buildHuffmanTree(frequencies, MAX_RANGE);

    // Generate the Huffman code table
    generateCodes(huffmanTree, codeTable, currentCode, 0);

    // Write (symbol, frequency, codeword) to the codebook
    FILE *codebook = fopen(argv[3], "w");
    writeCodebook(codebook, codeTable, frequencies);
    fclose(codebook);

    // ascii: write huffman data to hf_txt
    if (strcmp(argv[2], "ascii") == 0){
        FILE *hf_txt = fopen(argv[4], "w");
        write_dim(hf_txt, header); // write header information to hf_txt 
        fprintf(hf_txt, "\n\n"); // jump to row 3
        // write bitstream to hf_txt
        fprintf(hf_txt, "Bitstream: \n");
        long pos_Y = 0, pos_Cb = 0, pos_Cr = 0;
        int m, n;
        for (m = 0; m < M; m++) {
            for (n = 0; n < N; n++) {
                // Channel Y
                write_hf_txt(&rle_code_Y, hf_txt, codeTable, &pos_Y);

                // Channel Cb
                write_hf_txt(&rle_code_Cb, hf_txt, codeTable, &pos_Cb);

                // Channel Cr
                write_hf_txt(&rle_code_Cr, hf_txt, codeTable, &pos_Cr);
            }
        }
        fclose(hf_txt);
    }    

    // binary: write huffman data to hf_bin
    else if(strcmp(argv[2], "binary") == 0){
        FILE *hf_bin = fopen(argv[4], "wb");
        write_header(hf_bin, header); // write header information to hf_bin
        long pos_Y = 0, pos_Cb = 0, pos_Cr = 0, bitCount = 0, byteCount = 0;
        unsigned char byte = 0;
        int m, n;
        int temp = 0;
        for (m = 0; m < M; m++) {
            for (n = 0; n < N; n++) {
                // Channel Y
                write_hf_bin(&rle_code_Y, hf_bin, codeTable, &pos_Y, &bitCount, &byte, &byteCount, &temp);

                // Channel Cb
                write_hf_bin(&rle_code_Cb, hf_bin, codeTable, &pos_Cb, &bitCount, &byte, &byteCount, &temp);


                // Channel Cr
                write_hf_bin(&rle_code_Cr, hf_bin, codeTable, &pos_Cr, &bitCount, &byte, &byteCount, &temp);
            }
        }

        int k = 8 - bitCount % 8; // count of padding bits
        if (k != 0){
            int j;
            for (j = 8-k; j < 8; j++){
                byte |= 1 << (7-j);
                bitCount += 1;
            }
            fwrite(&byte, sizeof(unsigned char), 1, hf_bin); 
            byteCount += 1;                   
        }
        else{
            k = 0;
        }
        // printf("%d\n", k);
        // printf("%ld %ld\n", byteCount, bitCount);
        fwrite(&k, sizeof(int), 1, hf_bin);
        printf("Overall Compression Ratio(%s/huffman_code.bin): %f\n\n", argv[1], header.filesize/(float)byteCount);

        fclose(hf_bin);         
    } 
    freeHuffmanTree(huffmanTree);
    freeCodeTable(codeTable);                

    free(pix_rgb);

    // Free allocated memory for f_Y, f_Cb, f_Cr
    free_f(height+rM, f_Y);
    free_f(height+rM, f_Cb);
    free_f(height+rM, f_Cr);

    // Free allocated memory for F_Y, F_Cb, F_Cr, qF_Y, qF_Cb, qF_Cr
    free_F(M, N, F_Y);
    free_F(M, N, F_Cb);
    free_F(M, N, F_Cr);
    free_qF(M, N, qF_Y);
    free_qF(M, N, qF_Cb);
    free_qF(M, N, qF_Cr);

    // Free allocated memory for zz_Y, zz_Cb, zz_Cr
    free_zz(M, zz_Y);
    free_zz(M, zz_Cb);
    free_zz(M, zz_Cr);

    // Free allocated memory for rle_Y, rle_Cb, rle_Cr
    free(rle_code_Y); 
    free(rle_code_Cb); 
    free(rle_code_Cr); 

    fclose(bmp);        

    return 0;
}
