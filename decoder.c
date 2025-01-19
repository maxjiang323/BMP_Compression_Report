#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#define W 8 // dimension of basis vector (width)
#define H 8 // dimension of basis vector (height)
#define Pi 3.14159265359

typedef struct {
    char identifier[2];           // BMP identifier, typically "BM"
    unsigned int filesize;        // Total size of the BMP file in bytes
    unsigned short reserved;      // Reserved, must be set to 0
    unsigned short reserved2;     // Reserved, must be set to 0
    unsigned int bitmap_dataoffset; // Offset to the start of pixel data
    unsigned int bitmap_headersize; // Size of the DIB header (usually 40 bytes)
    unsigned int width;           // Width of the image in pixels
    unsigned int height;          // Height of the image in pixels
    unsigned short planes;        // Number of color planes, must be 1
    unsigned short bits_perpixel; // Bits per pixel (e.g., 1, 4, 8, 24)
    unsigned int compression;     // Compression type (0 = none, 1 = RLE-8, etc.)
    unsigned int bitmap_datasize; // Size of the bitmap data in bytes (may be 0 if uncompressed)
    unsigned int hresolution;     // Horizontal resolution in pixels per meter
    unsigned int vresolution;     // Vertical resolution in pixels per meter
    unsigned int usedcolors;      // Number of colors used in the color palette (0 = all colors)
    unsigned int importantcolors; // Number of important colors (0 = all colors are important)
} Bmpheader;


void read_dim(FILE *dim_txt, Bmpheader *header){
    fscanf(dim_txt, 
        "identifier: %c%c, filesize: %u, reserved: %hu, reserved2: %hu, "
        "bitmap_dataoffset: %u, bitmap_headersize: %u, width: %u, height: %u, "
        "planes: %hu, bits_perpixel: %hu, compression: %u, bitmap_datasize: %u, "
        "hresolution: %u, vresolution: %u, usedcolors: %u, importantcolors: %u",
        &header->identifier[0], &header->identifier[1],
        &header->filesize,
        &header->reserved, &header->reserved2,
        &header->bitmap_dataoffset,
        &header->bitmap_headersize,
        &header->width, &header->height,
        &header->planes, &header->bits_perpixel,
        &header->compression, &header->bitmap_datasize,
        &header->hresolution, &header->vresolution,
        &header->usedcolors, &header->importantcolors); 
}


void readheader(FILE* fp, Bmpheader *header) {
	fread(&header->identifier, sizeof(header->identifier), 1, fp);
	fread(&header->filesize, sizeof(header->filesize), 1, fp);
	fread(&header->reserved, sizeof(header->reserved), 1, fp);
	fread(&header->reserved2, sizeof(header->reserved2), 1, fp);
	fread(&header->bitmap_dataoffset, sizeof(header->bitmap_dataoffset), 1, fp);
	fread(&header->bitmap_headersize, sizeof(header->bitmap_headersize), 1, fp);
	fread(&header->width, sizeof(header->width), 1, fp);
	fread(&header->height, sizeof(header->height), 1, fp);
	fread(&header->planes, sizeof(header->planes), 1, fp);
	fread(&header->bits_perpixel, sizeof(header->bits_perpixel), 1, fp);
	fread(&header->compression, sizeof(header->compression), 1, fp);
	fread(&header->bitmap_datasize, sizeof(header->bitmap_datasize), 1, fp);
	fread(&header->hresolution, sizeof(header->hresolution), 1, fp);
	fread(&header->vresolution, sizeof(header->vresolution), 1, fp);
	fread(&header->usedcolors, sizeof(header->usedcolors), 1, fp);
	fread(&header->importantcolors, sizeof(header->importantcolors), 1, fp);
}

void write_header(FILE* outf, Bmpheader header1){
	fwrite(&header1.identifier, sizeof(char), 2, outf);
	fwrite(&header1.filesize, sizeof(unsigned int), 1, outf);
	fwrite(&header1.reserved, sizeof(unsigned short), 1, outf);
	fwrite(&header1.reserved2, sizeof(unsigned short), 1, outf);
	fwrite(&header1.bitmap_dataoffset, sizeof(unsigned int), 1, outf);
	fwrite(&header1.bitmap_headersize, sizeof(unsigned int), 1, outf);
	fwrite(&header1.width, sizeof(unsigned int), 1, outf);
	fwrite(&header1.height, sizeof(unsigned int), 1, outf);
	fwrite(&header1.planes, sizeof(unsigned short), 1, outf);
	fwrite(&header1.bits_perpixel, sizeof(unsigned short), 1, outf);
	fwrite(&header1.compression, sizeof(unsigned int), 1, outf);
	fwrite(&header1.bitmap_datasize, sizeof(unsigned int), 1, outf);
	fwrite(&header1.hresolution, sizeof(unsigned int), 1, outf);
	fwrite(&header1.vresolution, sizeof(unsigned int), 1, outf);
	fwrite(&header1.usedcolors, sizeof(unsigned int), 1, outf);
	fwrite(&header1.importantcolors, sizeof(unsigned int), 1, outf);
}


#define MAX_RANGE 512 // Total number of symbols including "0 0"

// Huffman tree node structure
typedef struct Node {
    int symbol; // Symbol represented by this node (-255 ~ 255 or 511 for "0 0")
    struct Node *left; // Left child
    struct Node *right; // Right child
} Node;

// Function to create a new node
Node *createNode(int symbol) {
    Node *newNode = (Node *)malloc(sizeof(Node));
    newNode->symbol = symbol;
    newNode->left = NULL;
    newNode->right = NULL;
    return newNode;
}

// Build Huffman tree from codebook
Node *buildHuffmanTreeFromCodebook(FILE *codebook) {
    char line[1024];
    Node *root = createNode(-1); // Root node

    while (fgets(line, sizeof(line), codebook)) {
        int symbol;
        char code[512];

        // Parse the codebook line (e.g., "Symbol: 74, Frequency: 10, Code: 010100")
        if (strstr(line, "Symbol: \"0 0\"")) { // strstr: find the string 
            sscanf(line, "Symbol: \"0 0\", Frequency: %*d, Code: %s", code);
            // line is char* not a file pointer, so can't replace sscanf with fscanf
            symbol = 511; // "0 0" symbol
        } else {
            sscanf(line, "Symbol: %d, Frequency: %*d, Code: %s", &symbol, code);
            symbol += 255; // Adjust symbol to match index
        }

        // Insert the code into the Huffman tree
        Node *current = root;
        int i;
        for (i = 0; code[i] != '\0'; i++) {
            if (code[i] == '0') {
                if (!current->left) current->left = createNode(-1);
                current = current->left;
            } else if (code[i] == '1') {
                if (!current->right) current->right = createNode(-1);
                current = current->right;
            }
        }
        current->symbol = symbol; // Set symbol at the leaf node
    }

    return root;
}

// Decode Huffman encoded bitstream
void decodeHuffmanTxt(Node *root, FILE *encodedFile, short **Y, short **Cb, short **Cr, long *pos_Y, long *pos_Cb, long *pos_Cr) {
    Node *current = root;
    char bit;
    int channel = 0; // 0: Y, 1: Cb, 2: Cr
    // printf("1\n");
    while ((bit = fgetc(encodedFile)) != EOF) {
        // printf("2\n");
        if (bit == '0') {
            current = current->left;
        } else if (bit == '1') {
            current = current->right;
        } else {
            continue; // Ignore invalid characters
        }

        // Check if we reached a leaf node
        if (current->left == NULL && current->right == NULL) {
            if (current->symbol == 511) { // "0 0" symbol
                if (channel == 0) {
                    (*Y)[(*pos_Y)] = 0;
                    (*Y)[(*pos_Y)+1] = 0;
                    (*pos_Y) += 2;
                } 
                else if (channel == 1) {
                    (*Cb)[(*pos_Cb)] = 0;
                    (*Cb)[(*pos_Cb)+1] = 0;
                    (*pos_Cb) += 2;
                } 
                else {
                    (*Cr)[(*pos_Cr)] = 0;
                    (*Cr)[(*pos_Cr)+1] = 0;
                    (*pos_Cr) += 2;
                }
                channel = (channel + 1) % 3; // Switch channel
            } 
            else {
                if (channel == 0) {
                    (*Y)[(*pos_Y)] = current->symbol -255;
                    // printf("%d ", (*Y)[(*pos_Y)]);
                    (*pos_Y) += 1;
                } else if (channel == 1) {
                    (*Cb)[(*pos_Cb)] = current->symbol -255;
                    (*pos_Cb) += 1;
                } else {
                    (*Cr)[(*pos_Cr)] = current->symbol -255;
                    (*pos_Cr) += 1;
                }
            }
            current = root; // Reset to root for next symbol
        }
    }
}

// Decode Huffman encoded bitstream from binary file
void decodeHuffmanBin(Node *root, FILE *binFile, short **Y, short **Cb, short **Cr, long *pos_Y, long *pos_Cb, long *pos_Cr) {
    Node *current = root;
    unsigned char byte;

    if (root == NULL) {
        fprintf(stderr, "Error: Huffman tree root is NULL.\n");
        exit(EXIT_FAILURE);
    }

    long pos0 = ftell(binFile); // the start position of the code data
    fseek(binFile, -sizeof(int), SEEK_END); // jump to count of padding bits
    long pos1 = ftell(binFile); // the end position of the code data
    // printf("%ld %ld\n", pos0, pos1);

    long byteCount = pos1 - pos0; // the count of data bytes 
    int k; // the remaining bits of data
    fread(&k, sizeof(int), 1, binFile);
    // printf("%ld\n", ftell(binFile));
    // printf("%d\n", k);
    // printf("%ld %ld\n", byteCount, byteCount*8); // bytes / bits

    fseek(binFile, pos0, SEEK_SET); // back to the start of the code data

    long i; 
    int bitIndex, channel = 0; // 0: Y, 1: Cb, 2: Cr
    char bit;
    for (i = 0; i < byteCount; i++) {
        fread(&byte, 1, 1, binFile);
        int x;
        if (i == (byteCount-1)){ // reach the last byte in bin file
            x = 8 - k;
        }
        else{
            x = 8;
        }
        for (bitIndex = 0; bitIndex < x; bitIndex++) {
            unsigned char bit = (byte >> (7 - bitIndex)) & 1; // Extract bit
            if (current == NULL || (bit && current->right == NULL) || (!bit && current->left == NULL)) {
                fprintf(stderr, "Invalid Huffman tree traversal.\n");
                exit(EXIT_FAILURE);
            }

            if (bit == 0) {
                current = current->left;
            } else if (bit == 1) {
                current = current->right;
            } else {
                continue; // Ignore invalid characters
            }

            // Check if we reached a leaf node
            if (current->left == NULL && current->right == NULL) {
                if (current->symbol == 511) { // "0 0" symbol
                    if (channel == 0) {
                        (*Y)[(*pos_Y)] = 0;
                        (*Y)[(*pos_Y)+1] = 0;
                        (*pos_Y) += 2;
                    } 
                    else if (channel == 1) {
                        (*Cb)[(*pos_Cb)] = 0;
                        (*Cb)[(*pos_Cb)+1] = 0;
                        (*pos_Cb) += 2;
                    } 
                    else {
                        (*Cr)[(*pos_Cr)] = 0;
                        (*Cr)[(*pos_Cr)+1] = 0;
                        (*pos_Cr) += 2;
                    }
                    channel = (channel + 1) % 3; // Switch channel
                } 
                else {
                    if (channel == 0) {
                        (*Y)[(*pos_Y)] = current->symbol -255;
                        // printf("%d ", (*Y)[(*pos_Y)]);
                        (*pos_Y) += 1;
                    } else if (channel == 1) {
                        (*Cb)[(*pos_Cb)] = current->symbol -255;
                        (*pos_Cb) += 1;
                    } else {
                        (*Cr)[(*pos_Cr)] = current->symbol -255;
                        (*pos_Cr) += 1;
                    }
                }
                current = root; // Reset to root for next symbol
            }
        }
    }
}


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


typedef struct {
    unsigned char R;
    unsigned char G;
    unsigned char B;
} RGB;

RGB YCbCr2RGB(float Y, float Cb, float Cr) { // Rec. 709
    RGB result;
    float R = Y + 1.5748*Cr;
    float G = Y - 0.1873*Cb - 0.4681*Cr;
    float B = Y + 1.8556*Cb;

    // Clamp values to 0-255 range
    if (R < 0){
        result.R = 0;
    } else if (R > 255){
        result.R = 255;
    } else{
        result.R = (unsigned char)(round(R));
    }

    if (G < 0){
        result.G = 0;
    } else if (G > 255){
        result.G = 255;
    } else{
        result.G = (unsigned char)(round(G));
    }

    if (B < 0){
        result.B = 0;
    } else if (B > 255){
        result.B = 255;
    } else{
        result.B = (unsigned char)(round(B));
    }

    return result;
}

void set_rgb_data(int width, int height, int *row_length, unsigned char **pix_rgb){
    int width_rgb = width * 3;  // rgb -> *3
    int padding = 4 - (width_rgb % 4);
    if ((width_rgb % 4) == 0){
        padding = 0;
    }
    *row_length = width_rgb + padding; 
    // add the padding to make (the modified width, the length of a row) % 4 == 0
    *pix_rgb = (unsigned char *)calloc((*row_length) * (height), sizeof(unsigned char));
}


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
