#ifndef HUFFMANENCODING_H
#define HUFFMANENCODING_H

#include "definitions.h"
#define MAX_RANGE 512 // -255 to 255 + one extra symbol for "0 0"

// Define the node structure for the Huffman tree
typedef struct Node {
    int value;         // Frequency value of the node
    int symbol;        // Corresponding symbol (-255 ~ 255 or 511 for "0 0")
    struct Node *left; // Left child node
    struct Node *right; // Right child node
} Node;

Node *createNode(int value, int symbol, Node *left, Node *right);
void countFrequency(long size, int (*frequencies)[512], short *data);
void findTwoSmallest(Node *nodes[], int size, int *first, int *second);
Node *buildHuffmanTree(int frequencies[], int size);
void generateCodes(Node *tree, char *codeTable[], char *currentCode, int depth);
void freeHuffmanTree(Node *tree);
void freeCodeTable(char *codeTable[]);
void writeCodebook(FILE *fp, char *codeTable[], int frequencies[]);
void write_hf_txt(short **rle_code, FILE *fp, char *codeTable[], long *pos);
void write_hf_bin(short **rle_code, FILE *fp, char *codeTable[], long *pos, long *bitCount, unsigned char *byte, long *byteCount, int *temp);
void bit_alignment(FILE *fp, char *code, unsigned char *byte, int *temp, long *bitCount, long *byteCount);

#endif // HUFFMANENCODING_H
