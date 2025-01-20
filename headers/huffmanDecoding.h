#ifndef HUFFMANDECODING_H
#define HUFFMANDECODING_H

#include "definitions.h"
#define MAX_RANGE 512 // -255 to 255 + one extra symbol for "0 0"

// Define the node structure for the Huffman tree
typedef struct Node {
    int symbol;        // Corresponding symbol (-255 ~ 255 or 511 for "0 0")
    struct Node *left; // Left child node
    struct Node *right; // Right child node
} Node;

Node *createNode(int symbol);
Node *buildHuffmanTreeFromCodebook(FILE *codebook);
void decodeHuffmanTxt(Node *root, FILE *encodedFile, short **Y, short **Cb, short **Cr, long *pos_Y, long *pos_Cb, long *pos_Cr);
void decodeHuffmanBin(Node *root, FILE *binFile, short **Y, short **Cb, short **Cr, long *pos_Y, long *pos_Cb, long *pos_Cr);
void freeHuffmanTree(Node *root);

#endif // HUFFMANDECODING_H
