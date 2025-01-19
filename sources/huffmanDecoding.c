#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../headers/huffmanDecoding.h"

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

// Function to free memory of the Huffman tree
void freeHuffmanTree(Node *root) {
    if (root == NULL) return;

    // Recursively free left and right subtrees
    freeHuffmanTree(root->left);
    freeHuffmanTree(root->right);

    // Free the current node
    free(root);
}
