#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "huffmanEncoding.h"

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
