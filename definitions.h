#ifndef DEFINITIONS_H
#define DEFINITIONS_H

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

#endif // DEFINITIONS_H