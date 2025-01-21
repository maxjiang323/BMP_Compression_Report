/*
Portions of this code are from ntpu-ce-mmsp-2023
Source: https://github.com/cychiang-ntpu/ntpu-ce-mmsp-2023
Copyright 2004 Chen-Yu Chiang
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Modifications made on 2025/01/19 by Ming Ju Chiang:
- Changed the structure name from `Bmpheader` to `Bmpheader` (for consistency).
- Updated field comments for better clarity and understanding.
- Removed the `palette` field as it is unused in this implementation.
*/


#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#define W 8 // dimension of basis vector (width)
#define H 8 // dimension of basis vector (height)
#define Pi 3.14159265359

// Modified version of the original `Bmpheader` structure from ntpu-ce-mmsp-2023
// header Info of BMP file
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


// pixeldata.h
typedef struct {
    unsigned char R;
    unsigned char G;
    unsigned char B;
} RGB;

typedef struct {
    float Y;
    float Cb;
    float Cr;
} YCbCr;

#endif // DEFINITIONS_H