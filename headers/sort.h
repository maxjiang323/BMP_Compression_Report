#ifndef SORT_H
#define SORT_H

// encoder.c
void dpcm(int height, int width, short *****qF);
void apply_zigzag(int height, int width, short ****dpcm, short ***zz);

// decoder.c
void un_zigzag(int height, int width, short ***zz, short *****dpcm);
void un_dpcm(int height, int width, short *****dpcm);

#endif // SORT_H
