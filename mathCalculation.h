#ifndef MATHCALCULATION_H
#define MATHCALCULATION_H

#define W 8 // dimension of basis vector (width)
#define H 8 // dimension of basis vector (height)
// encoder.c
void generate_basis_vector(float basis_vector[H][W][H][W]);
void compute_2D_DCT(int height, int width, float **f, float ****F, float basis_vector[H][W][H][W]);
void quantization(int height, int width, short *****qF, float ****F, char channel);

// decoder.c
void un_quantization(int height, int width, short ****qF, float *****F, char channel);
void compute_2D_IDCT(int height, int width, float **f, float ****F, float basis_vector[H][W][H][W]);

#endif // MATHCALCULATION_H
