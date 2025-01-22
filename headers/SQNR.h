#ifndef SQNR_H
#define SQNR_H

#include <stdio.h>
#include "definitions.h"

void calculate_Ps_Pe(int height, int width, float Ps[H][W], float Pe[H][W], float ****F, char channel);
void write_SQNR(FILE *fp, float sqnr[H][W], float Ps[H][W], float Pe[H][W]);

#endif // SQNR_H