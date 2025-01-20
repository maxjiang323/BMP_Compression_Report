#ifndef HEADERINFO_H
#define HEADERINFO_H

#include <stdio.h>
#include "definitions.h"

void read_dim(FILE *dim_txt, Bmpheader *header);
void write_dim(FILE *fp, Bmpheader header);
void readheader(FILE* fp, Bmpheader *header);
void write_header(FILE* outf, Bmpheader header1);

#endif // HEADERINFO_H