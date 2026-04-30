#ifndef OPS_DET_H
#define OPS_DET_H

#include "matrix.h"


int determinant_distributed(const Matrix *A, double *out_det);
int determinant_serial(const Matrix *A, double *out_det);
int determinant_omp(const Matrix *A, double *out_det);


#endif

