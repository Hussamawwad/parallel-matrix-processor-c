#ifndef OPS_MUL_H
#define OPS_MUL_H

#include "matrix.h"


int mul_matrices_distributed(const Matrix *A, const Matrix *B, Matrix *OUT);

int mul_matrices_serial(const Matrix *A, const Matrix *B, Matrix *OUT);

int mul_matrices_omp(const Matrix *A, const Matrix *B, Matrix *OUT);
#endif

