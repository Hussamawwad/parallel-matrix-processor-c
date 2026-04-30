#ifndef OPS_ADDSUB_H
#define OPS_ADDSUB_H

#include "matrix.h"


int add_matrices_distributed(const Matrix *A, const Matrix *B, Matrix *OUT);
int sub_matrices_distributed(const Matrix *A, const Matrix *B, Matrix *OUT);

int add_matrices_serial(const Matrix *A, const Matrix *B, Matrix *OUT);
int sub_matrices_serial(const Matrix *A, const Matrix *B, Matrix *OUT);

int add_matrices_omp(const Matrix *A, const Matrix *B, Matrix *OUT);
int sub_matrices_omp(const Matrix *A, const Matrix *B, Matrix *OUT);

#endif

