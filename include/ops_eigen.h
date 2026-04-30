#ifndef OPS_EIGEN_H
#define OPS_EIGEN_H

#include "matrix.h"


int eigen_power_distributed(const Matrix *A, int k, int max_iter, double tol,
                            Matrix *V_out, Matrix *L_out);
                            
                         
int eigen_power_serial(const Matrix *A, int k, int max_iter, double tol,
                       Matrix *V_out, Matrix *L_out);
                       
int eigen_power_omp(const Matrix *A, int k, int max_iter, double tol,
                    Matrix *V_out, Matrix *L_out);

#endif

