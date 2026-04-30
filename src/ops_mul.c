#include "ops_mul.h"
#include <unistd.h>
#include <sys/wait.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <signal.h>
#include <stdio.h> 
#include <sys/types.h>
#ifdef _OPENMP
#include <omp.h>
#endif



int kill(pid_t, int);


static volatile sig_atomic_t g_sig_done_mul = 0;
static void on_sigusr1_mul(int signo) {
    (void)signo;
    g_sig_done_mul++;
    
    signal(SIGUSR1, on_sigusr1_mul);
}


typedef struct {
    int i, j;   
} mul_meta_t;

typedef struct {
    int i, j;
    double v;
    int status; 
} mul_res_t;

/* pipes */
static ssize_t write_full(int fd, const void *buf, size_t n) {
    const char *p = (const char*)buf;
    size_t left = n;
    while (left > 0) {
        ssize_t w = write(fd, p, left);
        if (w < 0) {
            if (errno == EINTR) continue;
            return -1;
        }
        p += (size_t)w;
        left -= (size_t)w;
    }
    return (ssize_t)n;
}
static ssize_t read_full(int fd, void *buf, size_t n) {
    char *p = (char*)buf;
    size_t need = n;
    while (need > 0) {
        ssize_t r = read(fd, p, need);
        if (r == 0) return -1;               /* EOF */
        if (r < 0) {
            if (errno == EINTR) continue;
            return -1;
        }
        p += (size_t)r;
        need -= (size_t)r;
    }
    return (ssize_t)n;
}

int mul_matrices_distributed(const Matrix *A, const Matrix *B, Matrix *OUT) {
    if (!A || !B || !A->data || !B->data) return -1;
    if (A->cols != B->rows) return -2;

    if (mat_create(OUT, "PROD", A->rows, B->cols) != 0) return -1;

    const int R = A->rows;
    const int C = B->cols;
    const int K = A->cols; 

    
    signal(SIGPIPE, SIG_IGN);

    
    g_sig_done_mul = 0;
    signal(SIGUSR1, on_sigusr1_mul);

   
    for (int i = 0; i < R; ++i) {
        for (int j = 0; j < C; ++j) {

           
            double *row = (double*)malloc((size_t)K * sizeof(double));
            double *col = (double*)malloc((size_t)K * sizeof(double));
            if (!row || !col) {
                free(row); free(col);
                mat_destroy(OUT);
                return -1;
            }

            for (int k = 0; k < K; ++k) {
                row[k] = A->data[(size_t)i * A->cols + k];
                col[k] = B->data[(size_t)k * B->cols + j];
            }

            int p2c[2], c2p[2];
            if (pipe(p2c) != 0 || pipe(c2p) != 0) {
                free(row); free(col);
                mat_destroy(OUT);
                return -1;
            }

            pid_t pid = fork();
            if (pid < 0) {
                close(p2c[0]); close(p2c[1]);
                close(c2p[0]); close(c2p[1]);
                free(row); free(col);
                mat_destroy(OUT);
                return -1;
            }

            if (pid == 0) {
               
                close(p2c[1]);  
                close(c2p[0]);  

                mul_meta_t meta;
                int klen = 0;

                if (read_full(p2c[0], &meta, sizeof(meta)) != (ssize_t)sizeof(meta) ||
                    read_full(p2c[0], &klen, sizeof(klen)) != (ssize_t)sizeof(klen) ||
                    klen <= 0) {
                    close(p2c[0]); close(c2p[1]); _exit(1);
                }

                double *row_in = (double*)malloc((size_t)klen * sizeof(double));
                double *col_in = (double*)malloc((size_t)klen * sizeof(double));
                if (!row_in || !col_in) {
                    free(row_in); free(col_in);
                    close(p2c[0]); close(c2p[1]); _exit(1);
                }

                if (read_full(p2c[0], row_in, (size_t)klen * sizeof(double)) != (ssize_t)((size_t)klen * sizeof(double)) ||
                    read_full(p2c[0], col_in, (size_t)klen * sizeof(double)) != (ssize_t)((size_t)klen * sizeof(double))) {
                    free(row_in); free(col_in);
                    close(p2c[0]); close(c2p[1]); _exit(1);
                }
                close(p2c[0]);

                double sum = 0.0;
                for (int t = 0; t < klen; ++t)
                    sum += row_in[t] * col_in[t];

                mul_res_t res;
                res.i = meta.i;
                res.j = meta.j;
                res.v = sum;
                res.status = 0;

                (void)write_full(c2p[1], &res, sizeof(res));
                close(c2p[1]);
                free(row_in); free(col_in);

                
                kill(getppid(), SIGUSR1);

                _exit(0);
            }

            
            close(p2c[0]);  
            close(c2p[1]);  

            mul_meta_t meta = { .i = i, .j = j };
            int klen = K;

            if (write_full(p2c[1], &meta, sizeof(meta)) != (ssize_t)sizeof(meta) ||
                write_full(p2c[1], &klen, sizeof(klen)) != (ssize_t)sizeof(klen) ||
                write_full(p2c[1], row, (size_t)klen * sizeof(double)) != (ssize_t)((size_t)klen * sizeof(double)) ||
                write_full(p2c[1], col, (size_t)klen * sizeof(double)) != (ssize_t)((size_t)klen * sizeof(double))) {
                close(p2c[1]);
                close(c2p[0]);
                int st; (void)waitpid(pid, &st, 0);
                free(row); free(col);
                mat_destroy(OUT);
                return -1;
            }
            close(p2c[1]);

            mul_res_t rmsg;
            if (read_full(c2p[0], &rmsg, sizeof(rmsg)) != (ssize_t)sizeof(rmsg) || rmsg.status != 0) {
                close(c2p[0]);
                int st; (void)waitpid(pid, &st, 0);
                free(row); free(col);
                mat_destroy(OUT);
                return -1;
            }
            close(c2p[0]);

            *at(OUT, rmsg.i, rmsg.j) = rmsg.v;

            int st; (void)waitpid(pid, &st, 0);
            free(row); free(col);
        }
    }

    
       printf("[SIG] mul children finished: %d/%d\n",
              (int)g_sig_done_mul, R*C);
    

    return 0;
}

int mul_matrices_serial(const Matrix *A, const Matrix *B, Matrix *OUT) {
    if (!A || !B || !A->data || !B->data) return -1;
    if (A->cols != B->rows) return -2;

  
    if (mat_create(OUT, "PROD_SER", A->rows, B->cols) != 0) return -1;

    int R = A->rows;
    int C = B->cols;
    int K = A->cols;  

    for (int i = 0; i < R; ++i) {
        for (int j = 0; j < C; ++j) {
            double sum = 0.0;
            for (int k = 0; k < K; ++k) {
                double a = A->data[(size_t)i * A->cols + k];
                double b = B->data[(size_t)k * B->cols + j];
                sum += a * b;
            }
            *at(OUT, i, j) = sum;
        }
    }
    return 0;
}

int mul_matrices_omp(const Matrix *A, const Matrix *B, Matrix *OUT) {
    if (!A || !B || !A->data || !B->data) return -1;
    if (A->cols != B->rows) return -2;

    if (mat_create(OUT, "PROD_OMP", A->rows, B->cols) != 0) return -1;

    int R = A->rows;
    int C = B->cols;
    int K = A->cols;

   
    #pragma omp parallel for if ((long)R * C > 64)
    for (int i = 0; i < R; ++i) {
        for (int j = 0; j < C; ++j) {
            double sum = 0.0;
            for (int k = 0; k < K; ++k) {
                sum += A->data[(size_t)i * A->cols + k] *
                       B->data[(size_t)k * B->cols + j];
            }
            OUT->data[(size_t)i * C + j] = sum;
        }
    }

    return 0;
}



