#include "ops_addsub.h"
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


static volatile sig_atomic_t g_sig_done_addsub = 0;
static void on_sigusr1_addsub(int signo) {
    (void)signo;
    g_sig_done_addsub++;
   
    signal(SIGUSR1, on_sigusr1_addsub);
}


typedef struct {
    int i, j;
    double a, b;
} elem_msg_t;

typedef struct {
    int i, j;
    double v;
    int status; 
} res_msg_t;

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
        if (r == 0) return -1;              /* EOF */
        if (r < 0) {
            if (errno == EINTR) continue;
            return -1;
        }
        p += (size_t)r;
        need -= (size_t)r;
    }
    return (ssize_t)n;
}

static int do_op_distributed(const Matrix *A, const Matrix *B, Matrix *OUT, int is_add) {
    if (!A || !B || !A->data || !B->data) return -1;
    if (A->rows != B->rows || A->cols != B->cols) return -2;

    if (mat_create(OUT, is_add ? "SUM" : "DIFF", A->rows, A->cols) != 0) return -1;

    const int R = A->rows, C = A->cols;
    const int nchild = R * C;
    if (nchild <= 0) return 0;

  
    signal(SIGPIPE, SIG_IGN);

   
    g_sig_done_addsub = 0;
    signal(SIGUSR1, on_sigusr1_addsub);

    int k = 0;
    for (int i = 0; i < R; ++i) {
        for (int j = 0; j < C; ++j, ++k) {

            int p2c[2], c2p[2];
            if (pipe(p2c) != 0 || pipe(c2p) != 0) {
                mat_destroy(OUT);
                return -1;
            }

            pid_t pid = fork();
            if (pid < 0) {
                close(p2c[0]); close(p2c[1]);
                close(c2p[0]); close(c2p[1]);
                mat_destroy(OUT);
                return -1;
            }

            if (pid == 0) {
               
                close(p2c[1]);  /* read only */
                close(c2p[0]);  /* write only */

                elem_msg_t in;
                if (read_full(p2c[0], &in, sizeof(in)) != (ssize_t)sizeof(in)) {
                    close(p2c[0]); close(c2p[1]); _exit(1);
                }

                res_msg_t out;
                out.i = in.i;
                out.j = in.j;
                out.status = 0;
                out.v = is_add ? (in.a + in.b) : (in.a - in.b);

                (void)write_full(c2p[1], &out, sizeof(out));
                close(p2c[0]);
                close(c2p[1]);

               
                kill(getppid(), SIGUSR1);

                _exit(0);
            }

          
            close(p2c[0]);  
            close(c2p[1]);  

            elem_msg_t m;
            m.i = i;
            m.j = j;
            m.a = A->data[(size_t)i * A->cols + j];
            m.b = B->data[(size_t)i * B->cols + j];

            if (write_full(p2c[1], &m, sizeof(m)) != (ssize_t)sizeof(m)) {
                close(p2c[1]);
                close(c2p[0]);
                int st; (void)waitpid(pid, &st, 0);
                mat_destroy(OUT);
                return -1;
            }
            close(p2c[1]); 

            res_msg_t r;
            if (read_full(c2p[0], &r, sizeof(r)) != (ssize_t)sizeof(r) || r.status != 0) {
                close(c2p[0]);
                int st; (void)waitpid(pid, &st, 0);
                mat_destroy(OUT);
                return -1;
            }
            close(c2p[0]);

           
            *at(OUT, r.i, r.j) = r.v;

          
            int st; (void)waitpid(pid, &st, 0);
        }
    }

     printf("[SIG] add/sub children finished: %d/%d\n",
           (int)g_sig_done_addsub, nchild);

    return 0;
}

int add_matrices_distributed(const Matrix *A, const Matrix *B, Matrix *OUT) {
    return do_op_distributed(A, B, OUT, 1);
}

int sub_matrices_distributed(const Matrix *A, const Matrix *B, Matrix *OUT) {
    return do_op_distributed(A, B, OUT, 0);
}





int add_matrices_serial(const Matrix *A, const Matrix *B, Matrix *OUT) {
    if (!A || !B || !A->data || !B->data) return -1;
    if (A->rows != B->rows || A->cols != B->cols) return -2;

    if (mat_create(OUT, "SUM_SERIAL", A->rows, A->cols) != 0) return -1;

    for (int i = 0; i < A->rows; ++i) {
        for (int j = 0; j < A->cols; ++j) {
            *at(OUT, i, j) =
                A->data[(size_t)i * A->cols + j] +
                B->data[(size_t)i * B->cols + j];
        }
    }
    return 0;
}

int sub_matrices_serial(const Matrix *A, const Matrix *B, Matrix *OUT) {
    if (!A || !B || !A->data || !B->data) return -1;
    if (A->rows != B->rows || A->cols != B->cols) return -2;

    if (mat_create(OUT, "DIFF_SERIAL", A->rows, A->cols) != 0) return -1;

    for (int i = 0; i < A->rows; ++i) {
        for (int j = 0; j < A->cols; ++j) {
            *at(OUT, i, j) =
                A->data[(size_t)i * A->cols + j] -
                B->data[(size_t)i * B->cols + j];
        }
    }
    return 0;
}


int add_matrices_omp(const Matrix *A, const Matrix *B, Matrix *OUT) {
    if (!A || !B || !A->data || !B->data) return -1;
    if (A->rows != B->rows || A->cols != B->cols) return -2;

    if (mat_create(OUT, "SUM_OMP", A->rows, A->cols) != 0) return -1;

    int R = A->rows;
    int C = A->cols;

    /* parallel  */
    #pragma omp parallel for if (R * C > 256)
    for (int i = 0; i < R; ++i) {
        for (int j = 0; j < C; ++j) {
            OUT->data[(size_t)i * C + j] =
                A->data[(size_t)i * C + j] +
                B->data[(size_t)i * C + j];
        }
    }
    return 0;
}

int sub_matrices_omp(const Matrix *A, const Matrix *B, Matrix *OUT) {
    if (!A || !B || !A->data || !B->data) return -1;
    if (A->rows != B->rows || A->cols != B->cols) return -2;

    if (mat_create(OUT, "DIFF_OMP", A->rows, A->cols) != 0) return -1;

    int R = A->rows;
    int C = A->cols;

    #pragma omp parallel for if (R * C > 256)
    for (int i = 0; i < R; ++i) {
        for (int j = 0; j < C; ++j) {
            OUT->data[(size_t)i * C + j] =
                A->data[(size_t)i * C + j] -
                B->data[(size_t)i * C + j];
        }
    }
    return 0;
}


