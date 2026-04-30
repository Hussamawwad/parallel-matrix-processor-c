#include "ops_det.h"
#include <unistd.h>
#include <sys/wait.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <signal.h>
#include <stdio.h> 
#include <stdio.h>
#include <sys/types.h>


int kill(pid_t, int);


static volatile sig_atomic_t g_sig_done_det = 0;
static void on_sigusr1_det(int signo) {
    (void)signo;
    g_sig_done_det++;
    
    signal(SIGUSR1, on_sigusr1_det);
}

/*  pipes  */
static ssize_t write_full(int fd, const void *buf, size_t n) {
    const char *p = (const char*)buf;
    size_t left = n;
    while (left > 0) {
        ssize_t w = write(fd, p, left);
        if (w < 0) { if (errno == EINTR) continue; return -1; }
        p += (size_t)w; left -= (size_t)w;
    }
    return (ssize_t)n;
}
static ssize_t read_full(int fd, void *buf, size_t n) {
    char *p = (char*)buf;
    size_t need = n;
    while (need > 0) {
        ssize_t r = read(fd, p, need);
        if (r == 0) return -1; /* EOF */
        if (r < 0) { if (errno == EINTR) continue; return -1; }
        p += (size_t)r; need -= (size_t)r;
    }
    return (ssize_t)n;
}

/* ===== Helpers ===== */
static inline double M_at(const Matrix *A, int r, int c) {
    return A->data[(size_t)r * A->cols + c];
}


static double *build_minor_data(const Matrix *A, int row_skip, int col_skip, int *out_n) {
    int n = A->rows;
    int m = n - 1;
    *out_n = m;
    if (m <= 0) return NULL;
    double *buf = (double*)malloc((size_t)m * (size_t)m * sizeof(double));
    if (!buf) return NULL;

    int rr = 0;
    for (int r = 0; r < n; ++r) {
        if (r == row_skip) continue;
        int cc = 0;
        for (int c = 0; c < n; ++c) {
            if (c == col_skip) continue;
            buf[(size_t)rr * m + cc] = M_at(A, r, c);
            cc++;
        }
        rr++;
    }
    return buf;
}


static double det_serial_buf(const double *D, int n) {
    if (n == 1) return D[0];
    if (n == 2) return D[0]*D[3] - D[1]*D[2];

    double det = 0.0;
    for (int j = 0; j < n; ++j) {
        int m = n - 1;
        double *minor = (double*)malloc((size_t)m*(size_t)m*sizeof(double));
        if (!minor) return det; 
        int rr = 0;
        for (int r = 1; r < n; ++r) {
            int cc = 0;
            for (int c = 0; c < n; ++c) {
                if (c == j) continue;
                minor[(size_t)rr*m + cc] = D[(size_t)r*n + c];
                cc++;
            }
            rr++;
        }
        double sign = (j % 2) ? -1.0 : 1.0;
        double cof = sign * D[j] * det_serial_buf(minor, m);
        det += cof;
        free(minor);
    }
    return det;
}

/* Gaussian elimination with partial pivoting (returns det). */
static int det_gauss(double *D, int n, double *out_det) {
    double det = 1.0;
    for (int k = 0; k < n; ++k) {
        
        int piv = k;
        double v0 = D[(size_t)k*n + k];
        double maxabs = v0 >= 0 ? v0 : -v0;
        for (int r = k+1; r < n; ++r) {
            double v = D[(size_t)r*n + k]; if (v < 0) v = -v;
            if (v > maxabs) { maxabs = v; piv = r; }
        }
        if (maxabs == 0.0) { *out_det = 0.0; return 0; } /* singular */

        if (piv != k) {
            for (int c = 0; c < n; ++c) {
                double tmp = D[(size_t)k*n + c];
                D[(size_t)k*n + c] = D[(size_t)piv*n + c];
                D[(size_t)piv*n + c] = tmp;
            }
            det = -det; /* row swap flips sign */
        }

        double pivot = D[(size_t)k*n + k];
        det *= pivot;
        if (pivot == 0.0) { *out_det = 0.0; return 0; }

        for (int r = k+1; r < n; ++r) {
            double f = D[(size_t)r*n + k] / pivot;
            if (f == 0.0) continue;
            D[(size_t)r*n + k] = 0.0;
            for (int c = k+1; c < n; ++c) {
                D[(size_t)r*n + c] -= f * D[(size_t)k*n + c];
            }
        }
    }
    *out_det = det;
    return 0;
}



static int start_task_slot(const Matrix *A,
                           int (*p2c)[2], int (*c2p)[2], pid_t *pids,
                           int slot, int j) {
    (void)A; 
    if (pipe(p2c[slot]) != 0 || pipe(c2p[slot]) != 0) return -1;

    pid_t pid = fork();
    if (pid < 0) {
        close(p2c[slot][0]); close(p2c[slot][1]);
        close(c2p[slot][0]); close(c2p[slot][1]);
        return -1;
    }

    if (pid == 0) {
        /* ---- child ---- */
        
        close(p2c[slot][1]);  
        close(c2p[slot][0]);  

        int j_idx = -1;
        if (read_full(p2c[slot][0], &j_idx, sizeof(j_idx)) != (ssize_t)sizeof(j_idx)) {
            close(p2c[slot][0]); close(c2p[slot][1]); _exit(1);
        }
        close(p2c[slot][0]);

        int m_n = 0;
       
        double *minor = build_minor_data(A, 0, j_idx, &m_n);
        if (!minor) { 
            double bad = 0.0/0.0;
            (void)write_full(c2p[slot][1], &bad, sizeof(bad));
            close(c2p[slot][1]);
            _exit(1);
        }

        double d = det_serial_buf(minor, m_n);
        free(minor);

        double sign = (j_idx % 2) ? -1.0 : 1.0;
        double cof = sign * M_at(A,0,j_idx) * d;

        (void)write_full(c2p[slot][1], &cof, sizeof(cof));
        close(c2p[slot][1]);

        
        kill(getppid(), SIGUSR1);

        _exit(0);
    }

    /* ---- parent ---- */
    pids[slot] = pid;
    close(p2c[slot][0]);  
    close(c2p[slot][1]);  
    if (write_full(p2c[slot][1], &j, sizeof(j)) != (ssize_t)sizeof(j)) {
        close(p2c[slot][1]); close(c2p[slot][0]);
        int st; (void)waitpid(pid, &st, 0);
        return -1;
    }
    close(p2c[slot][1]); 
    return 0;
}


static int finish_task_slot(int (*c2p)[2], pid_t *pids, int slot, double *accum) {
    double contrib = 0.0;
    if (read_full(c2p[slot][0], &contrib, sizeof(contrib)) != (ssize_t)sizeof(contrib)) {
        close(c2p[slot][0]);
        int st; (void)waitpid(pids[slot], &st, 0);
        return -1;
    }
    close(c2p[slot][0]);
    {
        int st; (void)waitpid(pids[slot], &st, 0);
    }
    
    if (contrib != contrib) return -1;
    *accum += contrib;
    return 0;
}

/* ===== Main API ===== */
int determinant_distributed(const Matrix *A, double *out_det) {
    if (!A || !A->data || !out_det) return -1;
    if (A->rows != A->cols) return -2;

    int n = A->rows;

   
    if (n == 0) { *out_det = 0.0; return 0; }
    if (n == 1) { *out_det = M_at(A,0,0); return 0; }
    if (n == 2) {
        *out_det = M_at(A,0,0)*M_at(A,1,1) - M_at(A,0,1)*M_at(A,1,0);
        return 0;
    }

   
    long cores = sysconf(_SC_NPROCESSORS_ONLN);
    if (cores < 1) cores = 1;
    int POOL = (int)(cores * 2);
    if (POOL > 32) POOL = 32;
    if (POOL < 1)  POOL = 1;

    signal(SIGPIPE, SIG_IGN);

    /*  Gaussian  */
    if (n > 6) {
        double *buf = (double*)malloc((size_t)n*(size_t)n*sizeof(double));
        if (!buf) return -1;
        for (int r = 0; r < n; ++r)
            for (int c = 0; c < n; ++c)
                buf[(size_t)r*n + c] = M_at(A, r, c);

        int rc = det_gauss(buf, n, out_det);
        free(buf);
        return rc == 0 ? 0 : -1;
    }

    /* n <= 6 : Laplace , Process Pool  */
    double det_sum = 0.0;
    int next = 0;   
    int live = 0;   

   
    g_sig_done_det = 0;
    signal(SIGUSR1, on_sigusr1_det);

   
    int (*p2c)[2] = (int (*)[2])malloc((size_t)POOL * sizeof(int[2]));
    int (*c2p)[2] = (int (*)[2])malloc((size_t)POOL * sizeof(int[2]));
    pid_t *pids    = (pid_t*)malloc((size_t)POOL * sizeof(pid_t));
    if (!p2c || !c2p || !pids) {
        free(p2c); free(c2p); free(pids);
        return -1;
    }

    while (next < n || live > 0) {
     
        while (live < POOL && next < n) {
            int slot = live;
            if (start_task_slot(A, p2c, c2p, pids, slot, next) != 0) {
            
                for (int s = 0; s < live; ++s) {
                    close(c2p[s][0]);
                    int st; (void)waitpid(pids[s], &st, 0);
                }
                free(p2c); free(c2p); free(pids);
                return -1;
            }
            live++;
            next++;
        }

      
        if (live > 0) {
            if (finish_task_slot(c2p, pids, 0, &det_sum) != 0) {
                for (int s = 1; s < live; ++s) {
                    close(c2p[s][0]);
                    int st; (void)waitpid(pids[s], &st, 0);
                }
                free(p2c); free(c2p); free(pids);
                return -1;
            }
        
            for (int s = 0; s + 1 < live; ++s) {
                c2p[s][0] = c2p[s+1][0]; c2p[s][1] = c2p[s+1][1];
                pids[s]   = pids[s+1];
                p2c[s][0] = p2c[s+1][0]; p2c[s][1] = p2c[s+1][1];
            }
            live--;
        }
    }

    free(p2c); free(c2p); free(pids);
    *out_det = det_sum;

     
       printf("[SIG] det children finished: %d/%d\n",
              (int)g_sig_done_det, n);
    

    return 0;
}


int determinant_serial(const Matrix *A, double *out_det) {
    if (!A || !A->data || !out_det) return -1;
    if (A->rows != A->cols) return -2;

    int n = A->rows;
    if (n == 0) { *out_det = 0.0; return 0; }
    if (n == 1) { *out_det = A->data[0]; return 0; }
    if (n == 2) {
        *out_det = A->data[0] * A->data[3] - A->data[1] * A->data[2];
        return 0;
    }

    /*  Gaussian elimination */
    double *buf = (double*)malloc((size_t)n * (size_t)n * sizeof(double));
    if (!buf) return -1;

    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c)
            buf[(size_t)r * n + c] = A->data[(size_t)r * A->cols + c];

    double d = 0.0;
    int rc = det_gauss(buf, n, &d);
    free(buf);

    if (rc != 0) return -1;
    *out_det = d;
    return 0;
}
/* Gaussian elimination with partial pivoting (OpenMP) */
static int det_gauss_omp(double *D, int n, double *out_det) {
    double det = 1.0;
    for (int k = 0; k < n; ++k) {
       
        int piv = k;
        double v0 = D[(size_t)k*n + k];
        double maxabs = v0 >= 0 ? v0 : -v0;
        for (int r = k+1; r < n; ++r) {
            double v = D[(size_t)r*n + k]; if (v < 0) v = -v;
            if (v > maxabs) { maxabs = v; piv = r; }
        }
        if (maxabs == 0.0) { *out_det = 0.0; return 0; } /* singular */

        if (piv != k) {
            for (int c = 0; c < n; ++c) {
                double tmp = D[(size_t)k*n + c];
                D[(size_t)k*n + c] = D[(size_t)piv*n + c];
                D[(size_t)piv*n + c] = tmp;
            }
            det = -det; /* row swap flips sign */
        }

        double pivot = D[(size_t)k*n + k];
        det *= pivot;
        if (pivot == 0.0) { *out_det = 0.0; return 0; }

       
        #pragma omp parallel for
        for (int r = k+1; r < n; ++r) {
            double f = D[(size_t)r*n + k] / pivot;
            if (f == 0.0) continue;
            D[(size_t)r*n + k] = 0.0;
            for (int c = k+1; c < n; ++c) {
                D[(size_t)r*n + c] -= f * D[(size_t)k*n + c];
            }
        }
    }
    *out_det = det;
    return 0;
}

int determinant_omp(const Matrix *A, double *out_det) {
    if (!A || !A->data || !out_det) return -1;
    if (A->rows != A->cols) return -2;
    int n = A->rows;
    if (n == 0) { *out_det = 0.0; return 0; }
    if (n == 1) { *out_det = A->data[0]; return 0; }

    double *buf = (double*)malloc((size_t)n*(size_t)n*sizeof(double));
    if (!buf) return -1;

   
    #pragma omp parallel for
    for (int r=0; r<n; ++r) {
        for (int c=0; c<n; ++c) {
            buf[(size_t)r*n + c] = A->data[(size_t)r*A->cols + c];
        }
    }

    int rc = det_gauss_omp(buf, n, out_det);
    free(buf);
    return rc == 0 ? 0 : -1;
}


