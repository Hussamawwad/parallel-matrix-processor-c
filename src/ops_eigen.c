#include "ops_eigen.h"
#include <unistd.h>
#include <sys/wait.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <signal.h>
#include <math.h>
#include <stdio.h>

/* ===== I/O helpers ===== */
static ssize_t write_full(int fd, const void *buf, size_t n) {
    const char *p = (const char*)buf; size_t left = n;
    while (left > 0) {
        ssize_t w = write(fd, p, left);
        if (w < 0) { if (errno == EINTR) continue; return -1; }
        p += (size_t)w; left -= (size_t)w;
    }
    return (ssize_t)n;
}
static ssize_t read_full(int fd, void *buf, size_t n) {
    char *p = (char*)buf; size_t need = n;
    while (need > 0) {
        ssize_t r = read(fd, p, need);
        if (r == 0) return -1; /* EOF */
        if (r < 0) { if (errno == EINTR) continue; return -1; }
        p += (size_t)r; need -= (size_t)r;
    }
    return (ssize_t)n;
}

/* ===== Math helpers ===== */
static double vec_norm2(const double *v, int n) {
    double s = 0.0; for (int i=0;i<n;++i) s += v[i]*v[i]; return sqrt(s);
}
static double dot(const double *a, const double *b, int n) {
    double s = 0.0; for (int i=0;i<n;++i) s += a[i]*b[i]; return s;
}
static void vec_copy(double *dst, const double *src, int n) {
    memcpy(dst, src, (size_t)n*sizeof(double));
}
static void vec_scale(double *v, int n, double s) {
    for (int i=0;i<n;++i) v[i]*=s;
}

static int matvec_distributed(const double *W, int n,
                              const double *x, double *y,
                              int POOL) {
    if (n <= 0) return -1;
    if (POOL < 1) POOL = 1;
    if (POOL > n) POOL = n;

    int (*p2c)[2] = (int (*)[2])malloc((size_t)POOL * sizeof(int[2]));
    int (*c2p)[2] = (int (*)[2])malloc((size_t)POOL * sizeof(int[2]));
    pid_t *pids    = (pid_t*)malloc((size_t)POOL * sizeof(pid_t));
    if (!p2c || !c2p || !pids) { free(p2c); free(c2p); free(pids); return -1; }

    signal(SIGPIPE, SIG_IGN);

    int rows_per = (n + POOL - 1) / POOL; 
    int launched = 0;

    for (int t = 0; t < POOL; ++t) {
        int start = t * rows_per;
        int end   = start + rows_per;
        if (start >= n) break;
        if (end > n) end = n;

        if (pipe(p2c[t]) != 0 || pipe(c2p[t]) != 0) {
            launched = t;
            goto fail_launch;
        }
        pid_t pid = fork();
        if (pid < 0) {
            close(p2c[t][0]); close(p2c[t][1]);
            close(c2p[t][0]); close(c2p[t][1]);
            launched = t;
            goto fail_launch;
        }
        if (pid == 0) {
            /* child */
            close(p2c[t][1]); /* read only */
            close(c2p[t][0]); /* write only */

            int n_in=0, s=0, e=0;
            if (read_full(p2c[t][0], &n_in, sizeof(n_in)) != (ssize_t)sizeof(n_in) ||
                read_full(p2c[t][0], &s, sizeof(s)) != (ssize_t)sizeof(s) ||
                read_full(p2c[t][0], &e, sizeof(e)) != (ssize_t)sizeof(e)) {
                close(p2c[t][0]); close(c2p[t][1]); _exit(1);
            }
            double *x_in = (double*)malloc((size_t)n_in*sizeof(double));
            if (!x_in) { close(p2c[t][0]); close(c2p[t][1]); _exit(1); }
            if (read_full(p2c[t][0], x_in, (size_t)n_in*sizeof(double)) != (ssize_t)((size_t)n_in*sizeof(double))) {
                free(x_in); close(p2c[t][0]); close(c2p[t][1]); _exit(1);
            }
            close(p2c[t][0]);

            int len = e - s;
            double *y_part = (double*)malloc((size_t)len*sizeof(double));
            if (!y_part) { free(x_in); close(c2p[t][1]); _exit(1); }

          
            for (int i = s; i < e; ++i) {
                double acc = 0.0;
                const double *row = &W[(size_t)i*(size_t)n];
                for (int j = 0; j < n; ++j) acc += row[j] * x_in[j];
                y_part[i - s] = acc;
            }

            (void)write_full(c2p[t][1], &len, sizeof(len));
            (void)write_full(c2p[t][1], y_part, (size_t)len*sizeof(double));
            close(c2p[t][1]);
            free(x_in); free(y_part);
            _exit(0);
        }
        /* parent */
        pids[t] = pid;
        close(p2c[t][0]); /* write */
        close(c2p[t][1]); /* read */

        
        (void)write_full(p2c[t][1], &n, sizeof(n));
        (void)write_full(p2c[t][1], &start, sizeof(start));
        (void)write_full(p2c[t][1], &end, sizeof(end));
        (void)write_full(p2c[t][1], x, (size_t)n*sizeof(double));
        close(p2c[t][1]);

        launched = t+1;
    }

    /* collect */
    for (int t = 0; t < launched; ++t) {
        int len = 0;
        if (read_full(c2p[t][0], &len, sizeof(len)) != (ssize_t)sizeof(len)) {
            close(c2p[t][0]); int st; (void)waitpid(pids[t], &st, 0);
            goto fail_collect;
        }
        if (len < 0 || len > n) {
            close(c2p[t][0]); int st; (void)waitpid(pids[t], &st, 0);
            goto fail_collect;
        }
        double *tmp = (double*)malloc((size_t)len*sizeof(double));
        if (!tmp) {
            close(c2p[t][0]); int st; (void)waitpid(pids[t], &st, 0);
            goto fail_collect;
        }
        if (read_full(c2p[t][0], tmp, (size_t)len*sizeof(double)) != (ssize_t)((size_t)len*sizeof(double))) {
            free(tmp); close(c2p[t][0]); int st; (void)waitpid(pids[t], &st, 0);
            goto fail_collect;
        }
        close(c2p[t][0]);
        int start = t * rows_per;
        if (start >= n) start = n; 
        for (int i = 0; i < len && start+i < n; ++i) y[start+i] = tmp[i];
        free(tmp);
        int st; (void)waitpid(pids[t], &st, 0);
    }

    free(p2c); free(c2p); free(pids);
    return 0;

fail_collect:
    for (int t2 = 0; t2 < launched; ++t2) { 
        close(c2p[t2][0]);
        int st; (void)waitpid(pids[t2], &st, 0);
    }
    free(p2c); free(c2p); free(pids);
    return -1;

fail_launch:
    for (int t2 = 0; t2 < launched; ++t2) {
        close(p2c[t2][1]); close(c2p[t2][0]);
        int st; (void)waitpid(pids[t2], &st, 0);
    }
    free(p2c); free(c2p); free(pids);
    return -1;
}


static void deflate_rank1(double *W, int n, double lambda, const double *v) {
    for (int i=0;i<n;++i) {
        double li = lambda * v[i];
        double *row = &W[(size_t)i*(size_t)n];
        for (int j=0;j<n;++j) row[j] -= li * v[j];
    }
}

/* ===== API: eigen_power_distributed ===== */
int eigen_power_distributed(const Matrix *A, int k, int max_iter, double tol,
                            Matrix *V_out, Matrix *L_out) {
    if (!A || !A->data || !V_out || !L_out) return -1;
    if (A->rows != A->cols) return -2;
    int n = A->rows;
    if (k < 1 || k > n) return -3;

    
    long cores = sysconf(_SC_NPROCESSORS_ONLN);
    if (cores < 1) cores = 1;
    int POOL = (int)(cores * 2);
    if (POOL > 32) POOL = 32;
    if (POOL < 1)  POOL = 1;
    signal(SIGPIPE, SIG_IGN);

   
    double *W = (double*)malloc((size_t)n*(size_t)n*sizeof(double));
    if (!W) return -1;
    for (int r=0; r<n; ++r)
        for (int c=0; c<n; ++c)
            W[(size_t)r*(size_t)n + c] = A->data[(size_t)r*(size_t)A->cols + c];

    
    if (mat_create(V_out, "EIGVECS", n, k) != 0) { free(W); return -1; }
    if (mat_create(L_out, "EIGVALS", k, 1) != 0) { mat_destroy(V_out); free(W); return -1; }

    double *x   = (double*)malloc((size_t)n*sizeof(double));
    double *y   = (double*)malloc((size_t)n*sizeof(double));
    double *x_old = (double*)malloc((size_t)n*sizeof(double));
    if (!x || !y || !x_old) { free(x); free(y); free(x_old); mat_destroy(V_out); mat_destroy(L_out); free(W); return -1; }

    for (int eig=0; eig<k; ++eig) {
        
        for (int i=0;i<n;++i) x[i] = 1.0;
        double nrm = vec_norm2(x, n); if (nrm == 0.0) nrm = 1.0; vec_scale(x, n, 1.0/nrm);

        double lambda = 0.0;
        for (int it=0; it<max_iter; ++it) {
            vec_copy(x_old, x, n);

            
            if (matvec_distributed(W, n, x, y, POOL) != 0) {
                free(x); free(y); free(x_old); mat_destroy(V_out); mat_destroy(L_out); free(W);
                return -1;
            }

           
            double ny = vec_norm2(y, n);
            if (ny == 0.0) { 
               
                for (int i=0;i<n;++i) x[i] = 0.0;
                x[it % n] = 1.0;
                continue;
            }
            for (int i=0;i<n;++i) x[i] = y[i] / ny;

           
            if (matvec_distributed(W, n, x, y, POOL) != 0) {
                free(x); free(y); free(x_old); mat_destroy(V_out); mat_destroy(L_out); free(W);
                return -1;
            }
            lambda = dot(x, y, n);

         
            double diff = 0.0;
            for (int i=0;i<n;++i) {
                double d = x[i] - x_old[i];
                diff += d*d;
            }
            if (sqrt(diff) < tol) break;
        }

      
        for (int i=0;i<n;++i) V_out->data[(size_t)i*(size_t)V_out->cols + eig] = x[i];
        L_out->data[eig] = lambda;

        
        deflate_rank1(W, n, lambda, x);
    }

    free(x); free(y); free(x_old); free(W);
    return 0;
}





static void matvec_serial(const double *W, int n,
                          const double *x, double *y) {
    for (int i = 0; i < n; ++i) {
        double acc = 0.0;
        const double *row = &W[(size_t)i * (size_t)n];
        for (int j = 0; j < n; ++j) {
            acc += row[j] * x[j];
        }
        y[i] = acc;
    }
}
int eigen_power_serial(const Matrix *A, int k, int max_iter, double tol,
                       Matrix *V_out, Matrix *L_out) {
    if (!A || !A->data || !V_out || !L_out) return -1;
    if (A->rows != A->cols) return -2;
    int n = A->rows;
    if (k < 1 || k > n) return -3;

   
    double *W = (double*)malloc((size_t)n * (size_t)n * sizeof(double));
    if (!W) return -1;
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c)
            W[(size_t)r*(size_t)n + c] =
                A->data[(size_t)r*(size_t)A->cols + c];

    
    if (mat_create(V_out, "EIGVECS_SER", n, k) != 0) {
        free(W); 
        return -1;
    }
    if (mat_create(L_out, "EIGVALS_SER", k, 1) != 0) {
        mat_destroy(V_out);
        free(W);
        return -1;
    }

    double *x     = (double*)malloc((size_t)n * sizeof(double));
    double *y     = (double*)malloc((size_t)n * sizeof(double));
    double *x_old = (double*)malloc((size_t)n * sizeof(double));
    if (!x || !y || !x_old) {
        free(x); free(y); free(x_old);
        mat_destroy(V_out); mat_destroy(L_out);
        free(W);
        return -1;
    }

    for (int eig = 0; eig < k; ++eig) {
       
        for (int i = 0; i < n; ++i) x[i] = 1.0;
        double nrm = vec_norm2(x, n);
        if (nrm == 0.0) nrm = 1.0;
        vec_scale(x, n, 1.0/nrm);

        double lambda = 0.0;

        for (int it = 0; it < max_iter; ++it) {
            vec_copy(x_old, x, n);

            
            matvec_serial(W, n, x, y);

            double ny = vec_norm2(y, n);
            if (ny == 0.0) {
               
                for (int i = 0; i < n; ++i) x[i] = 0.0;
                x[it % n] = 1.0;
                continue;
            }
            for (int i = 0; i < n; ++i) x[i] = y[i] / ny;

           
            matvec_serial(W, n, x, y);
            lambda = dot(x, y, n);

            double diff = 0.0;
            for (int i = 0; i < n; ++i) {
                double d = x[i] - x_old[i];
                diff += d*d;
            }
            if (sqrt(diff) < tol) break;
        }

       
        for (int i = 0; i < n; ++i)
            V_out->data[(size_t)i*(size_t)V_out->cols + eig] = x[i];
        L_out->data[eig] = lambda;

        
        deflate_rank1(W, n, lambda, x);
    }

    free(x); free(y); free(x_old); free(W);
    return 0;
}

/*  (OpenMP, single-process) */
static void matvec_omp(const double *W, int n,
                       const double *x, double *y) {
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        double acc = 0.0;
        const double *row = &W[(size_t)i*(size_t)n];
        for (int j = 0; j < n; ++j) {
            acc += row[j] * x[j];
        }
        y[i] = acc;
    }
}

int eigen_power_omp(const Matrix *A, int k, int max_iter, double tol,
                    Matrix *V_out, Matrix *L_out) {
    if (!A || !A->data || !V_out || !L_out) return -1;
    if (A->rows != A->cols) return -2;
    int n = A->rows;
    if (k < 1 || k > n) return -3;

    
    double *W = (double*)malloc((size_t)n*(size_t)n*sizeof(double));
    if (!W) return -1;
    for (int r=0; r<n; ++r)
        for (int c=0; c<n; ++c)
            W[(size_t)r*(size_t)n + c] = A->data[(size_t)r*(size_t)A->cols + c];

    if (mat_create(V_out, "EIGVECS_OMP", n, k) != 0) { free(W); return -1; }
    if (mat_create(L_out, "EIGVALS_OMP", k, 1) != 0) {
        mat_destroy(V_out); free(W); return -1;
    }

    double *x     = (double*)malloc((size_t)n*sizeof(double));
    double *y     = (double*)malloc((size_t)n*sizeof(double));
    double *x_old = (double*)malloc((size_t)n*sizeof(double));
    if (!x || !y || !x_old) {
        free(x); free(y); free(x_old);
        mat_destroy(V_out); mat_destroy(L_out); free(W);
        return -1;
    }

    for (int eig=0; eig<k; ++eig) {
        
        for (int i=0;i<n;++i) x[i] = 1.0;
        double nrm = vec_norm2(x, n);
        if (nrm == 0.0) nrm = 1.0;
        vec_scale(x, n, 1.0/nrm);

        double lambda = 0.0;
        for (int it=0; it<max_iter; ++it) {
            vec_copy(x_old, x, n);

            /*  (OpenMP) */
            matvec_omp(W, n, x, y);

            double ny = vec_norm2(y, n);
            if (ny == 0.0) {
                for (int i=0;i<n;++i) x[i] = 0.0;
                x[it % n] = 1.0;
                continue;
            }
            for (int i=0;i<n;++i) x[i] = y[i] / ny;

           
            matvec_omp(W, n, x, y);
            lambda = dot(x, y, n);

            double diff = 0.0;
            for (int i=0;i<n;++i) {
                double d = x[i] - x_old[i];
                diff += d*d;
            }
            if (sqrt(diff) < tol) break;
        }

        
        for (int i=0;i<n;++i)
            V_out->data[(size_t)i*(size_t)V_out->cols + eig] = x[i];
        L_out->data[eig] = lambda;

        
        deflate_rank1(W, n, lambda, x);
    }

    free(x); free(y); free(x_old); free(W);
    return 0;
}


