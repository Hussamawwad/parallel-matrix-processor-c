#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>

/* ---------- list lifecycle ---------- */
void list_init(MatrixList *L){ L->n = 0; }
void list_clear(MatrixList *L){
    for (int i=0;i<L->n;i++) free(L->mats[i].data);
    L->n = 0;
}

/* ---------- matrix basics ---------- */
int mat_create(Matrix *M, const char *name, int r, int c){
    if (r<=0 || c<=0) return -1;
    strncpy(M->name, name && name[0] ? name : "M", sizeof(M->name)-1);
    M->name[sizeof(M->name)-1] = 0;
    M->rows = r; M->cols = c;
    M->data = (double*)calloc((size_t)r*c, sizeof(double));
    return M->data ? 0 : -1;
}

void mat_destroy(Matrix *M){
    free(M->data);
    M->data = NULL;
    M->rows = M->cols = 0;
    M->name[0] = 0;
}

void mat_print(const Matrix *M){
    if (!M || !M->data) { puts("(null)"); return; }
    printf("Matrix \"%s\" (%dx%d):\n", M->name, M->rows, M->cols);
    for (int i=0;i<M->rows;i++){
        for (int j=0;j<M->cols;j++){
            printf("%8.3f", getcM(M,i,j));
        }
        puts("");
    }
}

/* ---------- list ops ---------- */
int list_add(MatrixList *L, const Matrix *src){
    if (L->n >= MAX_MATS || !src || !src->data) return -1;
    Matrix *dst = &L->mats[L->n];
    if (mat_create(dst, src->name, src->rows, src->cols) != 0) return -1;
    memcpy(dst->data, src->data, (size_t)src->rows*src->cols*sizeof(double));
    return L->n++;
}

int list_delete(MatrixList *L, int idx){
    if (idx < 0 || idx >= L->n) return -1;
    mat_destroy(&L->mats[idx]);
    for (int i=idx+1;i<L->n;i++) L->mats[i-1] = L->mats[i];
    L->n--;
    return 0;
}

void list_print_all(const MatrixList *L){
    if (L->n == 0) { puts("(no matrices in memory)"); return; }
    for (int i=0;i<L->n;i++){
        printf("[%d] %s (%dx%d)\n", i, L->mats[i].name, L->mats[i].rows, L->mats[i].cols);
    }
}

/* ---------- modify ---------- */
int mat_set_value(Matrix *M, int r, int c, double v){
    if (!M || !M->data) return -1;
    if (r<0 || r>=M->rows || c<0 || c>=M->cols) return -1;
    *at(M,r,c) = v;
    return 0;
}

int mat_set_row(Matrix *M, int r, const double *row){
    if (!M || !M->data) return -1;
    if (r<0 || r>=M->rows) return -1;
    for (int j=0;j<M->cols;j++) *at(M,r,j) = row[j];
    return 0;
}

int mat_set_col(Matrix *M, int c, const double *col){
    if (!M || !M->data) return -1;
    if (c<0 || c>=M->cols) return -1;
    for (int i=0;i<M->rows;i++) *at(M,i,c) = col[i];
    return 0;
}



static int read_one_file(Matrix *out, const char *filepath){
    FILE *f = fopen(filepath, "r");
    if (!f) return -1;
    char nm[MAX_NAME]; int r,c;
    if (fscanf(f, " %31s %d %d", nm, &r, &c) != 3) { fclose(f); return -1; }
    if (mat_create(out, nm, r, c) != 0) { fclose(f); return -1; }
    for (int i=0;i<r;i++){
        for (int j=0;j<c;j++){
            if (fscanf(f, " %lf", at(out,i,j)) != 1){
                mat_destroy(out); fclose(f); return -1;
            }
        }
    }
    fclose(f);
    return 0;
}

int mat_read_file(MatrixList *L, const char *filepath){
    Matrix tmp = {0};
    if (read_one_file(&tmp, filepath) != 0) return -1;
    int idx = list_add(L, &tmp);
    mat_destroy(&tmp);
    return idx;
}

int list_read_folder(MatrixList *L, const char *folder){
    DIR *d = opendir(folder);
    if (!d) return -1;
    struct dirent *de;
    int added = 0;
    char path[1024];
    while ((de = readdir(d))){
        if (de->d_name[0] == '.') continue;
        /* accept .txt or .mtx */
        const char *name = de->d_name;
        size_t n = strlen(name);
        int ok = (n>4 && (strcmp(name+n-4, ".txt")==0 || strcmp(name+n-4, ".mtx")==0));
        if (!ok) continue;
        snprintf(path, sizeof(path), "%s/%s", folder, name);
        if (mat_read_file(L, path) >= 0) added++;
    }
    closedir(d);
    return added;
}

int mat_write_file(const Matrix *M, const char *filepath){
    FILE *f = fopen(filepath, "w");
    if (!f) return -1;
    fprintf(f, "%s %d %d\n", M->name, M->rows, M->cols);
    for (int i=0;i<M->rows;i++){
        for (int j=0;j<M->cols;j++){
            fprintf(f, (j+1==M->cols) ? "%.10g" : "%.10g ", getcM(M,i,j));
        }
        fputc('\n', f);
    }
    fclose(f);
    return 0;
}

int list_save_all(const MatrixList *L, const char *folder){
    /* ensure folder exists */
    struct stat st;
    if (stat(folder, &st) != 0){
        if (mkdir(folder, 0777) != 0 && errno != EEXIST) return -1;
    }
    char path[1024];
    int ok = 0;
    for (int i=0;i<L->n;i++){
        snprintf(path, sizeof(path), "%s/%s.txt", folder, L->mats[i].name[0]?L->mats[i].name:"M");
        if (mat_write_file(&L->mats[i], path) == 0) ok++;
    }
    return ok;
}
