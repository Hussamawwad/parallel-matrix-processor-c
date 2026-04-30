#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h>

#define MAX_MATS 32
#define MAX_NAME 32

typedef struct {
    char name[MAX_NAME];   
    int rows, cols;
    double *data;          
} Matrix;

typedef struct {
    Matrix mats[MAX_MATS];
    int n;
} MatrixList;

/* list lifecycle */
void list_init(MatrixList *L);
void list_clear(MatrixList *L);

/* matrix basics */
int  mat_create(Matrix *M, const char *name, int r, int c);
void mat_destroy(Matrix *M);
void mat_print(const Matrix *M);

/* helpers */
static inline double *at(Matrix *M, int r, int c) { return &M->data[(size_t)r * M->cols + c]; }
static inline double  getcM(const Matrix *M, int r, int c) { return M->data[(size_t)r * M->cols + c]; }

/* list ops (1–3,9) */
int  list_add(MatrixList *L, const Matrix *src);     // returns index or -1
int  list_delete(MatrixList *L, int idx);            // 0 ok, -1 bad idx
void list_print_all(const MatrixList *L);

/* modify (4) */
int  mat_set_value(Matrix *M, int r, int c, double v);
int  mat_set_row(Matrix *M, int r, const double *row);
int  mat_set_col(Matrix *M, int c, const double *col);

/* file I/O (5–8) */
int  mat_read_file(MatrixList *L, const char *filepath);           // adds into list
int  list_read_folder(MatrixList *L, const char *folder);          // adds many
int  mat_write_file(const Matrix *M, const char *filepath);        // save one
int  list_save_all(const MatrixList *L, const char *folder);       // save many

#endif
