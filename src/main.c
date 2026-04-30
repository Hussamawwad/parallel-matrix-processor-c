#include "matrix.h"
#include "util.h"
#include "ops_addsub.h"
#include "ops_mul.h"
#include "ops_det.h"
#include "ops_eigen.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static void print_menu(void){
    puts("=== Matrix Operations (Signals/Pipes/FIFOs Project) ===");
    puts("1. Enter a matrix");
    puts("2. Display a matrix");
    puts("3. Delete a matrix");
    puts("4. Modify a matrix (full row, full column, a particular value)");
    puts("5. Read a matrix from a file");
    puts("6. Read a set of matrices from a folder");
    puts("7. Save a matrix to a file");
    puts("8. Save all matrices in memory to a folder");
    puts("9. Display all matrices in memory");
    puts("10. Add 2 matrices");
    puts("11. Subtract 2 matrices");
    puts("12. Multiply 2 matrices");
    puts("13. Find the determinant");
    puts("14. Eigenvalues/vectors");
    puts("15. Exit");
}

static void op1_enter(MatrixList *L){
    if (L->n >= MAX_MATS){ puts("Memory full."); return; }
    char name[MAX_NAME];
    read_line("Name (e.g., A): ", name, sizeof(name));
    if (!name[0]) strcpy(name, "M");
    int r = read_int("Rows: ", 3);
    int c = read_int("Cols: ", 3);
    Matrix M = {0};
    if (mat_create(&M, name, r, c) != 0){ puts("Allocation failed."); return; }
    for (int i=0;i<r;i++){
        for (int j=0;j<c;j++){
            char prompt[64];
            snprintf(prompt, sizeof(prompt), "  [%d,%d] = ", i, j);
            *at(&M,i,j) = read_double(prompt, 0.0);
        }
    }
    int idx = list_add(L, &M);
    mat_destroy(&M);
    if (idx < 0) puts("Add failed.");
    else printf("Added at index %d.\n", idx);
}

static void op2_display(const MatrixList *L){
    if (L->n == 0){ puts("(no matrices)"); return; }
    list_print_all(L);
    int idx = read_int("Index to display: ", 0);
    if (idx < 0 || idx >= L->n){ puts("Bad index."); return; }
    mat_print(&L->mats[idx]);     
}

static void op3_delete(MatrixList *L){
    if (L->n == 0){ puts("(no matrices)"); return; }
    list_print_all(L);
    int idx = read_int("Index to delete: ", 0);
    if (list_delete(L, idx) == 0) puts("Deleted.");
    else puts("Bad index.");
}

static void op4_modify(MatrixList *L){
    if (L->n == 0){ puts("(no matrices)"); return; }
    list_print_all(L);
    int idx = read_int("Index to modify: ", 0);
    if (idx < 0 || idx >= L->n){ puts("Bad index."); return; }
    Matrix *M = &L->mats[idx];

    puts("Modify options:");
    puts(" 1) Set one value");
    puts(" 2) Set full row");
    puts(" 3) Set full column");
    int ch = read_int("Choose [1-3]: ", 1);
    if (ch == 1){
        int r = read_int("  row: ", 0);
        int c = read_int("  col: ", 0);
        double v = read_double("  value: ", 0);
        if (mat_set_value(M, r, c, v) == 0) puts("Updated.");
        else puts("Bad (r,c).");
    } else if (ch == 2){
        int r = read_int("  row index: ", 0);
        if (r < 0 || r >= M->rows){ puts("Bad row."); return; }
        double *row = (double*)malloc((size_t)M->cols * sizeof(double));
        for (int j=0;j<M->cols;j++){
            char pr[64]; snprintf(pr, sizeof(pr), "  [%d,%d] = ", r, j);
            row[j] = read_double(pr, 0);
        }
        mat_set_row(M, r, row); free(row); puts("Row updated.");
    } else if (ch == 3){
        int c = read_int("  column index: ", 0);
        if (c < 0 || c >= M->cols){ puts("Bad column."); return; }
        double *col = (double*)malloc((size_t)M->rows * sizeof(double));
        for (int i=0;i<M->rows;i++){
            char pr[64]; snprintf(pr, sizeof(pr), "  [%d,%d] = ", i, c);
            col[i] = read_double(pr, 0);
        }
        mat_set_col(M, c, col); free(col); puts("Column updated.");
    } else {
        puts("Unknown choice.");
    }
}

static void op5_read_file(MatrixList *L){
    char path[1024];
    read_line("File path to read: ", path, sizeof(path));
    int idx = mat_read_file(L, path);
    if (idx >= 0) printf("Loaded at index %d.\n", idx);
    else puts("Failed to read file.");
}

static void op6_read_folder(MatrixList *L){
    char dir[1024];
    read_line("Folder to scan: ", dir, sizeof(dir));
    int cnt = list_read_folder(L, dir);
    if (cnt >= 0) printf("Loaded %d matrices.\n", cnt);
    else puts("Failed to read folder.");
}

static void op7_save_one(const MatrixList *L){
    if (L->n == 0){ puts("(no matrices)"); return; }
    list_print_all(L);
    int idx = read_int("Index to save: ", 0);
    if (idx < 0 || idx >= L->n){ puts("Bad index."); return; }
    char path[1024];
    read_line("Output file path: ", path, sizeof(path));
    if (mat_write_file(&L->mats[idx], path) == 0) puts("Saved.");
    else puts("Write failed.");
}

static void op8_save_all(const MatrixList *L){
    char dir[1024];
    read_line("Output folder: ", dir, sizeof(dir));
    int cnt = list_save_all(L, dir);
    if (cnt >= 0) printf("Saved %d matrices.\n", cnt);
    else puts("Save failed.");
}

static void op9_display_all(const MatrixList *L){
    list_print_all(L);
}

static int choose_two(const MatrixList *L, int *iA, int *iB) {
    if (L->n < 2) { puts("Need at least 2 matrices."); return -1; }
    list_print_all(L);
    *iA = read_int("First index: ", 0);
    *iB = read_int("Second index: ", 1);
    if (*iA < 0 || *iA >= L->n || *iB < 0 || *iB >= L->n || *iA == *iB) {
        puts("Bad indices.");
        return -1;
    }
    return 0;
}

/* =========  Add (distributed + serial + OpenMP timing) ========= */
static void op10_add(MatrixList *L){
    int iA, iB;
    if (choose_two(L, &iA, &iB) != 0) return;

    Matrix OUTd = {0};  /* distributed */
    Matrix OUTs = {0};  /* serial */
    Matrix OUTo = {0};  /* serial + OpenMP */

    /* ---- distributed ---- */
    clock_t t1 = clock();
    int rc = add_matrices_distributed(&L->mats[iA], &L->mats[iB], &OUTd);
    clock_t t2 = clock();

    if (rc == -2) {
        puts("Dimension mismatch.");
        return;
    }
    if (rc != 0) {
        puts("Add (distributed) failed.");
        return;
    }
    double time_distributed = (double)(t2 - t1) / CLOCKS_PER_SEC;

    /* ---- serial ---- */
    t1 = clock();
    rc = add_matrices_serial(&L->mats[iA], &L->mats[iB], &OUTs);
    t2 = clock();

    if (rc == -2) {
        mat_destroy(&OUTd);
        puts("Dimension mismatch (serial).");
        return;
    }
    if (rc != 0) {
        mat_destroy(&OUTd);
        puts("Add (serial) failed.");
        return;
    }
    double time_serial = (double)(t2 - t1) / CLOCKS_PER_SEC;

    /* ---- serial + OpenMP ---- */
    t1 = clock();
    rc = add_matrices_omp(&L->mats[iA], &L->mats[iB], &OUTo);
    t2 = clock();

    if (rc == -2) {
        mat_destroy(&OUTd);
        mat_destroy(&OUTs);
        puts("Dimension mismatch (OpenMP).");
        return;
    }
    if (rc != 0) {
        mat_destroy(&OUTd);
        mat_destroy(&OUTs);
        puts("Add (OpenMP) failed.");
        return;
    }
    double time_omp = (double)(t2 - t1) / CLOCKS_PER_SEC;

   
    int idx_d = list_add(L, &OUTd);
    int idx_s = list_add(L, &OUTs);
    int idx_o = list_add(L, &OUTo);
    mat_destroy(&OUTd);
    mat_destroy(&OUTs);
    mat_destroy(&OUTo);

    if (idx_d >= 0 && idx_s >= 0 && idx_o >= 0) {
        printf("Distributed sum stored at index %d.\n", idx_d);
        printf("Serial sum        stored at index %d.\n", idx_s);
        printf("OpenMP sum        stored at index %d.\n", idx_o);
        printf("⏱ Multi-process (distributed) time: %.6f s\n", time_distributed);
        printf("⏱ Single-threaded (serial) time   : %.6f s\n", time_serial);
        printf("⏱ Serial + OpenMP time            : %.6f s\n", time_omp);
    } else {
        puts("Store failed.");
    }
}

/* =========  Sub (distributed + serial + OpenMP timing) ========= */
static void op11_sub(MatrixList *L){
    int iA, iB;
    if (choose_two(L, &iA, &iB) != 0) return;

    Matrix OUTd = {0};  /* distributed */
    Matrix OUTs = {0};  /* serial */
    Matrix OUTo = {0};  /* serial + OpenMP */

    /* ---- distributed ---- */
    clock_t t1 = clock();
    int rc = sub_matrices_distributed(&L->mats[iA], &L->mats[iB], &OUTd);
    clock_t t2 = clock();

    if (rc == -2) {
        puts("Dimension mismatch.");
        return;
    }
    if (rc != 0) {
        puts("Sub (distributed) failed.");
        return;
    }
    double time_distributed = (double)(t2 - t1) / CLOCKS_PER_SEC;

    /* ---- serial ---- */
    t1 = clock();
    rc = sub_matrices_serial(&L->mats[iA], &L->mats[iB], &OUTs);
    t2 = clock();

    if (rc == -2) {
        mat_destroy(&OUTd);
        puts("Dimension mismatch (serial).");
        return;
    }
    if (rc != 0) {
        mat_destroy(&OUTd);
        puts("Sub (serial) failed.");
        return;
    }
    double time_serial = (double)(t2 - t1) / CLOCKS_PER_SEC;

    /* ---- serial + OpenMP ---- */
    t1 = clock();
    rc = sub_matrices_omp(&L->mats[iA], &L->mats[iB], &OUTo);
    t2 = clock();

    if (rc == -2) {
        mat_destroy(&OUTd);
        mat_destroy(&OUTs);
        puts("Dimension mismatch (OpenMP).");
        return;
    }
    if (rc != 0) {
        mat_destroy(&OUTd);
        mat_destroy(&OUTs);
        puts("Sub (OpenMP) failed.");
        return;
    }
    double time_omp = (double)(t2 - t1) / CLOCKS_PER_SEC;

    
    int idx_d = list_add(L, &OUTd);
    int idx_s = list_add(L, &OUTs);
    int idx_o = list_add(L, &OUTo);
    mat_destroy(&OUTd);
    mat_destroy(&OUTs);
    mat_destroy(&OUTo);

    if (idx_d >= 0 && idx_s >= 0 && idx_o >= 0) {
        printf("Distributed diff stored at index %d.\n", idx_d);
        printf("Serial diff        stored at index %d.\n", idx_s);
        printf("OpenMP diff        stored at index %d.\n", idx_o);
        printf("⏱ Multi-process (distributed) time: %.6f s\n", time_distributed);
        printf("⏱ Single-threaded (serial) time   : %.6f s\n", time_serial);
        printf("⏱ Serial + OpenMP time            : %.6f s\n", time_omp);
    } else {
        puts("Store failed.");
    }
}


/* =========  Mul (distributed + serial + OpenMP timing) ========= */
static void op12_mul(MatrixList *L){
    if (L->n < 2) {
        puts("Need at least 2 matrices.");
        return;
    }
    list_print_all(L);
    int iA = read_int("Left matrix index (A): ", 0);
    int iB = read_int("Right matrix index (B): ", 1);
    if (iA < 0 || iA >= L->n || iB < 0 || iB >= L->n || iA == iB) {
        puts("Bad indices.");
        return;
    }

    Matrix OUTd = {0};  /* distributed */
    Matrix OUTs = {0};  /* serial      */
    Matrix OUTo = {0};  /* serial + OpenMP */

    /* ----------  (multi-process) ---------- */
    clock_t t1 = clock();
    int rc = mul_matrices_distributed(&L->mats[iA], &L->mats[iB], &OUTd);
    clock_t t2 = clock();

    if (rc == -2) {
        puts("Dimension mismatch: A.cols must equal B.rows.");
        return;
    }
    if (rc != 0) {
        puts("Multiply (distributed) failed.");
        return;
    }
    double time_distributed = (double)(t2 - t1) / CLOCKS_PER_SEC;

    /* ----------  (single-threaded) ---------- */
    t1 = clock();
    rc = mul_matrices_serial(&L->mats[iA], &L->mats[iB], &OUTs);
    t2 = clock();

    if (rc == -2) {
        mat_destroy(&OUTd);
        puts("Dimension mismatch (serial).");
        return;
    }
    if (rc != 0) {
        mat_destroy(&OUTd);
        puts("Multiply (serial) failed.");
        return;
    }
    double time_serial = (double)(t2 - t1) / CLOCKS_PER_SEC;

    /* ----------  OpenMP (serial + OpenMP) ---------- */
    t1 = clock();
    rc = mul_matrices_omp(&L->mats[iA], &L->mats[iB], &OUTo);
    t2 = clock();

    if (rc == -2) {
        mat_destroy(&OUTd);
        mat_destroy(&OUTs);
        puts("Dimension mismatch (OpenMP).");
        return;
    }
    if (rc != 0) {
        mat_destroy(&OUTd);
        mat_destroy(&OUTs);
        puts("Multiply (OpenMP) failed.");
        return;
    }
    double time_omp = (double)(t2 - t1) / CLOCKS_PER_SEC;

   
    int idx_d = list_add(L, &OUTd);
    int idx_s = list_add(L, &OUTs);
    int idx_o = list_add(L, &OUTo);
    mat_destroy(&OUTd);
    mat_destroy(&OUTs);
    mat_destroy(&OUTo);

    if (idx_d >= 0 && idx_s >= 0 && idx_o >= 0) {
        printf("Distributed product stored at index %d.\n", idx_d);
        printf("Single-threaded product stored at index %d.\n", idx_s);
        printf("OpenMP product        stored at index %d.\n", idx_o);
        printf("⏱ Multi-process (distributed) time: %.6f s\n", time_distributed);
        printf("⏱ Single-threaded (serial) time   : %.6f s\n", time_serial);
        printf("⏱ Serial + OpenMP time            : %.6f s\n", time_omp);
    } else {
        puts("Store failed.");
    }
}

/* =========  Determinant (distributed + serial + OpenMP) ========= */
static void op13_det(MatrixList *L){
    if (L->n < 1) { puts("No matrices."); return; }
    list_print_all(L);
    int idx = read_int("Matrix index: ", 0);
    if (idx < 0 || idx >= L->n) { puts("Bad index."); return; }

    const Matrix *A = &L->mats[idx];

    double d_dist = 0.0;
    double d_ser  = 0.0;
    double d_omp  = 0.0;

    /* ---- multi-process (distributed) ---- */
    clock_t t1 = clock();
    int rc = determinant_distributed(A, &d_dist);
    clock_t t2 = clock();

    if (rc == -2) { puts("Matrix must be square."); return; }
    if (rc != 0)  { puts("Determinant (distributed) failed."); return; }

    double time_distributed = (double)(t2 - t1) / CLOCKS_PER_SEC;

    /* ----  single-threaded (serial) ---- */
    t1 = clock();
    rc = determinant_serial(A, &d_ser);
    t2 = clock();

    if (rc == -2) { puts("Matrix must be square (serial)."); return; }
    if (rc != 0)  { puts("Determinant (serial) failed."); return; }

    double time_serial = (double)(t2 - t1) / CLOCKS_PER_SEC;

    /* ----  OpenMP ---- */
    t1 = clock();
    rc = determinant_omp(A, &d_omp);
    t2 = clock();

    if (rc == -2) { puts("Matrix must be square (OpenMP)."); return; }
    if (rc != 0)  { puts("Determinant (OpenMP) failed."); return; }

    double time_omp = (double)(t2 - t1) / CLOCKS_PER_SEC;

   
    const char *nm = (A->name && A->name[0]) ? A->name : "A";

    printf("det(%s) [multi-process]   = %.6f\n", nm, d_dist);
    printf("det(%s) [single-threaded] = %.6f\n", nm, d_ser);
    printf("det(%s) [OpenMP]          = %.6f\n", nm, d_omp);

    printf("⏱ Multi-process (distributed) time: %.6f s\n", time_distributed);
    printf("⏱ Single-threaded (serial) time  : %.6f s\n", time_serial);
    printf("⏱ Serial + OpenMP time           : %.6f s\n", time_omp);
}


static void op14_eigen(MatrixList *L){
    if (L->n < 1) { puts("No matrices."); return; }
    list_print_all(L);
    int idx = read_int("Matrix index: ", 0);
    if (idx < 0 || idx >= L->n) { puts("Bad index."); return; }
    const Matrix *A = &L->mats[idx];
    if (A->rows != A->cols) { puts("Matrix must be square."); return; }

    int n = A->rows;
    int k = read_int("How many eigenpairs (1..min(n,3))? ", 1);
    if (k < 1) k = 1;
    if (k > n) k = n;
    if (k > 3) k = 3;

    int max_iter = read_int("Max iterations per eigen (e.g., 1000): ", 1000);
    double tol   = read_double("Tolerance (e.g., 1e-8): ", 1e-8);

    Matrix Vd = {0}, Ld = {0};  /* distributed */
    Matrix Vs = {0}, Ls = {0};  /* serial      */
    Matrix Vo = {0}, Lo = {0};  /* OpenMP      */

    clock_t t1, t2;

    /* --------- distributed eigen --------- */
    t1 = clock();
    int rc = eigen_power_distributed(A, k, max_iter, tol, &Vd, &Ld);
    t2 = clock();

    if (rc == -2) { puts("Matrix must be square."); return; }
    if (rc == -3) { puts("Invalid k."); return; }
    if (rc != 0)  { puts("Eigen (distributed) computation failed."); return; }

    double time_distributed = (double)(t2 - t1) / CLOCKS_PER_SEC;

    int idxVd = list_add(L, &Vd);
    int idxLd = list_add(L, &Ld);
    mat_destroy(&Vd);
    mat_destroy(&Ld);
    if (idxVd < 0 || idxLd < 0) {
        puts("Store (distributed) failed.");
        return;
    }

    /* --------- serial eigen --------- */
    t1 = clock();
    rc = eigen_power_serial(A, k, max_iter, tol, &Vs, &Ls);
    t2 = clock();

    if (rc == -2) { puts("Matrix must be square (serial)."); return; }
    if (rc == -3) { puts("Invalid k (serial)."); return; }
    if (rc != 0)  { puts("Eigen (serial) computation failed."); return; }

    double time_serial = (double)(t2 - t1) / CLOCKS_PER_SEC;

    int idxVs = list_add(L, &Vs);
    int idxLs = list_add(L, &Ls);
    mat_destroy(&Vs);
    mat_destroy(&Ls);
    if (idxVs < 0 || idxLs < 0) {
        puts("Store (serial) failed.");
        return;
    }

    /* --------- OpenMP eigen --------- */
    t1 = clock();
    rc = eigen_power_omp(A, k, max_iter, tol, &Vo, &Lo);
    t2 = clock();

    if (rc == -2) { puts("Matrix must be square (OpenMP)."); return; }
    if (rc == -3) { puts("Invalid k (OpenMP)."); return; }
    if (rc != 0)  { puts("Eigen (OpenMP) computation failed."); return; }

    double time_omp = (double)(t2 - t1) / CLOCKS_PER_SEC;

    int idxVo = list_add(L, &Vo);
    int idxLo = list_add(L, &Lo);
    mat_destroy(&Vo);
    mat_destroy(&Lo);
    if (idxVo < 0 || idxLo < 0) {
        puts("Store (OpenMP) failed.");
        return;
    }

   
    printf("Distributed eigenvectors at index %d, eigenvalues at index %d.\n",
           idxVd, idxLd);
    printf("Serial eigenvectors      at index %d, eigenvalues at index %d.\n",
           idxVs, idxLs);
    printf("OpenMP eigenvectors      at index %d, eigenvalues at index %d.\n",
           idxVo, idxLo);

    printf("⏱ Eigen (multi-process) time: %.6f s\n", time_distributed);
    printf("⏱ Eigen (serial)        time: %.6f s\n", time_serial);
    printf("⏱ Eigen (OpenMP)        time: %.6f s\n", time_omp);

  
    printf("Distributed eigenvalues: ");
    for (int i = 0; i < L->mats[idxLd].rows; ++i) {
        printf("%.6f%s", L->mats[idxLd].data[i],
               (i+1 == L->mats[idxLd].rows ? "\n" : "  "));
    }

    printf("Serial eigenvalues     : ");
    for (int i = 0; i < L->mats[idxLs].rows; ++i) {
        printf("%.6f%s", L->mats[idxLs].data[i],
               (i+1 == L->mats[idxLs].rows ? "\n" : "  "));
    }

    printf("OpenMP eigenvalues     : ");
    for (int i = 0; i < L->mats[idxLo].rows; ++i) {
        printf("%.6f%s", L->mats[idxLo].data[i],
               (i+1 == L->mats[idxLo].rows ? "\n" : "  "));
    }
}



/* =======================  config file ======================= */
int main(int argc, char *argv[]){
    MatrixList L;
    list_init(&L);

  
    if (argc > 1) {
        const char *cfg_path = argv[1];
        FILE *cfg = fopen(cfg_path, "r");
        if (!cfg) {
            fprintf(stderr, "Warning: could not open config file '%s'\n", cfg_path);
        } else {
            char line[1024];
            while (fgets(line, sizeof(line), cfg)) {
               
                line[strcspn(line, "\r\n")] = 0;

                if (line[0] == '#' || line[0] == '\0')
                    continue;

             
                if (strncmp(line, "dir=", 4) == 0) {
                    char *dir = line + 4;
                    if (*dir) {
                        int cnt = list_read_folder(&L, dir);
                        if (cnt >= 0) {
                            printf("Loaded %d matrices from directory '%s' (from config).\n",
                                   cnt, dir);
                        } else {
                            fprintf(stderr,
                                    "Failed to load matrices from '%s' specified in config.\n",
                                    dir);
                        }
                    }
                }

               
            }
            fclose(cfg);
        }
    }

    for (;;){
        print_menu();
        int choice = read_int("Choose [1-15]: ", 15);
        switch (choice){
            case 1: op1_enter(&L); break;
            case 2: op2_display(&L); break;
            case 3: op3_delete(&L); break;
            case 4: op4_modify(&L); break;
            case 5: op5_read_file(&L); break;
            case 6: op6_read_folder(&L); break;
            case 7: op7_save_one(&L); break;
            case 8: op8_save_all(&L); break;
            case 9: op9_display_all(&L); break;
            case 10: op10_add(&L); break;
            case 11: op11_sub(&L); break;
            case 12: op12_mul(&L); break;
            case 13: op13_det(&L); break;
            case 14: op14_eigen(&L); break;
            case 15: list_clear(&L); return 0;
            default: puts("Option not implemented yet."); break;
        }
        puts("");
    }
}

