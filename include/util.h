#ifndef UTIL_H
#define UTIL_H

#include <stdbool.h>
#include "matrix.h"


int    read_int(const char *prompt, int def);
double read_double(const char *prompt, double def);
void   read_line(const char *prompt, char *buf, int max);


bool   file_exists(const char *path);
bool   ensure_dir(const char *path);  


int    load_all_matrices_from_dir(const char *dir, MatrixList *L);        
int    save_all_matrices_to_dir(const char *dir, const MatrixList *L);    

#endif
