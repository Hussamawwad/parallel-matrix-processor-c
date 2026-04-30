#include "util.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int read_int(const char *prompt, int def) {
    char buf[128];
    printf("%s", prompt);
    if (!fgets(buf, sizeof(buf), stdin)) return def;
    char *e = NULL;
    long v = strtol(buf, &e, 10);
    if (e == buf) return def;
    return (int)v;
}

double read_double(const char *prompt, double def) {
    char buf[128];
    printf("%s", prompt);
    if (!fgets(buf, sizeof(buf), stdin)) return def;
    char *e = NULL;
    double v = strtod(buf, &e);
    if (e == buf) return def;
    return v;
}

void read_line(const char *prompt, char *buf, int max) {
    printf("%s", prompt);
    if (!fgets(buf, max, stdin)) {
        if (max) buf[0] = 0;
        return;
    }
    size_t n = strlen(buf);
    if (n && buf[n-1] == '\n') buf[n-1] = 0;
}
