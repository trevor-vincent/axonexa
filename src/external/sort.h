#include <stdio.h>
#include <stdlib.h>

int compare_descend(const void *a, const void *b) {
    if (*(double *)a < *(double *)b) {
        return 1;
    } else if (*(double *)a > *(double *)b) {
        return -1;
    } else {
        return 0;
    }
}

int compare_ascend(const void *a, const void *b) {
    if (*(double *)a < *(double *)b) {
        return -1;
    } else if (*(double *)a > *(double *)b) {
        return 1;
    } else {
        return 0;
    }
}
