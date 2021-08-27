//
// Created by adamk on 22/08/2021.
//

#include <math.h>
#include <stdlib.h>
#include "utils.h"

/* find ideal k given eigenvalues*/
int get_elbow_k(double *eigenvalues, int n) {
    int i, k;
    double max = -1;
    bubbleSort(eigenvalues, n);
    double *gaps = calloc(n/2 + 1, sizeof (double)); // extra element for easy indexing
    for (i = 1; i < n/2 + 1 ; i++) {
        gaps[i] = fabs(eigenvalues[i] - eigenvalues[i+1]);
        if (gaps[i] > max) {
            k = i;
            max = gaps[i];
        }
    }
    return k;
}
