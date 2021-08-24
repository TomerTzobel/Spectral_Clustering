//
// Created by adamk on 24/08/2021.
//
#include <math.h>
#include "utils.c"

double **wam(double **datapoints, int n, int dimension) {
    double **wamMatrix;
    int i, j;
    double pointNorm, result;
    wamMatrix = init_matrix(n, n);
    for (i = 1; i < n - 1; i++) {
        for (j = i + 1; j < n; j++) {
            pointNorm = norm(datapoints[i], datapoints[j], dimension);
            result = (double) pow(M_E, ((-pointNorm) / 2)); //M_E is the mathematical constant e
            wamMatrix[i][j] = result;
            wamMatrix[j][i] = result; // symmetry of W
        }
    }
    return wamMatrix;
}
