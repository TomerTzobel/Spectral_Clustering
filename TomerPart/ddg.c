

#include "utils.c"
#include "wam.c"



double** eval_ddg(double **weighted_matrix, int n) {
    int i,j;
    double sum;
    double **diagonal_degree_matrix = init_matrix(n,n);
    for (i = 0; i < n; i++) {
        sum = 0;
        for (j = 0; j < n; j++) { // sum row
            sum += weighted_matrix[i][j];
        }
        diagonal_degree_matrix[i][i] = sum;
    }
    return diagonal_degree_matrix;
}

double** ddg(double **datapoints, int n, int dimension) {
    double **weighted_matrix = wam(datapoints, n, dimension);
    double **diag = eval_ddg(weighted_matrix, n);
    free_matrix(weighted_matrix, n);
    return diag;
}