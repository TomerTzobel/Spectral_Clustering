//
// Created by adamk on 22/08/2021.
//

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define ERR_MSG "An Error Has Occured"

double **init_matrix(int rows, int cols);

void print_matrix(double **matrix, int rows, int cols);

void free_matrix(double **matrix, int rows);

double **multiply_matrices(double **A, int rows_A, int cols_A, double **B, int rows_B, int cols_B);

void bubbleSort(double arr[], int n);

/* init matrix of zeros */
double **init_matrix(int rows, int cols) {
    double **matrix = malloc(rows * sizeof(double *));
    assert(matrix != NULL && ERR_MSG);
    for (int i = 0; i < rows; i++) {
        matrix[i] = calloc(cols, sizeof(double));
        assert(matrix[i] != NULL && ERR_MSG);
    }
    return matrix;
}

void print_matrix(double **matrix, int rows, int cols) {
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            printf("%.4f", matrix[i][j]);
            if (j < cols - 1)
                printf(",");
        }
        printf("\n");
    }
}

void free_matrix(double **matrix, int rows) {
    for (int i = 0; i < rows; i++)
        free(matrix[i]);
    free(matrix);
}

double **multiply_matrices(double **A, int rows_A, int cols_A, double **B, int rows_B, int cols_B) {
    assert(cols_A == rows_B && "invalid multiplication");
    double **res = init_matrix(rows_A, cols_B);
    for (int i = 0; i < rows_A; ++i) {
        for (int j = 0; j < cols_B; ++j) {
            for (int k = 0; k < cols_A; ++k) {
                res[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

/* code by GFG, with our minor optimization */
void swap(double *xp, double *yp) {
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}

/* code by GFG, with our minor optimization */
void bubbleSort(double arr[], int n) {
    int i, j, swapped;
    for (i = 0; i < n - 1; i++) {
        swapped = 0;
        // Last i elements are already in place
        for (j = 0; j < n - i - 1; j++) {
            if (arr[j] > arr[j + 1]) {
                swap(&arr[j], &arr[j + 1]);
                swapped = 1;
            }
        }
        if (!swapped) {
            break;
        }
    }
}
