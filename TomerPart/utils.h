//
// Created by adamk on 24/08/2021.
//

#ifndef MAIN_C_UTILS_H
#define MAIN_C_UTILS_H

double **init_matrix(int rows, int cols);

void print_matrix(double **matrix, int rows, int cols);

void free_matrix(double **matrix, int rows);

double **multiply_matrices(double **A, int rows_A, int cols_A, double **B, int rows_B, int cols_B);

double **multiply_matrices_same_dim(double **A, double **B, int n);

double norm(double *point1, double *point2, int dimension);

int *bubbleSort_index_tracked(double arr[], int n);

void power_matrix_elementwise(double **matrix, int n, double x);

void swap(double *xp, double *yp);

double **get_I_matrix(int n);

void copy_matrix(double **source, double **dest, int rows, int cols);

double **transpose_matrix(double **matrix, int rows, int cols);

void print_arr(double *arr, int n);

double **copy_columns_by_order(double **source, int rows, int cols, int *order);

void normalize_matrix(double **matrix, int rows, int cols);

#endif //MAIN_C_UTILS_H
