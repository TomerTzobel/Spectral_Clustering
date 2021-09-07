#ifndef SPKMEANS_H
#define SPKMEANS_H

/* The Weighted Adjacency Matrix Functions */
double **wam(double **datapoints, int n, int dimension);

/* The Diagonal Degree Matrix Functions */
double **ddg(double **datapoints, int n, int dimension);

double **eval_ddg(double **weighted_matrix, int n);

/* The Normalized Graph Laplacian Functions */
double **lnorm(double **datapoints, int n, int dimension);

/* Jacobi algorithm Functions */
double get_c(double **A, int i, int j);

double get_s(double **A, int i, int j);

double **get_rotation_matrix(double **A, int n, int i, int j);

double get_theta(double **A, int i, int j);

int sign(double theta);

double get_t(double theta);

double get_c_of_t(double t);

void update_Pivot(int *pivot, double **A, int n);

double frobenius_Norm_Pow(double **A, int n);

int quick_converged(double  **A, int i, int j);

double **transform_A(double **A, int n, int i, int j, double c, double s);

double **jacobi_eigenvectors(double **A, int n);

double **get_normalized_eigenvectors(int *k, double **points, int dimension, int n);

void do_wam(double **points, int pointsNumber, int dimension);

void do_ddg(double **points, int pointsNumber, int dimension);

void do_lnorm(double **points, int pointsNumber, int dimension);

void do_jacobi(double **points, int pointsNumber);

void do_spk_kmeans(double **points, int pointsNumber, int dimension, int k);

/* find ideal k given sorted eigenvalues*/
int get_elbow_k(double *eigenvalues, int n);

/* Kmeans algorithm functions */
double **kmeans(int k, int n, double **points, long *centroids_indices);
int findMinCent(double *point, double **centroids, int k, int dimension);

int UpdateAllAvg(double **centroids, int *clusters, double **points, int k, int pointsnumber, int dimension);

/* Helping methods */
void read_data(const char *filename, double ***points, int *dimension, int *pointsNumber);

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

void print_arr(double *arr, int n);

double **copy_columns_by_order(double **source, int rows, int cols, int *order);

void normalize_matrix(double **matrix, int rows, int cols);

double **transpose_matrix(double **matrix, int rows, int cols);

/* Output results depending on the Goal from the user */
void nsc(int k, char *goal, char *filename);

#endif
