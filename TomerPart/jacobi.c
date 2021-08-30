//
// Created by adamk on 24/08/2021.
//
#include <stdio.h> //delete this later, only for priting
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "utils.h"
#define Epsilon 1.0e-15
#define ERR_MSG "An Error Has Occured"

double get_c(double **A, int i, int j);
double get_s(double **A, int i, int j);


/* STEP 2: n is dim, i,j are A_ij of Pivot step, I is identity matrix */
void update_rotation_matrix(double **A, int n, int i, int j, double **I) {
    double c = get_c(A, i, j), s = get_s(A,i,j);
    I[i][i] = c, I[i][j] = s, I[j][i] = -s, I[j][j] = c;
}

void revert_P_to_Identity(double **P, int i, int j) {
    P[i][i] = 1, P[i][j] = 0, P[j][i] = 0, P[j][j] = 1;
}

/* i,j of Pivot step */
double get_theta(double **A, int i, int j){
    return (A[j][j]-A[i][i]) / (2 * A[i][j]);
}

int sign(double theta){
    if (theta < 0)
        return -1;
    return 1;
}

double get_t(double theta){
    return ((double) sign(theta)) / (fabs(theta) + sqrt(pow(theta, 2)+1));
}

double get_c_of_t(double t){
    return ((double) 1) / (sqrt(pow(t, 2)+1));
}

double get_c(double **A, int i, int j){
    double theta = get_theta(A, i, j);
    double t = get_t(theta);
    double c = get_c_of_t(t);
    return c;
}

double get_s(double **A, int i, int j){
    double theta = get_theta(A, i, j);
    double t = get_t(theta);
    double c = get_c_of_t(t);
    return (double)c*t;
}

/* STEP 3 PIVOT */
void update_Pivot(int* pivot, double **A, int n){
    double MaxValue = -1; // smaller than any abs val
    int i,j;
    for (i = 0; i < n - 1; i++) {
        for (j = i + 1; j < n; j++) {
            if (fabs(A[i][j]) > MaxValue){
                pivot[0] = i, pivot[1] = j;
                MaxValue = fabs(A[i][j]);
            }
        }
    }
}

double frobenius_Norm_Pow(double **A,int n){
    int i,j;
    double norm = 0;
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            norm += (double)pow(A[i][j],2);
        }
    }
    return norm;
}

double doubleOff(double **A,int n){
    double doubleOffA = 0;
    double sum_Digonal_Pow = 0;
    int i;
    for (i = 0; i < n; i++)
            sum_Digonal_Pow += (double)pow(A[i][i],2);
    doubleOffA = frobenius_Norm_Pow(A,n) - sum_Digonal_Pow;
    return doubleOffA;
}

int converged(double **A,double **ATag,int n) {
    if (doubleOff(A,n) - doubleOff(ATag,n) <= Epsilon)
        return 1;
    return 0;
}

/* transform A -> A' */
double ** transform_A(double **A, int n, int i, int j, double c, double s){
    double **ATag = init_matrix(n,n);
    copy_matrix(A, ATag, n, n);
    double tmp;
    for (int r = 0; r < n; r++){
        if (r == i || r == j)
            continue;
        tmp = c*A[r][i] - s*A[r][j];
        ATag[r][i] = tmp, ATag[i][r] = tmp;
        tmp = c*A[r][j] + s*A[r][i];
        ATag[r][j] = tmp, ATag[j][r] = tmp;
    }
    ATag[i][i] = pow(c,2)*A[i][i] + pow(s,2)*A[j][j]-2*s*c*A[i][j];
    ATag[j][j] = pow(s,2)*A[i][i] + pow(c,2)*A[j][j]+2*s*c*A[i][j];
    ATag[i][j] = 0, ATag[j][i] = 0;
    return ATag;
}

/*
 * STEP 1 - full jacobi
 * check if P_n needs to be calculated https://moodle.tau.ac.il/mod/forum/discuss.php?d=162730
 */
double **jacobi_eigenvectors(double **A, int n) {
    int i, j, ITERATIONS = 100;
    double **V = get_I_matrix(n), **P = get_I_matrix(n), **ATag;
    double c, s;
    int pivot[2];
    int is_converged = 0;
    while (!is_converged && ITERATIONS){
        update_Pivot(pivot, A, n);
        i = pivot[0];
        j = pivot[1];
        update_rotation_matrix(A, n, i, j, P); // step a
        c = P[i][i], s = P[i][j];
        ATag = transform_A(A, n, i, j, c, s); // step b
        is_converged = converged(A, ATag, n);
        free_matrix(A, n);
        A = ATag;
        double **tmp_multiply = V; // step e
        V = multiply_matrices_same_dim(V, P, n);
        revert_P_to_Identity(P, i, j);
        free_matrix(tmp_multiply, n);
        ITERATIONS -= 1;
    }
    free_matrix(P, n);
    double **output = init_matrix(n+1, n); // hack to include eigenvalues in matrix last row
    for (int i = 0; i < n; i++) { //ugly way
        for (int j = 0; j < n; j++) {
            output[i+1][j] = V[i][j];
        }
    }
    free_matrix(V, n);
    for (i = 0; i < n; i++) {
        output[0][i] = ATag[i][i];
    }
    return output;
}


//printf(">%d:\n", 5-ITERATIONS);
//printf("A':\n");
//print_matrix(ATag, n, n);
//
//printf("P:\n");
//print_matrix(V, n, n);
//printf("\n");
