//
// Created by adamk on 24/08/2021.
//

#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "utils.h"

double get_c(double **A, int i, int j);
double get_s(double c);


/* STEP 2: n is dim, i,j are A_ij Pivot step */
double **get_rotation_matrix(double **A, int n, int i, int j) {
    double **P = init_matrix(n, n);
    for (int k=0; k<n; k++){
        P[k][k] = 1;
    }
    double c = get_c(A, i, j);
    double s = get_s(c);
    //CHECK THIS PART, UNCLEAR
    P[i][i] = c;
    P[i][j] = s;
    P[j][i] = -s;
    P[j][j] = c;
    return P;
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
    return ((double) 1) / sqrt(pow(t, 2)+1);
}

double get_c(double **A, int i, int j){
    double theta = get_theta(A, i, j);
    double t = get_t(theta);
    double c = get_c_of_t(t);
    return c;
}

double get_s(double c){
    return sqrt(1 - pow(c, 2));
}

/*
 *
 * STEP 5
 *
 */
int is_converged() {return 0;};



/* transform A -> A' */
void transform_A(double **A, int n, int i, int j, double c, double s){
    double tmp_Ari, tmp_Arj;
    for (int r=0; r<n; r++){ // UNSURE THIS WORKS <<
        if (r == i || r == j)
            continue;
        tmp_Ari = c*A[r][i] - s*A[r][j];
        tmp_Arj = c*A[r][j] + s*A[r][i];
        A[r][i] = tmp_Ari, A[i][r] = tmp_Ari;
        A[r][j] = tmp_Arj, A[j][r] = tmp_Arj;
    }
    double tmp_Aii, tmp_Ajj;
    tmp_Aii = pow(c,2)*A[i][i] + pow(s,2)*A[j][j]-2*s*c*A[i][j];
    tmp_Ajj = pow(s,2)*A[i][i] + pow(c,2)*A[j][j]+2*s*c*A[i][j];
    A[i][j] = 0, A[j][i] = 0;
    A[i][i] = tmp_Aii, A[j][j] = tmp_Ajj;
}

/*
 * STEP 1 - full jacobi
 * check if P_n needs to be calculated https://moodle.tau.ac.il/mod/forum/discuss.php?d=162730
 */
double *jacobi_eigenvalues(double **A, int n){
    double *eigenvalues = calloc(n, sizeof(double));
    assert(eigenvalues != NULL && "calloc failed");
    int i, j, ITERATIONS = 100;
    double **V = get_I_matrix(n);
    double c, s;
    while((!is_converged) && ITERATIONS){
        /*
         * STEP 3 PIVOT - check this https://moodle.tau.ac.il/mod/forum/discuss.php?d=162299 and this https://moodle.tau.ac.il/mod/forum/discuss.php?d=159570
         * need i and j in scope of loop
         */

        double **P = get_rotation_matrix(A, n, i, j); // step a
        c = P[i][i], s = P[i][j];
        transform_A(A, n, i, j, c, s); // step b

        double **tmp_multiply = V; // step e
        V = multiply_matrices_same_dim(V, P, n);
        free_matrix(P, n);
        free_matrix(tmp_multiply, n);

        ITERATIONS -= 1;
    }
    for (i=0; i<n; i++){
        eigenvalues[i] = V[i][i];
    }
    return eigenvalues;
}

