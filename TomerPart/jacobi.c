//
// Created by adamk on 24/08/2021.
//

#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "utils.h"
#define Epsilon 1.0e-15

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

/* STEP 3 PIVOT */
void update_Pivot(int* pivot, double **A, int n){
    double MaxValue;
    pivots[0] = 0;
    pivots[1] = 0;
    MaxValue = A[0][0];
    int replace = 0;
    for (i = 0; i < n ; i++) {
        for (j = i; j < n; j++) {
            replace = 0;
            if (A[i][j] > MaxValue)
                replace = 1;
            if (A[i][j] == MaxValue){
                if (i < pivots[0])
                    replace = 1
                if (i == pivots[0] && j < pivots[1])
                    replace = 1;
            }
            if (replace){
                pivots[0] = i;
                pivots[1] = j;
                MaxValue = A[i][j];
            }
        }
    }

};
double frobenius_Norm_Pow(A,n){
    int i,j;
    double norm = 0;
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            norm += (double)pow(A[i][j],2);
        }
    }
    return norm;
}


double doubleOff(A,n){
    doubleOffA = 0;
    double sumDigonal = 0;
    for (i = 0; i < n; i++)
            sum_Digonal_Pow += (double)pow(A[i][i],2);
    doubleOffA = frobenius_Norm_Pow(A,n) - sum_Digonal_Pow;
    doubleOffA = (double)pow(doubleOffA,2);
    return doubleOffA;
}
/*
 *
 * STEP 5
 *
 */
int converged(A,ATag,n) {
    if (doubleOff(A,n) - doubleOff(ATag,n) <= Epsilon)
        return 1;
    return 0;
}

/* transform A -> A' */
double ** transform_A(double **A, int n, int i, int j, double c, double s){
    double **ATag
    Atag = init_matrix(n);
    double tmp_Ari, tmp_Arj;
    for (int r=0; r<n; r++){ // UNSURE THIS WORKS <<
        if (r == i || r == j)
            continue;
        tmp_Ari = c*A[r][i] - s*A[r][j];
        tmp_Arj = c*A[r][j] + s*A[r][i];
        Atag[r][i] = tmp_Ari, Atag[i][r] = tmp_Ari;
        Atag[r][j] = tmp_Arj, Atag[j][r] = tmp_Arj;
    }
    double tmp_Aii, tmp_Ajj;
    tmp_Aii = pow(c,2)*A[i][i] + pow(s,2)*A[j][j]-2*s*c*A[i][j];
    tmp_Ajj = pow(s,2)*A[i][i] + pow(c,2)*A[j][j]+2*s*c*A[i][j];
    Atag[i][j] = 0, Atag[j][i] = 0;
    Atag[i][i] = tmp_Aii, Atag[j][j] = tmp_Ajj;
    return ATag;
}

/*
 * STEP 1 - full jacobi
 * check if P_n needs to be calculated https://moodle.tau.ac.il/mod/forum/discuss.php?d=162730
 */
double *jacobi_eigenvalues(double **A, int n){
    double *eigenvalues = calloc(n, sizeof(double));
    assert(eigenvalues != NULL && "calloc failed");
    int i, j, ITERATIONS = 100;
    double **V = get_I_matrix(n), **ATag;
    double c, s;
    int pivot[2];
    int is_converged = 0;
    while(!is_converged && ITERATIONS){

//         * STEP 3 PIVOT - check this https://moodle.tau.ac.il/mod/forum/discuss.php?d=162299 and this https://moodle.tau.ac.il/mod/forum/discuss.php?d=159570
//         * need i and j in scope of loop

        update_Pivot(pivot,A,n);
        i = pivot[0];
        j = pivot[1];
        double **P = get_rotation_matrix(A, n, i, j); // step a
        c = P[i][i], s = P[i][j];
        ATag = transform_A(A, n, i, j, c, s); // step b
        is_converged = converged(A,ATag,n);
        A = ATag;

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

