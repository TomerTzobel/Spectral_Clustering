//
// Created by adamk on 24/08/2021.
//
#include <stdio.h> //delete this later, only for priting

#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "utils.h"
#define Epsilon 1.0e-15

double get_c(double **A, int i, int j);
double get_s(double **A, int i, int j);


/* STEP 2: n is dim, i,j are A_ij Pivot step */
double **get_rotation_matrix(double **A, int n, int i, int j) {
    double **P = init_matrix(n, n);
    for (int k=0; k<n; k++){
        P[k][k] = 1;
    }
    double c = get_c(A, i, j);
    double s = get_s(A,i,j);
    //CHECK THIS PART, UNCLEAR
    P[i][i] = c;
    P[i][j] = s;
    P[j][i] = -s;
    P[j][j] = c;
    printf("----\n");
    print_matrix(P,n,n);
    printf("----\n");
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
    //return sqrt(1 - pow(c, 2));
}

/* STEP 3 PIVOT */
void update_Pivot(int* pivot, double **A, int n){
    double MaxValue;
    int i,j;
    pivot[0] = 0;
    pivot[1] = 1;
    MaxValue = fabs(A[0][1]);
    int replace = 0;
    for (i = 0; i < n ; i++) {
        for (j = i; j < n; j++) {
            if (i != j){
                replace = 0;
                if (fabs(A[i][j]) > MaxValue)
                    replace = 1;
                if (fabs(A[i][j]) == MaxValue){
                    if (i < pivot[0])
                        replace = 1;
                    if (i == pivot[0] && j < pivot[1])
                        replace = 1;
                }
                if (replace) {
                    pivot[0] = i;
                    pivot[1] = j;
                    MaxValue = A[i][j];
                }
            }
        }
    }

};
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
    doubleOffA = (double)pow(doubleOffA,2);
    return doubleOffA;
}
/*
 *
 * STEP 5
 *
 */
int converged(double **A,double **ATag,int n) {
    if (doubleOff(A,n) - doubleOff(ATag,n) <= Epsilon)
        return 1;
    return 0;
}

/* transform A -> A' */
double ** transform_A(double **A, int n, int i, int j, double c, double s){
    double **ATag;
    int r;
    ATag = init_matrix(n,n);
    double tmp_Ari, tmp_Arj;
    for (r = 0; r < n; r++){ // UNSURE THIS WORKS <<
        if (r == i || r == j)
            continue;
        tmp_Ari = c*A[r][i] - s*A[r][j];
        tmp_Arj = c*A[r][j] + s*A[r][i];
        ATag[r][i] = tmp_Ari, ATag[i][r] = tmp_Ari;
        ATag[r][j] = tmp_Arj, ATag[j][r] = tmp_Arj;
    }
    double tmp_Aii, tmp_Ajj;
    tmp_Aii = pow(c,2)*A[i][i] + pow(s,2)*A[j][j]-2*s*c*A[i][j];
    tmp_Ajj = pow(s,2)*A[i][i] + pow(c,2)*A[j][j]+2*s*c*A[i][j];
//    ATag[i][j] = 0, ATag[j][i] = 0;
    ATag[i][i] = tmp_Aii, ATag[j][j] = tmp_Ajj;
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
//        printf("%d\n",i);
//        printf("%d\n",j);

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
        eigenvalues[i] =  ATag[i][i];
    }
    print_matrix(ATag,n,n);
    printf("-----\n");
    //delete later
    for (i = 0; i < n; i++) {
            printf("%.4f", eigenvalues[i]);
            if (j < n - 1)
                printf(",");
        }
        printf("\n");

return eigenvalues;
}

