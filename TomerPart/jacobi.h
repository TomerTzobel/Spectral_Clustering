//
// Created by adamk on 24/08/2021.
//

#ifndef MAIN_C_JACOBI_H
#define MAIN_C_JACOBI_H

double get_c(double **A, int i, int j);
double get_s(double c);
double **get_rotation_matrix(double **A, int n, int i, int j);
double get_theta(double **A, int i, int j);
int sign(double theta);
double get_t(double theta);
double get_c_of_t(double t);
void update_Pivot(int* pivot, double **A, int n);
double frobenius_Norm_Pow(double **A,int n);
double doubleOff(double **A,int n);
int converged(double **A,double **ATag,int n);
double ** transform_A(double **A, int n, int i, int j, double c, double s);
double *jacobi_eigenvectors(double **A, int n);


#endif //MAIN_C_JACOBI_H
