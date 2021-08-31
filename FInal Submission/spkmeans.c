#define ERR_MSG "An Error Has Occured"
#define Epsilon 1.0e-15
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define FLT_MAX 3.402823e+38
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "spkmeans.h"

/* The Weighted Adjacency Matrix Functions */
double **wam(double **datapoints, int n, int dimension) {
    double **wamMatrix;
    int i, j;
    double pointNorm, result;
    wamMatrix = init_matrix(n, n);
    for (i = 0; i < n - 1; i++) {
        for (j = i + 1; j < n; j++) {
            pointNorm = norm(datapoints[i], datapoints[j], dimension);
            result = (double) pow(M_E, ((-pointNorm) / 2)); //M_E is the mathematical constant e
            wamMatrix[i][j] = result;
            wamMatrix[j][i] = result; // symmetry of W
        }
    }
    free_matrix(datapoints, n);
    return wamMatrix;
}
/********************/

/* The Diagonal Degree Matrix Functions */
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
/********************/


/* The Normalized Graph Laplacian Functions */
double **lnorm(double **datapoints, int n, int dimension) {
    double **wamMatrix = wam(datapoints, n, dimension);
    double **DDGMatrix = eval_ddg(wamMatrix, n);
    double **tmp_multiply, **lnorm_matrix;
    int i, j;
    power_matrix_elementwise(DDGMatrix, n, ((double) -1 / (double) 2)); //D^(-1/2)
    tmp_multiply = multiply_matrices_same_dim(DDGMatrix, wamMatrix, n);
    free_matrix(wamMatrix, n);
    lnorm_matrix = multiply_matrices_same_dim(tmp_multiply, DDGMatrix, n); //D^(-1/2) * W * D^(-1/2)
    free_matrix(DDGMatrix, n);
    free_matrix(tmp_multiply, n);
    for (i = 0; i < n; i++) { //I - D^(-1/2) * W * D^(-1/2)
        for (j = 0; j < n; j++) {
            if (i == j)
                lnorm_matrix[i][j] = 1 - lnorm_matrix[i][j];
            else
                lnorm_matrix[i][j] = 0 - lnorm_matrix[i][j];
        }
    }
    return lnorm_matrix;
}
/********************/


/*
 * Jacobi algorithm Functions
 */
double get_c(double **A, int i, int j); //do we need this here?
double get_s(double **A, int i, int j); //do we need this here?


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
    int r;
    for (r = 0; r < n; r++){
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
 * Finding Eigenvalues and Eigenvectors
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
        i = pivot[0]; // ##################################### check ayelets question #################################################################
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
    double **output = init_matrix(n+1, n);
    copy_matrix(V, output+1, n, n);
    free_matrix(V, n);
    for (i = 0; i < n; i++) { // first line with eigenvalues
        output[0][i] = ATag[i][i];
    }
    return output;
}
/********************/

/* find ideal k given sorted eigenvalues*/
int get_elbow_k(double *eigenvalues, int n) {
    int i, k;
    double max = -1;
    double *gaps = calloc(n/2 + 1, sizeof (double)); // extra element for easy indexing
    for (i = 1; i < n/2 + 1 ; i++) {
        gaps[i] = fabs(eigenvalues[i] - eigenvalues[i+1]);
        if (gaps[i] > max) {
            k = i;
            max = gaps[i];
        }
    }
    free(gaps);
    return k;
}
/********************/

/* Kmeans algorithm functions
 * Based on H.W.1
 */
int findMinCent(double *point, double **centroids, int k, int dimension)
{
    double min_val = FLT_MAX;
    int min_idx;
    double curr_val;
    int i, j;
    double num;
    for (i = 0; i < k; i++)
    {
        curr_val = 0.0;
        for (j = 0; j < dimension; j++)
        {
            num = point[j] - centroids[i][j];
            curr_val = curr_val + (num * num);
        }
        if (curr_val < min_val)
        {
            min_val = curr_val;
            min_idx = i;
        }
    }
    return min_idx;
}

int UpdateAllAvg(double **centroids, int *clusters, double **points, int k, int pointsnumber, int dimension)
{
    int i, j, changed = 0;
    int curr;
    double **OriginCenter;
    int *CountCluster;

    OriginCenter = (double **)malloc(k * sizeof(double *));
    assert(OriginCenter != NULL && "malloc failed");
    for (i = 0; i < k; i++)
    {
        OriginCenter[i] = (double *)malloc(dimension * sizeof(double));
        assert(OriginCenter[i] != NULL);
        for (j = 0; j < dimension; j++)
            OriginCenter[i][j] = centroids[i][j];
    }

    assert(NULL != (CountCluster = calloc(k, sizeof(int))) && "calloc failed");
    /*RESET CENTROIDS*/
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < dimension; j++)
        {
            centroids[i][j] = 0.0;
        }
    }

    /*SUM ALL CORDINATES FOR EACH CLUSTERS VECTOR*/
    for (i = 0; i < pointsnumber; i++)
    {
        curr = clusters[i];
        CountCluster[curr]++;
        for (j = 0; j < dimension; j++)
        {
            centroids[curr][j] = centroids[curr][j] + points[i][j];
        }
    }

    /*CALCULATE AVG AND CHECK IF ANY CENTROID IS CHANGED*/
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < dimension; j++)
        {
            centroids[i][j] = (double)centroids[i][j] / (CountCluster[i]);
            if (centroids[i][j] != OriginCenter[i][j])
            {
                changed = 1;
            }
        }
    }

    for (i = 0; i < k; i++)
        free(OriginCenter[i]);
    free(OriginCenter);
    free(CountCluster);
    return changed;
}

double **kmeans(int k, int n, double **points)
{
    int *clusters;
    int max_iter = 300, dimension = k, pointsNumber = n, changed = 1;
    int i, j, ClusterNumber;


    double **centroids = init_matrix(k, k);
    for (i = 0; i < k; i++){
        for (j = 0; j < k; j++){
            centroids[i][j] = points[i][j];
        }
    }

    /*INIT CLUSTERS*/
    assert(NULL != (clusters = calloc(pointsNumber, sizeof(int))) && "calloc failed");
    for (i = 0; i < pointsNumber; i++)
    {
        if (i < k)
            clusters[i] = i;
        else
            clusters[i] = -1;
    }

    while (max_iter > 0 && changed)
    {
        changed = 0;
        max_iter--;
        for (i = 0; i < pointsNumber; i++)
        {
            ClusterNumber = findMinCent(points[i], centroids, k, dimension);
            clusters[i] = ClusterNumber;
        }
        changed = UpdateAllAvg(centroids, clusters, points, k, pointsNumber, dimension);
    }

    free(clusters);

    return centroids;
}
/********************/


/*
 * Helping methods
 */

/* init matrix of zeros */
double **init_matrix(int rows, int cols) {
    double **matrix = malloc(rows * sizeof(double *));
    assert(matrix != NULL && ERR_MSG);
    int i;
    for (i = 0; i < rows; i++) {
        matrix[i] = calloc(cols, sizeof(double));
        assert(matrix[i] != NULL && ERR_MSG);
    }
    return matrix;
}

void copy_matrix(double **source, double **dest, int rows, int cols) {
    int i,j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            dest[i][j] = source[i][j];
        }
    }
}

double **copy_columns_by_order(double **source, int rows, int cols, int* order){
    double **cpy = init_matrix(rows, cols);
    int curr_col;
    int j,i;
    for (j = 0; j < cols; j++) {
        curr_col = order[j];
        for (i = 0; i < rows; i++) {
            cpy[i][j] = source[i+1][curr_col]; //source has n+1 rows (eigenvaleus on the top)
        }
    }
    return cpy;
}

void normalize_matrix(double **matrix, int rows, int cols){
    //double **output = init_matrix(rows,cols);
    double sum_squared, denominator;
    int i,j;
    for (i = 0; i < rows; i++) {
        sum_squared = 0;
        for (j = 0; j < cols; j++) {
            sum_squared += (double)pow(matrix[i][j], 2);
        }
        denominator = (double)sqrt(sum_squared);
        for (j = 0; j < cols; j++) {
            matrix[i][j] = matrix[i][j] / denominator;
        }
    }
    //print_matrix(output,rows,cols);
}

double **transpose_matrix(double **matrix, int rows, int cols) {
    int i,j;
    double **transpose = init_matrix(cols, rows);
    for (i = 0; i < rows; ++i)
        for (j = 0; j < cols; ++j) {
            transpose[j][i] = matrix[i][j];
        }
    return transpose;
}

void print_matrix(double **matrix, int rows, int cols) {
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            if (matrix[i][j] < 0 && matrix[i][j] > -0.00005){
                printf("0.0000");
            }
            else {
                printf("%.4f", matrix[i][j]);
            }
            if (j < cols - 1)
                printf(",");
        }
        if (i != rows -1)
            printf("\n");
    }
}

void print_arr(double *arr, int n) {
    int i;
    for (i = 0; i < n; i++) {
        if ( i != n - 1)
            printf("%.4f,", arr[i]);
        else
            printf("%.4f", arr[i]);
    }
    printf("\n");
}

void free_matrix(double **matrix, int rows) {
    int i;
    for (i = 0; i < rows; i++)
        free(matrix[i]);
    free(matrix);
}

double **multiply_matrices(double **A, int rows_A, int cols_A, double **B, int rows_B, int cols_B) {
    assert(cols_A == rows_B && "invalid multiplication");
    double **res = init_matrix(rows_A, cols_B);
    int i,j,k;
    for (i = 0; i < rows_A; ++i) {
        for (j = 0; j < cols_B; ++j) {
            for (k = 0; k < cols_A; ++k) {
                res[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return res;
}

double **multiply_matrices_same_dim(double **A, double **B, int n){
    return multiply_matrices(A, n, n, B, n, n);
};

/* calculate norm of 2 vectors */
double norm(double *point1, double *point2, int dimension) {
    int i;
    double result = 0;
    double norm;
    for (i = 0; i < dimension; i++) {
        result += (double) pow((point1[i] - point2[i]), 2);
    }
    norm = (double) sqrt(result);
    return norm;
}

/* code by GFG, with our minor optimization */
void swap(double *xp, double *yp) {
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}
void swap_int(int *xp, int *yp) {
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}

/* code by GFG, with our minor optimization */
int* bubbleSort_index_tracked(double* arr, int n) {
    int i, j, swapped;
    int* indices = malloc(n * sizeof(int));
    assert(indices != NULL && ERR_MSG);
    for (i = 0; i < n; i++) {
        indices[i] = i;
    }
    for (i = 0; i < n - 1; i++) {
        swapped = 0;
        // Last i elements are already in place
        for (j = 0; j < n - i - 1; j++) {
            if (arr[j] > arr[j + 1]) {
                swap(&arr[j], &arr[j + 1]);
                swap_int(&indices[j], &indices[j + 1]);
                swapped = 1;
            }
        }
        if (!swapped) {
            break;
        }
    }
    return indices;
}

/* inplace power matrix elementwise by x */
void power_matrix_elementwise(double **matrix, int n, double x) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j)
                matrix[i][j] = (double) (pow(matrix[i][j], x));
        }
    }
}

double **get_I_matrix(int n) {
    int i;
    double **I = init_matrix(n, n);
    for (i = 0; i < n; i++) {
        I[i][i] = 1;
    }
    return I;
}
/********************/

/* Reads points from the file
 * Output results depending on the Goal from the user
*/
void nsc(int k, char *goal, char *filename) {
    double **points;
    int dimension = 1, pointsNumber = 1;
    int maxlinelen = 0, currlinelen = 0;
    int i, j, ch;
    FILE *fp;
    char *coordinate, *line;

    fp = fopen(filename, "r");
    assert(fp != NULL && ERR_MSG);
    while ((ch = fgetc(fp)) != 10) /*check the dimension of the vectors*/
    {
        if (ch == ',')
            dimension++;
        if (ch == EOF)
            break;
    }

    fp = fopen(filename, "r");
    while ((ch = fgetc(fp)) != EOF) /*check the number of the vectors*/
    {
        currlinelen++;
        if (ch == 10) {
            pointsNumber++;
            maxlinelen = MAX(maxlinelen, currlinelen);
            currlinelen = 0;
        }
    }

    /*malloc points matrix and read points */
    points = init_matrix(pointsNumber, dimension);
    fp = fopen(filename, "r");
    i = 0;
    line = malloc(maxlinelen * sizeof(char));
    while (fgets(line, maxlinelen + 1, fp) != NULL) {
        coordinate = strtok(line, ",");
        for (j = 0; j < dimension; j++) {
            points[i][j] = atof(coordinate);
            coordinate = strtok(NULL, ",");
        }
        i++;
    }
    fclose(fp);
    free(line);
    free(coordinate);

    if (strcmp(goal, "wam") == 0) {
        double **wamMatrix = wam(points, pointsNumber, dimension);
        print_matrix(wamMatrix, pointsNumber, pointsNumber);
        free_matrix(wamMatrix, pointsNumber);
    }

    if (strcmp(goal, "ddg") == 0) {
        double **ddgMatrix = ddg(points, pointsNumber, dimension);
        print_matrix(ddgMatrix, pointsNumber, pointsNumber);
        free_matrix(ddgMatrix, pointsNumber);
    }

    if (strcmp(goal, "lnorm") == 0) {
        double **lnormMatrix = lnorm(points, pointsNumber, dimension);
        print_matrix(lnormMatrix, pointsNumber, pointsNumber);
        free_matrix(lnormMatrix, pointsNumber);
    }

    if (strcmp(goal, "jacobi") == 0) {
        points = lnorm(points, pointsNumber, dimension); // remove later $$$
        double **eigenvectors = jacobi_eigenvectors(points, pointsNumber);
        print_matrix(eigenvectors, pointsNumber + 1, pointsNumber);
        free_matrix(eigenvectors, pointsNumber + 1);
    }

    if (strcmp(goal, "spk") == 0) {
        int n = pointsNumber;
        double **lnorm_matrix = lnorm(points, n, dimension);
        double **V = jacobi_eigenvectors(lnorm_matrix, n);
        double *eigenvalues = V[0];
        int *sorted_eigenvectors_indices = bubbleSort_index_tracked(eigenvalues, n);
        if (k == 0) { //we still need to test this
            k = get_elbow_k(eigenvalues, n);
        }
        double **U = copy_columns_by_order(V, n, k, sorted_eigenvectors_indices);
        free_matrix(V, n + 1);
        free(sorted_eigenvectors_indices);
        normalize_matrix(U, n, k); // U --> T
        double **centroids = kmeans(k, n, U);
        print_matrix(centroids, k, k);
        free_matrix(U, n);
        free_matrix(centroids, k);
    }
}
/********************/


/*
 * Reading user CMD arguments
 */
int main(int argc, char **argv) {
    int k;
    char *filename, *goal;
    assert(argc == 4 && ERR_MSG); //do we need this?
    k = atoi(argv[1]);
    goal = argv[2];
    filename = argv[3];
    nsc(k, goal, filename); //Output results depending on the Goal from the user
    //nsc(k, goal, filename, 0); //why zero in the end?
    return 0;
}
/********************/