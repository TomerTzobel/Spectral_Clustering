#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "utils.h"
#include "ddg.h"
#include "wam.h"
#include "jacobi.h"
#include "eigengap.h"
#include "kmeans.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define FLT_MAX 3.402823e+38
#define ERR_MSG "An Error Has Occured"

//#define LOCAL_PATH "C:\\Users\\adamk\\gitpractice\\spectral\\Spectral_Clustering\\TomerPart\\jacobi_input.txt"
//#define LOCAL_PATH "C:\\Users\\adamk\\gitpractice\\spectral\\Spectral_Clustering\\TomerPart\\test\\input_a_rami.txt"
//#define LOCAL_PATH "C:\\Users\\user\\Desktop\\input_2.txt"


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

int main(int argc, char **argv) {
    double **points;
    int dimension = 1, pointsNumber = 1;
    int maxlinelen = 0, currlinelen = 0;
    int k, i, j, ch;
    FILE *fp;
    char *coordinate, *line, *goal;

    assert(argc == 4 && ERR_MSG);
    k = atoi(argv[1]);
    goal = argv[2];
    fp = fopen(argv[3], "r");
    assert(fp != NULL && ERR_MSG);
    while ((ch = fgetc(fp)) != 10) /*check the dimension of the vectors*/
    {
        if (ch == ',')
            dimension++;
        if (ch == EOF)
            break;
    }

    fp = fopen(argv[3], "r");
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
    fp = fopen(argv[3], "r");
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

    if (strcmp(goal,"wam") == 0){
        double **wamMatrix = wam(points, pointsNumber, dimension);
        print_matrix(wamMatrix,pointsNumber,pointsNumber);
        free_matrix(wamMatrix,pointsNumber);
    }

    if (strcmp(goal,"ddg") == 0){
        double **ddgMatrix = ddg(points, pointsNumber, dimension);
        print_matrix(ddgMatrix, pointsNumber, pointsNumber);
        free_matrix(ddgMatrix, pointsNumber);
    }

    if (strcmp(goal,"lnorm") == 0){
        double **lnormMatrix = lnorm(points, pointsNumber, dimension);
        print_matrix(lnormMatrix, pointsNumber, pointsNumber);
        free_matrix(lnormMatrix, pointsNumber);
    }

    if (strcmp(goal,"jacobi") == 0){
        double **eigenvectors = jacobi_eigenvectors(points,pointsNumber);
        print_matrix(eigenvectors,pointsNumber,pointsNumber);
        free_matrix(eigenvectors,pointsNumber+1);
    }

    if (strcmp(goal,"spk") == 0){
      int n = pointsNumber;
      double **lnorm_matrix = lnorm(points, n, dimension); 
      double **V = jacobi_eigenvectors(lnorm_matrix, pointsNumber);
      double *eigenvalues = V[pointsNumber];
      int *sorted_eigenvectors_indices = bubbleSort_index_tracked(eigenvalues, n);
      if (k == 0){
          k = get_elbow_k(eigenvalues, pointsNumber);
      }
      double **U = copy_columns_by_order(V, n, k, sorted_eigenvectors_indices);
      free_matrix(V, n+1);
      free(sorted_eigenvectors_indices);
      normalize_matrix(U, n, k); // U --> T
      double **centroids = kmeans(k, pointsNumber, U);
      print_matrix(centroids, k, k);
    }


    return 0;
}
