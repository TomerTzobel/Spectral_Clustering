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
    double **tmp_multiply, **lnorm;
    int i, j;
    power_matrix_elementwise(DDGMatrix, n, ((double) -1 / (double) 2)); //D^(-1/2)
    tmp_multiply = multiply_matrices_same_dim(DDGMatrix, wamMatrix, n);
    free_matrix(wamMatrix, n);
    lnorm = multiply_matrices_same_dim(tmp_multiply, DDGMatrix, n); //D^(-1/2) * W * D^(-1/2)
    free_matrix(DDGMatrix, n);
    free_matrix(tmp_multiply, n);
    for (i = 0; i < n; i++) { //I - D^(-1/2) * W * D^(-1/2)
        for (j = 0; j < n; j++) {
            if (i == j)
                lnorm[i][j] = 1 - lnorm[i][j];
            else
                lnorm[i][j] = 0 - lnorm[i][j];
        }
    }
    return lnorm;
}

int main(int argc, char **argv) {
    double **points, **centroids, **wamMatrix;
    int *clusters;
    int max_itter = 200;
    int dimension = 1;
    int pointsNumber = 1;
    int maxlinelen = 0;
    int currlinelen = 0;
    int k, i, j, ch, ClusterNumber;
    FILE *fp;
    int changed = 1;
    char *cordinate;
    char *line;
    char *goal;

    assert(argc == 4 && ERR_MSG);
//    assert(isNumber(argv[1]) && "1st arg is not a number");
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
    points = malloc(pointsNumber * sizeof(double *));
    assert(points != NULL && ERR_MSG);

    for (i = 0; i < pointsNumber; i++) {
        points[i] = malloc(dimension * sizeof(double));
        assert(points[i] != NULL && ERR_MSG);
    }
    centroids = malloc(k * sizeof(double *));
    assert(centroids != NULL && ERR_MSG);
    for (i = 0; i < k; i++) {
        centroids[i] = malloc(dimension * sizeof(double));
        assert(centroids[i] != NULL && ERR_MSG);
    }
    fp = fopen(argv[3], "r");
    i = 0;
    line = malloc(maxlinelen * sizeof(char));
    while (fgets(line, maxlinelen + 1, fp) != NULL) {
        cordinate = strtok(line, ",");
        for (j = 0; j < dimension; j++) {
            points[i][j] = atof(cordinate); //problem
            if (i < k)
                centroids[i][j] = atof(cordinate);
            cordinate = strtok(NULL, ",");
        }
        i++;
    }
    fp = fopen(argv[3], "r");
    /*INIT CLUSTERS*/
    assert(NULL != (clusters = calloc(pointsNumber, sizeof(int))) && ERR_MSG);
    for (i = 0; i < pointsNumber; i++) {
        if (i < k)
            clusters[i] = i;
        else
            clusters[i] = -1;
    }
    fclose(fp);

    if (strcmp(goal,"wam") == 0){
        double **wamMatrix = wam(points, pointsNumber, dimension);
        print_matrix(wamMatrix,pointsNumber,pointsNumber);
        free_matrix(wamMatrix,pointsNumber);
    }

    if (strcmp(goal,"ddg") == 0){
        double **ddgMatirx = ddg(points, pointsNumber, dimension);
        print_matrix(ddgMatirx,pointsNumber,pointsNumber);
        free_matrix(ddgMatirx,pointsNumber);
    }

    if (strcmp(goal,"lnorm") == 0){
        double **lnormMatirx = lnorm(points, pointsNumber, dimension);
        print_matrix(lnormMatirx,pointsNumber,pointsNumber);
        free_matrix(lnormMatirx,pointsNumber);
    }

    if (strcmp(goal,"jacobi") == 0){
        double **eigenvectors = jacobi_eigenvectors(points,pointsNumber);
        print_matrix(eigenvectors,pointsNumber,pointsNumber);
        free_matrix(eigenvectors,pointsNumber);
    }

    if (strcmp(goal,"spk") == 0){
      //to fill
    }


    for (i = 0; i < k; i++)
        free(centroids[i]);
    free(centroids);
    for (i = 0; i < pointsNumber; i++)
        free(points[i]);
    free(points);
    free(clusters);
    free(line);

    return 0;
}
