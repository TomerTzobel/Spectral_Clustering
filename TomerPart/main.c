#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "utils.c"
#include "ddg.c"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define FLT_MAX 3.402823e+38
#define ERR_MSG "An Error Has Occured"

double **wam(double **datapoints, int n, int dimension) {
    double **wamMatrix;
    int i, j;
    double pointNorm, result;
    wamMatrix = init_matrix(n, n);
    for (i = 1; i < n - 1; i++) {
        for (j = i + 1; j < n; j++) {
            pointNorm = norm(datapoints[i], datapoints[j], dimension);
            result = (double) pow(M_E, ((-pointNorm) / 2)); //M_E is the mathematical constant e
            wamMatrix[i][j] = result;
            wamMatrix[j][i] = result; // symmetry of W
        }
    }
    return wamMatrix;
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

double **lnorm(double **datapoints, int n, int dimension) {
    double **wamMatrix = wam(datapoints, n, dimension);
    double **DDGMatrix = eval_ddg(datapoints, n);
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
    int pointsNumber = 0;
    int maxlinelen = 0;
    int currlinelen = 0;
    int k, i, j, ch, ClusterNumber;
    FILE *fp;
    int changed = 1;
    char *cordinate;
    char *line;

    assert(argc == 4 && "invalid number of args");
//    assert(isNumber(argv[1]) && "1st arg is not a number");
    k = atoi(argv[1]);

    fp = fopen("C:\\Users\\user\\Desktop\\input_1.txt", "r");
    assert(fp != NULL && "failed to open file");
    while ((ch = fgetc(fp)) != 10) /*check the dimension of the vectors*/
    {
        if (ch == ',')
            dimension++;
    }
    fp = fopen("C:\\Users\\user\\Desktop\\input_1.txt", "r");
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
    assert(points != NULL && "malloc failed");

    for (i = 0; i < pointsNumber; i++) {
        points[i] = malloc(dimension * sizeof(double));
        assert(points[i] != NULL && "malloc failed");
    }
    centroids = malloc(k * sizeof(double *));
    assert(centroids != NULL && "malloc failed");
    for (i = 0; i < k; i++) {
        centroids[i] = malloc(dimension * sizeof(double));
        assert(centroids[i] != NULL && "malloc failed");
    }
    fp = fopen("C:\\Users\\user\\Desktop\\input_1.txt", "r");
    i = 0;
    line = malloc(maxlinelen * sizeof(char));
    while (fgets(line, maxlinelen + 1, fp) != NULL) {
        cordinate = strtok(line, ",");
        for (j = 0; j < dimension; j++) {
            points[i][j] = atof(cordinate);
            if (i < k)
                centroids[i][j] = atof(cordinate);
            cordinate = strtok(NULL, ",");
        }
        i++;
    }
    fp = fopen("C:\\Users\\user\\Desktop\\input_1.txt", "r");
    /*INIT CLUSTERS*/
    assert(NULL != (clusters = calloc(pointsNumber, sizeof(int))) && "calloc failed");
    for (i = 0; i < pointsNumber; i++) {
        if (i < k)
            clusters[i] = i;
        else
            clusters[i] = -1;
    }
    fclose(fp);

    double **ddgmatirx;
    double **lnorm1;
    double **test;
    wamMatrix = wam(points, pointsNumber, dimension);
    ddgmatirx = ddg(points, pointsNumber, dimension);
    lnorm1 = lnorm(points, pointsNumber, dimension);

    print_matrix(lnorm1, pointsNumber, pointsNumber);

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
