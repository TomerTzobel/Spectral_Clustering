#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "ddg.h"
#include "wam.h"
#include "jacobi.h"
#include "eigengap.h"
#include "kmeans.h"
#include "lnorm.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define ERR_MSG "An Error Has Occured"

void nsc(int k, char *goal, char *filename, int kmeans_pp) {
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
