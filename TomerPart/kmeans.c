#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define FLT_MAX 3.402823e+38

int findMinCent(double *point, double **centroeids, int k, int dimension)
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
            num = point[j] - centroeids[i][j];
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
    int max_iter = 300;
    int dimension = k;
    int pointsNumber = n;
    int k, i, j, ClusterNumber;
    int changed = 1;


    double **centroids = malloc(k * sizeof(double *));
    assert(centroids != NULL && "malloc failed");
    for (i = 0; i < k; i++)
    {
        centroids[i] = malloc(dimension * sizeof(double));
        assert(centroids[i] != NULL && "malloc failed");
    }

    for (i = 0; i < k; i++){
        for (j = 0; j < k; j++){
            centroids[i][j] = points[i][j]
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

    while (max_itter > 0 && changed)
    {
        changed = 0;
        max_itter--;
        for (i = 0; i < pointsNumber; i++)
        {
            ClusterNumber = findMinCent(points[i], centroids, k, dimension);
            clusters[i] = ClusterNumber;
        }
        changed = UpdateAllAvg(centroids, clusters, points, k, pointsNumber, dimension);
    }

    for (i = 0; i < k; i++)
        free(centroids[i]);
    free(centroids);
    free(clusters);

    return 0;
}

int isNumber(char *stringNum)
{
    int length, i;
    length = strlen(stringNum);
    for (i = 0; i < length; i++)
        if (!isdigit(stringNum[i]))
        {
            return 0;
        }
    return 1;
}
