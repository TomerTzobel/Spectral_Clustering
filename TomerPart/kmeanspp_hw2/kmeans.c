#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define FLT_MAX 3.402823e+38

static double **Result(double **points, double **centroids, int max_iter, int k, int dimension, int pointsNumber, long *indexKmeans);
static PyObject *fit(PyObject *self, PyObject *args);

static PyMethodDef capiMethods[] = {
    {"fit",
     (PyCFunction)fit,
     METH_VARARGS,
     PyDoc_STR("Returns the final centroeids after Kmeans algorithm")},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    NULL,
    -1,
    capiMethods};

PyMODINIT_FUNC
PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m)
    {
        return NULL;
    }
    return m;
};

static int findMinCent(double *point, double **centroeids, int k, int dimension)
{
    double min_val = FLT_MAX;
    int min_idx = 0;
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

static int UpdateAllAvg(double **centroids, int *clusters, double **points, int k, int pointsnumber, int dimension)
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
    CountCluster = (int *)malloc(k * sizeof(int));
    assert(CountCluster != NULL && "malloc failed");
    for (j = 0; j < k; j++)
    {
        CountCluster[j] = 0;
    }

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

static double **pythonListToArrays(PyObject *pythonList, int pointsNumber, int dimension)
{
    double **points;
    int i;
    int j;
    PyObject *item;
    points = malloc(pointsNumber * sizeof(double *));
    assert(points != NULL && "malloc failed");
    for (i = 0; i < pointsNumber; i++)
    {
        points[i] = malloc(dimension * sizeof(double));
        assert(points[i] != NULL && "malloc failed");
    }
    for (i = 0; i < pointsNumber; i++)
    {
        for (j = 0; j < dimension; j++)
        {
            item = PyList_GetItem(pythonList, i * dimension + j);
            points[i][j] = PyFloat_AsDouble(item);
        }
    }
    return points;
}

static PyObject *fit(PyObject *self, PyObject *args)
{
    PyObject *pythonPoints;
    PyObject *pythonCentroids;
    PyObject *Kmeans;
    PyObject *indicies;
    PyObject *item;
    double **result;
    int pointsNumber, dimension, max_itter;
    double **points;
    double **centroids;
    long *indexKmeans;
    int i, k, j;

    if (!PyArg_ParseTuple(args, "OOOiiii", &pythonPoints, &pythonCentroids, &indicies, &dimension, &pointsNumber, &k, &max_itter))
        return NULL;

    indexKmeans = malloc(k * sizeof(long));
    assert(indexKmeans != NULL && "malloc failed");
    for (i = 0; i < k; i++)
    {
        item = PyList_GetItem(indicies, i);
        indexKmeans[i] = PyLong_AsLong(item);
    }

    points = pythonListToArrays(pythonPoints, pointsNumber, dimension);
    centroids = pythonListToArrays(pythonCentroids, k, dimension);

    result = Result(points, centroids, max_itter, k, dimension, pointsNumber, indexKmeans);
    Kmeans = PyList_New(k);
    assert(Kmeans != NULL && "PyList_New failed");
    for (i = 0; i < k; i++)
    {
        item = PyList_New(dimension);
        for (j = 0; j < dimension; j++)
        {
            PyList_SET_ITEM(item, j, Py_BuildValue("d", result[i][j]));
        }
        PyList_SET_ITEM(Kmeans, i, Py_BuildValue("O", item));
    }

    for (i = 0; i < k; i++)
        free(centroids[i]);
    free(centroids);
    for (i = 0; i < pointsNumber; i++)
        free(points[i]);
    free(points);
    free(indexKmeans);

    return Kmeans;
}

static double **Result(double **points, double **centroids, int max_iter, int k, int dimension, int pointsNumber, long *indexKmeans)
{
    int *clusters;
    int i, ClusterNumber;
    int changed = 1;
    int index;
    int counter;

    clusters = (int *)malloc(pointsNumber * sizeof(int));
    assert(clusters != NULL && "malloc failed");
    /*INIT CLUSTERS*/
    for (i = 0; i < pointsNumber; i++)
    {
        clusters[i] = -1;
    }
    counter = 0;
    for (i = 0; i < k; i++)
    {
        index = indexKmeans[i];
        clusters[index] = counter;
        counter++;
    }

    while (max_iter > 0 && changed)
    {
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
