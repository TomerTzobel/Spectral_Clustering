#define PY_SSIZE_T_CLEAN
#define ERR_MSG "An Error Has Occured"

#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "spkmeans.h"

static void fit(PyObject *self, PyObject *args);
static void kmeans_pp(PyObject *self, PyObject *args);
static PyObject *get_normalized_matrix_wrapper(PyObject *self, PyObject *args);
static double **pythonListToArrays(PyObject *pythonList, int pointsNumber, int dimension);

static PyMethodDef capiMethods[] = {
        {"get_normalized_matrix_wrapper", (PyCFunction)get_normalized_matrix_wrapper, METH_VARARGS, PyDoc_STR("returns normalized matrix")},
        {"kmeans_pp", (PyCFunction)kmeans_pp, METH_VARARGS, PyDoc_STR("run kmeans++")},
        {"fit", (PyCFunction)fit, METH_VARARGS, PyDoc_STR("runs any goal besides spk")},
        {NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "spkmeansmodule",
        NULL,
        -1,
        capiMethods};

PyMODINIT_FUNC
PyInit_spkmeansmodule(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m)
    {
        return NULL;
    }
    return m;
};

static void kmeans_pp(PyObject *self, PyObject *args){
    PyObject *pythonPoints , *python_centroids_indices, *item;
    int pointsNumber, dimension;
    double **points, **centroids;
    long *centroids_indices;
    int i, k;

    if (!PyArg_ParseTuple(args, "OOiii", &pythonPoints, &python_centroids_indices, &dimension, &pointsNumber, &k))
        assert(0 && ERR_MSG);

    centroids_indices = malloc(k * sizeof(long));
    assert(centroids_indices != NULL && ERR_MSG);
    for (i = 0; i < k; i++)
    {
        item = PyList_GetItem(python_centroids_indices, i);
        centroids_indices[i] = PyLong_AsLong(item);
    }

    points = pythonListToArrays(pythonPoints, pointsNumber, dimension);
    centroids = kmeans(k, pointsNumber, points, centroids_indices);
    free_matrix(points, pointsNumber);

    print_matrix(centroids, k, k);
    free_matrix(centroids, k);
}

static void fit(PyObject *self, PyObject *args)
{
    char *goal, *filename;
    if (!PyArg_ParseTuple(args, "ss", &goal, &filename))
        assert(0 && ERR_MSG);
    nsc(0, goal, filename); /* dummy value */
}

static PyObject *get_normalized_matrix_wrapper(PyObject *self, PyObject *args) {
    char *filename;
    double **T, **points;
    int k, i, j, dimension = 1, n = 1;
    PyObject *item, *normalized_eigenvectors;
    assert(flatten_T != NULL && ERR_MSG);

    if (!PyArg_ParseTuple(args, "is", &k, &filename))
        return NULL;

    read_data(filename, &points, &dimension, &n);
    T = get_normalized_eigenvectors(&k, points, dimension, n);
    normalized_eigenvectors = PyList_New(n);

    for (i = 0; i < n; i++) {
        item = PyList_New(k);
        assert(item  != NULL && ERR_MSG);
        for (j = 0; j < k; j++) {
            PyList_SET_ITEM(item, j, Py_BuildValue("d", T[i][j]));
        }
        PyList_SET_ITEM(normalized_eigenvectors, i, Py_BuildValue("O", item));
    }
    free_matrix(T, n);
    return normalized_eigenvectors;
}

static double **pythonListToArrays(PyObject *pythonList, int pointsNumber, int dimension)
{
    double **points;
    int i, j;
    PyObject *item;
    points = init_matrix(pointsNumber, dimension);
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
