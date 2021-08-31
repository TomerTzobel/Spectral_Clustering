#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include "nsc.h"


static void *fit(PyObject *self, PyObject *args);

static PyMethodDef capiMethods[] = {
        {"fit",
         (PyCFunction)fit,
         METH_VARARGS,
         PyDoc_STR("runs Normalized Spectral Clustering")},
         {NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "nsc_wrapper",
        NULL,
        -1,
        capiMethods};

PyMODINIT_FUNC
PyInit_nsc_wrapper(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m)
    {
        return NULL;
    }
    return m;
};


static void *fit(PyObject *self, PyObject *args)
{
    int k;
    char *goal, *filename;

    if (!PyArg_ParseTuple(args, "iss", &k, &goal, &filename))
        return NULL;

    nsc(k, goal, filename);
}
