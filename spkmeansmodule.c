#include "spkmeans.h"
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>

/* HELPER FUNCTIONS */

static int ArrayTo_PyList(PyObject *list, double *array, int length)
{
    int index;
    PyObject *CURR;
    int STATUS = 0;
    for (index = 0; index < length; index++){
        CURR = PyFloat_FromDouble(array[index]);
        if (PyList_SetItem(list, index, CURR) == -1){
            STATUS = 1;
        }
    }
    return STATUS; /* STATUS equals 0 means success, otherwise, failure */
}

static int PyListTo_Array(PyObject *list, double *array, int length)
{
    int index;
    PyObject *CURR;
    int STATUS = 0;
    for (index = 0; index < length; index++){
        CURR = PyList_GetItem(list, index);
        if (!PyFloat_Check(CURR)){
            STATUS = 1;
        }
        array[index] = PyFloat_AsDouble(CURR);  
    }
    return STATUS; /* STATUS equals 0 means success, otherwise, failure */
}

/* CALLING SPKMEANS.C ALGORITHM */

/* PO_POINTS and POINTO store the points
   T1 and target1 store the eigenvectors or NxN matrix
   T2 and target2 store the eigenvalues (if necessary)
   SPK and SPKtarget store the k largest vectors after normalization (if necessary) */
static PyObject* spkFit(PyObject* self, PyObject* args)
{
    PyObject *PO_POINTS, *T1, *T2, *SPK; 
    double *POINTO, *target1, *target2, *SPKtarget;
    PyObject* Output;
    int K;
    int D;
    int N;
    char *goal;
    enum mode m;

    if (!PyArg_ParseTuple(args, "siiiO", &goal, &K, &D, &N, &PO_POINTS))
    {
        return NULL;
    }

    POINTO = malloc(N * D * sizeof(double));

    if (NULL == POINTO) /*if allocation failed*/
        return NULL;

    if(PyListTo_Array(PO_POINTS, POINTO, N * D) != 0)
        return NULL;

    target1 = malloc(N * N *sizeof(double));
    target2 = malloc(N * sizeof(double));
    if (NULL == target1 || NULL == target2)
    {
        return NULL;
    }

    Main_args(goal, N, D, &K, POINTO, target1, target2); /* the algorithm */

    SPKtarget = malloc(N * K * sizeof(double));
    if (NULL == SPKtarget)
    {
        return NULL;
    }
    
    m = translation(goal);
    if(SPKmode == m)
    {
        spk_process2(target1, K, N, SPKtarget);
    }

    T1 = PyList_New(N * N);
    T2 = PyList_New(N);
    SPK = PyList_New(N * K);
    if (NULL == T1 || NULL == T2 || NULL == SPK){ /*if list creation failed*/
        return NULL;
    }

    ArrayTo_PyList(T1, target1, N * N);
    ArrayTo_PyList(T2, target2, N);
    ArrayTo_PyList(SPK, SPKtarget, N * K);

    free(POINTO);
    free(target1);
    free(target2);
    free(SPKtarget);
    
    Output = Py_BuildValue("OOOi", T1, T2, SPK, K);
    return Output;
}

/* C TO PYTHON MODULE HANDLING */

static PyMethodDef mymetoda[] = {
    {"spkFit", spkFit, METH_VARARGS, "All of SPK's processes"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef myspkmodule = {
    PyModuleDef_HEAD_INIT,
    "myspkmodule",
    "I hope it works",
    -1,
    mymetoda
};

PyMODINIT_FUNC PyInit_myspkmodule(void)
{
    PyObject *m;
    m = PyModule_Create(&myspkmodule);
    if (!m){
        return NULL;
    }
    return m;
}