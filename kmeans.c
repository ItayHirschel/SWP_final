#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Python.h>

/* FUNCTION DECLARATIONS */

void raiseError();
void raiseInputError();
double dist_sq(double *x, double *y, int VecDim);
int find_near_centroid(double *centroids, double *vec, int VecDim, int K);
double compConv(double *oldCent, double *newCent, int VecDim, int K);
void AddVecToCluster(double *Vec,double *sum_clusters, int *size_clusters, int cluster_index, int VecDim);
void newCentroids(double *sum_cluster, int *size_cluster,int VecDim, int K);
void setZeroArrayDouble(double *Mat, int arr_len);
void setZeroArrayInt(int *Mat, int arr_len);
void copyArr(double *source, double *target, int length);
void swapMatrices(double **A, double **B);


/* THE ALGORITHM USED */

int MAIN(int MAX_ITER, double *CENTO, int D, int N, double *POINTS, double EPS, int K)
{
    int MatDim[2] = {0,0};
    int iteration = 0;
    double *SUM_CLUSTER;
    int *SIZE_CLUSTER;
    double convy = 1000;
    int point_ind = 0;
    int closest_Cent = -1;
    double *VEC;
    double *CENTROIDS;

    MatDim[0] = D;
    MatDim[1] = N;


    
    SUM_CLUSTER = malloc(MatDim[0] * K * sizeof(double)); /* sum of the vectors assigned to each cluster */
    SIZE_CLUSTER = malloc(K * sizeof(double)); /* size of each cluster */
    CENTROIDS = malloc(sizeof(double) * MatDim[0] * K); 
    copyArr(CENTO, CENTROIDS, MatDim[0] * K);

    if (SUM_CLUSTER == NULL || SIZE_CLUSTER == NULL || CENTROIDS == NULL){ /* if allocation failed, raise error */
        raiseError();
        
    }

    while (convy >= EPS && iteration < MAX_ITER){ 
        setZeroArrayDouble(SUM_CLUSTER, MatDim[0] * K); /* resets SUM_CLUSTER array */
        setZeroArrayInt(SIZE_CLUSTER, K); /* resets SIZE_CLUSTER array */
        for (point_ind = 0; point_ind < MatDim[1]; point_ind++){ /* loop which assigns each vector to the appropriate cluster */
            VEC = POINTS + point_ind * MatDim[0]; 
            closest_Cent = find_near_centroid(CENTROIDS, VEC, MatDim[0], K);
            AddVecToCluster(VEC, SUM_CLUSTER, SIZE_CLUSTER, closest_Cent, MatDim[0]);
        }
        newCentroids(SUM_CLUSTER, SIZE_CLUSTER, MatDim[0], K); /* calculates new centroids and writes them to SUM_CLUSTER */
        convy = compConv(CENTROIDS, SUM_CLUSTER, MatDim[0], K); /* calculates the difference between the old and the new centroids */
        swapMatrices(&CENTROIDS, &SUM_CLUSTER); /* updates the centroids */
        iteration ++;
    }
    copyArr(CENTROIDS, CENTO, MatDim[0] * K); /* copies the final result into CENTO */

    free(SUM_CLUSTER);
    free(SIZE_CLUSTER);
    free(CENTROIDS);

    return 0;
}

/* FUNCTION DEFINING */

void raiseError()
{
    printf("An Error Has Occured!");
    exit(1);
}

void raiseInputError()
{
    printf("Invalid Input!");
    exit(1);
}

/* calculates the squared distance between the vectors x and y of VecDim dimension */
double dist_sq(double *x, double *y, int VecDim)
{
    double sum = 0;
    int i = 0;
    while (i < VecDim){
        sum += pow(x[i] - y[i], 2);
        i ++;
    }
    return sum;
}

/* returns the index of the nearest centroid to a given vector vec */
int find_near_centroid(double *centroids, double *vec, int VecDim, int K)
{
    int ind = 0 ;
    double argmin = dist_sq(centroids, vec, VecDim);
    int j = 1;
    double d_2;
    while (j < K) /* loops over all the centroids */
    {
        d_2 = dist_sq(centroids + j * VecDim, vec, VecDim);
        if (d_2 < argmin)
        {
            ind = j;
            argmin = d_2;
        }
        j++;
    }
    return ind;
}

/* calculates the difference between the old and the new centroids */
double compConv(double *oldCent, double *newCent, int VecDim, int K)
{                             
    double delta = dist_sq(oldCent, newCent, VecDim * K);
    return sqrt(delta);
}

/* adds a given vector Vec to a cluster with given cluster index */
void AddVecToCluster(double *Vec,double *sum_clusters, int *size_clusters, int cluster_index, int VecDim)
{
    int i = 0;
    while (i < VecDim){ /* adds Vec to the sum of vectors in the given cluster */
        sum_clusters[cluster_index * VecDim + i] += Vec[i]; 
        i ++;
    }
    size_clusters[cluster_index] += 1; /* increases the size of the given cluster by 1 */
}

/* calculates new centroids and writes them into sum_cluster */
void newCentroids(double *sum_cluster, int *size_cluster,int VecDim, int K)
{
    int i = 0;
    int j = 0;
    while (i < K){
        while (j < VecDim){ /* divides each cluster's sum by its size */
            if (0 != size_cluster[i]) /* prevented division by 0 */
            {
                sum_cluster[i * VecDim + j] = sum_cluster[i * VecDim + j] / size_cluster[i];
            }
            j++;
        }
        j = 0;
        i++;
    }
}

void setZeroArrayDouble(double *Mat, int arr_len)
{
    int i = 0;
    for (i = 0; i < arr_len; i++)
    {
        Mat[i] = 0.0;
    }
}

void setZeroArrayInt(int *Mat, int arr_len)
{
    int i = 0;
    for (i = 0; i < arr_len; i++)
    {
        Mat[i] = 0;
    }
}

/* copies source array into target array */
void copyArr(double *source, double *target, int length)
{
    int i = 0;
    for (i = 0; i < length; i++){
        target[i] = source[i];
    }
}

void swapMatrices(double **A, double **B)
{
    double *temp = *A;
    *A = *B;
    *B = temp;
}


/* C TO PYTHON MODULE HANDLING */

static int ArrayTo_PyList(PyObject *list, double *array, int length)
{
    int index;
    PyObject *CURR;
    int STATUS = 0;
    for (index = 0; index < length; index++)
    {
        CURR = PyFloat_FromDouble(array[index]);
        if (PyList_SetItem(list, index, CURR) == -1) /* if setting item failed */
        {
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
    for (index = 0; index < length; index++)
    {
        CURR = PyList_GetItem(list, index);
        if (!PyFloat_Check(CURR))
        {
            STATUS = 1;
        }
        array[index] = PyFloat_AsDouble(CURR);  
    }
    return STATUS; /* STATUS equals 0 means success, otherwise, failure */
}

static PyObject* fit(PyObject* self, PyObject* args)
{
    PyObject *PO_POINTS; /* PyList of points */
    PyObject *PO_CENTROIDS; /* PyList of centroids */
    double *POINTO; /* array of points */
    double *CENTO; /* array of centroids */
    int K;
    int D; /* vector dimension */
    int N; /* number of points */
    double EPS;
    int MAX_ITER;
    PyObject *OUTPUT; /* PyList storing the result */

    if (!PyArg_ParseTuple(args, "OOiiiid", &PO_POINTS, &PO_CENTROIDS, &N, &D, &K, &MAX_ITER, &EPS)) /*try Parsing*/
        return NULL;
    POINTO = malloc(N * D * sizeof(double));
    CENTO = malloc(K * D * sizeof(double));

    if (NULL == POINTO || NULL == CENTO){ /*if allocation failed*/
        raiseError();
    }

    if(PyListTo_Array(PO_POINTS, POINTO, N * D) != 0)
        raiseError();
    
    if(PyListTo_Array(PO_CENTROIDS, CENTO, K * D) != 0)
        raiseError();
    

    MAIN(MAX_ITER, CENTO, D, N, POINTO, EPS, K); /* the algorithm */
    OUTPUT = PyList_New(D * K); 
    if (NULL == OUTPUT){ /*if list creation failed*/
        raiseError();
    }
    ArrayTo_PyList(OUTPUT, CENTO, D * K);
    
    free(CENTO);
    free(POINTO);
    
    return OUTPUT;
}

static PyMethodDef mymetoda[] = {
    {"fit", fit, METH_VARARGS, "fits K centroids according to K means algorithm"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef mykmeanssp = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    "I hope it works",
    -1,
    mymetoda
};

PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&mykmeanssp);
    if (!m){
        return NULL;
    }
    return m;
}



