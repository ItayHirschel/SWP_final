
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>


#define EPS_JACOBI 0.00001
#define MAX_ITER_JACOBI 100

enum mode { /* WORK MODES ENUM */
    ERROR = -1,
    WAMmode,
    DDGmode,
    LNORMmode,
    JACOBImode,
    SPKmode
};

enum ARGUMENT { /* Argument placement */
    PROGRAM_NAME = 0,
    GOAL_PLACE,
    FILE_NAME_PLACE
};

/* FUNCTION DECLARATIONS */

void raiseError();
void raiseInputError();
void copyArr(double *source, double *target, int length);
void PrintMatrix(double *Mat, int rows, int columns);
double dist_sq(double *x, double *y, int VecDim);
double Norma(double *x, int VecDim);
double dist(double *x, double *y, int VecDim);
void Matrix_Mult(double *LeftMat, double *RightMat, double *TargetMat, int i_dim, int j_dim, int k_dim);
void AddArray(double *left, double *right, double *target, int length);
void SubArray(double *left, double *right, double *target, int length);
void Identity(double *target, int n);
void WAM_create(double *points, double *target, int N, int VecDim);
void DDM_create(double *wam, double *target, int N);
void sqrt_Diag(double *matrix, int N);
void Lnorm_create(double *D_sqrt, double *W, double *target, int N);
void find_pivot(double *matrix, int *i_coor, int *j_coor, int N);
void obtain_CS(double *matrix, int i_coor, int j_coor, int N, double *C_value, double *S_value);
void Rotation_Mat_create(double *target, int N, double c, double s, int i_coor, int j_coor);
void Rotate_Sym(double *S_mat, int N, double c, double s, int i_coor, int j_coor);
double sq_sum_off_diag(double *Mat, int N);
void Jacobi(double *Sym, int N, double *target_vec, double *target_val);
void Swap_helper(double *Vals, double *Vecs, int i, int j, int N);
void Bubble_sort(double *Vals, double *Vecs, int N);
int Eigengap(double *Vals, int N);
void Largest_k_vecs(double *Sorted_mat, int K, int N, double *target);
void Renorma(double *U, int N, int K, double *target);
void ddg_process(double *Points, int N, int VecDim, double *target);
void Lnorm_process(double *Points, int N, int VecDim, double *target);
void spk_process1(double *Points, int N, int VecDim, int *K, double *target);
void spk_process2(double *Vecs, int K, int N, double *target);
enum mode translation(char* goal);
void Main_args(char *goal, int N, int VecDim, int *K, double *Points, double *target1, double *target2);
void get_MatDim(char *in_file, int *K, int *N);
void read_IN_FILE(int N, int D, char *in_file, double *Points);

/* MAIN FUNCTION */

int main(int argc, char **argv)
{
    char *goal, *file_name;
    double *Points, *target1, *target2;
    int D, N, k = -1;
    enum mode MODE;
    if (3 != argc) /* check for correct input format */
    {
        raiseInputError();
    }
    /* extract arguments */
    goal = argv[GOAL_PLACE];
    file_name = argv[FILE_NAME_PLACE];
    /* get dimensions into D, N */
    get_MatDim(file_name, &D, &N);
    Points = malloc(D * N * sizeof(double));
    if (NULL == Points) /* if allocation failed raise error */
    {
        raiseError();
    }
    /* extract data into Points */
    read_IN_FILE(N, D, file_name, Points);

    MODE = translation(goal);

    if (ERROR == MODE || SPKmode == MODE) /* check for legal work mode */
    {
        raiseInputError();
    }

    target1 = malloc(N * N * sizeof(double));
    target2 = malloc(N * sizeof(double));
    if (NULL == target1 || NULL == target2) /* if allocation failed raise error */
    {
        raiseError();
    }

    if (D != N && JACOBImode == MODE)  /* Make sure that Points is a symmetric matrix for Jacobi */
    {
        raiseError();
    }
    /* apply the main algorithm */
    Main_args(goal, N, D, &k, Points, target1, target2);

    /* print output */
    if (JACOBImode == MODE)
    {
        PrintMatrix(target2, 1, N);
    }
    PrintMatrix(target1, N, N);

    free(Points);
    free(target1);
    free(target2);
    return 0;
}

/* FUNCTION DEFINING */

/* raises a runtime error and exits */
void raiseError()
{
    printf("An Error Has Occured!");
    exit(1);
}

/* raises an Input Error for wrong input formats */
void raiseInputError()
{
    printf("Invalid Input!");
    exit(1);
}

/* Copies source array into target array */
void copyArr(double *source, double *target, int length)
{
    int i = 0;
    for (i = 0; i < length; i++){
        target[i] = source[i];
    }
}

/* prints an array as matrix of dimensions rows x columns */
void PrintMatrix(double *Mat, int rows, int columns)
{
    int i = 0;
    int j = 0;
    for (i = 0; i < rows; i++){
        for (j = 0; j < columns; j++){
            printf("%.4f", Mat[i * columns + j]);
            if (columns - 1 != j){
                printf(",");
            }
        }
        printf("\n");
    }
}

/* calculate the squared distance between vectors x, y of length VecDim */
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

/* Calculates the Size (distance from 0) of vector x */
double Norma(double *x, int VecDim)
{
    double sum = 0;
    int i;
    for(i = 0; i < VecDim; i++)
    {
        sum += pow(x[i], 2);
    }
    return sqrt(sum);
}

/* calculate the distance between vectors x, y of length VecDim */
double dist(double *x, double *y, int VecDim)
{
    return sqrt(dist_sq(x, y, VecDim));
}

/* Calculate the matrix multlipication of RightMat (i_dim x j_dim) and LeftMat (j_dim x k_dim)
   and writes result into TargetMat (i_dim x k_dim) 

   !!! WARNING !!! putting LeftMat or RightMat as TargetMat may give wrong results
*/
void Matrix_Mult(double *LeftMat, double *RightMat, double *TargetMat, int i_dim, int j_dim, int k_dim) /*Matrix Multlipication*/
{
    int i = 0;
    int j = 0;
    int k = 0;
    double sum;
    for (i = 0; i < i_dim; i++)
    {
        for (k = 0; k < k_dim; k++)
        {
            sum = 0;
            for (j = 0; j < j_dim; j++)
            {
                sum += LeftMat[i * j_dim + j] * RightMat[j * k_dim + k];
            }
            TargetMat[i * k_dim + k] = sum;
        }
    }
}

/*adding arrays*/
void AddArray(double *left, double *right, double *target, int length)
{ 
    int i;
    for (i = 0; length > i; i++)
    {
        target[i] = left[i] + right[i];
    }
}

/*Subtracting arrays*/
void SubArray(double *left, double *right, double *target, int length)
{ 
    int i;
    for (i = 0; length > i; i++){
        target[i] = left[i] - right[i];
    }
}

/*Create an Identity Matrix (N x N) in target*/
void Identity(double *target, int n)
{ 
    int i;
    int j;
    for (i = 0; n > i; i++){
        for (j = 0; n > j; j++){
            target[i * n + j] = (i == j ? 1 : 0);
        }
    }
}

/* Creates the  WAM matrix (N x N) appropriate to the given points (N x D) and writes into target*/
void WAM_create(double *points, double *target, int N, int VecDim)
{
    int i = 0;
    int j = 0;
    double val = 0;
    for (i = 0; N > i; i++){
        for (j = 0; N > j; j++){
            if (i == j)
            {
                val = 0; /*Diagonal is 0 */
            }
            else
            {
                val = exp(-0.5 * dist(points + i * VecDim, points + j * VecDim, VecDim)); /*set correct value to val */
            }
            target[i * N + j] = val; /* write val to appropriate location in target matrix */
        }
    }
}

/* Creates the Diagonal Degree Matrix (N x N) according to the given wam matrix (N x N) and writes into target*/
void DDM_create(double *wam, double *target, int N)
{
    int i = 0;
    int j = 0;
    int z = 0;
    double sum = 0;

    for (i = 0; N > i; i++){
        for (j = 0; N > j; j++){
            sum = 0;
            if (i == j)
            {
                for (z = 0; z < N; z++) /* summing row i in WAM*/
                { 
                    sum += wam[i * N + z];
                }
            }
            target[i * N + j] = sum; /*Diagonal values setting */
        }
    }
}

/* divides 1 by the square roots the diagonal of given N x N matrix */
void sqrt_Diag(double *matrix, int N)
{
    int i;
    for (i = 0; N > i; i++){
        if (0 != matrix[i * N + i])
        {
            matrix[i * N + i] = 1.0 / sqrt(matrix[i * N + i]);
        }
    }
}

/* Calculate the Lnorm matrix (N x N) for given WAM matrix (wam, N x N) and D^(-1/2) (D_sqrt, N x N) and writes into target*/
void Lnorm_create(double *D_sqrt, double *W, double *target, int N)
{
    double *eye;
    double *DW;
    eye = malloc(N * N * sizeof(double)); 
    DW = malloc(N * N * sizeof(double));
    if (NULL == eye || NULL == DW){ /* if allocation failed */
        raiseError();
    }
    Identity(eye, N); /* eye = I */
    Matrix_Mult(D_sqrt, W, DW, N, N, N); /* DW = D_sqrt * W  */
    Matrix_Mult(DW, D_sqrt, target, N, N, N); /*  target = D_sqrt * W * D_sqrt  */
    SubArray(eye, target, target, N * N); /*target = I - D_sqrt * W * D_sqrt  */
    free(DW);
    free(eye);
}

/* find the off diagonal maximum value, AKA pivot of rotation, for symmetrical
matrix (N x N) and writes the found i,j indices into i_coor, j_coor */
void find_pivot(double *matrix, int *i_coor, int *j_coor, int N) 
{
    int i, j;
    double m = fabs(matrix[1]);
    *i_coor = 0;
    *j_coor = 1;
    for (i = 0; i < N; i++)
    {
        for (j = i + 1; j < N; j++)
        {
            if (fabs(matrix[i * N + j]) > m)
            {
                *i_coor = i;
                *j_coor = j;
                m = fabs(matrix[i * N + j]);
            }
        }
    }
}

/* calculate the C, S values necessary for the rotation for given matrix (N x N), i_coor, j_coor 
and writes into C_value, S_value accordingly */
void obtain_CS(double *matrix, int i_coor, int j_coor, int N, double *C_value, double *S_value)
{
    double theta, t;
    if (0 != matrix[i_coor * N + j_coor])
        theta = (matrix[j_coor * N + j_coor] - matrix[i_coor * N + i_coor]) / (2 * matrix[i_coor * N + j_coor]);
    else
    {
        theta = 0;
    }
    t = (theta < 0.0 ? -1.0 : 1.0) / (fabs(theta) + sqrt(theta * theta + 1));
    *C_value = 1 / sqrt(t * t + 1);
    *S_value = t * *C_value;
}

/* writes Rotation MAtrix based on i_coor, j_coor, c, s into target */
void Rotation_Mat_create(double *target, int N, double c, double s, int i_coor, int j_coor)
{
    Identity(target, N);
    target[i_coor * N + i_coor] = c;
    target[j_coor * N + j_coor] = c;
    target[i_coor * N + j_coor] = s;
    target[j_coor * N + i_coor] = -s;
}

/* Rotates (Multiplies by R^T * S_mat * R where R is the rotation matrix) S_mat (N x N) in place */
void Rotate_Sym(double *S_mat, int N, double c, double s, int i_coor, int j_coor)
{
    int r = 0;
    double *Copy_S;
    double val = 0;
    Copy_S = malloc(N * N * sizeof(double));
    if (NULL == Copy_S){
        raiseError();
    }
    copyArr(S_mat, Copy_S, N * N);
    /*determine A[r][i] and A[i][r],  Working on row and column I*/
    for (r = 0; r < N; r++){
        if (r == i_coor) /* in case of Aii*/
        { 
            val = c * c * Copy_S[i_coor * N + i_coor]
            + s * s * Copy_S[j_coor * N + j_coor]
            - 2 * c * s * Copy_S[i_coor * N + j_coor];
        }
        else
        {
            if (r == j_coor) /* in case of Aij*/
            { 
                val = 0;
            }
            else /* in case of r different from i,j*/
            {   
                val = c * Copy_S[r * N + i_coor] - s * Copy_S[r * N + j_coor];
            }
        }
        S_mat[i_coor * N + r] = val; 
        S_mat[r * N + i_coor] = val;
    }

    /*determine A[r][j] and A[j][r],  Working on row and column I*/

    for (r = 0; r < N; r++){
        if (r == j_coor) /* in case of Ajj*/
        { 
            val = s * s * Copy_S[i_coor * N + i_coor]
            + c * c * Copy_S[j_coor * N + j_coor]
            + 2 * c * s * Copy_S[i_coor * N + j_coor];
        }
        else
        {
            if (r == i_coor) /* in case of Aji*/
            { 
                val = 0;
            }
            else /* in case of r different from i,j*/
            { 
                val = c * Copy_S[r * N + j_coor] + s * Copy_S[r * N + i_coor];
            }
        }
        S_mat[j_coor * N + r] = val;
        S_mat[r * N + j_coor] = val;
    }


    free(Copy_S);
}

/* calculate the squared sum of all off diagonal elements in matrix Mat (N x N) */
double sq_sum_off_diag(double *Mat, int N) /*compute the convergence*/
{
    int i = 0;
    int j = 0;
    double sum = 0;
    for(i = 0; i < N; i++)
    {
        for(j = 0; j < N; j++)
        {
            if (i != j)
            {
                sum += pow(Mat[i * N + j], 2);
            }
        }
    }
    return sum;
}

/* Apply the Jacobi Algorithm to a Symmetric Matrix Sym (N x N)
and write eigenvalues and eigen vector into target_val (N vector), target_vec (N x N) accordingly */
void Jacobi(double *Sym, int N, double *target_vec, double *target_val)
{
    int iter = 0, i_coor, j_coor, i;
    double c, s, *P, *V, conv = 1.0 , sum, prev_sum;

    P = malloc(N * N * sizeof(double));
    V = malloc(N * N * sizeof(double));

    if (NULL == P || NULL == V){
        raiseError();
    }

    Identity(target_vec, N);
    Identity(V, N);
    /* through V and target_vec we multiply all rotation matrices in order */
    prev_sum = sq_sum_off_diag(Sym, N);

    while(iter < MAX_ITER_JACOBI && conv > EPS_JACOBI && 0.0 != prev_sum)
    {
        find_pivot(Sym, &i_coor, &j_coor, N);
        obtain_CS(Sym, i_coor, j_coor, N, &c, &s);
        Rotation_Mat_create(P, N, c, s, i_coor, j_coor); /* P is Rotation Matrix */
        Matrix_Mult(target_vec, P, V, N, N, N); /* V = target_vec * P */
        copyArr(V, target_vec, N*N); /* target_vec = V = (previous target_vec) * P */
        Rotate_Sym(Sym, N, c, s, i_coor, j_coor);
        sum = sq_sum_off_diag(Sym, N);
        conv = fabs(sum - prev_sum);
        iter += 1;
        prev_sum = sum;
    }
    /* write the diagonal into an N-array of eigenvalues */
    for(i = 0; i < N; i++)
    {
        target_val[i] = Sym[i * N + i];
    }
    free(P);
    free(V);
}

/* swaps elements indexed i, j in Vals (N array) and appropriate columns (i, j) in Vec (N x N)  */
void Swap_helper(double *Vals, double *Vecs, int i, int j, int N)
{
    double temp;
    int ind;

    /* swap elements i,j in array Vals */
    temp = Vals[i];
    Vals[i]= Vals[j];
    Vals[j] = temp;

    /*swap i,j columns */
    for(ind = 0; ind < N; ind++)
    {
        temp = Vecs[ind * N + i];
        Vecs[ind * N + i] = Vecs[ind * N + j];
        Vecs[ind * N + j] = temp;
    }
}

/* Applies descending bubble sort algorithm to Vals (N array) and Vecs (N x N) */
void Bubble_sort(double *Vals, double *Vecs, int N)
{
    int iter, ind;

    for (iter = 0; iter < N; iter++)
    {
        for (ind = 0; ind < N-1; ind++)
        {
            if(Vals[ind] < Vals[ind + 1])
            {
                Swap_helper(Vals, Vecs, ind, ind + 1, N);
            }
        }
    }
}

/* Calculate the EigenGap heuristic given sorted eigenvalues Vals */
int Eigengap(double *Vals, int N)
{
    int ind = 1, i;
    double max = 0.0 ,delta;

    for (i = 0; i < N / 2; i++)
    {
        delta = Vals[i] - Vals[i + 1];
        if (delta > max)
        {
            max = delta;
            ind = i + 1;
        }
    }
    return ind;
}

/* extracts k leftmost columns of Sorted_Mat (N x N) into target (N x K) */
void Largest_k_vecs(double *Sorted_mat, int K, int N, double *target)
{
    int i, j;
    for(i = 0; i < N; i++)
    {
        for(j = 0; j < K; j++){
            target[i * K + j]= Sorted_mat[i * N + j];
        }
    }
}

/* Renormalize each row to norm of 1 and write into target */
void Renorma(double *U, int N, int K, double *target)
{
    int i, j;
    double norm;

    for(i = 0; i < N; i++) 
    {
        norm = Norma(U + i * K, K); /* Calculate the size of row i */
        for (j = 0; j < K; j++) /* normalize */
        {
            if (0 != norm)
            {
                target[i * K + j] = U[i * K + j] / norm; 
            }
            else
            {
                target[i * K + j] = 0;
            }
        }
    }
}

/*  Process needed for ddg goal, calculate DDM (N x N) according to given points (N x VecDim) and write result
into target */
void ddg_process(double *Points, int N, int VecDim, double *target)
{
    double *wam = malloc(N * N * sizeof(double));
    if (NULL == wam) /* if allocation failed raise error */
    {
        raiseError();
    }
    WAM_create(Points, wam, N, VecDim);
    DDM_create(wam, target, N);
    free(wam);
}

/*  Process needed for lnorm goal, calculate Normalized Graph Laplacian (N x N) 
according to given points (N x VecDim) and write result into target */
void Lnorm_process(double *Points, int N, int VecDim, double *target)
{
    double *wam = malloc(N * N * sizeof(double));
    double *ddm = malloc(N * N * sizeof(double));
    if (NULL == wam || NULL == ddm){ /* if allocation failed raise error */
        raiseError();
    }
    WAM_create(Points, wam, N, VecDim);
    DDM_create(wam, ddm, N);
    sqrt_Diag(ddm, N);
    Lnorm_create(ddm, wam, target, N);
    free(wam);
    free(ddm);
}

/* First Half of the process needed for spk goal, calculate Eigenvectors (N x N) and eigenvalue (N array)
of the Normalized Graph Laplacian according to given points (N x VecDim) and write result into target 
In addition calculates K if necessary and writes in K*/
void spk_process1(double *Points, int N, int VecDim, int *K, double *target)
{
    double *L = malloc(N * N * sizeof(double));
    double *vals = malloc(N * sizeof(double));

    if (NULL == L || NULL == vals){ /* if allocation failed raise error */
        raiseError();
    }

    Lnorm_process(Points, N, VecDim, L);
    Jacobi(L, N, target, vals);
    Bubble_sort(vals, target, N);
    if(0 == *K) /* if K is 0 then it needs to be calculated */
    {
        *K = Eigengap(vals, N);
    }
    free(L);
    free(vals);
}

/* Second Half of the process needed for spk goal, get K Eigenvectors (N x K) with larget eigenvalues
from Vecs and renormalize them - write into target*/
void spk_process2(double *Vecs, int K, int N, double *target)
{
    double *largestK = malloc(K * N * sizeof(double));
    if (NULL == largestK){
        raiseError();
    }
    Largest_k_vecs(Vecs, K, N, largestK);
    Renorma(largestK, N, K, target);
    free(largestK);
}

/* translate goal string into appropriate enum mode */
enum mode translation(char* goal)
{
    enum mode m = ERROR;
    if(strcmp("wam", goal) == 0)
    {
        m = WAMmode;
    }
    if(strcmp("ddg", goal) == 0)
    {
        m = DDGmode;
    }
    if(strcmp("lnorm", goal) == 0)
    {
        m = LNORMmode;
    }
    if(strcmp("jacobi", goal) == 0)
    {
        m = JACOBImode;
    }
    if(strcmp("spk", goal) == 0)
    {
        m = SPKmode;
    }
    return m;
}

/* target1 is a N*N matrix for matrices or eigen vectors, target2 is an N-sized array for eigenvalues*/
/* The Main algorithm used in spkmeans.c
gets points (Points), goal (goal), DImensions (N, K, D)
write results into target1 (N x N), target2(N array) */
void Main_args(char *goal, int N, int VecDim, int *K, double *Points, double *target1, double *target2)
{
   enum mode m = translation(goal);
   switch (m)
   {
    case WAMmode:
     WAM_create(Points, target1, N, VecDim);
     break;
   
    case DDGmode:
     ddg_process(Points, N, VecDim,target1);
     break;
    
    case LNORMmode:
     Lnorm_process(Points, N, VecDim, target1);
     break;
    
    case JACOBImode:
     Jacobi(Points, N, target1, target2);
     break;
    
    case SPKmode:
     spk_process1(Points, N, VecDim, K, target1);
     break;
   
    default: 
     break;
   }
}


/* Gets the dimensions of matrix in file named in_file, writes into K,N */
void get_MatDim(char *in_file, int *K, int *N)
{
    char *c;
    int row_len = 0;
    FILE *ftpr;
    FILE *ftpr2;
    char current;
    char *str_num;
    char *delim = ",";
    int k = 0;
    int lineNotEmpty = 0;
    int i = 0;

    if ((ftpr = fopen(in_file,"r")) == NULL){               /*if cannot open exit*/
        raiseError();
    }
    if ((ftpr2 = fopen(in_file,"r")) == NULL){               /*if cannot open exit*/
        raiseError();
    }

    while (current != '\n' && !feof(ftpr)){                 /*measure length of first row*/
        current = fgetc(ftpr);
        row_len ++;
    }

    fseek(ftpr, 0, 0);                                      /*return to start*/
    c = calloc(row_len, sizeof(char));                      /*allocate necessary space for reading row*/
    if (c == NULL){ /* if allocation failed raise error */
        raiseError();
    }
    fscanf(ftpr, "%[^\n]", c);                              /*scan until '\n' is encountered*/

    str_num = strtok(c, delim);                             /*split string by ',' */
    while (str_num != NULL){          /*measure dimension of vector by the amount of tokens str_num is split to*/
        k ++;
        str_num = strtok(NULL, delim);
    }
    *K = k;         /* write number of columns to K */

    fseek(ftpr, 0, 0);                                      /*return to start*/
    row_len += k;                                           /*each number might be negative thus adding a - sign*/
    free(c);                                                /* free c for further allocation */
    c = calloc(row_len, sizeof(char));                      /*allocate neccessary space for reading next row*/
    if (c == NULL){ /* if allocation failed raise error */
        raiseError();
    }
    k = 0;

    current = ' ';
    while (current != '\n' && !feof(ftpr2)){                 /*measure length of first row*/
        current = fgetc(ftpr2);
    }
    while (fgets(c, row_len, ftpr)){                        /*while there is a next line*/
        lineNotEmpty = 0;                                   /* "boolean" 1 iff line is not empty, otherwise 0*/
        for (i = 0; i < (int) strlen(c); i++){              /*checking if line is empty*/
            if (! isspace(*(c + i))){
                lineNotEmpty = 1;
            }
        }
        if (lineNotEmpty){      
            k ++; /* count number of not empty lines in in_file */
        }
        current = ' ';
        row_len = *K;
        while (current != '\n' && !feof(ftpr2)){                 /*measure length of row*/
            current = fgetc(ftpr2);
            row_len ++;
        }
        free(c);
        c = calloc(row_len, sizeof(char));                      /*allocate neccessary space for reading next row*/
        if (NULL == c){
            raiseError();
        }
    }
    *N = k; /* write number of lines to N */
    fclose(ftpr);                                           /*close file*/
    fclose(ftpr2);
    free(c);                                                /*free space allocation*/
}

/* read the input file and extract data into Points (N x D) */
void read_IN_FILE(int N, int D, char *in_file, double *Points)
{
    int i = 0;
    int j = 0;
    FILE *file;
    double temp;

    if ((file = fopen(in_file,"r")) == NULL){               /*if cannot open exit*/
        raiseError();
    }

    while(i < N){
        while(j < D){
            fscanf(file, "%lf,", &temp);
            Points[i * D + j] = temp;
            j++;
        }
        j = 0;
        i++;
    }
    fclose(file);
}
