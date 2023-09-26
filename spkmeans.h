#ifndef SPKMEANS_H
#define SPKMEANS_H

enum mode {
    WAMmode = 0,
    DDGmode = 1,
    LNORMmode = 2,
    JACOBImode = 3,
    SPKmode = 4,
    ERROR = -1
};
void Main_args(char *goal, int N, int VecDim, int *K, double *Points, double *target1, double *target2);
enum mode translation(char* goal);
void spk_process2(double *Vecs, int K, int N, double *target);

#endif