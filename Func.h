#ifndef FUNC_H_INCLUDED
#define FUNC_H_INCLUDED


#include <time.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include <time.h>
#include <sys/time.h>
#include <pthread.h>
double mach_eps(void);
double NormSol(double * Matrix, double* InvMatrix, int n, double* X);
int Gauss(double* Matrix, double* InvMatrix, int n, int* swap, double* X, double* Y);
int Input(double *A,char *file_name, int n, double* Norm, int p, int number);
double generate(int k, int n, int i, int j);
void POut2(double *InvMatrix,int l,int n,int m);
void POut(double *Matrix,int l,int n,int m, double* X);
#endif // FUNC_H_INCLUDED
