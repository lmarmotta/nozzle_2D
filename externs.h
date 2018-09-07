#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/* Public variables. */

double max_rhs_rho;
double max_rhs_rhou;
double max_rhs_rhov;
double max_rhs_e;

/* LAPACK externals. */

extern void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
extern void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
