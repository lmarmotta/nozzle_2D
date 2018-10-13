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

/* Allocation... */

extern double * alloc_dvector(int imax);
extern double ** alloc_dmatrix(int imax, int jmax);
extern double *** alloc_dcube(int imax, int jmax, int kmax);
extern void free_vector(double * vector);
extern void free_dmatrix(double ** matrix, int imax);
extern void free_dcube(double *** cube, int imax, int jmax);
