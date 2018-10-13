#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "cgnslib.h"
#include "structs.h"
#include "externs.h"

/* Uses LAPACK library to invert a matrix. */

void inv(double ** matrix, int N){

    int * IPIV = calloc((N+1),sizeof(int));
    int LWORK = N*N;
    double * WORK = calloc(LWORK,sizeof(double));
    int INFO;

    double A[N*N];

    for (int i = 0; i<N; i++)
        for (int j = 0; j<N; j++)
            A[(i)*(N)+j] = matrix[i][j];

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    for (int i = 0; i<N; i++)
        for (int j = 0; j<N; j++)
            matrix[i][j] = A[(i)*(N)+j];

    free(IPIV);
    free(WORK);
}

/* Multiplies squared matrices. */

void dmuls(double ** mA, double ** mB, double ** mC, int size){

    for (int i=0; i<size; i++)
        for (int j=0; j<size; j++) mC[i][j] = 0.0;

    for (int i=0; i<size; i++)
        for (int j=0; j<size; j++) 
            for (int k=0; k<size; k++)
                mC[i][j] = mC[i][j] + mA[i][k]*mB[k][j];
}

/* Matrix multiplication any size. Remember what I always forget:
 * [m x n] X [n x p] = [m x p] Where [a,b] are lines and columns. 
 *    A    x    B    =    C */ 

void dmgss(double ** mA, double ** mB, double ** mC, int m, int n, int p){

    for (int i=0; i<m; i++)
        for (int j=0; j<p; j++) mC[i][j] = 0.0;

    for (int i = 0; i < m; i++)
        for (int j = 0; j < p; j++)
            for (int k = 0; k < n; k++)
                mC[i][j] = mC[i][j] + mA[i][k]*mB[k][j];
}

/* Solves a system composed of block tridiagonal matrices. */

void blk_tri(double *** main, double *** lower, double *** upper, int size_m, int num_m, double ** XB, double ** X){

    /* Define auxiliar variables. */

    double *** gamm = alloc_dcube(size_m,size_m,num_m);

    double **  beta = alloc_dmatrix(size_m,num_m);

    double ** aux_copy = alloc_dmatrix(size_m,size_m);

    /* Get the first gamma. */

    for (int i = 0; i<size_m; i++)
        for (int j = 0; j<size_m; j++) aux_copy[i][j] = main[i][j][0];

    /* Invert matrix. */
            
    inv(aux_copy, size_m);

    double ** auxm1 = alloc_dmatrix(size_m,size_m);
    double ** auxm2 = alloc_dmatrix(size_m,size_m);

    for (int i = 0; i<size_m; i++)
        for (int j = 0; j<size_m; j++)
            auxm1[i][j] = upper[i][j][0];

    /* Multiply to obtain the first gamma. */

     dmuls(aux_copy, auxm1, auxm2, size_m);

    /* Get the gamma, aleluia !! */

    for (int i = 0; i<size_m; i++)
        for (int j = 0; j<size_m; j++)
            gamm[i][j][0] = auxm2[i][j];

    /* Get the rest of the gammas. */

    double ** aux_mult = alloc_dmatrix(size_m,size_m);
    double ** aux_summ = alloc_dmatrix(size_m,size_m);

    for (int m = 1; m<num_m-1; m++){

        for (int i = 0; i<size_m; i++){
            for (int j = 0; j<size_m; j++){ 
                auxm1[i][j] = lower[i][j][m];
                auxm2[i][j] = gamm[i][j][m-1];
            }
        }

        dmuls(auxm1, auxm2, aux_mult, size_m);

        for (int i = 0; i<size_m; i++){
            for (int j = 0; j<size_m; j++){
                aux_summ[i][j] = main[i][j][m] - aux_mult[i][j];
            }
        }

        inv(aux_summ, size_m);

        for (int i = 0; i<size_m; i++){
            for (int j = 0; j<size_m; j++){ 
                auxm1[i][j] = upper[i][j][m];
            }
        }

        dmuls(aux_summ, auxm1, auxm2, size_m);

        for (int i = 0; i<size_m; i++){
            for (int j = 0; j<size_m; j++){ 
                gamm[i][j][m] = auxm2[i][j];
            }
        }
    }

    /* Zero out our arrays. */

    for (int i = 0; i<size_m; i++)
        for (int j = 0; j<size_m; j++) aux_copy[i][j] = main[i][j][0];

    inv(aux_copy,size_m);

    double ** auxv1 = alloc_dmatrix(size_m,1);
    double ** auxv2 = alloc_dmatrix(size_m,1);

    for (int i = 0; i<size_m; i++) auxv1[i][0] = XB[i][0];

    dmgss(aux_copy, auxv1, auxv2, size_m, size_m, 1);

    double * aux_dumm = alloc_dvector(size_m);

    /* Form the first beta. */

    for (int i = 0; i<size_m; i++) beta[i][0] = auxv2[i][0]; 

    /* Go and grab the rest. */

    for (int m = 1; m<num_m; m++){

        for (int i = 0; i<size_m; i++)
            for (int j = 0; j<size_m; j++){
                aux_mult[i][j] = 0.0;
                aux_summ[i][j] = 0.0;
        }

        for (int i = 0; i<size_m; i++){
            for (int j = 0; j<size_m; j++){ 
                auxm1[i][j] = lower[i][j][m];
                auxm2[i][j] = gamm[i][j][m-1];
            }
        }

        dmuls(auxm1, auxm2, aux_mult, size_m);

        for (int i = 0; i<size_m; i++){
            for (int j = 0; j<size_m; j++){
                aux_summ[i][j] = main[i][j][m] - aux_mult[i][j];
            }
        }

        inv(aux_summ, size_m);

        for (int i = 0; i<size_m; i++){
            for (int j = 0; j<size_m; j++){
                auxm1[i][j] = lower[i][j][m];
            }
        }

        for (int i = 0; i<size_m; i++) auxv1[i][0] = beta[i][m-1];

        dmgss(auxm1, auxv1, auxv2, size_m, size_m, 1);

        for (int i = 0; i<size_m; i++)
            aux_dumm[i] = XB[i][m] - auxv2[i][0];

        for (int i = 0; i<size_m; i++) auxv1[i][0] = aux_dumm[i];

        dmgss(aux_summ, auxv1, auxv2, size_m, size_m, 1);

        for (int i = 0; i<size_m; i++) beta[i][m] = auxv2[i][0];

    }

    /* Start backward sweep. */

    for (int i = 0; i<size_m; i++)
        for (int j = 0; j<num_m; j++) X[i][j] = 0.0;

    for (int i = 0; i<size_m; i++) X[i][num_m-1] = beta[i][num_m-1];

    for (int m = num_m-2; m >= 0; m--){

        for (int i = 0; i<size_m; i++) 
            for (int j = 0; j<size_m; j++) auxm1[i][j] = gamm[i][j][m];

        for (int i = 0; i<size_m; i++) auxv1[i][0] = X[i][m+1];

        dmgss(auxm1, auxv1, auxv2, size_m, size_m, 1);

        for (int i = 0; i<size_m; i++) 
            X[i][m] = beta[i][m] - auxv2[i][0];
    }

    /* Avoid leaks baby. */

    free_dcube(gamm,size_m,size_m);

    free_dmatrix(beta,size_m);

    free_dmatrix(aux_copy,size_m);

    free_dmatrix(auxm1,size_m);

    free_dmatrix(auxm2,size_m);

    free_dmatrix(aux_mult,size_m);

    free_dmatrix(aux_summ,size_m);

    free_dmatrix(auxv1,size_m);

    free_dmatrix(auxv2,size_m);

    free_vector(aux_dumm);
}
