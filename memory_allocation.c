#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "cgnslib.h"
#include "structs.h"
#include "externs.h"

/* 
 * Allocate the main struct of the code.
 */

void alloc_struct_matrix(t_points *** pntss, int imax, int jmax){

    /* Malloc the first row of the memory. */

    t_points ** pnts = (t_points**)calloc(imax,sizeof(t_points*));

    /* Check if the allocation of the first row is successfull. */

    if (pnts == NULL){ 
        printf("ERROR: Memory Allocation\n"); exit(1); 
    }

    /* Now, loop through pointers initializing the columns. */

    for (int i = 0; i<imax; i++){

        pnts[i] = NULL;

        pnts[i] = (t_points*)calloc(jmax,sizeof(t_points)); 

        if (pnts[i] == NULL){
            printf("ERROR: Memory Allocation\n"); exit(1); 
        }
    }

    *pntss = pnts;
}

/* 
 * Free the main struct of the code.
 */

void free_struct_matrix(t_points ** pnts, int imax){

    for (int i = 0; i<imax; i++) free(pnts[i]);
        
    free(pnts);
}

/* Alloc double vector. */

double * alloc_dvector(int imax){

    double * out = NULL;

    out = (double*)calloc(imax,sizeof(double));

    if (out == NULL) {printf("ERROR: double matrix allocation.\n"); exit(1);}

    return out;
}

/* Allocate a double matrix. */

double ** alloc_dmatrix(int imax, int jmax){

    double ** out = NULL;

    out = (double**)calloc(imax,sizeof(double));

    if (out == NULL) {printf("ERROR: double matrix allocation.\n"); exit(1);}

    for (int i = 0; i<imax; i++){
        out[i] = (double*)calloc(jmax,sizeof(double));
        if (out[i] == NULL) {printf("ERROR: double matrix allocation.\n"); exit(1);}
    }

    /* Initialize the matrix. */

    for (int i = 0; i<imax; i++)
        for (int j = 0; j<jmax; j++) out[i][j] = 0.0;

    return out;
}

/* Free vector */

void free_vector(double * vector){
    free(vector);
}

/* Free double matrix. */

void free_dmatrix(double ** matrix, int imax){

    for (int i = 0; i<imax; i++)
        free(matrix[i]);

    free(matrix);
}

/* Free cube. */

void free_dcube(double *** cube, int imax, int jmax){

    for (int i = 0; i<imax; i++)
        for (int j = 0; j<jmax; j++) free(cube[i][j]);
            
    for (int i = 0; i<imax; i++) free(cube[i]);

    free(cube);

}

/* Alloc a double cube. */

double *** alloc_dcube(int imax, int jmax, int kmax){

    double *** out = NULL;

    /* Allocate the first line. */

    out = (double***)calloc(imax,sizeof(double**));

    /* Check the allocation pointer position. */

    if (out == NULL) {printf("ERROR: double matrix allocation.\n"); exit(1);}

    /* Allocate the matrix. */

    for (int i = 0; i<imax; i++){

        out[i] = (double**)calloc(jmax,sizeof(double*));
        if (out[i] == NULL) {printf("ERROR: double matrix allocation.\n"); exit(1);}

        /* Allocate the block. */

        for (int j = 0; j<jmax; j++){
            out[i][j] = (double*)calloc(kmax,sizeof(double));
            if (out[i][j] == NULL) {printf("ERROR: double matrix allocation.\n"); exit(1);}
        }
    }

    /* Initialize the cube. */

    for (int i = 0; i<imax; i++)
        for (int j = 0; j<jmax; j++) 
            for (int k = 0; k<kmax; k++) out[i][j][k] = 0.0;

    return out;
}

/*
 * Free a double matrix.
 */

void free_double_matrix(double ** matrix, int imax){

    for (int i = 0; i<imax; i++){
        free(matrix[i]);
    }

    free(matrix);
}

/* 
 * Initialize all values of the struct. 
 */

void initialize_structs(t_define p_setup, t_points ** pnts){

    /* Separate the limits of the mesh. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* initialize all variables. */

    for (int i = 0; i < imax; i++){
        for (int j = 0; j < jmax; j++){

            pnts[i][j].x = 0.0;
            pnts[i][j].y = 0.0;

            pnts[i][j].J  = 0.0;
            pnts[i][j].J1 = 0.0;

            pnts[i][j].x_ksi = 0.0;
            pnts[i][j].y_ksi = 0.0;

            pnts[i][j].x_eta = 0.0;
            pnts[i][j].y_eta = 0.0;

            pnts[i][j].ksi_x = 0.0;
            pnts[i][j].ksi_y = 0.0;

            pnts[i][j].eta_x = 0.0;
            pnts[i][j].eta_y = 0.0;

            pnts[i][j].cov_u = 0.0;
            pnts[i][j].cov_v = 0.0;

            pnts[i][j].q_hat[0] = 0.0;
            pnts[i][j].q_hat[1] = 0.0;
            pnts[i][j].q_hat[2] = 0.0;
            pnts[i][j].q_hat[3] = 0.0;

            pnts[i][j].e_hat[0] = 0.0;
            pnts[i][j].e_hat[1] = 0.0;
            pnts[i][j].e_hat[2] = 0.0;
            pnts[i][j].e_hat[3] = 0.0;

            pnts[i][j].f_hat[0] = 0.0;
            pnts[i][j].f_hat[1] = 0.0;
            pnts[i][j].f_hat[2] = 0.0;
            pnts[i][j].f_hat[3] = 0.0;

            pnts[i][j].RHS[0] = 0.0;
            pnts[i][j].RHS[1] = 0.0;
            pnts[i][j].RHS[2] = 0.0;
            pnts[i][j].RHS[3] = 0.0;

            pnts[i][j].a = 0.0;
            pnts[i][j].m = 0.0;
            pnts[i][j].p = 0.0;
            pnts[i][j].t = 0.0;

            pnts[i][j].dt = 0.0;
        }
    }
}

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
