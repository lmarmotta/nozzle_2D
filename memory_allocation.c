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

    out = (double**)calloc(imax,sizeof(double*));

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
