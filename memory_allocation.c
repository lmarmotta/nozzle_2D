#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "cgnslib.h"
#include "structs.h"

/* Allocate the main struct of the code. */

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

/* Free the main struct of the code. */

void free_struct_matrix(t_points ** pnts, int imax){

    for (int i = 0; i<imax; i++) free(pnts[i]);
        
    free(pnts);
}

/* Allocate matrix. */

double ** alloc_double_matrix(int imax, int jmax){

    double ** matrix = NULL;

    matrix = (double**)calloc(imax, sizeof(double*));

    if (matrix == NULL) {
        printf("ERROR: Memory Allocation\n"); exit(1); 
    }

    for (int i = 0; i<imax; i++){

        matrix[i] = (double*)calloc(jmax, sizeof(double));

        if (matrix[i] == NULL) {
            printf("ERROR: Memory Allocation\n"); exit(1); 
        }
    }

    return matrix;
}

void free_double_matrix(double ** matrix, int imax){

    for (int i = 0; i<imax; i++){
        free(matrix[i]);
    }

    free(matrix);
}

/* Initialize all values of the struct. */

void initialize_structs(t_define p_setup, t_points ** pnts){

    /* Separate the limits of the mesh. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* initialize all variables. */

    for (int i = 0; i < imax; i++){
        for (int j = 0; j < jmax; j++){

            pnts[i][j].x = -9999.0;
            pnts[i][j].y = -9999.0;

            pnts[i][j].J  = -9999.0;
            pnts[i][j].J1 = -9999.0;

            pnts[i][j].x_ksi = -9999.0;
            pnts[i][j].y_ksi = -9999.0;

            pnts[i][j].x_eta = -9999.0;
            pnts[i][j].y_eta = -9999.0;

            pnts[i][j].ksi_x = -9999.0;
            pnts[i][j].ksi_y = -9999.0;

            pnts[i][j].eta_x = -9999.0;
            pnts[i][j].eta_y = -9999.0;

            pnts[i][j].cov_u = -9999.0;
            pnts[i][j].cov_v = -9999.0;

            pnts[i][j].q_hat[0] = -9999.0;
            pnts[i][j].q_hat[1] = -9999.0;
            pnts[i][j].q_hat[2] = -9999.0;
            pnts[i][j].q_hat[3] = -9999.0;

            pnts[i][j].e_hat[0] = -9999.0;
            pnts[i][j].e_hat[1] = -9999.0;
            pnts[i][j].e_hat[2] = -9999.0;
            pnts[i][j].e_hat[3] = -9999.0;

            pnts[i][j].f_hat[0] = -9999.0;
            pnts[i][j].f_hat[1] = -9999.0;
            pnts[i][j].f_hat[2] = -9999.0;
            pnts[i][j].f_hat[3] = -9999.0;

            pnts[i][j].RHS[0] = -9999.0;
            pnts[i][j].RHS[1] = -9999.0;
            pnts[i][j].RHS[2] = -9999.0;
            pnts[i][j].RHS[3] = -9999.0;

            pnts[i][j].a = -9999.0;
            pnts[i][j].m = -9999.0;
            pnts[i][j].p = -9999.0;
            pnts[i][j].t = -9999.0;

            pnts[i][j].dt = -9999.0;
        }
    }
}
