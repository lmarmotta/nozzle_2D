#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "structs.h"
#include "externs.h"
#include "prototypes.h"

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))


/* Computes the pullian and Chausee matrices. */

void compute_pc_matrices(t_define p_setup, t_points ** pnts){

    /* Separate bounds of the field. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Loop through internal points and compute the needed stuff.*/

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){

            /* Compute the matrices direction. */

            double mu = 1.0/pow(2.0,0.5);

            double m1 = pnts[i][j].ksi_x*pnts[i][j].eta_x + pnts[i][j].ksi_y*pnts[i][j].eta_y;

            double m2 = pnts[i][j].ksi_x*pnts[i][j].eta_y + pnts[i][j].ksi_y*pnts[i][j].eta_x;

            /* Compute the matrix. */

            pnts[i][j].nh[0][0] = 1.0;
            pnts[i][j].nh[0][1] = 0.0; 
            pnts[i][j].nh[0][2] = 0.0; 
            pnts[i][j].nh[0][3] = 0.0; 

            pnts[i][j].nh[1][0] = 0.0;
            pnts[i][j].nh[1][1] = m1;
            pnts[i][j].nh[1][2] = -mu*m2;
            pnts[i][j].nh[1][3] = mu*m2;

            pnts[i][j].nh[2][0] = 0.0;
            pnts[i][j].nh[2][1] = mu*m2;
            pnts[i][j].nh[2][2] = pow(mu,2.0)*(1.0 + m1);
            pnts[i][j].nh[2][3] = pow(mu,2.0)*(1.0 - m1);

            pnts[i][j].nh[3][0] = 0.0;
            pnts[i][j].nh[3][1] = -mu*m2;
            pnts[i][j].nh[3][2] = pow(mu,2.0)*(1.0 - m1);
            pnts[i][j].nh[3][3] = pow(mu,2.0)*(1.0 + m1);

            /* Compute the inverse of matrix. */

            pnts[i][j].nhi[0][0] = 1.0;
            pnts[i][j].nhi[0][1] = 0.0; 
            pnts[i][j].nhi[0][2] = 0.0; 
            pnts[i][j].nhi[0][3] = 0.0; 

            pnts[i][j].nhi[1][0] = 0.0;
            pnts[i][j].nhi[1][1] = m1;
            pnts[i][j].nhi[1][2] = mu*m2;
            pnts[i][j].nhi[1][3] = -mu*m2;

            pnts[i][j].nhi[2][0] = 0.0;
            pnts[i][j].nhi[2][1] = -mu*m2;
            pnts[i][j].nhi[2][2] = pow(mu,2.0)*(1.0 + m1);
            pnts[i][j].nhi[2][3] = pow(mu,2.0)*(1.0 - m1);

            pnts[i][j].nhi[3][0] = 0.0;
            pnts[i][j].nhi[3][1] = mu*m2;
            pnts[i][j].nhi[3][2] = pow(mu,2.0)*(1.0 - m1);
            pnts[i][j].nhi[3][3] = pow(mu,2.0)*(1.0 + m1);


        }
    }

    /* Compute the diagonal eigenvalue matrices. */

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){

            /* Bizu vars. */

            U = pnts[i][j].cov_u;
            V = pnts[i][j].cov_u;

            a = pnts[i][j].a;

            ksix = pnts[i][j].ksi_x;
            ksiy = pnts[i][j].ksi_y;

            etax = pnts[i][j].eta_x;
            etay = pnts[i][j].eta_y;

            /* Ksi direction. */

            pnts[i][j].lambda_ksi[0][0] = U;
            pnts[i][j].lambda_ksi[0][1] = 0.0; 
            pnts[i][j].lambda_ksi[0][2] = 0.0; 
            pnts[i][j].lambda_ksi[0][3] = 0.0; 

            pnts[i][j].lambda_ksi[1][0] = 0.0; 
            pnts[i][j].lambda_ksi[1][1] = U;
            pnts[i][j].lambda_ksi[1][2] = 0.0; 
            pnts[i][j].lambda_ksi[1][3] = 0.0; 

            pnts[i][j].lambda_ksi[2][0] = 0.0; 
            pnts[i][j].lambda_ksi[2][1] = 0.0; 
            pnts[i][j].lambda_ksi[2][2] = U + a*pow(pow(ksix,2) + pow(ksiy,2.0),0.5);
            pnts[i][j].lambda_ksi[2][3] = 0.0; 

            pnts[i][j].lambda_ksi[3][0] = 0.0; 
            pnts[i][j].lambda_ksi[3][1] = 0.0; 
            pnts[i][j].lambda_ksi[3][2] = 0.0; 
            pnts[i][j].lambda_ksi[3][3] = U - a*pow(pow(ksix,2) + pow(ksiy,2.0),0.5);

            /* Eta direction. */

            pnts[i][j].lambda_eta[0][0] = V;
            pnts[i][j].lambda_eta[0][1] = 0.0; 
            pnts[i][j].lambda_eta[0][2] = 0.0; 
            pnts[i][j].lambda_eta[0][3] = 0.0; 

            pnts[i][j].lambda_eta[1][0] = 0.0; 
            pnts[i][j].lambda_eta[1][1] = V;
            pnts[i][j].lambda_eta[1][2] = 0.0; 
            pnts[i][j].lambda_eta[1][3] = 0.0; 

            pnts[i][j].lambda_eta[2][0] = 0.0; 
            pnts[i][j].lambda_eta[2][1] = 0.0; 
            pnts[i][j].lambda_eta[2][2] = V + a*pow(pow(etax,2) + pow(etay,2.0),0.5);
            pnts[i][j].lambda_eta[2][3] = 0.0; 

            pnts[i][j].lambda_eta[3][0] = 0.0; 
            pnts[i][j].lambda_eta[3][1] = 0.0; 
            pnts[i][j].lambda_eta[3][2] = 0.0; 
            pnts[i][j].lambda_eta[3][3] = V - a*pow(pow(etax,2) + pow(etay,2.0),0.5);

        }
    }
}

/* Compute the pullian chausse implicit operator. */

void imp_pullian_chaussee(t_define p_setup, t_points ** pnts){

    /* Separate bounds of the field. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Auxiliar matrices. */

    double ** aux1 = alloc_dmatrix(4, 4);
    double ** aux2 = alloc_dmatrix(4, 4);
    double ** aux3 = alloc_dmatrix(4, 4);

    /* Multiply the operators and obtain the Y1. */

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){

            /* Copy the residue. */

            for (ii = 0; ii<4; ii++)
            for (jj = 0; jj<4; jj++) aux1[ii][jj] = pnts[i][j].RHS[ii];

            /* Copy the T1 matrix. */

            for (ii = 0; ii<4; ii++)
            for (jj = 0; jj<4; jj++) aux2[ii][jj] = T1_ksi[4][4];

            /* Multiply the matrix. */

            dmuls(aux2, aux1, aux3, 4);

            /* Store Y1 in ksi direction. */

            for (ii = 0; ii<4; ii++)
            for (jj = 0; jj<4; jj++) pnts[i][j].Y1[ii][jj] = aux3[ii][jj];

        }
    }

    /* Identity matrix. */

    double I[4][4];

    I[0][0] = 1.0; I[0][1] = 0.0; I[0][2] = 0.0; I[0][3] = 0.0; 
    I[1][0] = 0.0; I[1][1] = 1.0; I[1][2] = 0.0; I[1][3] = 0.0; 
    I[2][0] = 0.0; I[2][1] = 0.0; I[2][2] = 1.0; I[2][3] = 0.0; 
    I[3][0] = 0.0; I[3][1] = 0.0; I[3][2] = 0.0; I[3][3] = 1.0; 

    /* Now, solve the system for the first direction. */

    double *** X1  = alloc_dcube(imax-1, jmax-1, 4);
    double *** Y1c = alloc_dcube(imax-1, jmax-1, 4);

    /* Auxiliar. */

    double *** upper_j = alloc_dcube(4,4,jmax-1);
    double *** maind_j = alloc_dcube(4,4,jmax-1);
    double *** lower_j = alloc_dcube(4,4,jmax-1);

    for (int j = 0; j<jmax; j++){
        for (int i = 0; i<imax; i++){

            /* Estimate the constant. */

            double cte_i = pnts[i][j].dt;

            /* Build the three diagonals. */

            for (int ii = 0; ii<nim; ii++){
                for (int jj = 0; jj<nim; jj++){

                    upper_j[ii][jj][i] = cte_i * pnts[i+1][j].lambda_ksi[ii][jj];

                    main_j[ii][jj][i] = I[ii][jj];

                    lower_j[ii][jj][i] = cte_i * pnts[i-1][j].lambda_ksi[ii][jj];

                }
            }

            /* Build the B vector using the RHS. */

            s_rhsi[0][i] = pnts[i][j].RHS[0];
            s_rhsi[1][i] = pnts[i][j].RHS[1];
            s_rhsi[2][i] = pnts[i][j].RHS[2];
            s_rhsi[3][i] = pnts[i][j].RHS[3];
        }
    }


    /* Free the matrices. */

    free_dmatrix(aux1, 4);
    free_dmatrix(aux2, 4);
    free_dmatrix(aux3, 4);
    free_dmatrix(X1,  imax-1, jmax-1);
    free_dmatrix(Y1c, imax-1, jmax-1);
}
