#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"
#include "externs.h"
#include "prototypes.h"

/* Advance the solution in time using the Beam-Warming Scheme. */

/* The definition of the artificial dissipation is not well stablished in the
 * code. I am sure in need of some internet right now. To check the proper
 * implementation, of course. */

void beam_warming(t_define p_setup, t_points ** pnts){

    /* Select bounds. */
    
    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Define inner matrices size. */

    int nim = 4;

    /* Dissipation. */

    double diss = 0.0;

    /* Define the identity matrix. */

    double I[4][4];

    I[0][0] = 1.0; I[0][1] = 0.0; I[0][2] = 0.0; I[0][3] = 0.0; 
    I[1][0] = 0.0; I[1][1] = 1.0; I[1][2] = 0.0; I[1][3] = 0.0; 
    I[2][0] = 0.0; I[2][1] = 0.0; I[2][2] = 1.0; I[2][3] = 0.0; 
    I[3][0] = 0.0; I[3][1] = 0.0; I[3][2] = 0.0; I[3][3] = 1.0; 
    
    /* Allocate the du_s vector. It has to be the same size of the residue. */

    double *** DU_S  = alloc_dcube(imax-1, jmax-1, 4);

    /* Transport vector for the i and j passes. */

    double **  du_si = alloc_dmatrix(4,imax-1);
    double **  du_sj = alloc_dmatrix(4,jmax-1);

    /* Transport resudue vector for system solve. */

    double **  s_rhsi = alloc_dmatrix(4,imax-1);
    double **  s_rhsj = alloc_dmatrix(4,jmax-1);

    /* Allocate the three diagonals for the i pass */

    double *** upper_i = alloc_dcube(4,4,imax-1);
    double *** maind_i = alloc_dcube(4,4,imax-1);
    double *** lower_i = alloc_dcube(4,4,imax-1);

    /* Allocate the three diagonals for the j pass */

    double *** upper_j = alloc_dcube(4,4,jmax-1);
    double *** maind_j = alloc_dcube(4,4,jmax-1);
    double *** lower_j = alloc_dcube(4,4,jmax-1);

    /* Now loop through every line of the mesh. */

    for (int j = 1; j<jmax-1; j++){
        for (int i = 1; i<imax-1; i++){

            /* Estimate the constant. */

            double cte_i = pnts[i][j].dt;

            /* Build the upper diagonal. */

            for (int ii = 0; ii<nim; ii++)
                for (int jj = 0; jj<nim; jj++){

                    /* Dissipation free operator. */

                    upper_i[ii][jj][i] = (cte_i * pnts[i+1][j].A_hat[ii][jj])/2.0;

                    /* Compose the dissipative term. */

                    diss = - (2.0*p_setup.dissp2) * pnts[i+1][j].J1*pnts[i+1][j].J;

                    /* Add the dissipative term. */

                    if (ii == jj) upper_i[ii][jj][i] = upper_i[ii][jj][i] + diss; 
                }

            /* Build the main diagonal. */

            for (int ii = 0; ii<nim; ii++)
                for (int jj = 0; jj<nim; jj++){

                    /* Dissipation free operator. */

                    maind_i[ii][jj][i] = I[ii][jj];

                    /* Compose the dissipative term. */

                    diss = - (2.0*p_setup.dissp2) * pnts[i][j].J1*(-2.0*pnts[i][j].J);

                    /* Add the dissipative term. */

                    if (ii == jj) maind_i[ii][jj][i] = maind_i[ii][jj][i] + diss;

                }

            /* Build the lower diagonal. */

            for (int ii = 0; ii<nim; ii++)
                for (int jj = 0; jj<nim; jj++){

                    /* Dissipation free operator. */

                    lower_i[ii][jj][i] = - (cte_i * pnts[i-1][j].A_hat[ii][jj])/2.0;

                    /* Compose the dissipative term */

                    diss = - (2.0*p_setup.dissp2) * pnts[i-1][j].J1*pnts[i-1][j].J;

                    /* Add the dissipative term. */

                    if (ii == jj) lower_i[ii][jj][i] = lower_i[ii][jj][i] + diss;
                }

            /* Build the B vector using the RHS. */

            s_rhsi[0][i] = pnts[i][j].dt * pnts[i][j].RHS[0];
            s_rhsi[1][i] = pnts[i][j].dt * pnts[i][j].RHS[1];
            s_rhsi[2][i] = pnts[i][j].dt * pnts[i][j].RHS[2];
            s_rhsi[3][i] = pnts[i][j].dt * pnts[i][j].RHS[3];

        }

        /* Now, solve the system and store the du_s. */

        blk_tri(maind_i, lower_i, upper_i, nim, imax-1, s_rhsi, du_si);

        /* Store the final DU_S. */

        for (int i = 1; i<imax-1; i++){

            DU_S[i][j][0] = du_si[0][i];
            DU_S[i][j][1] = du_si[1][i];
            DU_S[i][j][2] = du_si[2][i];
            DU_S[i][j][3] = du_si[3][i];

        }
    }

    /* Now lets do the same in the other direction. */

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            /* Estimate the constant. */

            double cte_i = pnts[i][j].dt;

            /* Build the upper diagonal. */

            for (int ii = 0; ii<nim; ii++)
                for (int jj = 0; jj<nim; jj++){

                    /* Dissipation free operator. */

                    upper_j[ii][jj][j] = (cte_i * pnts[i][j+1].B_hat[ii][jj])/2.0;

                    /* Compose the dissipative term. */

                    diss = - (2.0*p_setup.dissp2) * pnts[i][j+1].J1*pnts[i][j+1].J;

                    /* Add the dissipative term. */

                    if (ii == jj) upper_j[ii][jj][j] = upper_j[ii][jj][j] + diss;
                }

            /* Build the main diagonal. */

            for (int ii = 0; ii<nim; ii++)
                for (int jj = 0; jj<nim; jj++){

                    /* Dissipation free operator. */

                    maind_j[ii][jj][j] = I[ii][jj];

                    /* Compose the dissipative term. */

                    diss = - (2.0*p_setup.dissp2) * pnts[i][j].J1*(-2.0*pnts[i][j].J);

                    /* Add the dissipative term. */

                    if (ii == jj) maind_j[ii][jj][j] = maind_j[ii][jj][j] + diss;

                }

            /* Build the lower diagonal. */

            for (int ii = 0; ii<nim; ii++)
                for (int jj = 0; jj<nim; jj++){

                    /* Dissipation free operator. */

                    lower_j[ii][jj][j] = - (cte_i * pnts[i][j-1].B_hat[ii][jj])/2.0;

                    /* Compose the dissipative operator. */

                    diss = - (2.0*p_setup.dissp2) * pnts[i][j-1].J1*pnts[i][j-1].J;

                    /* Add the dissipative term. */

                    if (ii == jj) lower_j[ii][jj][j] = lower_j[ii][jj][j] + diss;
                }

            /* Build the B vector using the RHS. */

            s_rhsj[0][j] = - pnts[i][j].dt * DU_S[i][j][0];
            s_rhsj[1][j] = - pnts[i][j].dt * DU_S[i][j][1];
            s_rhsj[2][j] = - pnts[i][j].dt * DU_S[i][j][2];
            s_rhsj[3][j] = - pnts[i][j].dt * DU_S[i][j][3];

        }

        /* Now, solve the system and store the du_s. */

        blk_tri(maind_j, lower_j, upper_j, nim, jmax-1, s_rhsj, du_sj);

        /* Update the solution. */

        for (int j = 1; j<jmax-1; j++){

           pnts[i][j].q_hat[0] = pnts[i][j].q_hat[0] + du_sj[0][j];  
           pnts[i][j].q_hat[1] = pnts[i][j].q_hat[1] + du_sj[1][j];  
           pnts[i][j].q_hat[2] = pnts[i][j].q_hat[2] + du_sj[2][j];  
           pnts[i][j].q_hat[3] = pnts[i][j].q_hat[3] + du_sj[3][j];  

        }
    }

    /* Free everyone. */

    free_dcube(upper_i,4,4);
    free_dcube(maind_i,4,4);
    free_dcube(lower_i,4,4);
    free_dcube(upper_j,4,4);
    free_dcube(maind_j,4,4);
    free_dcube(lower_j,4,4);
    free_dmatrix(du_si,4);
    free_dmatrix(du_sj,4);
    free_dmatrix(s_rhsi,4);
    free_dmatrix(s_rhsj,4);
    free_dcube(DU_S,imax-1, jmax-1);
}
