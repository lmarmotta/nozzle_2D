#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"
#include "externs.h"
#include "prototypes.h"

/* Advance the solution in time using the Beam-Warming Scheme. */

void beam_warming(t_define p_setup, t_points ** pnts){

    /* Select bounds. */
    
    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Define inner matrices size. */

    int nim = 4;

    /* Set the constants of the algorithm. */

    double theta = 1.0;
    double c_ksi = 0.5;

    /* Define the identity matrix. */

    double I[4][4];

    I[0][0] = 1.0; I[0][1] = 0.0; I[0][2] = 0.0; I[0][3] = 0.0; 
    I[1][0] = 0.0; I[1][1] = 1.0; I[1][2] = 0.0; I[1][3] = 0.0; 
    I[2][0] = 0.0; I[2][1] = 0.0; I[2][2] = 1.0; I[2][3] = 0.0; 
    I[3][0] = 0.0; I[3][1] = 0.0; I[3][2] = 0.0; I[3][3] = 1.0; 
    
    /* Allocate the du_s vector. It has to be the same size of the residue. */

    double *** DU_S  = alloc_dcube(imax, jmax, 4);
    double **  du_si = alloc_dmatrix(imax-1,4);
    double **  s_rhs = alloc_dmatrix(imax-1,4);

    /* Allocate the three diagonals. */

    double *** upper = alloc_dcube(4,4,imax-1);
    double *** maind = alloc_dcube(4,4,imax-1);
    double *** lower = alloc_dcube(4,4,imax-1);

    /* Now loop through every line of the mesh. */

    for (int j = 1; j<jmax-1; j++){
        for (int i = 1; i<imax-1; i++){

            /* Estimate the constant. */

            double cte_i = (theta*pnts[i][j].dt)/(1.0 + c_ksi);

            /* Build the upper diagonal. */

            for (int ii = 0; ii<nim; ii++)
                for (int jj = 0; jj<nim; jj++){

                    /* Dissipation free operator. */

                    upper[ii][jj][i] = cte_i * (pnts[i+1][j].A_hat[ii][jj])/2.0;

                    /* Add the dissipative term. */

                    upper[ii][jj][i] = upper[ii][jj][i] + (p_setup.dissp2) * (pnts[i+1][j].J*pnts[i+1][j].q_hat[0]);
                }

            /* Build the main diagonal. */

            for (int ii = 0; ii<nim; ii++)
                for (int jj = 0; jj<nim; jj++){

                    /* Dissipation free operator. */

                    maind[ii][jj][i] = cte_i * (I[ii][jj])/2.0;

                    /* Add the dissipative term. */

                    maind[ii][jj][i] = maind[ii][jj][i] + (p_setup.dissp2) * - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[0]);
                }

            /* Build the lower diagonal. */

            for (int ii = 0; ii<nim; ii++)
                for (int jj = 0; jj<nim; jj++){

                    /* Dissipation free operator. */

                    lower[ii][jj][i] = cte_i * (pnts[i-1][j].A_hat[ii][jj])/2.0;

                    /* Add the dissipative term. */

                    lower[ii][jj][i] = lower[ii][jj][i] + (p_setup.dissp2) * (pnts[i-1][j].J*pnts[i-1][j].q_hat[0]);
                }

            /* Build the B vector using the RHS. */

            s_rhs[i][0] = pnts[i][j].RHS[0];
            s_rhs[i][1] = pnts[i][j].RHS[1];
            s_rhs[i][2] = pnts[i][j].RHS[2];
            s_rhs[i][3] = pnts[i][j].RHS[3];

        }

        /* Now, solve the system and store the du_s. */

        blk_tri(maind, lower, upper, nim, imax-1, s_rhs, du_si);

        /* Store the final DU_S. */

        for (int i = 1; i<imax-1; i++){

            DU_S[i][j][0] = du_si[i][0];
            DU_S[i][j][1] = du_si[i][1];
            DU_S[i][j][2] = du_si[i][2];
            DU_S[i][j][3] = du_si[i][3];

        }
    }

    /* Now, free the space used in the i pass. */

    free_dmatrix(du_si, imax-1);
    free_dmatrix(s_rhs, imax-1);
    free_dcube(upper, 4, 4);
    free_dcube(maind, 4, 4);
    free_dcube(lower, 4, 4);

    /* Allocate again every needed guy with proper dimension. */

    double **  du_sj = alloc_dmatrix(jmax-1,4);
    double **  s_rhj = alloc_dmatrix(jmax-1,4);

    /* Allocate the three diagonals. */

    double *** jupper = alloc_dcube(4,4,jmax-1);
    double *** jmaind = alloc_dcube(4,4,jmax-1);
    double *** jlower = alloc_dcube(4,4,jmax-1);

    /* Now loop through every column of the mesh. */

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            /* Estimate the constant. */

            double cte_i = (theta*pnts[i][j].dt)/(1.0 + c_ksi);

            /* Build the upper diagonal. */

            for (int ii = 0; ii<nim; ii++)
                for (int jj = 0; jj<nim; jj++){

                    /* Dissipation free operator. */

                    jupper[ii][jj][j] = cte_i * (pnts[i][j+1].A_hat[ii][jj])/2.0;

                    /* Add the dissipative term. */

                    jupper[ii][jj][j] = jupper[ii][jj][i] + (p_setup.dissp2) * (pnts[i][j+1].J*pnts[i][j+1].q_hat[0]);
                }

            /* Build the main diagonal. */

            for (int ii = 0; ii<nim; ii++)
                for (int jj = 0; jj<nim; jj++){

                    /* Dissipation free operator. */

                    jmaind[ii][jj][j] = cte_i * (I[ii][jj])/2.0;

                    /* Add the dissipative term. */

                    jmaind[ii][jj][j] = jmaind[ii][jj][j] + (p_setup.dissp2) * - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[0]);
                }

            /* Build the lower diagonal. */

            for (int ii = 0; ii<nim; ii++)
                for (int jj = 0; jj<nim; jj++){

                    /* Dissipation free operator. */

                    jlower[ii][jj][j] = cte_i * (pnts[i][j-1].A_hat[ii][jj])/2.0;

                    /* Add the dissipative term. */

                    jlower[ii][jj][j] = lower[ii][jj][i] + (p_setup.dissp2) * (pnts[i][j-1].J*pnts[i][j-1].q_hat[0]);
                }

            /* Build the B vector using the RHS. */

            s_rhj[j][0] = pnts[i][j].RHS[0];
            s_rhj[j][1] = pnts[i][j].RHS[1];
            s_rhj[j][2] = pnts[i][j].RHS[2];
            s_rhj[j][3] = pnts[i][j].RHS[3];

            }

        /* Now, solve the system and store the du_s. */

        blk_tri(jmaind, jlower, jupper, nim, jmax-1, s_rhj, du_sj);

        /* Store the final DU_S. */

        for (int j = 1; j<jmax-1; j++){

            DU_S[i][j][0] = du_sj[j][0];
            DU_S[i][j][1] = du_sj[j][1];
            DU_S[i][j][2] = du_sj[j][2];
            DU_S[i][j][3] = du_sj[j][3];
        }
    }

    /* Do not forget the last vector allocated. */

    free_dcube(DU_S, imax, jmax);

    /* Free every auxiliary array. */

    free_dmatrix(du_sj, jmax-1);
    free_dmatrix(s_rhs, jmax-1);
    free_dcube(upper, 4, 4);
    free_dcube(maind, 4, 4);
    free_dcube(lower, 4, 4);

}
