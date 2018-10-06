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

    /* Separate the diagonals. */

    double *** m_lower = alloc_dcube(4,4,jmax-1);
    double *** m_main  = alloc_dcube(4,4,jmax-1);
    double *** m_upper = alloc_dcube(4,4,jmax-1);

    /* Separate solution vectors. */

    double ** du_s = alloc_dmatrix(4,jmax-1);
    double ** srhs = alloc_dmatrix(4,jmax-1);

    /* Define the identity matrix. */

    double I[4][4];

    I[0][0] = 1.0; I[0][1] = 0.0; I[0][2] = 0.0; I[0][3] = 0.0; 
    I[1][0] = 0.0; I[1][1] = 1.0; I[1][2] = 0.0; I[1][3] = 0.0; 
    I[2][0] = 0.0; I[2][1] = 0.0; I[2][2] = 1.0; I[2][3] = 0.0; 
    I[3][0] = 0.0; I[3][1] = 0.0; I[3][2] = 0.0; I[3][3] = 1.0; 

    /* Loop through the field and solve the system for all points in eta. */

    for (int i = 1; i<imax-1; i++){

        /* Do the first direction. */

        for (int j = 1; j<jmax-1; j++){

            /* Separate the scalar of the system. */

            double cnts = (pnts[i][j].dt * theta) / (1.0 + c_ksi);

            /* Separate the lower diagonal. */

            for (int ii = 0; ii<nim; ii++)
                for (int jj = 0; jj<nim; jj++)
                    m_lower[ii][jj][j] = - (cnts*pnts[i-1][j].A_hat[ii][jj]) / 2.0;

            /* Separate the main diagonal. */

            for (int ii = 0; ii<nim; ii++)
                for (int jj = 0; jj<nim; jj++)
                    m_main[ii][jj][j] = I[ii][jj];

            /* Separate the lower diagonal. */

            for (int ii = 0; ii<nim; ii++)
                for (int jj = 0; jj<nim; jj++)
                    m_upper[ii][jj][j] = - (cnts*pnts[i+1][j].A_hat[ii][jj]) / 2.0;

            /* Separate the srhs. */

            srhs[0][j] = pnts[i][j].RHS[0];
            srhs[1][j] = pnts[i][j].RHS[1];
            srhs[2][j] = pnts[i][j].RHS[2];
            srhs[3][j] = pnts[i][j].RHS[3];
        }

        /* Solve the system. */

        blk_tri(m_main, m_lower, m_upper, 4, jmax-1, srhs, du_s);

        /* Do the other direction. */

        for (int j = 1; j<jmax-1; j++){

            /* Separate the scalar of the system. */

            double cnts = (pnts[i][j].dt * theta) / (1.0 + c_ksi);

            /* Separate the lower diagonal. */

            for (int ii = 0; ii<nim; ii++)
                for (int jj = 0; jj<nim; jj++)
                    m_lower[ii][jj][j] = - (cnts*pnts[i-1][j].B_hat[ii][jj]) / 2.0;

            /* Separate the main diagonal. */

            for (int ii = 0; ii<nim; ii++)
                for (int jj = 0; jj<nim; jj++)
                    m_main[ii][jj][j] = I[ii][jj];

            /* Separate the lower diagonal. */

            for (int ii = 0; ii<nim; ii++)
                for (int jj = 0; jj<nim; jj++)
                    m_upper[ii][jj][j] = - (cnts*pnts[i+1][j].B_hat[ii][jj]) / 2.0;

            /* Separate the srhs. */

            srhs[0][j] = du_s[0][j];
            srhs[1][j] = du_s[1][j];
            srhs[2][j] = du_s[2][j];
            srhs[3][j] = du_s[3][j];
        }

        /* Solve the system. */

        blk_tri(m_main, m_lower, m_upper, 4, jmax-1, srhs, du_s);

        /* Update the solution. */

        for (int j = 1; j<jmax-1; j++){

            pnts[i][j].q_hat[0] = pnts[i][j].q_hat[0] + du_s[j][0]; 
            pnts[i][j].q_hat[1] = pnts[i][j].q_hat[1] + du_s[j][1]; 
            pnts[i][j].q_hat[2] = pnts[i][j].q_hat[2] + du_s[j][2]; 
            pnts[i][j].q_hat[3] = pnts[i][j].q_hat[3] + du_s[j][3]; 

        }

        /* Clean the solution vectors. */

        for (int j = 1; j<jmax-1; j++){

            /* Clean the du_s */

            du_s[0][j] = 0.0;
            du_s[1][j] = 0.0;
            du_s[2][j] = 0.0;
            du_s[3][j] = 0.0;

            /* Clean the srhs */

            srhs[0][j] = 0.0;
            srhs[1][j] = 0.0;
            srhs[2][j] = 0.0;
            srhs[3][j] = 0.0;
        }
    }

    /* Free the vectors. */

    free_dcube(m_lower,4,4);
    free_dcube(m_main,4,4);
    free_dcube(m_upper,4,4);

    /* Separate solution vectors. */

    free_dmatrix(du_s,4);
    free_dmatrix(srhs,4);
}
