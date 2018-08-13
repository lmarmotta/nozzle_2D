#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"
#include "externs.h"

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

/* Use the explicit euler scheme. */

void explicitEuler(t_define p_setup, t_points ** pnts){

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    for (int i = 1; i < imax-1; i++){
        for (int j = 1; j < jmax-1; j++){

            pnts[i][j].q[0] = pnts[i][j].q[0] - pnts[i][j].dt*pnts[i][j].RHS[0];
            pnts[i][j].q[1] = pnts[i][j].q[1] - pnts[i][j].dt*pnts[i][j].RHS[1];
            pnts[i][j].q[2] = pnts[i][j].q[2] - pnts[i][j].dt*pnts[i][j].RHS[2];
            pnts[i][j].q[3] = pnts[i][j].q[3] - pnts[i][j].dt*pnts[i][j].RHS[3];

        }
    }

}

/* This function marches the solution in time using RK5. */

void rungeKuttaJST(t_define p_setup, t_points ** pnts){

    /* Separate bounds of the field. */

    int neq = 4;

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Declare the constants of the scheme. */

    double alpha[5];

    alpha[0] = 1.0/4.0; 
    alpha[1] = 1.0/6.0; 
    alpha[2] = 3.0/8.0; 
    alpha[3] = 1.0/2.0;
    alpha[4] = 1.0;

    /* Separate the proper aux vector. */

    double q_0[imax][jmax][neq];

    /* Loop through the field and copy the q vector. */

    for (int i = 0; i < imax; i++){
        for (int j = 0; j < jmax; j++){

            q_0[i][j][0] = pnts[i][j].q[0];
            q_0[i][j][1] = pnts[i][j].q[1];
            q_0[i][j][2] = pnts[i][j].q[2];
            q_0[i][j][3] = pnts[i][j].q[3];
        }
    }

    /* ---------------------------------------------------------------------- */
    /* Do the FIRST stage of the RK scheme. */

    for (int i = 1; i < imax-1; i++){
        for (int j = 1; j < jmax-1; j++){

            pnts[i][j].q[0] = q_0[i][j][0] - alpha[0]*pnts[i][j].dt*pnts[i][j].RHS[0];
            pnts[i][j].q[1] = q_0[i][j][1] - alpha[0]*pnts[i][j].dt*pnts[i][j].RHS[1];
            pnts[i][j].q[2] = q_0[i][j][2] - alpha[0]*pnts[i][j].dt*pnts[i][j].RHS[2];
            pnts[i][j].q[3] = q_0[i][j][3] - alpha[0]*pnts[i][j].dt*pnts[i][j].RHS[3];

        }
    }

    /* Re-Build the fluxes with the new Q. */

    build_fluxes(p_setup, pnts);

    /* Re-Compute the RHS.*/

    compute_rhs(p_setup, pnts);

    /* Artificial dissipation.*/

    jst_art_dissip(p_setup, pnts);

    /* Update Boundary conditions. */

    boundary_condition_euler(p_setup, pnts);

    /* Compute local time step. */

    local_time(p_setup, pnts);

    /* ---------------------------------------------------------------------- */
    /* Do the SECOND stage of the RK scheme. */

    for (int i = 1; i < imax-1; i++){
        for (int j = 1; j < jmax-1; j++){

            pnts[i][j].q[0] = q_0[i][j][0] - alpha[1]*pnts[i][j].dt*pnts[i][j].RHS[0];
            pnts[i][j].q[1] = q_0[i][j][1] - alpha[1]*pnts[i][j].dt*pnts[i][j].RHS[1];
            pnts[i][j].q[2] = q_0[i][j][2] - alpha[1]*pnts[i][j].dt*pnts[i][j].RHS[2];
            pnts[i][j].q[3] = q_0[i][j][3] - alpha[1]*pnts[i][j].dt*pnts[i][j].RHS[3];

        }
    }

    /* Re-Build the fluxes with the new Q. */

    build_fluxes(p_setup, pnts);

    /* Re-Compute the RHS.*/

    compute_rhs(p_setup, pnts);

    /* Artificial dissipation.*/

    jst_art_dissip(p_setup, pnts);

    /* Update Boundary conditions. */

    boundary_condition_euler(p_setup, pnts);

    /* Compute local time step. */

    local_time(p_setup, pnts);

    /* ---------------------------------------------------------------------- */
    /* Do the THIRD stage of the RK scheme. */

    for (int i = 1; i < imax-1; i++){
        for (int j = 1; j < jmax-1; j++){

            pnts[i][j].q[0] = q_0[i][j][0] - alpha[2]*pnts[i][j].dt*pnts[i][j].RHS[0];
            pnts[i][j].q[1] = q_0[i][j][1] - alpha[2]*pnts[i][j].dt*pnts[i][j].RHS[1];
            pnts[i][j].q[2] = q_0[i][j][2] - alpha[2]*pnts[i][j].dt*pnts[i][j].RHS[2];
            pnts[i][j].q[3] = q_0[i][j][3] - alpha[2]*pnts[i][j].dt*pnts[i][j].RHS[3];

        }
    }

    /* Re-Build the fluxes with the new Q. */

    build_fluxes(p_setup, pnts);

    /* Re-Compute the RHS.*/

    compute_rhs(p_setup, pnts);

    /* Artificial dissipation.*/

    jst_art_dissip(p_setup, pnts);

    /* Update Boundary conditions. */

    boundary_condition_euler(p_setup, pnts);

    /* Compute local time step. */

    local_time(p_setup, pnts);

    /* ---------------------------------------------------------------------- */
    /* Do the FOURTH stage of the RK scheme. */

    for (int i = 1; i < imax-1; i++){
        for (int j = 1; j < jmax-1; j++){

            pnts[i][j].q[0] = q_0[i][j][0] - alpha[3]*pnts[i][j].dt*pnts[i][j].RHS[0];
            pnts[i][j].q[1] = q_0[i][j][1] - alpha[3]*pnts[i][j].dt*pnts[i][j].RHS[1];
            pnts[i][j].q[2] = q_0[i][j][2] - alpha[3]*pnts[i][j].dt*pnts[i][j].RHS[2];
            pnts[i][j].q[3] = q_0[i][j][3] - alpha[3]*pnts[i][j].dt*pnts[i][j].RHS[3];

        }
    }

    /* Re-Build the fluxes with the new Q. */

    build_fluxes(p_setup, pnts);

    /* Re-Compute the RHS.*/

    compute_rhs(p_setup, pnts);

    /* Artificial dissipation.*/

    jst_art_dissip(p_setup, pnts);

    /* Update Boundary conditions. */

    boundary_condition_euler(p_setup, pnts);

    /* Compute local time step. */

    local_time(p_setup, pnts);

    /* ---------------------------------------------------------------------- */
    /* Do the FIFTH stage of the RK scheme. */

    for (int i = 1; i < imax-1; i++){
        for (int j = 1; j < jmax-1; j++){

            pnts[i][j].q[0] = q_0[i][j][0] - alpha[4]*pnts[i][j].dt*pnts[i][j].RHS[0];
            pnts[i][j].q[1] = q_0[i][j][1] - alpha[4]*pnts[i][j].dt*pnts[i][j].RHS[1];
            pnts[i][j].q[2] = q_0[i][j][2] - alpha[4]*pnts[i][j].dt*pnts[i][j].RHS[2];
            pnts[i][j].q[3] = q_0[i][j][3] - alpha[4]*pnts[i][j].dt*pnts[i][j].RHS[3];

        }
    }

    /* Re-Build the fluxes with the new Q. */

    build_fluxes(p_setup, pnts);

    /* Re-Compute the RHS.*/

    compute_rhs(p_setup, pnts);

    /* Artificial dissipation.*/

    jst_art_dissip(p_setup, pnts);

    /* Update Boundary conditions. */

    boundary_condition_euler(p_setup, pnts);

    /* Compute local time step. */

    local_time(p_setup, pnts);
}

/* This function computes the local time step. */

void local_time(t_define p_setup, t_points ** pnts){

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Loop through the points computing the dt. */

    for (int i = 0; i < imax; i++){
        for (int j = 0; j < jmax; j++){

            /* Compute the sound velocity. */

            double rho = pnts[i][j].q[0];
            double u   = pnts[i][j].q[1] / pnts[i][j].q[0];
            double v   = pnts[i][j].q[2] / pnts[i][j].q[0];

            double e = pnts[i][j].q[3];
            double p = (p_setup.gamma - 1.0)*(e - 0.5*rho*( pow(u,2.0) + pow(v,2.0)) );
            double a =  sqrt(p_setup.gamma*p/pnts[i][j].q[0]);

            /* Compute the characteristic speed. */

            double term_1 = fabs(pnts[i][j].cov_u) + a*sqrt( pow(pnts[i][j].ksi_x,2.0) + pow(pnts[i][j].ksi_y,2.0) );
            double term_2 = fabs(pnts[i][j].cov_v) + a*sqrt( pow(pnts[i][j].eta_x,2.0) + pow(pnts[i][j].eta_y,2.0) );

            double c_ij = max(term_1, term_2);

            /* Now, compute and store the local dt. */

            pnts[i][j].dt = p_setup.CFL / c_ij;

        }
    }
}
