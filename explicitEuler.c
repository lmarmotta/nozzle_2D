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

            /* March in time. */

            pnts[i][j].q[0] = pnts[i][j].q[0] - pnts[i][j].dt*pnts[i][j].RHS[0];
            pnts[i][j].q[1] = pnts[i][j].q[1] - pnts[i][j].dt*pnts[i][j].RHS[1];
            pnts[i][j].q[2] = pnts[i][j].q[2] - pnts[i][j].dt*pnts[i][j].RHS[2];
            pnts[i][j].q[3] = pnts[i][j].q[3] - pnts[i][j].dt*pnts[i][j].RHS[3];

            /* Re-Build the transformed fluxes. */

            pnts[i][j].q_hat[0] = pnts[i][j].jm1 * pnts[i][j].q[0]; 
            pnts[i][j].q_hat[1] = pnts[i][j].jm1 * pnts[i][j].q[1];
            pnts[i][j].q_hat[2] = pnts[i][j].jm1 * pnts[i][j].q[2];
            pnts[i][j].q_hat[3] = pnts[i][j].jm1 * pnts[i][j].q[3];

        }
    }
}

/* This function computes the local time step. */

void local_time(t_define p_setup, t_points ** pnts){

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Loop through the points computing the dt. */

    for (int i = 1; i < imax-1; i++){
        for (int j = 1; j < jmax-1; j++){

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

            pnts[i][j].dt = 0.0;

            pnts[i][j].dt = p_setup.CFL / c_ij;

        }
    }
}
