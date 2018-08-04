#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "cgnslib.h"
#include "structs.h"

/*
 * Apply initial condition.
 */

void apply_initial_condition(t_define p_setup, t_points ** pnts){

    /* Separate the limits of the mesh. */

    int i, j;

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Apply a basic initial condition. */

    for (i = 0; i<imax; i++){
        for (j = 0; j<jmax; j++){

            pnts[i][j].q[0] = p_setup.i_rho;
            pnts[i][j].q[1] = p_setup.i_rhou;
            pnts[i][j].q[2] = p_setup.i_rhov;
            pnts[i][j].q[3] = p_setup.i_e;
        }
    }
}

/*
 * Build transformed fluxes.
 */

void build_fluxes(t_define p_setup, t_points ** pnts){

    /* Separate the limits of the mesh. */

    int i, j;

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Build the basic fluxes without transformation. */

    double rho, u, v, e, p;

    for (i = 0; i<imax; i++){
        for (j = 0; j<jmax; j++){

            /* Separate the properties we need. */

            rho = pnts[i][j].q[0];
            u   = pnts[i][j].q[1] / pnts[i][j].q[0];
            v   = pnts[i][j].q[2] / pnts[i][j].q[0];
            e   = pnts[i][j].q[3];
            p   = (p_setup.gamma - 1.0)*(e - 0.5*rho*( pow(u,2.0) + pow(v,2.0)) );

            /* Compute the covariant velocity components. */

            pnts[i][j].cov_u =  pnts[i][j].ksi_x*u + pnts[i][j].ksi_y*v;
            pnts[i][j].cov_v =  pnts[i][j].eta_x*u + pnts[i][j].eta_y*v;

            /* Build the Q fluxes in curvilinear coordinates. */

            pnts[i][j].q_hat[0] = pnts[i][j].jm1*rho;
            pnts[i][j].q_hat[1] = pnts[i][j].jm1*rho*u;
            pnts[i][j].q_hat[2] = pnts[i][j].jm1*rho*v;
            pnts[i][j].q_hat[3] = pnts[i][j].jm1*e;

            /* Build the E fluxes in curvilinear coordinates. */

            pnts[i][j].e_hat[0] = pnts[i][j].jm1 * (rho*pnts[i][j].cov_u);
            pnts[i][j].e_hat[1] = pnts[i][j].jm1 * (rho*u*pnts[i][j].cov_u + pnts[i][j].ksi_x*p);
            pnts[i][j].e_hat[2] = pnts[i][j].jm1 * (rho*v*pnts[i][j].cov_u + pnts[i][j].ksi_y*p);
            pnts[i][j].e_hat[3] = pnts[i][j].jm1 * (pnts[i][j].cov_u*(e + p) - p);

            /* Build the F fluxes in curvilinear coordinates. */

            pnts[i][j].f_hat[0] = pnts[i][j].jm1 * (rho*pnts[i][j].cov_v);
            pnts[i][j].f_hat[1] = pnts[i][j].jm1 * (rho*u*pnts[i][j].cov_v + pnts[i][j].eta_x*p);
            pnts[i][j].f_hat[2] = pnts[i][j].jm1 * (rho*v*pnts[i][j].cov_v + pnts[i][j].eta_y*p);
            pnts[i][j].f_hat[3] = pnts[i][j].jm1 * (pnts[i][j].cov_v*(e + p) - p);

        }
    }


}
