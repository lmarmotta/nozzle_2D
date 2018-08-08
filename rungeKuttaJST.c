#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

/* This function marches the solution in time. */

void rungeKuttaJST(t_define p_setup, t_points ** pnts){

    /* Declare the constants of the scheme. */

    double alpha_1, alpha_2, alpha_3, alpha_4, alpha_5;

    alpha_1 = 1.0/4.0; alpha_2 = 1.0/6.0; alpha_3 = 3.0/8.0; alpha_4 = 1.0/2.0, alpha_5 = 1.0;

    /* Allocate a simple vector to store everyone. */

}

/* This function computes the local time step. */

void local_time(t_define p_setup, t_points ** pnts){

    int i, j;

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Loop through the points computing the dt. */

    for (i = 0; i < imax; i++){
        for (j = 0; j < jmax; j++){

            /* Compute the sound velocity. */

            double rho = pnts[i][j].q[0];
            double u   = pnts[i][j].q[1] / pnts[i][j].q[0];
            double v   = pnts[i][j].q[2] / pnts[i][j].q[0];

            double e = pnts[i][j].q[4];
            double p = (p_setup.gamma - 1.0)*(e - 0.5*rho*( pow(u,2.0) + pow(v,2.0)) );
            double a =  sqrt(p_setup.gamma*p/pnts[i][j].q[0]);

            /* Compute the characteristic speed. */

            double term_1 = abs(pnts[i][j].cov_u) + a*sqrt( pow(pnts[i][j].ksi_x,2.0) + pow(pnts[i][j].ksi_y,2.0) );
            double term_2 = abs(pnts[i][j].cov_v) + a*sqrt( pow(pnts[i][j].eta_x,2.0) + pow(pnts[i][j].eta_y,2.0) );

            double c_ij = max(term_1, term_2);

            /* Now, compute and store the local dt. */

            pnts[i][j].dt = p_setup.CFL / c_ij;


        }
    }

}
