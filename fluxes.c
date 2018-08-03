#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "cgnslib.h"
#include "structs.h"

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

            rho = pnts[i][j].q[1];
            u   = pnts[i][j].q[2] / pnts[i][j].q[1];
            v   = pnts[i][j].q[3] / pnts[i][j].q[1];
            e   = pnts[i][j].q[4];
            p   = (p_setup.gamma - 1.0)*(e - 0.5*rho*( pow(u,2.0) + pow(v,2.0)) );

        }
    }

}
