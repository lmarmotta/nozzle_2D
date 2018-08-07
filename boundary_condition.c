#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"

void boundary_condition(t_define p_setup, t_points ** pnts){

    /* Separate bounds of the field. */

    int i, j;

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Now, apply the symmetry B.C. */

    for (i = 1; i < imax-1; i++){

        pnts[i][jmax-1].q[0] = pnts[i][jmax-3].q[0];
        pnts[i][jmax-1].q[1] = pnts[i][jmax-3].q[1];
        pnts[i][jmax-1].q[2] = pnts[i][jmax-3].q[2];
        pnts[i][jmax-1].q[3] = pnts[i][jmax-3].q[3];

    }

    /* Now, apply the inlet boundary condition. */

    for (j = 0; j < jmax; j++){

        pnts[0][j].q[1] =  p_setup.BCIN_u*pnts[0][j].q[0]; 
        pnts[0][j].q[2] =  p_setup.BCIN_v*pnts[0][j].q[0]; 
        pnts[0][j].q[3] = (p_setup.BCIN_p/p_setup.gamma) + 0.5*pnts[0][j].q[0]*(pow(p_setup.BCIN_u,2.0) + pow(p_setup.BCIN_v,2.0));

    }

    /* Now, apply the wall boundary condition. */

    for (i = 0; i < imax; i++){

        double u, v;

        /* Tangent boundary condition. */

    }

    /* Apply the outlet boundary condition. */

    for (j = 1; j < jmax-1; j++){

        /* Here the boundary condition of the exit goes. */

    }

}
