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

        j = jmax-1;

        pnts[i][j].q[0] = pnts[i][j-2].q[0];
        pnts[i][j].q[1] = pnts[i][j-2].q[1];
        pnts[i][j].q[2] = pnts[i][j-2].q[2];
        pnts[i][j].q[3] = pnts[i][j-2].q[3];

    }

    /* Now, apply the inlet boundary condition. */

    for (j = 0; j < jmax; j++){

        i = 0;

        /* Separate the primitives. */

        double u = pnts[i][j].q[1]/pnts[i][j].q[0];
        double v = pnts[i][j].q[2]/pnts[i][j].q[0];

        double u_mag = sqrt( pow(u,2.0) + pow(v,2.0) );

        /* Separate the critical sound speed. */

        double a_ssqr = 2.0 * p_setup.gamma * ( (p_setup.gamma - 1.0)/(p_setup.gamma + 1.0) ) * p_setup.F_Cv * p_setup.BCIN_tt;

        /* Compute the pressure based on the total pressure definition. */
        
        double p = p_setup.BCIN_pt * ( 1.0 - ((p_setup.gamma - 1.0)/(p_setup.gamma + 1.0)) * u_mag/a_ssqr); 

        /* Almost forgot ! */

        p = pow(p,(p_setup.gamma/p_setup.gamma - 1.0));

        /* Update the cartesian Q. */

        pnts[i][j].q[0] = pnts[i][j].q[0];
        pnts[i][j].q[1] = pnts[i][j].q[0]*u;
        pnts[i][j].q[2] = pnts[i][j].q[0]*v;
        pnts[i][j].q[3] = ( p / (p_setup.gamma - 1.0) ) + 0.5 * pnts[i][j].q[0] * ( pow(u,2.0) + pow(v,2.0) );

        /* Update the transformed Q_hat. */

        pnts[i][j].q_hat[0] = pnts[i][j].jm1 * pnts[i][j].q[0];
        pnts[i][j].q_hat[1] = pnts[i][j].jm1 * pnts[i][j].q[0]*u;
        pnts[i][j].q_hat[2] = pnts[i][j].jm1 * pnts[i][j].q[0]*v;
        pnts[i][j].q_hat[3] = pnts[i][j].jm1 * ( p / (p_setup.gamma - 1.0) ) + 0.5 * pnts[i][j].q[0] * ( pow(u,2.0) + pow(v,2.0) );

    }

    /* Now, apply the wall boundary condition. Look... This boundary condition
     * came from a lot of talks with Azevedo.*/

    for (i = 0; i < imax; i++){

        j = 0; 

        /* Separate the primitive variables. */

        double rho = pnts[i][j].q[0]; 
        double e   = pnts[i][j].q[3]; 

        /* Zero out the V boundary condition. */

        pnts[i][j].cov_v = 0.0;

        /* Now extrapolate the covariant velocity components. */

        double u_jp1 = pnts[i][j+1].q[1]/pnts[i][j+1].q[1];
        double v_jp1 = pnts[i][j+1].q[2]/pnts[i][j+1].q[1];

        pnts[i][j].cov_u = pnts[i][j+1].ksi_x*u_jp1 + pnts[i][j+1].ksi_y*v_jp1;

        /* Now, compute the pressure in order to reconstruct. */

        double p_jp1 = (p_setup.gamma - 1.0)*(e - 0.5*rho*( pow(u_jp1,2.0) + pow(v_jp1,2.0)) );

        /* Update the cartesian Q. */

        pnts[i][j].q[0] = pnts[i][j+1].q[0];
        pnts[i][j].q[1] = pnts[i][j+1].q[0]*u_jp1;
        pnts[i][j].q[2] = pnts[i][j+1].q[0]*v_jp1;
        pnts[i][j].q[3] = ( p_jp1 / (p_setup.gamma - 1.0) ) + 0.5 * pnts[i][j+1].q[0] * ( pow(u_jp1,2.0) + pow(v_jp1,2.0) );

        /* Update the transformed Q_hat. */

        pnts[i][j].q_hat[0] = pnts[i][j].jm1 * pnts[i][j+1].q[0];
        pnts[i][j].q_hat[1] = pnts[i][j].jm1 * pnts[i][j+1].q[0]*u_jp1;
        pnts[i][j].q_hat[2] = pnts[i][j].jm1 * pnts[i][j+1].q[0]*v_jp1;
        pnts[i][j].q_hat[3] = pnts[i][j].jm1 * ( p_jp1 / (p_setup.gamma - 1.0) ) + 0.5 * pnts[i][j+1].q[0] * ( pow(u_jp1,2.0) + pow(v_jp1,2.0) );

    }

    /* Apply the outlet boundary condition. */

    for (j = 0; j < jmax; j++){

        i = imax-1;

        /* Separate the primitives. */

        double u = pnts[i][j].q[1]/pnts[i][j].q[0];
        double v = pnts[i][j].q[0]/pnts[i][j].q[0];

        /* Compute the exit boundary condition. Here I am considering, for now,
         * that the flow is exiting subsonic. */
        
        double p = p_setup.BCIN_pt/3.0;

        /* Update the cartesian Q. */

        pnts[i][j].q[0] = pnts[i][j].q[0];
        pnts[i][j].q[1] = pnts[i][j].q[0]*u;
        pnts[i][j].q[2] = pnts[i][j].q[0]*v;
        pnts[i][j].q[3] = ( p / (p_setup.gamma - 1.0) ) + 0.5 * pnts[i][j].q[0] * ( pow(u,2.0) + pow(v,2.0) );

        /* Update the transformed Q_hat. */

        pnts[i][j].q_hat[0] = pnts[i][j].jm1 * pnts[i][j].q[0];
        pnts[i][j].q_hat[1] = pnts[i][j].jm1 * pnts[i][j].q[0]*u;
        pnts[i][j].q_hat[2] = pnts[i][j].jm1 * pnts[i][j].q[0]*v;
        pnts[i][j].q_hat[3] = pnts[i][j].jm1 * ( p / (p_setup.gamma - 1.0) ) + 0.5 * pnts[i][j].q[0] * ( pow(u,2.0) + pow(v,2.0) );

        /* Here the boundary condition of the exit goes. */

    }

}
