#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"

void boundary_condition(t_define p_setup, t_points ** pnts){

    /* Separate bounds of the field. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Now, apply the symmetry B.C. */

    for (int i = 1; i < imax-1; i++){

        int j = jmax-1;

        pnts[i][j].q[0] = pnts[i][j-2].q[0];
        pnts[i][j].q[1] = pnts[i][j-2].q[1];
        pnts[i][j].q[2] = pnts[i][j-2].q[2];
        pnts[i][j].q[3] = pnts[i][j-2].q[3];

    }

    /* Now, apply the inlet boundary condition. */

    for (int j = 0; j < jmax; j++){

        int i = 0;

        /* Sepaate the primitives. */

        double u = pnts[i][j].q[1]/pnts[i][j].q[0];
        double v = pnts[i][j].q[2]/pnts[i][j].q[0];

        /* Separate bizu boundary conditions. */

        double gamm1 = p_setup.gamma - 1.0;

        /* Compute the velocity magnitude. */

        double velf = sqrt( pow(u,2.0) + pow(v,2.0) );

        /* Compute the squared critical speed of sound. */

        double acrsq = 2.0*(p_setup.gamma * p_setup.F_R * p_setup.BCIN_pt)/(p_setup.gamma + 1.0);

        /* Compute the static pressure. */

        double auxcte1 = gamm1/(p_setup.gamma + 1.0);
        double auxcte2 = p_setup.gamma/gamm1;

        double aux1  = 1.0 - (auxcte1 * ( pow(velf,2.0) )/acrsq);
        double eintf = p_setup.F_Cv * p_setup.BCIN_pt * aux1;
        double pf    = p_setup.BCIN_pt * (pow(aux1,auxcte2));

        double rhof = pf / (gamm1*eintf);
        double ef   = rhof*( eintf + 0.5*( pow(velf,2.0) ) );

        /* Update the Q. */

        pnts[i][j].q[0] = rhof;
        pnts[i][j].q[1] = rhof*u*p_setup.BCIN_udir;
        pnts[i][j].q[2] = rhof*v*p_setup.BCIN_vdir;
        pnts[i][j].q[3] = ef;

        /* Update the transformed Q_hat. */

        pnts[i][j].q_hat[0] = pnts[i][j].jm1 * (rhof);
        pnts[i][j].q_hat[1] = pnts[i][j].jm1 * (rhof*u*p_setup.BCIN_udir);
        pnts[i][j].q_hat[2] = pnts[i][j].jm1 * (rhof*v*p_setup.BCIN_vdir);
        pnts[i][j].q_hat[3] = pnts[i][j].jm1 * (ef);

    }

    /* Now, apply the wall boundary condition. Look... This boundary condition
     * came from a lot of talks with Azevedo.*/

    for (int i = 0; i < imax; i++){

        int j = 0; 

        /* Zero out the V boundary condition. */

        pnts[i][j].cov_v = 0.0;

        /* Now extrapolate the covariant velocity components. */

        double u_jp1 = pnts[i][j+1].q[1]/pnts[i][j+1].q[0];
        double v_jp1 = pnts[i][j+1].q[2]/pnts[i][j+1].q[0];

        pnts[i][j].cov_u = pnts[i][j+1].ksi_x*u_jp1 + pnts[i][j+1].ksi_y*v_jp1;

        /* Now, compute the pressure in order to reconstruct. */

        double e_jp1   = pnts[i][j+1].q[3]; 

        double rho_jp1 = pnts[i][j].q[0]; 

        double p_jp1 = (p_setup.gamma - 1.0)*(e_jp1 - 0.5*rho_jp1*( pow(u_jp1,2.0) + pow(v_jp1,2.0)) );

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

    for (int j = 0; j < jmax; j++){

        int i = imax-1;

        /* Separate the primitives. */

        double u = pnts[i][j].q[1]/pnts[i][j].q[0];
        double v = pnts[i][j].q[2]/pnts[i][j].q[0];

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

    }
}
