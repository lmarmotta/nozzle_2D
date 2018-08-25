#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"

/*
 * This function computes the boundary conditions. Note that the properties
 * inside the boundary calculations shall be performed with the real properties
 * which are re-transoformed after during the q_hat reconstruction. 
 */

void boundary_condition_euler(t_define p_setup, t_points ** pnts){

    /* Separate bounds of the field. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Now, apply the inlet boundary condition. */

    for (int j = 0; j < jmax-1; j++){

        int i = 0;

        double theta = 0.0;

        /* Bizu properties. */

        double P_t = p_setup.BCIN_pt;
        double T_t = p_setup.BCIN_tt;

        /* Compute the stagnation properties. */

        double u = pnts[i+1][j].q_hat[1]/pnts[i+1][j].q_hat[0];

        double v = u*tan(theta); // pi

        double a_star = sqrt(2.0*p_setup.gamma*( (p_setup.gamma - 1.0)/(p_setup.gamma + 1.0) )*p_setup.F_Cv*T_t);
        
        double aux = 1.0 - ( ( (p_setup.gamma - 1.0)/(p_setup.gamma + 1.0) ) * ( (pow(u,2.0) + pow(v,2.0))/pow(a_star,2.0) ) );

        double T = T_t*aux;

        double p = P_t*pow(aux,(p_setup.gamma/(p_setup.gamma-1.0)));

        /* Update the Q. */

        pnts[i][j].q_hat[0] = pnts[i][j].J1 * (p/(p_setup.F_R*T));
        pnts[i][j].q_hat[1] = pnts[i][j].q_hat[0]*u;
        pnts[i][j].q_hat[2] = pnts[i][j].q_hat[0]*v;
        pnts[i][j].q_hat[3] = pnts[i][j].q_hat[0]*( (p_setup.F_Cv*T) + 0.5*( pow(u,2.0) + pow(v,2.0) ) );

    }

    /* Now, apply the wall boundary condition. */

    for (int i = 1; i < imax-1; i++){

        int j = 0; 

        /* Obtain the Cartesian velocity components. */

        double U_i1 = pnts[i][j+1].cov_u;
        double V_i1 = 0.0;

        double u = pnts[i][j].x_ksi*U_i1 + pnts[i][j].x_eta*V_i1;
        double v = pnts[i][j].y_ksi*U_i1 + pnts[i][j].y_eta*V_i1;

        /* Now, compute the energy and pressure in order to reconstruct. */

        double e_jp1   = pnts[i][j+1].J * pnts[i][j+1].q_hat[3]; 

        double rho_jp1 = pnts[i][j+1].J * pnts[i][j+1].q_hat[0]; 

        double u_jp1 = pnts[i][j+1].q_hat[1]/pnts[i][j+1].q_hat[0];
        double v_jp1 = pnts[i][j+1].q_hat[2]/pnts[i][j+1].q_hat[0];

        double p_jp1   = (p_setup.gamma - 1.0)*(e_jp1 - 0.5*rho_jp1*( pow(u_jp1,2.0) + pow(v_jp1,2.0)) );

        /* Update the cartesian Q. */

        pnts[i][j].q_hat[0] = pnts[i][j].J1 * rho_jp1;
        pnts[i][j].q_hat[1] = pnts[i][j].q_hat[0]*u;
        pnts[i][j].q_hat[2] = pnts[i][j].q_hat[0]*v;
        pnts[i][j].q_hat[3] = pnts[i][j].J1 * ( ( p_jp1 / (p_setup.gamma - 1.0) ) + 0.5 * rho_jp1 * ( pow(u,2.0) + pow(v,2.0) ) );

    }

    /* Apply the outlet boundary condition. */

    for (int j = 0; j < jmax-1; j++){

        int i = imax-1;

        /* Separate the primitives. */

        double u = pnts[i-1][j].q_hat[1]/pnts[i-1][j].q_hat[0];
        double v = pnts[i-1][j].q_hat[2]/pnts[i-1][j].q_hat[0];

        /* Compute the magnitude velocity. */

        double u_mag = sqrt( pow(u,2.0) + pow(v,2.0) );

        double e     = pnts[i-1][j].J * pnts[i-1][j].q_hat[3];

        double rho   = pnts[i-1][j].J * pnts[i-1][j].q_hat[0];

        double p     = (p_setup.gamma - 1.0)*(e - 0.5*rho*( pow(u,2.0) + pow(v,2.0)) );

        double a     = sqrt( (p_setup.gamma*p)/ rho);

        double mach  = u_mag/a;

        /* Separate the subsonic from supersonic. */

        if (mach < 1.0){

            /* Compute the exit boundary condition. Here I am considering, for now,
             * that the flow is exiting subsonic. */
            
            double pp = p_setup.BCIN_pt/10.0;

            /* Update the cartesian Q. */

            pnts[i][j].q_hat[0] = pnts[i][j].J1 * (pnts[i-1][j].J * pnts[i-1][j].q_hat[0]);
            pnts[i][j].q_hat[1] = pnts[i][j].q_hat[0]*u;
            pnts[i][j].q_hat[2] = pnts[i][j].q_hat[0]*v;
            pnts[i][j].q_hat[3] = pnts[i][j].J1 * ( ( pp / (p_setup.gamma - 1.0) ) + 0.5 * rho * ( pow(u,2.0) + pow(v,2.0) ) );

        /* Now, deal with supersonic case. */

        } else {

            /* Update the cartesian Q. */

            pnts[i][j].q_hat[0] = pnts[i][j].J1 * (pnts[i-1][j].J * pnts[i-1][j].q_hat[0]); 
            pnts[i][j].q_hat[1] = pnts[i][j].J1 * (pnts[i-1][j].J * pnts[i-1][j].q_hat[1]);
            pnts[i][j].q_hat[2] = pnts[i][j].J1 * (pnts[i-1][j].J * pnts[i-1][j].q_hat[2]);
            pnts[i][j].q_hat[3] = pnts[i][j].J1 * (pnts[i-1][j].J * pnts[i-1][j].q_hat[3]);

        }
    }

    /* Now, apply the symmetry B.C. */

    for (int i = 0; i < imax; i++){

        int j = jmax-1;

        pnts[i][j].q_hat[0] =   pnts[i][j].J1 * (pnts[i][j-2].J * pnts[i][j-2].q_hat[0]);
        pnts[i][j].q_hat[1] =   pnts[i][j].J1 * (pnts[i][j-2].J * pnts[i][j-2].q_hat[1]);
        pnts[i][j].q_hat[2] = - pnts[i][j].J1 * (pnts[i][j-2].J * pnts[i][j-2].q_hat[2]);
        pnts[i][j].q_hat[3] =   pnts[i][j].J1 * (pnts[i][j-2].J * pnts[i][j-2].q_hat[3]);

    }
}
