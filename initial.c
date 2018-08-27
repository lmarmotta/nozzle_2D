#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"

/*
 * Apply initial condition.
 */

void apply_initial_condition(t_define p_setup, t_points ** pnts){

    /* Separate the limits of the mesh. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    double P_t = p_setup.BCIN_pt;
    double T_t = p_setup.BCIN_tt;
    double e_i = p_setup.F_Cv*T_t;

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){

            pnts[i][j].q_hat[0] = pnts[i][j].J1 * ( (P_t*pow(1.0,(p_setup.gamma/(p_setup.gamma-1.0))))/(p_setup.F_R*T_t) );
            pnts[i][j].q_hat[1] = 0.0;
            pnts[i][j].q_hat[2] = 0.0;
            pnts[i][j].q_hat[3] = pnts[i][j].q_hat[0]*e_i;

        }
    }
}

/*
 * Make the calculations non-dimensional.
 */

void dim2nondim(t_define p_setup, t_points ** pnts){

    /* Separate the limits of the mesh. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    double P_t = p_setup.BCIN_pt;
    double T_t = p_setup.BCIN_tt;

    double rho_0 = (P_t*pow(1.0,(p_setup.gamma/(p_setup.gamma-1.0))))/(p_setup.F_R*T_t);
    double a_star = sqrt(2.0*p_setup.gamma*( (p_setup.gamma - 1.0)/(p_setup.gamma + 1.0) )*p_setup.F_Cv*T_t);

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){

            /* Compute the primitives. */

            double rho = (pnts[i][j].J * pnts[i][j].q_hat[0])/rho_0;
            double u = (pnts[i][j].q_hat[1]/pnts[i][j].q_hat[0])/a_star;
            double v = (pnts[i][j].q_hat[2]/pnts[i][j].q_hat[0])/a_star;
            double e = (pnts[i][j].J * pnts[i][j].q_hat[3])/(rho*pow(a_star,2.0));

            /* Feed the structs with proper transformation. */

            pnts[i][j].q_hat[0] = pnts[i][j].J1 * rho;
            pnts[i][j].q_hat[1] = pnts[i][j].q_hat[0]*u; 
            pnts[i][j].q_hat[2] = pnts[i][j].q_hat[0]*v; 
            pnts[i][j].q_hat[3] = pnts[i][j].J1 * e;

        }
    }
}

/*
 * Make the calculations dimensional.
 */

void nondim2dim(t_define p_setup, t_points ** pnts){

    /* Separate the limits of the mesh. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    double P_t = p_setup.BCIN_pt;
    double T_t = p_setup.BCIN_tt;

    double rho_0 = (P_t*pow(1.0,(p_setup.gamma/(p_setup.gamma-1.0))))/(p_setup.F_R*T_t);
    double a_star = sqrt(2.0*p_setup.gamma*( (p_setup.gamma - 1.0)/(p_setup.gamma + 1.0) )*p_setup.F_Cv*T_t);

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){

            /* Compute the primitives. */

            double rho = (pnts[i][j].J * pnts[i][j].q_hat[0])*rho_0;
            double u   = (pnts[i][j].q_hat[1]/pnts[i][j].q_hat[0])*a_star;
            double v   = (pnts[i][j].q_hat[2]/pnts[i][j].q_hat[0])*a_star;
            double e   = (pnts[i][j].J * pnts[i][j].q_hat[3])*(rho*pow(a_star,2.0));

            /* Feed the structs with proper transformation. */

            pnts[i][j].q_hat[0] = pnts[i][j].J1 * rho;
            pnts[i][j].q_hat[1] = pnts[i][j].q_hat[0]*u; 
            pnts[i][j].q_hat[2] = pnts[i][j].q_hat[0]*v; 
            pnts[i][j].q_hat[3] = pnts[i][j].J1 * e;

        }
    }
}
