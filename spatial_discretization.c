#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "cgnslib.h"
#include "structs.h"
#include "externs.h"

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

/*
 * Apply initial condition.
 */

void apply_initial_condition(t_define p_setup, t_points ** pnts){

    /* Separate the limits of the mesh. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Apply a basic initial condition. */

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){

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

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Build the basic fluxes without transformation. */

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){

            /* Limit the scope of these variables. */

            double rho, u, v, e;

            /* Separate the properties we need. */

            rho = pnts[i][j].q[0];
            u   = pnts[i][j].q[1] / pnts[i][j].q[0];
            v   = pnts[i][j].q[2] / pnts[i][j].q[0];
            e   = pnts[i][j].q[3];

            /* Compute pressure. */

            pnts[i][j].p = (p_setup.gamma - 1.0)*(e - 0.5*rho*( pow(u,2.0) + pow(v,2.0)) );

            /* Compute the speed of sound. */

            pnts[i][j].a = sqrt( (p_setup.gamma*pnts[i][j].p)/rho );

            /* Compute the Mach number. */

            pnts[i][j].m = sqrt( pow(u,2.0) + pow(v,2.0) )/pnts[i][j].a;

            /* Compute the Mach number. */

            pnts[i][j].t = (e/rho - 0.5*(pow(u,2.0) + pow(v,2)))/p_setup.F_Cv;

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
            pnts[i][j].e_hat[1] = pnts[i][j].jm1 * (rho*u*pnts[i][j].cov_u + pnts[i][j].ksi_x*pnts[i][j].p);
            pnts[i][j].e_hat[2] = pnts[i][j].jm1 * (rho*v*pnts[i][j].cov_u + pnts[i][j].ksi_y*pnts[i][j].p);
            pnts[i][j].e_hat[3] = pnts[i][j].jm1 * (pnts[i][j].cov_u*(e + pnts[i][j].p));

            /* Build the F fluxes in curvilinear coordinates. */

            pnts[i][j].f_hat[0] = pnts[i][j].jm1 * (rho*pnts[i][j].cov_v);
            pnts[i][j].f_hat[1] = pnts[i][j].jm1 * (rho*u*pnts[i][j].cov_v + pnts[i][j].eta_x*pnts[i][j].p);
            pnts[i][j].f_hat[2] = pnts[i][j].jm1 * (rho*v*pnts[i][j].cov_v + pnts[i][j].eta_y*pnts[i][j].p);
            pnts[i][j].f_hat[3] = pnts[i][j].jm1 * (pnts[i][j].cov_v*(e + pnts[i][j].p));

        }
    }
}

void compute_rhs(t_define p_setup, t_points ** pnts){

    /* Separate bounds of the field. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Initialize the residue. */

    max_rhs_rho  = -99999999.0;
    max_rhs_rhou = -99999999.0;
    max_rhs_rhov = -99999999.0;
    max_rhs_e    = -99999999.0;

    /* Compute the RHS for the internal points. Remenber that the boundaries
     * are computed by the boundary conditions. */

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            double d_Eh0, d_Eh1, d_Eh2, d_Eh3;
            double d_Fh0, d_Fh1, d_Fh2, d_Fh3;

            /* Compute the E fluxes. */

            d_Eh0 = ( pnts[i+1][j].e_hat[0] - pnts[i-1][j].e_hat[0] ) / 2.0;
            d_Eh1 = ( pnts[i+1][j].e_hat[1] - pnts[i-1][j].e_hat[1] ) / 2.0;
            d_Eh2 = ( pnts[i+1][j].e_hat[2] - pnts[i-1][j].e_hat[2] ) / 2.0;
            d_Eh3 = ( pnts[i+1][j].e_hat[3] - pnts[i-1][j].e_hat[3] ) / 2.0;

            /* Compute the F fluxes. */

            d_Fh0 = ( pnts[i][j+1].f_hat[0] - pnts[i][j-1].f_hat[0] ) / 2.0; 
            d_Fh1 = ( pnts[i][j+1].f_hat[1] - pnts[i][j-1].f_hat[1] ) / 2.0; 
            d_Fh2 = ( pnts[i][j+1].f_hat[2] - pnts[i][j-1].f_hat[2] ) / 2.0; 
            d_Fh3 = ( pnts[i][j+1].f_hat[3] - pnts[i][j-1].f_hat[3] ) / 2.0; 

            /* Store the RHS properly. */

            pnts[i][j].RHS[0] = d_Eh0 + d_Fh0;
            pnts[i][j].RHS[1] = d_Eh1 + d_Fh1;
            pnts[i][j].RHS[2] = d_Eh2 + d_Fh2;
            pnts[i][j].RHS[3] = d_Eh3 + d_Fh3;

            /* Store the max residue. */

            if ( fabs( pnts[i][j].RHS[0]) > max_rhs_rho )  max_rhs_rho  = log10(fabs(pnts[i][j].RHS[0]));
            if ( fabs( pnts[i][j].RHS[1]) > max_rhs_rhou)  max_rhs_rhou = log10(fabs(pnts[i][j].RHS[1]));
            if ( fabs( pnts[i][j].RHS[2]) > max_rhs_rhov)  max_rhs_rhov = log10(fabs(pnts[i][j].RHS[2]));
            if ( fabs( pnts[i][j].RHS[3]) > max_rhs_e   )  max_rhs_e    = log10(fabs(pnts[i][j].RHS[3]));

        }
    }
}

void art_dissip_2nd(t_define p_setup, t_points ** pnts){

    /* Separate bounds. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Allocate the dissip vector. */

    double diss_ksi[imax][jmax][4]; 
    double diss_eta[imax][jmax][4]; 

    /* Compute the internal points. */

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            /* Ksi direction. */

            diss_ksi[i][j][0] = pnts[i+1][j].q_hat[0] - 2.0 * pnts[i][j].q_hat[0] + pnts[i-1][j].q_hat[0];
            diss_ksi[i][j][1] = pnts[i+1][j].q_hat[1] - 2.0 * pnts[i][j].q_hat[1] + pnts[i-1][j].q_hat[1];
            diss_ksi[i][j][2] = pnts[i+1][j].q_hat[2] - 2.0 * pnts[i][j].q_hat[2] + pnts[i-1][j].q_hat[2];
            diss_ksi[i][j][3] = pnts[i+1][j].q_hat[3] - 2.0 * pnts[i][j].q_hat[3] + pnts[i-1][j].q_hat[3];

            /* Eta direction. */

            diss_eta[i][j][0] = pnts[i][j+1].q_hat[0] - 2.0 * pnts[i][j].q_hat[0] + pnts[i][j-1].q_hat[0];
            diss_eta[i][j][1] = pnts[i][j+1].q_hat[1] - 2.0 * pnts[i][j].q_hat[1] + pnts[i][j-1].q_hat[1];
            diss_eta[i][j][2] = pnts[i][j+1].q_hat[2] - 2.0 * pnts[i][j].q_hat[2] + pnts[i][j-1].q_hat[2];
            diss_eta[i][j][3] = pnts[i][j+1].q_hat[3] - 2.0 * pnts[i][j].q_hat[3] + pnts[i][j-1].q_hat[3];

        }
    }

    /* Dow, add to the residue. */

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            pnts[i][j].RHS[0] = pnts[i][j].RHS[0] - (diss_ksi[i][j][0]*pnts[i][j].RHS[0] + diss_eta[i][j][0]*pnts[i][j].RHS[0]);
            pnts[i][j].RHS[1] = pnts[i][j].RHS[1] - (diss_ksi[i][j][1]*pnts[i][j].RHS[1] + diss_eta[i][j][1]*pnts[i][j].RHS[1]);
            pnts[i][j].RHS[2] = pnts[i][j].RHS[2] - (diss_ksi[i][j][2]*pnts[i][j].RHS[2] + diss_eta[i][j][2]*pnts[i][j].RHS[2]);
            pnts[i][j].RHS[3] = pnts[i][j].RHS[3] - (diss_ksi[i][j][3]*pnts[i][j].RHS[3] + diss_eta[i][j][3]*pnts[i][j].RHS[3]);

        }
    }

}
