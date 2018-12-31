#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "structs.h"
#include "externs.h"
#include "prototypes.h"

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

/*
 * Build transformed fluxes.
 */

void compute_fluxes(t_define p_setup, t_points ** pnts){

    /* Separate the limits of the mesh. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Build the basic fluxes without transformation. */

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){

            /* Separate the properties we need. */

            double rho = pnts[i][j].J * pnts[i][j].q_hat[0];
            double u   = pnts[i][j].q_hat[1] / pnts[i][j].q_hat[0];
            double v   = pnts[i][j].q_hat[2] / pnts[i][j].q_hat[0];
            double e   = pnts[i][j].J * pnts[i][j].q_hat[3];

            /* Compute pressure. */

            pnts[i][j].p = (p_setup.gamma - 1.0)*(e - 0.5*rho*( pow(u,2.0) + pow(v,2.0)) );

            /* Compute the speed of sound. */

            pnts[i][j].a = sqrt( (p_setup.gamma*pnts[i][j].p)/rho );

            /* Compute the Mach number. */

            pnts[i][j].m = sqrt( pow(u,2.0) + pow(v,2.0) )/pnts[i][j].a;

            /* Compute the Temperature. */

            pnts[i][j].t = ( (e/rho) - 0.5*(pow(u,2.0) + pow(v,2.0)) )/p_setup.F_Cv;

            /* Compute the covariant velocity components. */

            pnts[i][j].cov_u =  pnts[i][j].ksi_x*u + pnts[i][j].ksi_y*v;
            pnts[i][j].cov_v =  pnts[i][j].eta_x*u + pnts[i][j].eta_y*v;

            if (j == 0) pnts[i][j].cov_v = 0.0;

            /* Build the Q fluxes in curvilinear coordinates. */

            pnts[i][j].q_hat[0] = pnts[i][j].J1*rho;
            pnts[i][j].q_hat[1] = pnts[i][j].J1*rho*u;
            pnts[i][j].q_hat[2] = pnts[i][j].J1*rho*v;
            pnts[i][j].q_hat[3] = pnts[i][j].J1*e;

            /* Build the E fluxes in curvilinear coordinates. */

            pnts[i][j].e_hat[0] = pnts[i][j].J1 * (rho*pnts[i][j].cov_u);
            pnts[i][j].e_hat[1] = pnts[i][j].J1 * (rho*u*pnts[i][j].cov_u + pnts[i][j].ksi_x*pnts[i][j].p);
            pnts[i][j].e_hat[2] = pnts[i][j].J1 * (rho*v*pnts[i][j].cov_u + pnts[i][j].ksi_y*pnts[i][j].p);
            pnts[i][j].e_hat[3] = pnts[i][j].J1 * (pnts[i][j].cov_u*(e + pnts[i][j].p));

            /* Build the F fluxes in curvilinear coordinates. */

            pnts[i][j].f_hat[0] = pnts[i][j].J1 * (rho*pnts[i][j].cov_v);
            pnts[i][j].f_hat[1] = pnts[i][j].J1 * (rho*u*pnts[i][j].cov_v + pnts[i][j].eta_x*pnts[i][j].p);
            pnts[i][j].f_hat[2] = pnts[i][j].J1 * (rho*v*pnts[i][j].cov_v + pnts[i][j].eta_y*pnts[i][j].p);
            pnts[i][j].f_hat[3] = pnts[i][j].J1 * (pnts[i][j].cov_v*(e + pnts[i][j].p));

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

    /* Compute the residue for the internal points. */ 

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

            if (fabs( pnts[i][j].RHS[0]) > max_rhs_rho ) max_rhs_rho  = fabs(pnts[i][j].RHS[0]+DBL_EPSILON);
            if (fabs( pnts[i][j].RHS[1]) > max_rhs_rhou) max_rhs_rhou = fabs(pnts[i][j].RHS[1]+DBL_EPSILON);
            if (fabs( pnts[i][j].RHS[2]) > max_rhs_rhov) max_rhs_rhov = fabs(pnts[i][j].RHS[2]+DBL_EPSILON);
            if (fabs( pnts[i][j].RHS[3]) > max_rhs_e   ) max_rhs_e    = fabs(pnts[i][j].RHS[3]+DBL_EPSILON);

        }
    }
}

void art_dissip_2nd(t_define p_setup, t_points ** pnts){

    /* Separate bounds. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Compute the internal points. */

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            /* Ksi direction. */

            pnts[i][j].diss_ksi[0] = (pnts[i+1][j].J*pnts[i+1][j].q_hat[0]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[0]) + (pnts[i-1][j].J*pnts[i-1][j].q_hat[0]);
            pnts[i][j].diss_ksi[1] = (pnts[i+1][j].J*pnts[i+1][j].q_hat[1]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[1]) + (pnts[i-1][j].J*pnts[i-1][j].q_hat[1]);
            pnts[i][j].diss_ksi[2] = (pnts[i+1][j].J*pnts[i+1][j].q_hat[2]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[2]) + (pnts[i-1][j].J*pnts[i-1][j].q_hat[2]);
            pnts[i][j].diss_ksi[3] = (pnts[i+1][j].J*pnts[i+1][j].q_hat[3]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[3]) + (pnts[i-1][j].J*pnts[i-1][j].q_hat[3]);

            /* Eta direction. */

            pnts[i][j].diss_eta[0] = (pnts[i][j+1].J*pnts[i][j+1].q_hat[0]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[0]) + (pnts[i][j-1].J*pnts[i][j-1].q_hat[0]);
            pnts[i][j].diss_eta[1] = (pnts[i][j+1].J*pnts[i][j+1].q_hat[1]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[1]) + (pnts[i][j-1].J*pnts[i][j-1].q_hat[1]);
            pnts[i][j].diss_eta[2] = (pnts[i][j+1].J*pnts[i][j+1].q_hat[2]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[2]) + (pnts[i][j-1].J*pnts[i][j-1].q_hat[2]);
            pnts[i][j].diss_eta[3] = (pnts[i][j+1].J*pnts[i][j+1].q_hat[3]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[3]) + (pnts[i][j-1].J*pnts[i][j-1].q_hat[3]);

        }
    }

    /* Now, add the second difference dissipation to the residue. */

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            pnts[i][j].RHS[0] = pnts[i][j].RHS[0] - pnts[i][j].J1 * ( (p_setup.dissp2)*(pnts[i][j].diss_ksi[0] + pnts[i][j].diss_eta[0]) );
            pnts[i][j].RHS[1] = pnts[i][j].RHS[1] - pnts[i][j].J1 * ( (p_setup.dissp2)*(pnts[i][j].diss_ksi[1] + pnts[i][j].diss_eta[1]) );
            pnts[i][j].RHS[2] = pnts[i][j].RHS[2] - pnts[i][j].J1 * ( (p_setup.dissp2)*(pnts[i][j].diss_ksi[2] + pnts[i][j].diss_eta[2]) );
            pnts[i][j].RHS[3] = pnts[i][j].RHS[3] - pnts[i][j].J1 * ( (p_setup.dissp2)*(pnts[i][j].diss_ksi[3] + pnts[i][j].diss_eta[3]) );

        }
    }
}

void art_dissip_4th(t_define p_setup, t_points ** pnts){

    /* Separate bounds. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Temporalely takes away the coordinate transformation. */

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){

            pnts[i][j].q_hat[0] = pnts[i][j].J * pnts[i][j].q_hat[0];   
            pnts[i][j].q_hat[1] = pnts[i][j].J * pnts[i][j].q_hat[1];
            pnts[i][j].q_hat[2] = pnts[i][j].J * pnts[i][j].q_hat[2];
            pnts[i][j].q_hat[3] = pnts[i][j].J * pnts[i][j].q_hat[3];

        }
    }

    /* Compute the fourth difference dissipation for the internal points. */

    for (int i = 2; i<imax-2; i++){
        for (int j = 2; j<jmax-2; j++){

            /* Compute in ksi direction. */

            pnts[i][j].diss_ksi[0] = pnts[i+2][j].q_hat[0] - 4.0*pnts[i+1][j].q_hat[0] + 6.0*pnts[i][j].q_hat[0] - 4.0*pnts[i-1][j].q_hat[0] + pnts[i-2][j].q_hat[0];
            pnts[i][j].diss_ksi[1] = pnts[i+2][j].q_hat[1] - 4.0*pnts[i+1][j].q_hat[1] + 6.0*pnts[i][j].q_hat[1] - 4.0*pnts[i-1][j].q_hat[1] + pnts[i-2][j].q_hat[1];
            pnts[i][j].diss_ksi[2] = pnts[i+2][j].q_hat[2] - 4.0*pnts[i+1][j].q_hat[2] + 6.0*pnts[i][j].q_hat[2] - 4.0*pnts[i-1][j].q_hat[2] + pnts[i-2][j].q_hat[2];
            pnts[i][j].diss_ksi[3] = pnts[i+2][j].q_hat[3] - 4.0*pnts[i+1][j].q_hat[3] + 6.0*pnts[i][j].q_hat[3] - 4.0*pnts[i-1][j].q_hat[3] + pnts[i-2][j].q_hat[3];

            /* Compute the eta direction. */

            pnts[i][j].diss_eta[0] = pnts[i][j+2].q_hat[0] - 4.0*pnts[i][j+1].q_hat[0] + 6.0*pnts[i][j].q_hat[0] - 4.0*pnts[i][j-1].q_hat[0] + pnts[i][j-2].q_hat[0];
            pnts[i][j].diss_eta[1] = pnts[i][j+2].q_hat[1] - 4.0*pnts[i][j+1].q_hat[1] + 6.0*pnts[i][j].q_hat[1] - 4.0*pnts[i][j-1].q_hat[1] + pnts[i][j-2].q_hat[1];
            pnts[i][j].diss_eta[2] = pnts[i][j+2].q_hat[2] - 4.0*pnts[i][j+1].q_hat[2] + 6.0*pnts[i][j].q_hat[2] - 4.0*pnts[i][j-1].q_hat[2] + pnts[i][j-2].q_hat[2];
            pnts[i][j].diss_eta[3] = pnts[i][j+2].q_hat[3] - 4.0*pnts[i][j+1].q_hat[3] + 6.0*pnts[i][j].q_hat[3] - 4.0*pnts[i][j-1].q_hat[3] + pnts[i][j-2].q_hat[3];

        }
    }

    /* Compute the symmetry line. */

    for (int i = 2; i<imax-2; i++){

        int j = jmax - 2;

        /* Compute in ksi direction. */

        pnts[i][j].diss_ksi[0] = pnts[i+2][j].q_hat[0] - 4.0*pnts[i+1][j].q_hat[0] + 6.0*pnts[i][j].q_hat[0] - 4.0*pnts[i-1][j].q_hat[0] + pnts[i-2][j].q_hat[0];
        pnts[i][j].diss_ksi[1] = pnts[i+2][j].q_hat[1] - 4.0*pnts[i+1][j].q_hat[1] + 6.0*pnts[i][j].q_hat[1] - 4.0*pnts[i-1][j].q_hat[1] + pnts[i-2][j].q_hat[1];
        pnts[i][j].diss_ksi[2] = pnts[i+2][j].q_hat[2] - 4.0*pnts[i+1][j].q_hat[2] + 6.0*pnts[i][j].q_hat[2] - 4.0*pnts[i-1][j].q_hat[2] + pnts[i-2][j].q_hat[2];
        pnts[i][j].diss_ksi[3] = pnts[i+2][j].q_hat[3] - 4.0*pnts[i+1][j].q_hat[3] + 6.0*pnts[i][j].q_hat[3] - 4.0*pnts[i-1][j].q_hat[3] + pnts[i-2][j].q_hat[3];

        /* Compute the eta direction. */

        pnts[i][j].diss_eta[0] = pnts[i][j].q_hat[0] - 4.0*pnts[i][j-1].q_hat[0] + 6.0*pnts[i][j-2].q_hat[0] - 4.0*pnts[i][j-3].q_hat[0] + pnts[i][j-4].q_hat[0];
        pnts[i][j].diss_eta[1] = pnts[i][j].q_hat[1] - 4.0*pnts[i][j-1].q_hat[1] + 6.0*pnts[i][j-2].q_hat[1] - 4.0*pnts[i][j-3].q_hat[1] + pnts[i][j-4].q_hat[1];
        pnts[i][j].diss_eta[2] = pnts[i][j].q_hat[2] - 4.0*pnts[i][j-1].q_hat[2] + 6.0*pnts[i][j-2].q_hat[2] - 4.0*pnts[i][j-3].q_hat[2] + pnts[i][j-4].q_hat[2];
        pnts[i][j].diss_eta[3] = pnts[i][j].q_hat[3] - 4.0*pnts[i][j-1].q_hat[3] + 6.0*pnts[i][j-2].q_hat[3] - 4.0*pnts[i][j-3].q_hat[3] + pnts[i][j-4].q_hat[3];

    }


    /* Compute the lower surface. */

    for (int i = 2; i<imax-2; i++){

        int j = 1;

        /* Compute in ksi direction. */

        pnts[i][j].diss_ksi[0] = pnts[i+2][j].q_hat[0] - 4.0*pnts[i+1][j].q_hat[0] + 6.0*pnts[i][j].q_hat[0] - 4.0*pnts[i-1][j].q_hat[0] + pnts[i-2][j].q_hat[0];
        pnts[i][j].diss_ksi[1] = pnts[i+2][j].q_hat[1] - 4.0*pnts[i+1][j].q_hat[1] + 6.0*pnts[i][j].q_hat[1] - 4.0*pnts[i-1][j].q_hat[1] + pnts[i-2][j].q_hat[1];
        pnts[i][j].diss_ksi[2] = pnts[i+2][j].q_hat[2] - 4.0*pnts[i+1][j].q_hat[2] + 6.0*pnts[i][j].q_hat[2] - 4.0*pnts[i-1][j].q_hat[2] + pnts[i-2][j].q_hat[2];
        pnts[i][j].diss_ksi[3] = pnts[i+2][j].q_hat[3] - 4.0*pnts[i+1][j].q_hat[3] + 6.0*pnts[i][j].q_hat[3] - 4.0*pnts[i-1][j].q_hat[3] + pnts[i-2][j].q_hat[3];

        /* Compute the eta direction. */

        pnts[i][j].diss_eta[0] = pnts[i][j+4].q_hat[0] - 4.0*pnts[i][j+3].q_hat[0] + 6.0*pnts[i][j+2].q_hat[0] - 4.0*pnts[i][j+1].q_hat[0] + pnts[i][j].q_hat[0];
        pnts[i][j].diss_eta[1] = pnts[i][j+4].q_hat[1] - 4.0*pnts[i][j+3].q_hat[1] + 6.0*pnts[i][j+2].q_hat[1] - 4.0*pnts[i][j+1].q_hat[1] + pnts[i][j].q_hat[1];
        pnts[i][j].diss_eta[2] = pnts[i][j+4].q_hat[2] - 4.0*pnts[i][j+3].q_hat[2] + 6.0*pnts[i][j+2].q_hat[2] - 4.0*pnts[i][j+1].q_hat[2] + pnts[i][j].q_hat[2];
        pnts[i][j].diss_eta[3] = pnts[i][j+4].q_hat[3] - 4.0*pnts[i][j+3].q_hat[3] + 6.0*pnts[i][j+2].q_hat[3] - 4.0*pnts[i][j+1].q_hat[3] + pnts[i][j].q_hat[3];

    }

    /* Compute the inlet surface. */

    for (int j = 2; j<jmax-2; j++){

        int i = 1;

        /* Compute in ksi direction. */

        pnts[i][j].diss_ksi[0] = pnts[i+4][j].q_hat[0] - 4.0*pnts[i+3][j].q_hat[0] + 6.0*pnts[i+2][j].q_hat[0] - 4.0*pnts[i+1][j].q_hat[0] + pnts[i][j].q_hat[0];
        pnts[i][j].diss_ksi[1] = pnts[i+4][j].q_hat[1] - 4.0*pnts[i+3][j].q_hat[1] + 6.0*pnts[i+2][j].q_hat[1] - 4.0*pnts[i+1][j].q_hat[1] + pnts[i][j].q_hat[1];
        pnts[i][j].diss_ksi[2] = pnts[i+4][j].q_hat[2] - 4.0*pnts[i+3][j].q_hat[2] + 6.0*pnts[i+2][j].q_hat[2] - 4.0*pnts[i+1][j].q_hat[2] + pnts[i][j].q_hat[2];
        pnts[i][j].diss_ksi[3] = pnts[i+4][j].q_hat[3] - 4.0*pnts[i+3][j].q_hat[3] + 6.0*pnts[i+2][j].q_hat[3] - 4.0*pnts[i+1][j].q_hat[3] + pnts[i][j].q_hat[3];

        /* Compute the eta direction. */

        pnts[i][j].diss_eta[0] = pnts[i][j+2].q_hat[0] - 4.0*pnts[i][j+1].q_hat[0] + 6.0*pnts[i][j].q_hat[0] - 4.0*pnts[i][j-1].q_hat[0] + pnts[i][j-2].q_hat[0];
        pnts[i][j].diss_eta[1] = pnts[i][j+2].q_hat[1] - 4.0*pnts[i][j+1].q_hat[1] + 6.0*pnts[i][j].q_hat[1] - 4.0*pnts[i][j-1].q_hat[1] + pnts[i][j-2].q_hat[1];
        pnts[i][j].diss_eta[2] = pnts[i][j+2].q_hat[2] - 4.0*pnts[i][j+1].q_hat[2] + 6.0*pnts[i][j].q_hat[2] - 4.0*pnts[i][j-1].q_hat[2] + pnts[i][j-2].q_hat[2];
        pnts[i][j].diss_eta[3] = pnts[i][j+2].q_hat[3] - 4.0*pnts[i][j+1].q_hat[3] + 6.0*pnts[i][j].q_hat[3] - 4.0*pnts[i][j-1].q_hat[3] + pnts[i][j-2].q_hat[3];

    }

    /* Compute the outlet surface. */

    for (int j = 2; j<jmax-2; j++){

        int i = imax - 2;

        /* Compute in ksi direction. */

        pnts[i][j].diss_ksi[0] = pnts[i][j].q_hat[0] - 4.0*pnts[i-1][j].q_hat[0] + 6.0*pnts[i-2][j].q_hat[0] - 4.0*pnts[i-3][j].q_hat[0] + pnts[i-4][j].q_hat[0];
        pnts[i][j].diss_ksi[1] = pnts[i][j].q_hat[1] - 4.0*pnts[i-1][j].q_hat[1] + 6.0*pnts[i-2][j].q_hat[1] - 4.0*pnts[i-3][j].q_hat[1] + pnts[i-4][j].q_hat[1];
        pnts[i][j].diss_ksi[2] = pnts[i][j].q_hat[2] - 4.0*pnts[i-1][j].q_hat[2] + 6.0*pnts[i-2][j].q_hat[2] - 4.0*pnts[i-3][j].q_hat[2] + pnts[i-4][j].q_hat[2];
        pnts[i][j].diss_ksi[3] = pnts[i][j].q_hat[3] - 4.0*pnts[i-1][j].q_hat[3] + 6.0*pnts[i-2][j].q_hat[3] - 4.0*pnts[i-3][j].q_hat[3] + pnts[i-4][j].q_hat[3];

        /* Compute the eta direction. */

        pnts[i][j].diss_eta[0] = pnts[i][j+2].q_hat[0] - 4.0*pnts[i][j+1].q_hat[0] + 6.0*pnts[i][j].q_hat[0] - 4.0*pnts[i][j-1].q_hat[0] + pnts[i][j-2].q_hat[0];
        pnts[i][j].diss_eta[1] = pnts[i][j+2].q_hat[1] - 4.0*pnts[i][j+1].q_hat[1] + 6.0*pnts[i][j].q_hat[1] - 4.0*pnts[i][j-1].q_hat[1] + pnts[i][j-2].q_hat[1];
        pnts[i][j].diss_eta[2] = pnts[i][j+2].q_hat[2] - 4.0*pnts[i][j+1].q_hat[2] + 6.0*pnts[i][j].q_hat[2] - 4.0*pnts[i][j-1].q_hat[2] + pnts[i][j-2].q_hat[2];
        pnts[i][j].diss_eta[3] = pnts[i][j+2].q_hat[3] - 4.0*pnts[i][j+1].q_hat[3] + 6.0*pnts[i][j].q_hat[3] - 4.0*pnts[i][j-1].q_hat[3] + pnts[i][j-2].q_hat[3];

    }

    /* Now, lower left corner. */

    int i = 1; int j = 1;

    pnts[i][j].diss_ksi[0] = pnts[i+4][j].q_hat[0] - 4.0*pnts[i+3][j].q_hat[0] + 6.0*pnts[i+2][j].q_hat[0] - 4.0*pnts[i+1][j].q_hat[0] + pnts[i][j].q_hat[0];
    pnts[i][j].diss_ksi[1] = pnts[i+4][j].q_hat[1] - 4.0*pnts[i+3][j].q_hat[1] + 6.0*pnts[i+2][j].q_hat[1] - 4.0*pnts[i+1][j].q_hat[1] + pnts[i][j].q_hat[1];
    pnts[i][j].diss_ksi[2] = pnts[i+4][j].q_hat[2] - 4.0*pnts[i+3][j].q_hat[2] + 6.0*pnts[i+2][j].q_hat[2] - 4.0*pnts[i+1][j].q_hat[2] + pnts[i][j].q_hat[2];
    pnts[i][j].diss_ksi[3] = pnts[i+4][j].q_hat[3] - 4.0*pnts[i+3][j].q_hat[3] + 6.0*pnts[i+2][j].q_hat[3] - 4.0*pnts[i+1][j].q_hat[3] + pnts[i][j].q_hat[3];

    pnts[i][j].diss_eta[0] = pnts[i][j+4].q_hat[0] - 4.0*pnts[i][j+3].q_hat[0] + 6.0*pnts[i][j+2].q_hat[0] - 4.0*pnts[i][j+1].q_hat[0] + pnts[i][j].q_hat[0];
    pnts[i][j].diss_eta[1] = pnts[i][j+4].q_hat[1] - 4.0*pnts[i][j+3].q_hat[1] + 6.0*pnts[i][j+2].q_hat[1] - 4.0*pnts[i][j+1].q_hat[1] + pnts[i][j].q_hat[1];
    pnts[i][j].diss_eta[2] = pnts[i][j+4].q_hat[2] - 4.0*pnts[i][j+3].q_hat[2] + 6.0*pnts[i][j+2].q_hat[2] - 4.0*pnts[i][j+1].q_hat[2] + pnts[i][j].q_hat[2];
    pnts[i][j].diss_eta[3] = pnts[i][j+4].q_hat[3] - 4.0*pnts[i][j+3].q_hat[3] + 6.0*pnts[i][j+2].q_hat[3] - 4.0*pnts[i][j+1].q_hat[3] + pnts[i][j].q_hat[3];

    /* Upper left corner. */

    i = 1; j = jmax - 2;

    pnts[i][j].diss_ksi[0] = pnts[i+4][j].q_hat[0] - 4.0*pnts[i+3][j].q_hat[0] + 6.0*pnts[i+2][j].q_hat[0] - 4.0*pnts[i+1][j].q_hat[0] + pnts[i][j].q_hat[0];
    pnts[i][j].diss_ksi[1] = pnts[i+4][j].q_hat[1] - 4.0*pnts[i+3][j].q_hat[1] + 6.0*pnts[i+2][j].q_hat[1] - 4.0*pnts[i+1][j].q_hat[1] + pnts[i][j].q_hat[1];
    pnts[i][j].diss_ksi[2] = pnts[i+4][j].q_hat[2] - 4.0*pnts[i+3][j].q_hat[2] + 6.0*pnts[i+2][j].q_hat[2] - 4.0*pnts[i+1][j].q_hat[2] + pnts[i][j].q_hat[2];
    pnts[i][j].diss_ksi[3] = pnts[i+4][j].q_hat[3] - 4.0*pnts[i+3][j].q_hat[3] + 6.0*pnts[i+2][j].q_hat[3] - 4.0*pnts[i+1][j].q_hat[3] + pnts[i][j].q_hat[3];

    pnts[i][j].diss_eta[0] = pnts[i][j].q_hat[0] - 4.0*pnts[i][j-1].q_hat[0] + 6.0*pnts[i][j-2].q_hat[0] - 4.0*pnts[i][j-3].q_hat[0] + pnts[i][j-4].q_hat[0];
    pnts[i][j].diss_eta[1] = pnts[i][j].q_hat[1] - 4.0*pnts[i][j-1].q_hat[1] + 6.0*pnts[i][j-2].q_hat[1] - 4.0*pnts[i][j-3].q_hat[1] + pnts[i][j-4].q_hat[1];
    pnts[i][j].diss_eta[2] = pnts[i][j].q_hat[2] - 4.0*pnts[i][j-1].q_hat[2] + 6.0*pnts[i][j-2].q_hat[2] - 4.0*pnts[i][j-3].q_hat[2] + pnts[i][j-4].q_hat[2];
    pnts[i][j].diss_eta[3] = pnts[i][j].q_hat[3] - 4.0*pnts[i][j-1].q_hat[3] + 6.0*pnts[i][j-2].q_hat[3] - 4.0*pnts[i][j-3].q_hat[3] + pnts[i][j-4].q_hat[3];

    /* Lower right corner. */

    i = imax - 2; j = 1;

    pnts[i][j].diss_ksi[0] = pnts[i][j].q_hat[0] - 4.0*pnts[i-1][j].q_hat[0] + 6.0*pnts[i-2][j].q_hat[0] - 4.0*pnts[i-3][j].q_hat[0] + pnts[i-4][j].q_hat[0];
    pnts[i][j].diss_ksi[1] = pnts[i][j].q_hat[1] - 4.0*pnts[i-1][j].q_hat[1] + 6.0*pnts[i-2][j].q_hat[1] - 4.0*pnts[i-3][j].q_hat[1] + pnts[i-4][j].q_hat[1];
    pnts[i][j].diss_ksi[2] = pnts[i][j].q_hat[2] - 4.0*pnts[i-1][j].q_hat[2] + 6.0*pnts[i-2][j].q_hat[2] - 4.0*pnts[i-3][j].q_hat[2] + pnts[i-4][j].q_hat[2];
    pnts[i][j].diss_ksi[3] = pnts[i][j].q_hat[3] - 4.0*pnts[i-1][j].q_hat[3] + 6.0*pnts[i-2][j].q_hat[3] - 4.0*pnts[i-3][j].q_hat[3] + pnts[i-4][j].q_hat[3];

    pnts[i][j].diss_eta[0] = pnts[i][j+4].q_hat[0] - 4.0*pnts[i][j+3].q_hat[0] + 6.0*pnts[i][j+2].q_hat[0] - 4.0*pnts[i][j+1].q_hat[0] + pnts[i][j].q_hat[0];
    pnts[i][j].diss_eta[1] = pnts[i][j+4].q_hat[1] - 4.0*pnts[i][j+3].q_hat[1] + 6.0*pnts[i][j+2].q_hat[1] - 4.0*pnts[i][j+1].q_hat[1] + pnts[i][j].q_hat[1];
    pnts[i][j].diss_eta[2] = pnts[i][j+4].q_hat[2] - 4.0*pnts[i][j+3].q_hat[2] + 6.0*pnts[i][j+2].q_hat[2] - 4.0*pnts[i][j+1].q_hat[2] + pnts[i][j].q_hat[2];
    pnts[i][j].diss_eta[3] = pnts[i][j+4].q_hat[3] - 4.0*pnts[i][j+3].q_hat[3] + 6.0*pnts[i][j+2].q_hat[3] - 4.0*pnts[i][j+1].q_hat[3] + pnts[i][j].q_hat[3];

    /* Upper right corner. */

    i = imax - 2; j = jmax -2;

    pnts[i][j].diss_ksi[0] = pnts[i][j].q_hat[0] - 4.0*pnts[i-1][j].q_hat[0] + 6.0*pnts[i-2][j].q_hat[0] - 4.0*pnts[i-3][j].q_hat[0] + pnts[i-4][j].q_hat[0];
    pnts[i][j].diss_ksi[1] = pnts[i][j].q_hat[1] - 4.0*pnts[i-1][j].q_hat[1] + 6.0*pnts[i-2][j].q_hat[1] - 4.0*pnts[i-3][j].q_hat[1] + pnts[i-4][j].q_hat[1];
    pnts[i][j].diss_ksi[2] = pnts[i][j].q_hat[2] - 4.0*pnts[i-1][j].q_hat[2] + 6.0*pnts[i-2][j].q_hat[2] - 4.0*pnts[i-3][j].q_hat[2] + pnts[i-4][j].q_hat[2];
    pnts[i][j].diss_ksi[3] = pnts[i][j].q_hat[3] - 4.0*pnts[i-1][j].q_hat[3] + 6.0*pnts[i-2][j].q_hat[3] - 4.0*pnts[i-3][j].q_hat[3] + pnts[i-4][j].q_hat[3];

    pnts[i][j].diss_eta[0] = pnts[i][j].q_hat[0] - 4.0*pnts[i][j-1].q_hat[0] + 6.0*pnts[i][j-2].q_hat[0] - 4.0*pnts[i][j-3].q_hat[0] + pnts[i][j-4].q_hat[0];
    pnts[i][j].diss_eta[1] = pnts[i][j].q_hat[1] - 4.0*pnts[i][j-1].q_hat[1] + 6.0*pnts[i][j-2].q_hat[1] - 4.0*pnts[i][j-3].q_hat[1] + pnts[i][j-4].q_hat[1];
    pnts[i][j].diss_eta[2] = pnts[i][j].q_hat[2] - 4.0*pnts[i][j-1].q_hat[2] + 6.0*pnts[i][j-2].q_hat[2] - 4.0*pnts[i][j-3].q_hat[2] + pnts[i][j-4].q_hat[2];
    pnts[i][j].diss_eta[3] = pnts[i][j].q_hat[3] - 4.0*pnts[i][j-1].q_hat[3] + 6.0*pnts[i][j-2].q_hat[3] - 4.0*pnts[i][j-3].q_hat[3] + pnts[i][j-4].q_hat[3];

    /* Makes it coordinate transformed. */

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){

            pnts[i][j].q_hat[0] = pnts[i][j].J1 * pnts[i][j].q_hat[0];  
            pnts[i][j].q_hat[1] = pnts[i][j].J1 * pnts[i][j].q_hat[1];
            pnts[i][j].q_hat[2] = pnts[i][j].J1 * pnts[i][j].q_hat[2];
            pnts[i][j].q_hat[3] = pnts[i][j].J1 * pnts[i][j].q_hat[3];

        }
    }

    /* Now, add to the residue. */

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            pnts[i][j].RHS[0] = pnts[i][j].RHS[0] + pnts[i][j].J1 * p_setup.dissp2 * (pnts[i][j].diss_ksi[0] + pnts[i][j].diss_eta[0]);
            pnts[i][j].RHS[1] = pnts[i][j].RHS[1] + pnts[i][j].J1 * p_setup.dissp2 * (pnts[i][j].diss_ksi[1] + pnts[i][j].diss_eta[1]);
            pnts[i][j].RHS[2] = pnts[i][j].RHS[2] + pnts[i][j].J1 * p_setup.dissp2 * (pnts[i][j].diss_ksi[2] + pnts[i][j].diss_eta[2]);
            pnts[i][j].RHS[3] = pnts[i][j].RHS[3] + pnts[i][j].J1 * p_setup.dissp2 * (pnts[i][j].diss_ksi[3] + pnts[i][j].diss_eta[3]);

        }
    }
}
