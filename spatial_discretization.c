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

            double rho, u, v, e, p;

            /* Separate the properties we need. */

            rho = pnts[i][j].q[0];
            u   = pnts[i][j].q[1] / pnts[i][j].q[0];
            v   = pnts[i][j].q[2] / pnts[i][j].q[0];
            e   = pnts[i][j].q[3];

            /* Compute pressure. */

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

void compute_rhs(t_define p_setup, t_points ** pnts){

    /* Separate bounds of the field. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Initialize the residue. */

    max_rhs_rho  = -1.0;
    max_rhs_rhou = -1.0;
    max_rhs_rhov = -1.0;
    max_rhs_e    = -1.0;

    /* Compute the RHS for the internal points. */

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

            if ( fabs( pnts[i][j].RHS[0]) > max_rhs_rho )  max_rhs_rho  = fabs(pnts[i][j].RHS[0]);
            if ( fabs( pnts[i][j].RHS[1]) > max_rhs_rhou)  max_rhs_rhou = fabs(pnts[i][j].RHS[1]);
            if ( fabs( pnts[i][j].RHS[2]) > max_rhs_rhov)  max_rhs_rhov = fabs(pnts[i][j].RHS[2]);
            if ( fabs( pnts[i][j].RHS[3]) > max_rhs_e   )  max_rhs_e    = fabs(pnts[i][j].RHS[3]);

        }
    }
}


/* This function computes the artificial dissipation. */

void jst_art_dissip(t_define p_setup, t_points ** pnts){

    /* Separate bounds and constants.*/

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    double k2 = 1.0/4.0; double k4 = 1.0/256.0;

    /* Allocate some aux vectors. */

    double ** p =  alloc_double_matrix(imax, jmax);

    /* First, compute the pressure in all points. */

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){

            double rho = pnts[i][j].q[0];
            double u   = pnts[i][j].q[1]/pnts[i][j].q[0];
            double v   = pnts[i][j].q[2]/pnts[i][j].q[0];
            double e   = pnts[i][j].q[3];

            p[i][j] = (p_setup.gamma - 1.0) * (e - 0.5*rho*(pow(u,2.0) + pow(v,2.0)));

        }
    }

    /* Compute the v_ij in all internal points. */

    double ** v_ij = alloc_double_matrix(imax, jmax);

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            v_ij[i][j] = abs(p[i+1][j] - 2.0*p[i][j] + p[i-1][j])/
                         abs(p[i+1][j]) + 2.0*abs(p[i][j]) + abs(p[i-1][j]);
        }
    }

    /* Compute the eps2 in ksi direction. */

    double ** eps2_ksi   = alloc_double_matrix(imax, jmax);
    double ** aux_eps = alloc_double_matrix(imax, jmax);

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            aux_eps[i][j] = k2*max(v_ij[i+1][j],v_ij[i][j]);
            
        }
    }

    /* Now store in the half points. */

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            eps2_ksi[i][j] = (aux_eps[i+1][j] - aux_eps[i-1][j])/2.0;
        }
    }

    /* Now compute eps 2 in eta direction. */

    double ** eps2_eta = alloc_double_matrix(imax, jmax);

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            aux_eps[i][j] = k2*max(v_ij[i][j+1],v_ij[i][j]);

        }
    }

    /* Now store in the half points. */

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            eps2_eta[i][j] = (aux_eps[i][j+1] - aux_eps[i][j-1])/2.0;

        }
    }

    /* Now, compute the eps4 in ksi direction. */

    double ** eps4_ksi = alloc_double_matrix(imax, jmax);

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            aux_eps[i][j] = max(0.0, k4 - eps2_ksi[i][j]);

        }
    }

    /* Now store in the half points. */

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            eps4_ksi[i][j] = (aux_eps[i+1][j] - aux_eps[i+1][j])/2.0;

        }
    }

    /* Now, compute the eps4 in eta direction. */

    double ** eps4_eta = alloc_double_matrix(imax, jmax);

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            aux_eps[i][j] = max(0.0, k4 - eps2_ksi[i][j]);

        }
    }

    /* Now store in the half points. */

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            eps4_eta[i][j] = (aux_eps[i][j+1] - aux_eps[i][j+1])/2.0;

        }
    }

    double ** d_iphj = alloc_double_matrix(imax, jmax);
    double ** d_ijph = alloc_double_matrix(imax, jmax);

    /* Needs to understand the d_ definition in eta. */
    /* Remember to ask how to deal with boundaries. */

    for (int i = 2; i<imax-2; i++){
        for (int j = 2; j<jmax-2; j++){

            double dt = pnts[i][j].dt;

            d_iphj[i][j] = (1.0/dt)*( (eps2_ksi[i][j]*pnts[i][j].q[0]) - 
                    eps4_ksi[i][j]*(pnts[i+2][j].q[0] - 3.0*pnts[i+1][j].q[0] + 3.0*pnts[i][j].q[0] - 3.0*pnts[i-1][j].q[0]) );

            d_ijph[i][j] = (1.0/dt)*( (eps2_eta[i][j]*pnts[i][j].q[0]) - 
                    eps4_eta[i][j]*(pnts[i][j+2].q[0] - 3.0*pnts[i][j+1].q[0] + 3.0*pnts[i][j].q[0] - 3.0*pnts[i][j-1].q[0]) );

        }
    }

    /* Apply the artificial dissipation to the residue. */

    for (int i = 2; i<imax-2; i++){
        for (int j = 2; j<jmax-2; j++){

            double d_ksi = d_iphj[i][j] + d_iphj[i-1][j];
            double d_eta = d_ijph[i][j] - d_ijph[i][j-1];

            double d[4];

            /* Apply store the dissipation term. */

            d[0] = pnts[i][j].q_hat[0]*d_ksi + pnts[i][j].q_hat[0]*d_eta;
            d[1] = pnts[i][j].q_hat[1]*d_ksi + pnts[i][j].q_hat[1]*d_eta;
            d[2] = pnts[i][j].q_hat[2]*d_ksi + pnts[i][j].q_hat[2]*d_eta;
            d[3] = pnts[i][j].q_hat[3]*d_ksi + pnts[i][j].q_hat[3]*d_eta;

            /* Apply to residue. */

            pnts[i][j].RHS[0] = pnts[i][j].RHS[0] - d[0]; 
            pnts[i][j].RHS[1] = pnts[i][j].RHS[1] - d[1]; 
            pnts[i][j].RHS[2] = pnts[i][j].RHS[2] - d[2]; 
            pnts[i][j].RHS[3] = pnts[i][j].RHS[3] - d[3]; 


        }
    }

    /* Free everyone. */

    free_double_matrix(p, imax);
    free_double_matrix(v_ij, imax);
    free_double_matrix(eps2_ksi, imax);
    free_double_matrix(aux_eps, imax);
    free_double_matrix(eps2_eta, imax);
    free_double_matrix(eps4_ksi, imax);
    free_double_matrix(eps4_eta, imax);
    free_double_matrix(d_iphj, imax);
    free_double_matrix(d_ijph, imax);

}

