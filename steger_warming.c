#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "structs.h"
#include "externs.h"
#include "prototypes.h"

/* Here all the routines needed in order to compute the flow solution 
 * using the steger-warming flux vector splitting. */

void compute_sw_fluxes(t_define p_setup, t_points ** pnts){

    /* Separate bounds of the field. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Allocate the fluxes. */

    double eig[4], eig_p[4], eig_m[4], aux[4];

    /* Loop through internal points and compute the needed stuff.*/

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){

            /* Separate the properties we need. */

            double rho = pnts[i][j].J * pnts[i][j].q_hat[0];
            double u   = pnts[i][j].q_hat[1] / pnts[i][j].q_hat[0];
            double v   = pnts[i][j].q_hat[2] / pnts[i][j].q_hat[0];

            /* Separate useful variables. */

            double a   = pnts[i][j].a;

            double k1 = pnts[i][j].ksi_x;
            double k2 = pnts[i][j].ksi_y;

            double kt1 = k1/(pow(k1*k1 + k2*k2,0.5));
            double kt2 = k2/(pow(k1*k1 + k2*k2,0.5));

            /* Compute the auxiliar and the needed variables to compute the fluxes. */

            /* I tried to use Eq. B5 of SW original paper to compute the
             * eigenvalues using the Kt variables. However, in general
             * curvilinear coordinates they are not working properly while
             * building the eigenvalues. */

            eig[0] = kt1*pnts[i][j].cov_u;
            eig[1] = kt1*pnts[i][j].cov_u;
            eig[2] = kt1*pnts[i][j].cov_u + a*pow(k1*k1 + k2*k2,0.5);
            eig[3] = kt1*pnts[i][j].cov_u - a*pow(k1*k1 + k2*k2,0.5);

            /* Eq. 4.4 of SW original paper. */

            eig_p[0] = (eig[0] + fabs(eig[0])) / 2.0;
            eig_p[1] = (eig[1] + fabs(eig[1])) / 2.0;
            eig_p[2] = (eig[2] + fabs(eig[2])) / 2.0;
            eig_p[3] = (eig[3] + fabs(eig[3])) / 2.0;

            eig_m[0] = (eig[0] - fabs(eig[0])) / 2.0;
            eig_m[1] = (eig[1] - fabs(eig[1])) / 2.0;
            eig_m[2] = (eig[2] - fabs(eig[2])) / 2.0;
            eig_m[3] = (eig[3] - fabs(eig[3])) / 2.0;

            /* Get the positive F fluxes. */

            double w2 = ((3.0 - p_setup.gamma) * (eig_p[2] + eig_p[3])*pow(a,2.0))/(2.0*(p_setup.gamma - 1.0));

            aux[0] = 2.0*(p_setup.gamma - 1.0)*eig_p[0]   + eig_p[2]             + eig_p[3];
            aux[1] = 2.0*(p_setup.gamma - 1.0)*eig_p[0]*u + eig_p[2]*(u + a*kt1) + eig_p[3]*(u - a*kt1);
            aux[2] = 2.0*(p_setup.gamma - 1.0)*eig_p[0]*v + eig_p[2]*(v + a*kt2) + eig_p[3]*(v - a*kt2);
            aux[3] = (p_setup.gamma - 1.0)*eig_p[0]*(u*u + v*v)                  + 
                        (eig_p[2]/2.0)*(pow(u + a*kt1,2.0) + pow(v + a*kt2,2.0)) + 
                        (eig_p[3]/2.0)*(pow(u - a*kt1,2.0) + pow(v - a*kt2,2.0)) + w2;

            /* Do the proper scaling. */

            pnts[i][j].f_plus[0] = pnts[i][j].J1 * (rho/(2.0*p_setup.gamma))*aux[0];
            pnts[i][j].f_plus[1] = pnts[i][j].J1 * (rho/(2.0*p_setup.gamma))*aux[1];
            pnts[i][j].f_plus[2] = pnts[i][j].J1 * (rho/(2.0*p_setup.gamma))*aux[2];
            pnts[i][j].f_plus[3] = pnts[i][j].J1 * (rho/(2.0*p_setup.gamma))*aux[3];

            /* Get the negative F fluxes. */

            w2 = ((3.0 - p_setup.gamma) * (eig_m[2] + eig_m[3])*pow(a,2.0))/(2.0*(p_setup.gamma - 1.0));

            aux[0] = 2.0*(p_setup.gamma - 1.0)*eig_m[0]   + eig_m[2]             + eig_m[3];
            aux[1] = 2.0*(p_setup.gamma - 1.0)*eig_m[0]*u + eig_m[2]*(u + a*kt1) + eig_m[3]*(u - a*kt1);
            aux[2] = 2.0*(p_setup.gamma - 1.0)*eig_m[0]*v + eig_m[2]*(v + a*kt2) + eig_m[3]*(v - a*kt2);
            aux[3] = (p_setup.gamma - 1.0)*eig_m[0]*(u*u + v*v)                  + 
                        (eig_m[2]/2.0)*(pow(u + a*kt1,2.0) + pow(v + a*kt2,2.0)) + 
                        (eig_m[3]/2.0)*(pow(u - a*kt1,2.0) + pow(v - a*kt2,2.0)) + w2;

            /* Do the proper scaling. */

            pnts[i][j].f_minu[0] = pnts[i][j].J1 * (rho/(2.0*p_setup.gamma))*aux[0];
            pnts[i][j].f_minu[1] = pnts[i][j].J1 * (rho/(2.0*p_setup.gamma))*aux[1];
            pnts[i][j].f_minu[2] = pnts[i][j].J1 * (rho/(2.0*p_setup.gamma))*aux[2];
            pnts[i][j].f_minu[3] = pnts[i][j].J1 * (rho/(2.0*p_setup.gamma))*aux[3];

            /* Now compute the G positive and negative fluxes. */

            k1 = pnts[i][j].eta_x;
            k2 = pnts[i][j].eta_y;

            kt1 = k1/(pow(k1*k1 + k2*k2,0.5));
            kt2 = k2/(pow(k1*k1 + k2*k2,0.5));

            /* Compute the auxiliar and the needed variables to compute the fluxes. */

            eig[0] = kt2*pnts[i][j].cov_v;
            eig[1] = kt2*pnts[i][j].cov_v;
            eig[2] = kt2*pnts[i][j].cov_v + a*pow(k1*k1 + k2*k2,0.5);
            eig[3] = kt2*pnts[i][j].cov_v - a*pow(k1*k1 + k2*k2,0.5);

            eig_p[0] = (eig[0] + fabs(eig[0])) / 2.0;
            eig_p[1] = (eig[1] + fabs(eig[1])) / 2.0;
            eig_p[2] = (eig[2] + fabs(eig[2])) / 2.0;
            eig_p[3] = (eig[3] + fabs(eig[3])) / 2.0;

            eig_m[0] = (eig[0] - fabs(eig[0])) / 2.0;
            eig_m[1] = (eig[1] - fabs(eig[1])) / 2.0;
            eig_m[2] = (eig[2] - fabs(eig[2])) / 2.0;
            eig_m[3] = (eig[3] - fabs(eig[3])) / 2.0;

            /* Get the positive F fluxes. */

            w2 = ((3.0 - p_setup.gamma) * (eig_p[2] + eig_p[3])*pow(a,2.0))/(2.0*(p_setup.gamma - 1.0));

            aux[0] = 2.0*(p_setup.gamma - 1.0)*eig_p[0]   + eig_p[2]             + eig_p[3];
            aux[1] = 2.0*(p_setup.gamma - 1.0)*eig_p[0]*u + eig_p[2]*(u + a*kt1) + eig_p[3]*(u - a*kt1);
            aux[2] = 2.0*(p_setup.gamma - 1.0)*eig_p[0]*v + eig_p[2]*(v + a*kt2) + eig_p[3]*(v - a*kt2);
            aux[3] = (p_setup.gamma - 1.0)*eig_p[0]*(u*u + v*v)                  + 
                        (eig_p[2]/2.0)*(pow(u + a*kt1,2.0) + pow(v + a*kt2,2.0)) + 
                        (eig_p[3]/2.0)*(pow(u - a*kt1,2.0) + pow(v - a*kt2,2.0)) + w2;

            /* Do the proper scaling. */

            pnts[i][j].g_plus[0] = pnts[i][j].J1 * (rho/(2.0*p_setup.gamma))*aux[0];
            pnts[i][j].g_plus[1] = pnts[i][j].J1 * (rho/(2.0*p_setup.gamma))*aux[1];
            pnts[i][j].g_plus[2] = pnts[i][j].J1 * (rho/(2.0*p_setup.gamma))*aux[2];
            pnts[i][j].g_plus[3] = pnts[i][j].J1 * (rho/(2.0*p_setup.gamma))*aux[3];

            /* Get the negative F fluxes. */

            w2 = ((3.0 - p_setup.gamma) * (eig_m[2] + eig_m[3])*pow(a,2.0))/(2.0*(p_setup.gamma - 1.0));

            aux[0] = 2.0*(p_setup.gamma - 1.0)*eig_m[0]   + eig_m[2]             + eig_m[3];
            aux[1] = 2.0*(p_setup.gamma - 1.0)*eig_m[0]*u + eig_m[2]*(u + a*kt1) + eig_m[3]*(u - a*kt1);
            aux[2] = 2.0*(p_setup.gamma - 1.0)*eig_m[0]*v + eig_m[2]*(v + a*kt2) + eig_m[3]*(v - a*kt2);
            aux[3] = (p_setup.gamma - 1.0)*eig_m[0]*(u*u + v*v)                  + 
                        (eig_m[2]/2.0)*(pow(u + a*kt1,2.0) + pow(v + a*kt2,2.0)) + 
                        (eig_m[3]/2.0)*(pow(u - a*kt1,2.0) + pow(v - a*kt2,2.0)) + w2;

            /* Do the proper scaling. */

            pnts[i][j].g_minu[0] = pnts[i][j].J1 * (rho/(2.0*p_setup.gamma))*aux[0];
            pnts[i][j].g_minu[1] = pnts[i][j].J1 * (rho/(2.0*p_setup.gamma))*aux[1];
            pnts[i][j].g_minu[2] = pnts[i][j].J1 * (rho/(2.0*p_setup.gamma))*aux[2];
            pnts[i][j].g_minu[3] = pnts[i][j].J1 * (rho/(2.0*p_setup.gamma))*aux[3];

        }
    }
}

/* Computes the steger warming residues. */

void compute_sw_residue_1sto(t_define p_setup, t_points ** pnts){

    /* Separate bounds of the field. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Initialize the residue. */

    max_rhs_rho  = -99999999.0;
    max_rhs_rhou = -99999999.0;
    max_rhs_rhov = -99999999.0;
    max_rhs_e    = -99999999.0;

    /* Now, with the splited fluxes computed, lets join everything in our residue. For internal points */

    double dbFp_ksi[4];
    double dfFm_ksi[4];

    double dbGp_eta[4];
    double dfGm_eta[4];

    /* Initialize the residues to avoid problems. */

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){

            pnts[i][j].RHS[0] = -10.0;
            pnts[i][j].RHS[1] = -10.0;
            pnts[i][j].RHS[2] = -10.0;
            pnts[i][j].RHS[3] = -10.0;
        }
    }

    /* Build the fluxes for real. */

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            /* F flux in ksi direction */

            dbFp_ksi[0] = pnts[i][j].f_plus[0] - pnts[i-1][j].f_plus[0];
            dbFp_ksi[1] = pnts[i][j].f_plus[1] - pnts[i-1][j].f_plus[1];
            dbFp_ksi[2] = pnts[i][j].f_plus[2] - pnts[i-1][j].f_plus[2];
            dbFp_ksi[3] = pnts[i][j].f_plus[3] - pnts[i-1][j].f_plus[3];

            dfFm_ksi[0] = pnts[i+1][j].f_minu[0] - pnts[i][j].f_minu[0];
            dfFm_ksi[1] = pnts[i+1][j].f_minu[1] - pnts[i][j].f_minu[1];
            dfFm_ksi[2] = pnts[i+1][j].f_minu[2] - pnts[i][j].f_minu[2];
            dfFm_ksi[3] = pnts[i+1][j].f_minu[3] - pnts[i][j].f_minu[3];

            /* G flux in eta direction */

            dbGp_eta[0] = pnts[i][j].g_plus[0] - pnts[i][j-1].g_plus[0];
            dbGp_eta[1] = pnts[i][j].g_plus[1] - pnts[i][j-1].g_plus[1];
            dbGp_eta[2] = pnts[i][j].g_plus[2] - pnts[i][j-1].g_plus[2];
            dbGp_eta[3] = pnts[i][j].g_plus[3] - pnts[i][j-1].g_plus[3];

            dfGm_eta[0] = pnts[i][j+1].g_minu[0] - pnts[i][j].g_minu[0];
            dfGm_eta[1] = pnts[i][j+1].g_minu[1] - pnts[i][j].g_minu[1];
            dfGm_eta[2] = pnts[i][j+1].g_minu[2] - pnts[i][j].g_minu[2];
            dfGm_eta[3] = pnts[i][j+1].g_minu[3] - pnts[i][j].g_minu[3];

            /* Build the residues. */

            pnts[i][j].RHS[0] = dbFp_ksi[0] + dfFm_ksi[0] + dbGp_eta[0] + dfGm_eta[0];
            pnts[i][j].RHS[1] = dbFp_ksi[1] + dfFm_ksi[1] + dbGp_eta[1] + dfGm_eta[1];
            pnts[i][j].RHS[2] = dbFp_ksi[2] + dfFm_ksi[2] + dbGp_eta[2] + dfGm_eta[2];
            pnts[i][j].RHS[3] = dbFp_ksi[3] + dfFm_ksi[3] + dbGp_eta[3] + dfGm_eta[3];

            
        }
    }

    /* Collect the maximun residue. */

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            /* Store the max residue. */

            if (fabs( pnts[i][j].RHS[0]) > max_rhs_rho ) max_rhs_rho  = fabs(pnts[i][j].RHS[0]+DBL_EPSILON);
            if (fabs( pnts[i][j].RHS[1]) > max_rhs_rhou) max_rhs_rhou = fabs(pnts[i][j].RHS[1]+DBL_EPSILON);
            if (fabs( pnts[i][j].RHS[2]) > max_rhs_rhov) max_rhs_rhov = fabs(pnts[i][j].RHS[2]+DBL_EPSILON);
            if (fabs( pnts[i][j].RHS[3]) > max_rhs_e   ) max_rhs_e    = fabs(pnts[i][j].RHS[3]+DBL_EPSILON);
        }
    }
}

/* Computes the steger warming residues. */

void compute_sw_residue_2ndo(t_define p_setup, t_points ** pnts){

    /* Separate bounds of the field. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Initialize the residue. */

    max_rhs_rho  = -99999999.0;
    max_rhs_rhou = -99999999.0;
    max_rhs_rhov = -99999999.0;
    max_rhs_e    = -99999999.0;

    /* Now, with the splited fluxes computed, lets join everything in our residue. For internal points */

    double dbFp_ksi[4];
    double dfFm_ksi[4];

    double dbGp_eta[4];
    double dfGm_eta[4];

    /* Initialize the residues to avoid problems. */

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){

            pnts[i][j].RHS[0] = -10.0;
            pnts[i][j].RHS[1] = -10.0;
            pnts[i][j].RHS[2] = -10.0;
            pnts[i][j].RHS[3] = -10.0;
        }
    }

    /* As we will use first order operators over the borders, initialize the
     * field with then. */

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            /* F flux in ksi direction */

            dbFp_ksi[0] = pnts[i][j].f_plus[0] - pnts[i-1][j].f_plus[0];
            dbFp_ksi[1] = pnts[i][j].f_plus[1] - pnts[i-1][j].f_plus[1];
            dbFp_ksi[2] = pnts[i][j].f_plus[2] - pnts[i-1][j].f_plus[2];
            dbFp_ksi[3] = pnts[i][j].f_plus[3] - pnts[i-1][j].f_plus[3];

            dfFm_ksi[0] = pnts[i+1][j].f_minu[0] - pnts[i][j].f_minu[0];
            dfFm_ksi[1] = pnts[i+1][j].f_minu[1] - pnts[i][j].f_minu[1];
            dfFm_ksi[2] = pnts[i+1][j].f_minu[2] - pnts[i][j].f_minu[2];
            dfFm_ksi[3] = pnts[i+1][j].f_minu[3] - pnts[i][j].f_minu[3];

            /* G flux in eta direction */

            dbGp_eta[0] = pnts[i][j].g_plus[0] - pnts[i][j-1].g_plus[0];
            dbGp_eta[1] = pnts[i][j].g_plus[1] - pnts[i][j-1].g_plus[1];
            dbGp_eta[2] = pnts[i][j].g_plus[2] - pnts[i][j-1].g_plus[2];
            dbGp_eta[3] = pnts[i][j].g_plus[3] - pnts[i][j-1].g_plus[3];

            dfGm_eta[0] = pnts[i][j+1].g_minu[0] - pnts[i][j].g_minu[0];
            dfGm_eta[1] = pnts[i][j+1].g_minu[1] - pnts[i][j].g_minu[1];
            dfGm_eta[2] = pnts[i][j+1].g_minu[2] - pnts[i][j].g_minu[2];
            dfGm_eta[3] = pnts[i][j+1].g_minu[3] - pnts[i][j].g_minu[3];

            /* Build the residues. */

            pnts[i][j].RHS[0] = dbFp_ksi[0] + dfFm_ksi[0] + dbGp_eta[0] + dfGm_eta[0];
            pnts[i][j].RHS[1] = dbFp_ksi[1] + dfFm_ksi[1] + dbGp_eta[1] + dfGm_eta[1];
            pnts[i][j].RHS[2] = dbFp_ksi[2] + dfFm_ksi[2] + dbGp_eta[2] + dfGm_eta[2];
            pnts[i][j].RHS[3] = dbFp_ksi[3] + dfFm_ksi[3] + dbGp_eta[3] + dfGm_eta[3];
            
        }
    }

    /* Build the fluxes using second order operators. */

    for (int i = 2; i<imax-2; i++){
        for (int j = 2; j<jmax-2; j++){

            /* F flux in ksi direction */

            dbFp_ksi[0] = (3.0*pnts[i][j].f_plus[0] - 4.0*pnts[i-1][j].f_plus[0] + pnts[i-2][j].f_plus[0])/2.0;
            dbFp_ksi[1] = (3.0*pnts[i][j].f_plus[1] - 4.0*pnts[i-1][j].f_plus[1] + pnts[i-2][j].f_plus[1])/2.0;
            dbFp_ksi[2] = (3.0*pnts[i][j].f_plus[2] - 4.0*pnts[i-1][j].f_plus[2] + pnts[i-2][j].f_plus[2])/2.0;
            dbFp_ksi[3] = (3.0*pnts[i][j].f_plus[3] - 4.0*pnts[i-1][j].f_plus[3] + pnts[i-2][j].f_plus[3])/2.0;

            dfFm_ksi[0] = (-3.0*pnts[i][j].f_minu[0] + 4.0*pnts[i+1][j].f_minu[0] - pnts[i+2][j].f_minu[0])/2.0; 
            dfFm_ksi[1] = (-3.0*pnts[i][j].f_minu[1] + 4.0*pnts[i+1][j].f_minu[1] - pnts[i+2][j].f_minu[1])/2.0; 
            dfFm_ksi[2] = (-3.0*pnts[i][j].f_minu[2] + 4.0*pnts[i+1][j].f_minu[2] - pnts[i+2][j].f_minu[2])/2.0; 
            dfFm_ksi[3] = (-3.0*pnts[i][j].f_minu[3] + 4.0*pnts[i+1][j].f_minu[3] - pnts[i+2][j].f_minu[3])/2.0; 

            /* G flux in eta direction */

            dbGp_eta[0] = (3.0*pnts[i][j].g_plus[0] - 4.0*pnts[i][j-1].g_plus[0] + pnts[i][j-2].g_plus[0])/2.0;
            dbGp_eta[1] = (3.0*pnts[i][j].g_plus[1] - 4.0*pnts[i][j-1].g_plus[1] + pnts[i][j-2].g_plus[1])/2.0;
            dbGp_eta[2] = (3.0*pnts[i][j].g_plus[2] - 4.0*pnts[i][j-1].g_plus[2] + pnts[i][j-2].g_plus[2])/2.0;
            dbGp_eta[3] = (3.0*pnts[i][j].g_plus[3] - 4.0*pnts[i][j-1].g_plus[3] + pnts[i][j-2].g_plus[3])/2.0;

            dfGm_eta[0] = (-3.0*pnts[i][j].g_minu[0] + 4.0*pnts[i][j+1].g_minu[0] - pnts[i][j+2].g_minu[0])/2.0; 
            dfGm_eta[1] = (-3.0*pnts[i][j].g_minu[1] + 4.0*pnts[i][j+1].g_minu[1] - pnts[i][j+2].g_minu[1])/2.0; 
            dfGm_eta[2] = (-3.0*pnts[i][j].g_minu[2] + 4.0*pnts[i][j+1].g_minu[2] - pnts[i][j+2].g_minu[2])/2.0; 
            dfGm_eta[3] = (-3.0*pnts[i][j].g_minu[3] + 4.0*pnts[i][j+1].g_minu[3] - pnts[i][j+2].g_minu[3])/2.0; 

            /* Build the residues. */

            pnts[i][j].RHS[0] = dbFp_ksi[0] + dfFm_ksi[0] + dbGp_eta[0] + dfGm_eta[0];
            pnts[i][j].RHS[1] = dbFp_ksi[1] + dfFm_ksi[1] + dbGp_eta[1] + dfGm_eta[1];
            pnts[i][j].RHS[2] = dbFp_ksi[2] + dfFm_ksi[2] + dbGp_eta[2] + dfGm_eta[2];
            pnts[i][j].RHS[3] = dbFp_ksi[3] + dfFm_ksi[3] + dbGp_eta[3] + dfGm_eta[3];

            
        }
    }

    /* Collect the maximun residue. */

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            /* Store the max residue. */

            if (fabs( pnts[i][j].RHS[0]) > max_rhs_rho ) max_rhs_rho  = fabs(pnts[i][j].RHS[0]+DBL_EPSILON);
            if (fabs( pnts[i][j].RHS[1]) > max_rhs_rhou) max_rhs_rhou = fabs(pnts[i][j].RHS[1]+DBL_EPSILON);
            if (fabs( pnts[i][j].RHS[2]) > max_rhs_rhov) max_rhs_rhov = fabs(pnts[i][j].RHS[2]+DBL_EPSILON);
            if (fabs( pnts[i][j].RHS[3]) > max_rhs_e   ) max_rhs_e    = fabs(pnts[i][j].RHS[3]+DBL_EPSILON);
        }
    }
}

/* Build the implicit operator of the SW scheme. */

void compute_sw_impicit_operator(t_define p_setup, t_points ** pnts){

    /* Select bounds. */
    
    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Define inner matrices size. */

    int nim = 4;

    /* Define the identity matrix. */

    double I[4][4];

    I[0][0] = 1.0; I[0][1] = 0.0; I[0][2] = 0.0; I[0][3] = 0.0; 
    I[1][0] = 0.0; I[1][1] = 1.0; I[1][2] = 0.0; I[1][3] = 0.0; 
    I[2][0] = 0.0; I[2][1] = 0.0; I[2][2] = 1.0; I[2][3] = 0.0; 
    I[3][0] = 0.0; I[3][1] = 0.0; I[3][2] = 0.0; I[3][3] = 1.0; 
    
    /* Allocate the du_s vector. It has to be the same size of the residue. */

    double *** DU_S  = alloc_dcube(imax-1, jmax-1, 4);

    /* Transport vector for the i and j passes. */

    double **  du_si = alloc_dmatrix(4,imax-1);
    double **  du_sj = alloc_dmatrix(4,jmax-1);

    /* Transport resudue vector for system solve. */

    double **  s_rhsi = alloc_dmatrix(4,imax-1);
    double **  s_rhsj = alloc_dmatrix(4,jmax-1);

    /* Allocate the three diagonals for the i pass */

    double *** upper_i = alloc_dcube(4,4,imax-1);
    double *** maind_i = alloc_dcube(4,4,imax-1);
    double *** lower_i = alloc_dcube(4,4,imax-1);

    /* Allocate the three diagonals for the j pass */

    double *** upper_j = alloc_dcube(4,4,jmax-1);
    double *** maind_j = alloc_dcube(4,4,jmax-1);
    double *** lower_j = alloc_dcube(4,4,jmax-1);

    /* Now loop through every line of the mesh. */

    for (int j = 1; j<jmax-1; j++){
        for (int i = 1; i<imax-1; i++){

            /* Estimate the constant. */

            double cte_i = pnts[i][j].dt;

            /* Build the upper diagonal. */

            for (int ii = 0; ii<nim; ii++)
            for (int jj = 0; jj<nim; jj++) upper_i[ii][jj][i] = cte_i * pnts[i+1][j].A_minu[ii][jj];

            /* Build the main diagonal. */

            for (int ii = 0; ii<nim; ii++)
            for (int jj = 0; jj<nim; jj++) maind_i[ii][jj][i] = I[ii][jj] + cte_i*pnts[i][j].A_plus[ii][jj] - cte_i*pnts[i][j].A_minu[ii][jj];

            /* Build the lower diagonal. */

            for (int ii = 0; ii<nim; ii++)
            for (int jj = 0; jj<nim; jj++) lower_i[ii][jj][i] = - cte_i * pnts[i-1][j].A_plus[ii][jj];

            /* Build the B vector using the RHS. */

            s_rhsi[0][i] = - pnts[i][j].dt * pnts[i][j].RHS[0];
            s_rhsi[1][i] = - pnts[i][j].dt * pnts[i][j].RHS[1];
            s_rhsi[2][i] = - pnts[i][j].dt * pnts[i][j].RHS[2];
            s_rhsi[3][i] = - pnts[i][j].dt * pnts[i][j].RHS[3];

        }

        /* Now, solve the system and store the du_s. */

        blk_tri(maind_i, lower_i, upper_i, nim, imax-1, s_rhsi, du_si);

        /* Store the final DU_S. */

        for (int i = 1; i<imax-1; i++){

            DU_S[i][j][0] = du_si[0][i];
            DU_S[i][j][1] = du_si[1][i];
            DU_S[i][j][2] = du_si[2][i];
            DU_S[i][j][3] = du_si[3][i];

        }
    }

    /* Now lets do the same in the other direction. */

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            /* Estimate the constant. */

            double cte_i = pnts[i][j].dt;

            /* Build the upper diagonal. */

            for (int ii = 0; ii<nim; ii++)
            for (int jj = 0; jj<nim; jj++) upper_j[ii][jj][j] = cte_i * pnts[i][j+1].B_minu[ii][jj];

            /* Build the main diagonal. */

            for (int ii = 0; ii<nim; ii++)
            for (int jj = 0; jj<nim; jj++) maind_j[ii][jj][j] = I[ii][jj] + cte_i*pnts[i][j].B_plus[ii][jj] - cte_i*pnts[i][j].B_minu[ii][jj];

            /* Build the lower diagonal. */

            for (int ii = 0; ii<nim; ii++)
            for (int jj = 0; jj<nim; jj++) lower_j[ii][jj][j] = - cte_i * pnts[i][j-1].B_plus[ii][jj];

            /* Build the B vector using the RHS. */

            s_rhsj[0][j] = pnts[i][j].dt * DU_S[i][j][0];
            s_rhsj[1][j] = pnts[i][j].dt * DU_S[i][j][1];
            s_rhsj[2][j] = pnts[i][j].dt * DU_S[i][j][2];
            s_rhsj[3][j] = pnts[i][j].dt * DU_S[i][j][3];

        }

        /* Now, solve the system and store the du_s. */

        blk_tri(maind_j, lower_j, upper_j, nim, jmax-1, s_rhsj, du_sj);

        /* Update the solution. */

        for (int j = 1; j<jmax-1; j++){

           pnts[i][j].q_hat[0] = pnts[i][j].q_hat[0] + du_sj[0][j];  
           pnts[i][j].q_hat[1] = pnts[i][j].q_hat[1] + du_sj[1][j];  
           pnts[i][j].q_hat[2] = pnts[i][j].q_hat[2] + du_sj[2][j];  
           pnts[i][j].q_hat[3] = pnts[i][j].q_hat[3] + du_sj[3][j];  

        }
    }

    /* Free everyone. */

    free_dcube(upper_i,4,4);
    free_dcube(maind_i,4,4);
    free_dcube(lower_i,4,4);
    free_dcube(upper_j,4,4);
    free_dcube(maind_j,4,4);
    free_dcube(lower_j,4,4);
    free_dmatrix(du_si,4);
    free_dmatrix(du_sj,4);
    free_dmatrix(s_rhsi,4);
    free_dmatrix(s_rhsj,4);
    free_dcube(DU_S,imax-1, jmax-1);

}
