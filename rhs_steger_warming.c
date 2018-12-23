#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"
#include "externs.h"
#include "prototypes.h"

/* Here all the routines needed in order to compute the flow solution 
 * using the steger-warming flux vector splitting. */

void sw_residue(t_define p_setup, t_points ** pnts){

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
            double e   = pnts[i][j].J * pnts[i][j].q_hat[3];

            /* Separate useful variables. */

            double a   = pnts[i][j].a;

            double k1 = pnts[i][j].ksi_x;
            double k2 = pnts[i][j].ksi_y;

            double kt1 = k1/(pow(k1*k1 + k2*k2,0.5));
            double kt2 = k2/(pow(k1*k1 + k2*k2,0.5));

            /* Compute the auxiliar and the needed variables to compute the fluxes. */

            eig[0] = kt1*u + kt2*v;
            eig[1] = kt1*u + kt2*v;
            eig[2] = kt1*u + kt2*v + a*pow(kt1*kt1 + kt2*kt2,0.5);
            eig[3] = kt1*u + kt2*v - a*pow(kt1*kt1 + kt2*kt2,0.5);

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

            pnts[i][j].f_plus[0] = (rho/(2.0*p_setup.gamma))*aux[0];
            pnts[i][j].f_plus[1] = (rho/(2.0*p_setup.gamma))*aux[1];
            pnts[i][j].f_plus[2] = (rho/(2.0*p_setup.gamma))*aux[2];
            pnts[i][j].f_plus[3] = (rho/(2.0*p_setup.gamma))*aux[3];

            /* Get the negative F fluxes. */

            w2 = ((3.0 - p_setup.gamma) * (eig_m[2] + eig_m[3])*pow(a,2.0))/(2.0*(p_setup.gamma - 1.0));

            aux[0] = 2.0*(p_setup.gamma - 1.0)*eig_m[0]   + eig_m[2]             + eig_m[3];
            aux[1] = 2.0*(p_setup.gamma - 1.0)*eig_m[0]*u + eig_m[2]*(u + a*kt1) + eig_m[3]*(u - a*kt1);
            aux[2] = 2.0*(p_setup.gamma - 1.0)*eig_m[0]*v + eig_m[2]*(v + a*kt2) + eig_m[3]*(v - a*kt2);
            aux[3] = (p_setup.gamma - 1.0)*eig_m[0]*(u*u + v*v)                  + 
                        (eig_m[2]/2.0)*(pow(u + a*kt1,2.0) + pow(v + a*kt2,2.0)) + 
                        (eig_m[3]/2.0)*(pow(u - a*kt1,2.0) + pow(v - a*kt2,2.0)) + w2;

            /* Do the proper scaling. */

            pnts[i][j].f_minu[0] = (rho/(2.0*p_setup.gamma))*aux[0];
            pnts[i][j].f_minu[1] = (rho/(2.0*p_setup.gamma))*aux[1];
            pnts[i][j].f_minu[2] = (rho/(2.0*p_setup.gamma))*aux[2];
            pnts[i][j].f_minu[3] = (rho/(2.0*p_setup.gamma))*aux[3];

            /* Now compute the G positive and negative fluxes. */

            k1 = pnts[i][j].eta_x;
            k2 = pnts[i][j].eta_y;

            a = pnts[i][j].a;

            kt1 = k1/(pow(k1*k1 + k2*k2,0.5));
            kt2 = k2/(pow(k1*k1 + k2*k2,0.5));

            /* Compute the auxiliar and the needed variables to compute the fluxes. */

            eig[0] = kt1*u + kt2*v;
            eig[1] = kt1*u + kt2*v;
            eig[2] = kt1*u + kt2*v + a*pow(kt1*kt1 + kt2*kt2,0.5);
            eig[3] = kt1*u + kt2*v - a*pow(kt1*kt1 + kt2*kt2,0.5);

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

            pnts[i][j].g_plus[0] = (rho/(2.0*p_setup.gamma))*aux[0];
            pnts[i][j].g_plus[1] = (rho/(2.0*p_setup.gamma))*aux[1];
            pnts[i][j].g_plus[2] = (rho/(2.0*p_setup.gamma))*aux[2];
            pnts[i][j].g_plus[3] = (rho/(2.0*p_setup.gamma))*aux[3];

            /* Get the negative F fluxes. */

            w2 = ((3.0 - p_setup.gamma) * (eig_m[2] + eig_m[3])*pow(a,2.0))/(2.0*(p_setup.gamma - 1.0));

            aux[0] = 2.0*(p_setup.gamma - 1.0)*eig_m[0]   + eig_m[2]             + eig_m[3];
            aux[1] = 2.0*(p_setup.gamma - 1.0)*eig_m[0]*u + eig_m[2]*(u + a*kt1) + eig_m[3]*(u - a*kt1);
            aux[2] = 2.0*(p_setup.gamma - 1.0)*eig_m[0]*v + eig_m[2]*(v + a*kt2) + eig_m[3]*(v - a*kt2);
            aux[3] = (p_setup.gamma - 1.0)*eig_m[0]*(u*u + v*v)                  + 
                        (eig_m[2]/2.0)*(pow(u + a*kt1,2.0) + pow(v + a*kt2,2.0)) + 
                        (eig_m[3]/2.0)*(pow(u - a*kt1,2.0) + pow(v - a*kt2,2.0)) + w2;

            /* Do the proper scaling. */

            pnts[i][j].g_minu[0] = (rho/(2.0*p_setup.gamma))*aux[0];
            pnts[i][j].g_minu[1] = (rho/(2.0*p_setup.gamma))*aux[1];
            pnts[i][j].g_minu[2] = (rho/(2.0*p_setup.gamma))*aux[2];
            pnts[i][j].g_minu[3] = (rho/(2.0*p_setup.gamma))*aux[3];

        }
    }

    /* Now, with the splited fluxes computed, lets join everything in our residue. For internal points */

    double dbFp_ksi[4];
    double dfFm_ksi[4];

    double dbGp_eta[4];
    double dfGm_eta[4];

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

    /* Do the same for the second inlet internal points. */

    for (int j = 1; j < jmax-1; j++){

        int i = 1;

        dbFp_ksi[0] = pnts[i][j].f_plus[0] - pnts[i-1][j].f_plus[0]; 
        dbFp_ksi[1] = pnts[i][j].f_plus[1] - pnts[i-1][j].f_plus[1]; 
        dbFp_ksi[2] = pnts[i][j].f_plus[2] - pnts[i-1][j].f_plus[2]; 
        dbFp_ksi[3] = pnts[i][j].f_plus[3] - pnts[i-1][j].f_plus[3]; 

        dfFm_ksi[0] = pnts[i+1][j].f_plus[0] - pnts[i][j].f_plus[0]; 
        dfFm_ksi[1] = pnts[i+1][j].f_plus[1] - pnts[i][j].f_plus[1]; 
        dfFm_ksi[2] = pnts[i+1][j].f_plus[2] - pnts[i][j].f_plus[2]; 
        dfFm_ksi[3] = pnts[i+1][j].f_plus[3] - pnts[i][j].f_plus[3]; 

        dbGp_eta[0] = pnts[i][j].g_plus[0] - pnts[i][j-1].g_plus[0];
        dbGp_eta[1] = pnts[i][j].g_plus[1] - pnts[i][j-1].g_plus[1];
        dbGp_eta[2] = pnts[i][j].g_plus[2] - pnts[i][j-1].g_plus[2];
        dbGp_eta[3] = pnts[i][j].g_plus[3] - pnts[i][j-1].g_plus[3];

        dfGm_eta[0] = pnts[i][j+1].g_plus[0] - pnts[i][j].g_plus[0];
        dfGm_eta[1] = pnts[i][j+1].g_plus[1] - pnts[i][j].g_plus[1];
        dfGm_eta[2] = pnts[i][j+1].g_plus[2] - pnts[i][j].g_plus[2];
        dfGm_eta[3] = pnts[i][j+1].g_plus[3] - pnts[i][j].g_plus[3];

        /* Build the residues. */

        pnts[i][j].RHS[0] = dbFp_ksi[0] + dfFm_ksi[0] + dbGp_eta[0] + dfGm_eta[0];
        pnts[i][j].RHS[1] = dbFp_ksi[1] + dfFm_ksi[1] + dbGp_eta[1] + dfGm_eta[1];
        pnts[i][j].RHS[2] = dbFp_ksi[2] + dfFm_ksi[2] + dbGp_eta[2] + dfGm_eta[2];
        pnts[i][j].RHS[3] = dbFp_ksi[3] + dfFm_ksi[3] + dbGp_eta[3] + dfGm_eta[3];

    }

    /* Symmetry internal points. */

    for (int i = 1; i < imax-1; i++){

        int j = jmax-2;

        dbFp_ksi[0] = pnts[i][j].f_plus[0] - pnts[i-1][j].f_plus[0]; 
        dbFp_ksi[1] = pnts[i][j].f_plus[1] - pnts[i-1][j].f_plus[1]; 
        dbFp_ksi[2] = pnts[i][j].f_plus[2] - pnts[i-1][j].f_plus[2]; 
        dbFp_ksi[3] = pnts[i][j].f_plus[3] - pnts[i-1][j].f_plus[3]; 

        dfFm_ksi[0] = pnts[i+1][j].f_plus[0] - pnts[i][j].f_plus[0]; 
        dfFm_ksi[1] = pnts[i+1][j].f_plus[1] - pnts[i][j].f_plus[1]; 
        dfFm_ksi[2] = pnts[i+1][j].f_plus[2] - pnts[i][j].f_plus[2]; 
        dfFm_ksi[3] = pnts[i+1][j].f_plus[3] - pnts[i][j].f_plus[3]; 

        dbGp_eta[0] = pnts[i][j].g_plus[0] - pnts[i][j-1].g_plus[0];
        dbGp_eta[1] = pnts[i][j].g_plus[1] - pnts[i][j-1].g_plus[1];
        dbGp_eta[2] = pnts[i][j].g_plus[2] - pnts[i][j-1].g_plus[2];
        dbGp_eta[3] = pnts[i][j].g_plus[3] - pnts[i][j-1].g_plus[3];

        dfGm_eta[0] = pnts[i][j+1].g_plus[0] - pnts[i][j].g_plus[0];
        dfGm_eta[1] = pnts[i][j+1].g_plus[1] - pnts[i][j].g_plus[1];
        dfGm_eta[2] = pnts[i][j+1].g_plus[2] - pnts[i][j].g_plus[2];
        dfGm_eta[3] = pnts[i][j+1].g_plus[3] - pnts[i][j].g_plus[3];

        /* Build the residues. */

        pnts[i][j].RHS[0] = dbFp_ksi[0] + dfFm_ksi[0] + dbGp_eta[0] + dfGm_eta[0];
        pnts[i][j].RHS[1] = dbFp_ksi[1] + dfFm_ksi[1] + dbGp_eta[1] + dfGm_eta[1];
        pnts[i][j].RHS[2] = dbFp_ksi[2] + dfFm_ksi[2] + dbGp_eta[2] + dfGm_eta[2];
        pnts[i][j].RHS[3] = dbFp_ksi[3] + dfFm_ksi[3] + dbGp_eta[3] + dfGm_eta[3];

    }

    /* Outlet internal points. */

    for (int j = 1; j < jmax-1; j++){

        int i = imax-2;

        dbFp_ksi[0] = pnts[i][j].f_plus[0] - pnts[i-1][j].f_plus[0]; 
        dbFp_ksi[1] = pnts[i][j].f_plus[1] - pnts[i-1][j].f_plus[1]; 
        dbFp_ksi[2] = pnts[i][j].f_plus[2] - pnts[i-1][j].f_plus[2]; 
        dbFp_ksi[3] = pnts[i][j].f_plus[3] - pnts[i-1][j].f_plus[3]; 

        dfFm_ksi[0] = pnts[i+1][j].f_plus[0] - pnts[i][j].f_plus[0]; 
        dfFm_ksi[1] = pnts[i+1][j].f_plus[1] - pnts[i][j].f_plus[1]; 
        dfFm_ksi[2] = pnts[i+1][j].f_plus[2] - pnts[i][j].f_plus[2]; 
        dfFm_ksi[3] = pnts[i+1][j].f_plus[3] - pnts[i][j].f_plus[3]; 

        dbGp_eta[0] = pnts[i][j].g_plus[0] - pnts[i][j-1].g_plus[0];
        dbGp_eta[1] = pnts[i][j].g_plus[1] - pnts[i][j-1].g_plus[1];
        dbGp_eta[2] = pnts[i][j].g_plus[2] - pnts[i][j-1].g_plus[2];
        dbGp_eta[3] = pnts[i][j].g_plus[3] - pnts[i][j-1].g_plus[3];

        dfGm_eta[0] = pnts[i][j+1].g_plus[0] - pnts[i][j].g_plus[0];
        dfGm_eta[1] = pnts[i][j+1].g_plus[1] - pnts[i][j].g_plus[1];
        dfGm_eta[2] = pnts[i][j+1].g_plus[2] - pnts[i][j].g_plus[2];
        dfGm_eta[3] = pnts[i][j+1].g_plus[3] - pnts[i][j].g_plus[3];

        /* Build the residues. */

        pnts[i][j].RHS[0] = dbFp_ksi[0] + dfFm_ksi[0] + dbGp_eta[0] + dfGm_eta[0];
        pnts[i][j].RHS[1] = dbFp_ksi[1] + dfFm_ksi[1] + dbGp_eta[1] + dfGm_eta[1];
        pnts[i][j].RHS[2] = dbFp_ksi[2] + dfFm_ksi[2] + dbGp_eta[2] + dfGm_eta[2];
        pnts[i][j].RHS[3] = dbFp_ksi[3] + dfFm_ksi[3] + dbGp_eta[3] + dfGm_eta[3];

    }

    /* Wall internal points. */

    for (int i = 1; i < imax-1; i++){

        int j = 1;

        dbFp_ksi[0] = pnts[i][j].f_plus[0] - pnts[i-1][j].f_plus[0]; 
        dbFp_ksi[1] = pnts[i][j].f_plus[1] - pnts[i-1][j].f_plus[1]; 
        dbFp_ksi[2] = pnts[i][j].f_plus[2] - pnts[i-1][j].f_plus[2]; 
        dbFp_ksi[3] = pnts[i][j].f_plus[3] - pnts[i-1][j].f_plus[3]; 

        dfFm_ksi[0] = pnts[i+1][j].f_plus[0] - pnts[i][j].f_plus[0]; 
        dfFm_ksi[1] = pnts[i+1][j].f_plus[1] - pnts[i][j].f_plus[1]; 
        dfFm_ksi[2] = pnts[i+1][j].f_plus[2] - pnts[i][j].f_plus[2]; 
        dfFm_ksi[3] = pnts[i+1][j].f_plus[3] - pnts[i][j].f_plus[3]; 

        dbGp_eta[0] = pnts[i][j].g_plus[0] - pnts[i][j-1].g_plus[0];
        dbGp_eta[1] = pnts[i][j].g_plus[1] - pnts[i][j-1].g_plus[1];
        dbGp_eta[2] = pnts[i][j].g_plus[2] - pnts[i][j-1].g_plus[2];
        dbGp_eta[3] = pnts[i][j].g_plus[3] - pnts[i][j-1].g_plus[3];

        dfGm_eta[0] = pnts[i][j+1].g_plus[0] - pnts[i][j].g_plus[0];
        dfGm_eta[1] = pnts[i][j+1].g_plus[1] - pnts[i][j].g_plus[1];
        dfGm_eta[2] = pnts[i][j+1].g_plus[2] - pnts[i][j].g_plus[2];
        dfGm_eta[3] = pnts[i][j+1].g_plus[3] - pnts[i][j].g_plus[3];

        /* Build the residues. */

        pnts[i][j].RHS[0] = dbFp_ksi[0] + dfFm_ksi[0] + dbGp_eta[0] + dfGm_eta[0];
        pnts[i][j].RHS[1] = dbFp_ksi[1] + dfFm_ksi[1] + dbGp_eta[1] + dfGm_eta[1];
        pnts[i][j].RHS[2] = dbFp_ksi[2] + dfFm_ksi[2] + dbGp_eta[2] + dfGm_eta[2];
        pnts[i][j].RHS[3] = dbFp_ksi[3] + dfFm_ksi[3] + dbGp_eta[3] + dfGm_eta[3];

    }

    /* Lower left corner point. */

    int i = 1; int j = 1;

    dbFp_ksi[0] = pnts[i][j].f_plus[0] - pnts[i-1][j].f_plus[0]; 
    dbFp_ksi[1] = pnts[i][j].f_plus[1] - pnts[i-1][j].f_plus[1]; 
    dbFp_ksi[2] = pnts[i][j].f_plus[2] - pnts[i-1][j].f_plus[2]; 
    dbFp_ksi[3] = pnts[i][j].f_plus[3] - pnts[i-1][j].f_plus[3]; 

    dfFm_ksi[0] = pnts[i+1][j].f_plus[0] - pnts[i][j].f_plus[0]; 
    dfFm_ksi[1] = pnts[i+1][j].f_plus[1] - pnts[i][j].f_plus[1]; 
    dfFm_ksi[2] = pnts[i+1][j].f_plus[2] - pnts[i][j].f_plus[2]; 
    dfFm_ksi[3] = pnts[i+1][j].f_plus[3] - pnts[i][j].f_plus[3]; 

    dbGp_eta[0] = pnts[i][j].g_plus[0] - pnts[i][j-1].g_plus[0];
    dbGp_eta[1] = pnts[i][j].g_plus[1] - pnts[i][j-1].g_plus[1];
    dbGp_eta[2] = pnts[i][j].g_plus[2] - pnts[i][j-1].g_plus[2];
    dbGp_eta[3] = pnts[i][j].g_plus[3] - pnts[i][j-1].g_plus[3];

    dfGm_eta[0] = pnts[i][j+1].g_plus[0] - pnts[i][j].g_plus[0];
    dfGm_eta[1] = pnts[i][j+1].g_plus[1] - pnts[i][j].g_plus[1];
    dfGm_eta[2] = pnts[i][j+1].g_plus[2] - pnts[i][j].g_plus[2];
    dfGm_eta[3] = pnts[i][j+1].g_plus[3] - pnts[i][j].g_plus[3];

    /* Build the residues. */

    pnts[i][j].RHS[0] = dbFp_ksi[0] + dfFm_ksi[0] + dbGp_eta[0] + dfGm_eta[0];
    pnts[i][j].RHS[1] = dbFp_ksi[1] + dfFm_ksi[1] + dbGp_eta[1] + dfGm_eta[1];
    pnts[i][j].RHS[2] = dbFp_ksi[2] + dfFm_ksi[2] + dbGp_eta[2] + dfGm_eta[2];
    pnts[i][j].RHS[3] = dbFp_ksi[3] + dfFm_ksi[3] + dbGp_eta[3] + dfGm_eta[3];

    /* Lower right point. */

    i = imax-2; j = 1;

    dbFp_ksi[0] = pnts[i][j].f_plus[0] - pnts[i-1][j].f_plus[0]; 
    dbFp_ksi[1] = pnts[i][j].f_plus[1] - pnts[i-1][j].f_plus[1]; 
    dbFp_ksi[2] = pnts[i][j].f_plus[2] - pnts[i-1][j].f_plus[2]; 
    dbFp_ksi[3] = pnts[i][j].f_plus[3] - pnts[i-1][j].f_plus[3]; 

    dfFm_ksi[0] = pnts[i+1][j].f_plus[0] - pnts[i][j].f_plus[0]; 
    dfFm_ksi[1] = pnts[i+1][j].f_plus[1] - pnts[i][j].f_plus[1]; 
    dfFm_ksi[2] = pnts[i+1][j].f_plus[2] - pnts[i][j].f_plus[2]; 
    dfFm_ksi[3] = pnts[i+1][j].f_plus[3] - pnts[i][j].f_plus[3]; 

    dbGp_eta[0] = pnts[i][j].g_plus[0] - pnts[i][j-1].g_plus[0];
    dbGp_eta[1] = pnts[i][j].g_plus[1] - pnts[i][j-1].g_plus[1];
    dbGp_eta[2] = pnts[i][j].g_plus[2] - pnts[i][j-1].g_plus[2];
    dbGp_eta[3] = pnts[i][j].g_plus[3] - pnts[i][j-1].g_plus[3];

    dfGm_eta[0] = pnts[i][j+1].g_plus[0] - pnts[i][j].g_plus[0];
    dfGm_eta[1] = pnts[i][j+1].g_plus[1] - pnts[i][j].g_plus[1];
    dfGm_eta[2] = pnts[i][j+1].g_plus[2] - pnts[i][j].g_plus[2];
    dfGm_eta[3] = pnts[i][j+1].g_plus[3] - pnts[i][j].g_plus[3];

    /* Build the residues. */

    pnts[i][j].RHS[0] = dbFp_ksi[0] + dfFm_ksi[0] + dbGp_eta[0] + dfGm_eta[0];
    pnts[i][j].RHS[1] = dbFp_ksi[1] + dfFm_ksi[1] + dbGp_eta[1] + dfGm_eta[1];
    pnts[i][j].RHS[2] = dbFp_ksi[2] + dfFm_ksi[2] + dbGp_eta[2] + dfGm_eta[2];
    pnts[i][j].RHS[3] = dbFp_ksi[3] + dfFm_ksi[3] + dbGp_eta[3] + dfGm_eta[3];

    /* Upper right point. */

    i = imax-2; j = jmax-2;

    dbFp_ksi[0] = pnts[i][j].f_plus[0] - pnts[i-1][j].f_plus[0]; 
    dbFp_ksi[1] = pnts[i][j].f_plus[1] - pnts[i-1][j].f_plus[1]; 
    dbFp_ksi[2] = pnts[i][j].f_plus[2] - pnts[i-1][j].f_plus[2]; 
    dbFp_ksi[3] = pnts[i][j].f_plus[3] - pnts[i-1][j].f_plus[3]; 

    dfFm_ksi[0] = pnts[i+1][j].f_plus[0] - pnts[i][j].f_plus[0]; 
    dfFm_ksi[1] = pnts[i+1][j].f_plus[1] - pnts[i][j].f_plus[1]; 
    dfFm_ksi[2] = pnts[i+1][j].f_plus[2] - pnts[i][j].f_plus[2]; 
    dfFm_ksi[3] = pnts[i+1][j].f_plus[3] - pnts[i][j].f_plus[3]; 

    dbGp_eta[0] = pnts[i][j].g_plus[0] - pnts[i][j-1].g_plus[0];
    dbGp_eta[1] = pnts[i][j].g_plus[1] - pnts[i][j-1].g_plus[1];
    dbGp_eta[2] = pnts[i][j].g_plus[2] - pnts[i][j-1].g_plus[2];
    dbGp_eta[3] = pnts[i][j].g_plus[3] - pnts[i][j-1].g_plus[3];

    dfGm_eta[0] = pnts[i][j+1].g_plus[0] - pnts[i][j].g_plus[0];
    dfGm_eta[1] = pnts[i][j+1].g_plus[1] - pnts[i][j].g_plus[1];
    dfGm_eta[2] = pnts[i][j+1].g_plus[2] - pnts[i][j].g_plus[2];
    dfGm_eta[3] = pnts[i][j+1].g_plus[3] - pnts[i][j].g_plus[3];

    /* Build the residues. */

    pnts[i][j].RHS[0] = dbFp_ksi[0] + dfFm_ksi[0] + dbGp_eta[0] + dfGm_eta[0];
    pnts[i][j].RHS[1] = dbFp_ksi[1] + dfFm_ksi[1] + dbGp_eta[1] + dfGm_eta[1];
    pnts[i][j].RHS[2] = dbFp_ksi[2] + dfFm_ksi[2] + dbGp_eta[2] + dfGm_eta[2];
    pnts[i][j].RHS[3] = dbFp_ksi[3] + dfFm_ksi[3] + dbGp_eta[3] + dfGm_eta[3];

    /* Upper left point. */

    i = 1; j = jmax-2;

    dbFp_ksi[0] = pnts[i][j].f_plus[0] - pnts[i-1][j].f_plus[0]; 
    dbFp_ksi[1] = pnts[i][j].f_plus[1] - pnts[i-1][j].f_plus[1]; 
    dbFp_ksi[2] = pnts[i][j].f_plus[2] - pnts[i-1][j].f_plus[2]; 
    dbFp_ksi[3] = pnts[i][j].f_plus[3] - pnts[i-1][j].f_plus[3]; 

    dfFm_ksi[0] = pnts[i+1][j].f_plus[0] - pnts[i][j].f_plus[0]; 
    dfFm_ksi[1] = pnts[i+1][j].f_plus[1] - pnts[i][j].f_plus[1]; 
    dfFm_ksi[2] = pnts[i+1][j].f_plus[2] - pnts[i][j].f_plus[2]; 
    dfFm_ksi[3] = pnts[i+1][j].f_plus[3] - pnts[i][j].f_plus[3]; 

    dbGp_eta[0] = pnts[i][j].g_plus[0] - pnts[i][j-1].g_plus[0];
    dbGp_eta[1] = pnts[i][j].g_plus[1] - pnts[i][j-1].g_plus[1];
    dbGp_eta[2] = pnts[i][j].g_plus[2] - pnts[i][j-1].g_plus[2];
    dbGp_eta[3] = pnts[i][j].g_plus[3] - pnts[i][j-1].g_plus[3];

    dfGm_eta[0] = pnts[i][j+1].g_plus[0] - pnts[i][j].g_plus[0];
    dfGm_eta[1] = pnts[i][j+1].g_plus[1] - pnts[i][j].g_plus[1];
    dfGm_eta[2] = pnts[i][j+1].g_plus[2] - pnts[i][j].g_plus[2];
    dfGm_eta[3] = pnts[i][j+1].g_plus[3] - pnts[i][j].g_plus[3];

    /* Build the residues. */

    pnts[i][j].RHS[0] = dbFp_ksi[0] + dfFm_ksi[0] + dbGp_eta[0] + dfGm_eta[0];
    pnts[i][j].RHS[1] = dbFp_ksi[1] + dfFm_ksi[1] + dbGp_eta[1] + dfGm_eta[1];
    pnts[i][j].RHS[2] = dbFp_ksi[2] + dfFm_ksi[2] + dbGp_eta[2] + dfGm_eta[2];
    pnts[i][j].RHS[3] = dbFp_ksi[3] + dfFm_ksi[3] + dbGp_eta[3] + dfGm_eta[3];

}
