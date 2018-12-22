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

    double * f_plus = alloc_dvector(4);
    double * f_minu = alloc_dvector(4);
    double * eig    = alloc_dvector(3);

    /* Loop through internal points and compute the needed stuff.*/

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            /* Separate the properties we need. */

            double rho = pnts[i][j].J * pnts[i][j].q_hat[0];
            double u   = pnts[i][j].q_hat[1] / pnts[i][j].q_hat[0];
            double v   = pnts[i][j].q_hat[2] / pnts[i][j].q_hat[0];
            double e   = pnts[i][j].J * pnts[i][j].q_hat[3];

            /* Separate useful variables. */

            double a   = pnts[i][j].a;
            double kx  = pnts[i][j].ksi_x;
            double ky  = pnts[i][j].ksi_y;
            double ktx = kx/(pow(kx*kx + ky*ky,0.5));
            double kty = ky/(pow(kx*kx + ky*ky,0.5));

            /* Store the eigenvalues for x direction. */

            double k1 = 1.0;
            double k1 = 0.0;

            eig[0] = k1*u + k2*v;
            eig[1] = k1*u + k2*v + a*pow(k1*k1 + k2*k2,0.5);
            eig[2] = k1*u + k2*v - a*pow(k1*k1 + k2*k2,0.5);

            double w2 = ( (3.0 - p_setup.gamma)*(eig[1] + eig[3])*pow(a,2.0) )/(2.0*(p_setup.gamma - 1.0));

            /* Now, build the fluxes in X direction. */

            f_plus[0] = 2.0*(p_setup.gamma)*eig[0] + eig[1] + eig[2];
            f_plus[1] = 2.0*(p_setup.gamma)*eig[0]*u + eig[1]*(u + a*ktx) + eig[2]*(u - a*ktx);
            f_plus[2] = 2.0*(p_setup.gamma)*eig[0]*v + eig[1]*(v + a*kty) + eig[2]*(v - a*kty);
            f_plus[3] = (p_setup.gamma)*eig[0]*(u*u + v*v) + 
                        (eig[1]/2.0)*(pow((u + a*ktx),2.0) + pow((v + a*kty),2.0)) + 
                        (eig[3]/2.0)*(pow((u - a*ktx),2.0) + pow((v - a*kty),2.0)) + w2;

        }
    }

    free_vector(f_plus);
    free_vector(f_minu);
    free_vector(eig);
}

void sw_implicit(t_define p_setup, t_points ** pnts){

    /* Separate bounds of the field. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;


}

