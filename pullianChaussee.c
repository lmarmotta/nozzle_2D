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


/* Computes the pullian and Chausee matrices. */

void compute_pc_matrices(t_define p_setup, t_points ** pnts){

    /* Separate bounds of the field. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Loop through internal points and compute the needed stuff.*/

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){

            /* Compute the matrices direction. */

            double mu = 1.0/pow(2.0,0.5);

            double m1 = pnts[i][j].ksi_x*pnts[i][j].eta_x + pnts[i][j].ksi_y*pnts[i][j].eta_y;

            double m2 = pnts[i][j].ksi_x*pnts[i][j].eta_y + pnts[i][j].ksi_y*pnts[i][j].eta_x;

            /* Compute the matrix. */

            pnts[i][j].nh[0][0] = 1.0;
            pnts[i][j].nh[0][1] = 0.0; 
            pnts[i][j].nh[0][2] = 0.0; 
            pnts[i][j].nh[0][3] = 0.0; 

            pnts[i][j].nh[1][0] = 0.0;
            pnts[i][j].nh[1][1] = m1;
            pnts[i][j].nh[1][2] = -mu*m2;
            pnts[i][j].nh[1][3] = mu*m2;

            pnts[i][j].nh[2][0] = 0.0;
            pnts[i][j].nh[2][1] = mu*m2;
            pnts[i][j].nh[2][2] = pow(mu,2.0)*(1.0 + m1);
            pnts[i][j].nh[2][3] = pow(mu,2.0)*(1.0 - m1);

            pnts[i][j].nh[3][0] = 0.0;
            pnts[i][j].nh[3][1] = -mu*m2;
            pnts[i][j].nh[3][2] = pow(mu,2.0)*(1.0 - m1);
            pnts[i][j].nh[3][3] = pow(mu,2.0)*(1.0 + m1);

            /* Compute the inverse of matrix. */

            pnts[i][j].nhi[0][0] = 1.0;
            pnts[i][j].nhi[0][1] = 0.0; 
            pnts[i][j].nhi[0][2] = 0.0; 
            pnts[i][j].nhi[0][3] = 0.0; 

            pnts[i][j].nhi[1][0] = 0.0;
            pnts[i][j].nhi[1][1] = m1;
            pnts[i][j].nhi[1][2] = mu*m2;
            pnts[i][j].nhi[1][3] = -mu*m2;

            pnts[i][j].nhi[2][0] = 0.0;
            pnts[i][j].nhi[2][1] = -mu*m2;
            pnts[i][j].nhi[2][2] = pow(mu,2.0)*(1.0 + m1);
            pnts[i][j].nhi[2][3] = pow(mu,2.0)*(1.0 - m1);

            pnts[i][j].nhi[3][0] = 0.0;
            pnts[i][j].nhi[3][1] = mu*m2;
            pnts[i][j].nhi[3][2] = pow(mu,2.0)*(1.0 - m1);
            pnts[i][j].nhi[3][3] = pow(mu,2.0)*(1.0 + m1);


        }
    }

    /* Compute the diagonal eigenvalue matrices. */

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){

            /* Bizu vars. */

            U = pnts[i][j].cov_u;
            V = pnts[i][j].cov_u;

            a = pnts[i][j].a;

            ksix = pnts[i][j].ksi_x;
            ksiy = pnts[i][j].ksi_y;

            etax = pnts[i][j].eta_x;
            etay = pnts[i][j].eta_y;

            /* Ksi direction. */

            pnts[i][j].lambda_ksi[0][0] = U;
            pnts[i][j].lambda_ksi[0][1] = 0.0; 
            pnts[i][j].lambda_ksi[0][2] = 0.0; 
            pnts[i][j].lambda_ksi[0][3] = 0.0; 

            pnts[i][j].lambda_ksi[1][0] = 0.0; 
            pnts[i][j].lambda_ksi[1][1] = U;
            pnts[i][j].lambda_ksi[1][2] = 0.0; 
            pnts[i][j].lambda_ksi[1][3] = 0.0; 

            pnts[i][j].lambda_ksi[2][0] = 0.0; 
            pnts[i][j].lambda_ksi[2][1] = 0.0; 
            pnts[i][j].lambda_ksi[2][2] = U + a*pow(pow(ksix,2) + pow(ksiy,2.0),0.5);
            pnts[i][j].lambda_ksi[2][3] = 0.0; 

            pnts[i][j].lambda_ksi[3][0] = 0.0; 
            pnts[i][j].lambda_ksi[3][1] = 0.0; 
            pnts[i][j].lambda_ksi[3][2] = 0.0; 
            pnts[i][j].lambda_ksi[3][3] = U - a*pow(pow(ksix,2) + pow(ksiy,2.0),0.5);

            /* Eta direction. */

            pnts[i][j].lambda_eta[0][0] = V;
            pnts[i][j].lambda_eta[0][1] = 0.0; 
            pnts[i][j].lambda_eta[0][2] = 0.0; 
            pnts[i][j].lambda_eta[0][3] = 0.0; 

            pnts[i][j].lambda_eta[1][0] = 0.0; 
            pnts[i][j].lambda_eta[1][1] = V;
            pnts[i][j].lambda_eta[1][2] = 0.0; 
            pnts[i][j].lambda_eta[1][3] = 0.0; 

            pnts[i][j].lambda_eta[2][0] = 0.0; 
            pnts[i][j].lambda_eta[2][1] = 0.0; 
            pnts[i][j].lambda_eta[2][2] = V + a*pow(pow(etax,2) + pow(etay,2.0),0.5);
            pnts[i][j].lambda_eta[2][3] = 0.0; 

            pnts[i][j].lambda_eta[3][0] = 0.0; 
            pnts[i][j].lambda_eta[3][1] = 0.0; 
            pnts[i][j].lambda_eta[3][2] = 0.0; 
            pnts[i][j].lambda_eta[3][3] = V - a*pow(pow(etax,2) + pow(etay,2.0),0.5);


        }
    }

}

/* Compute the pullian chausse implicit operator. */

void imp_pullian_chaussee(t_define p_setup, t_points ** pnts){

    /* Separate bounds of the field. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Loop through internal points and compute the needed stuff.*/

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){
        }
    }
}
