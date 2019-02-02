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
}
