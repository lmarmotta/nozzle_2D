#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "cgnslib.h"
#include "structs.h"


/* 
 * This function is concerned with the calculation and storage of the metric
 * terms which are used to pass the cartesian space to the computational space. 
 */

void calc_metric_relations(t_define p_setup, t_points ** pnts){

    /* Separate the limits of the mesh. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* In the computational space, the grid is always uniform and cartesian. */

    double dksi = 1.0;
    double deta = 1.0;

    /* Now, compute for each direction the proper derivatives. Start with the
     * internal points respecting the symmetry boundary conditions where jmax-1
     * is the centerline and the values of jmax = jmax - 2.*/

    /* Internal points. */

    for (int i = 1; i < imax-1; i++){
        for (int j = 1; j < jmax-1; j++){

            pnts[i][j].x_ksi = 0.5 * ( (pnts[i+1][j].x - pnts[i-1][j].x) / dksi );
            pnts[i][j].y_ksi = 0.5 * ( (pnts[i+1][j].y - pnts[i-1][j].y) / dksi );

            pnts[i][j].x_eta = 0.5 * ( (pnts[i][j+1].x - pnts[i][j-1].x) / deta );
            pnts[i][j].y_eta = 0.5 * ( (pnts[i][j+1].y - pnts[i][j-1].y) / deta );

        }
    }

    /* Inlet internal points. */

    for (int j = 1; j < jmax-1; j++){

        int i = 0;

        pnts[i][j].x_ksi = ( pnts[i+1][j].x - pnts[i][j].x ) / dksi;
        pnts[i][j].y_ksi = ( pnts[i+1][j].y - pnts[i][j].y ) / dksi;
        
        pnts[i][j].x_eta = 0.5 * ( (pnts[i][j+1].x - pnts[i][j-1].x) / deta );
        pnts[i][j].y_eta = 0.5 * ( (pnts[i][j+1].y - pnts[i][j-1].y) / deta );

    }

    /* Symmetry internal points. */

    for (int i = 1; i < imax-1; i++){

        int j = jmax-1;

        pnts[i][j].x_ksi = 0.5 * ( ( pnts[i+1][j].x - pnts[i-1][j].x ) / dksi );
        pnts[i][j].y_ksi = 0.5 * ( ( pnts[i+1][j].y - pnts[i-1][j].y ) / dksi );

        pnts[i][j].x_eta = ( pnts[i][j].x - pnts[i][j-1].x ) / dksi;
        pnts[i][j].y_eta = ( pnts[i][j].y - pnts[i][j-1].y ) / dksi;

    }

    /* Outlet internal points. */

    for (int j = 1; j < jmax-1; j++){

        int i = imax-1;

        pnts[i][j].x_ksi = ( pnts[i][j].x - pnts[i-1][j].x ) / dksi;
        pnts[i][j].y_ksi = ( pnts[i][j].y - pnts[i-1][j].y ) / dksi;

        pnts[i][j].x_eta = 0.5 * ( (pnts[i][j+1].x - pnts[i][j-1].x) / deta );
        pnts[i][j].y_eta = 0.5 * ( (pnts[i][j+1].y - pnts[i][j-1].y) / deta );

    }

    /* Wall internal points. */

    for (int i = 1; i < imax-1; i++){

        int j = 0;

        pnts[i][j].x_ksi = 0.5 * ( (pnts[i+1][j].x - pnts[i-1][j].x) / dksi );
        pnts[i][j].y_ksi = 0.5 * ( (pnts[i+1][j].y - pnts[i-1][j].y) / dksi );

        pnts[i][j].x_eta = ( pnts[i][j+1].x - pnts[i][j].x ) /deta;
        pnts[i][j].y_eta = ( pnts[i][j+1].y - pnts[i][j].y ) /deta;
    }

    /* Lower left corner point. */

    int i = 0; int j = 0;

    pnts[i][j].x_ksi = ( pnts[i+1][j].x - pnts[i][j].x ) / dksi;
    pnts[i][j].y_ksi = ( pnts[i+1][j].y - pnts[i][j].y ) / dksi;

    pnts[i][j].x_eta = ( pnts[i][j+1].x - pnts[i][j].x ) / deta;
    pnts[i][j].y_eta = ( pnts[i][j+1].y - pnts[i][j].y ) / deta;

    /* Lower right point. */

    i = imax-1; j = 0;

    pnts[i][j].x_ksi = ( pnts[i][j].x - pnts[i-1][j].x ) / dksi;
    pnts[i][j].y_ksi = ( pnts[i][j].y - pnts[i-1][j].y ) / dksi;

    pnts[i][j].x_eta = ( pnts[i][j+1].x - pnts[i][j].x ) / deta;
    pnts[i][j].y_eta = ( pnts[i][j+1].y - pnts[i][j].y ) / deta;

    /* Upper right point. */

    i = imax-1; j = jmax-1;

    pnts[i][j].x_ksi = ( pnts[i][j].x - pnts[i-1][j].x ) / dksi;
    pnts[i][j].y_ksi = ( pnts[i][j].y - pnts[i-1][j].y ) / dksi;

    pnts[i][j].x_eta = ( pnts[i][j].x - pnts[i][j-1].x ) / deta;
    pnts[i][j].y_eta = ( pnts[i][j].y - pnts[i][j-1].y ) / deta;

    /* Upper left point. */

    i = 0; j = jmax-1;

    pnts[i][j].x_ksi = ( pnts[i+1][j].x - pnts[i][j].x ) / dksi;
    pnts[i][j].y_ksi = ( pnts[i+1][j].y - pnts[i][j].y ) / dksi;

    pnts[i][j].x_eta = ( pnts[i][j].x - pnts[i][j-1].x ) / deta;
    pnts[i][j].y_eta = ( pnts[i][j].y - pnts[i][j-1].y ) / deta;

    /* Now compute the transformation matrix. */

    for (int i = 0; i < imax; i++){
        for (int j = 0; j < jmax; j++){

            pnts[i][j].j   = pow( (pnts[i][j].x_ksi*pnts[i][j].y_eta 
                                   - pnts[i][j].x_eta*pnts[i][j].y_ksi), -1.0);

            pnts[i][j].jm1 =  pnts[i][j].x_ksi*pnts[i][j].y_eta 
                            - pnts[i][j].x_eta*pnts[i][j].y_ksi;

        }
    }

    /* Now, compute the proper derivatives which we need to include the
     * formulation itself. */

    for (int i = 0; i < imax; i++){
        for (int j = 0; j < jmax; j++){

            pnts[i][j].ksi_x =   pnts[i][j].j*pnts[i][j].y_eta;
            pnts[i][j].ksi_y = - pnts[i][j].j*pnts[i][j].x_eta;

            pnts[i][j].eta_x = - pnts[i][j].j*pnts[i][j].y_ksi;
            pnts[i][j].eta_y =   pnts[i][j].j*pnts[i][j].x_ksi;

        }
    }
}
