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

void calc_metric_relations(T_DEFINE p_setup, T_POINTS ** pnts){

    /* Separate the limits of the mesh. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* In the computational space, the grid is always uniform and cartesian. */

    double dksi = 1.0;
    double deta = 1.0;

    /* Now, compute for each direction the proper derivatives. Start with the
     * internal points respecting the symmetry boundary conditions where jmax-1
     * is the centerline and the values of jmax = jmax - 2.*/

    int i, j; 

    /* Internal points. */

    for (i = 1; i < imax-1; i++){
        for (j = 1; j < jmax-1; j++){

            pnts[i][j].x_ksi = 0.5 * ( (pnts[i+1][j].x - pnts[i-1][j].x) / dksi );
            pnts[i][j].y_ksi = 0.5 * ( (pnts[i+1][j].y - pnts[i-1][j].y) / dksi );

            pnts[i][j].x_eta = 0.5 * ( (pnts[i][j+1].x - pnts[i][j-1].x) / deta );
            pnts[i][j].y_eta = 0.5 * ( (pnts[i][j+1].y - pnts[i][j-1].y) / deta );

        }
    }

    /* Inlet internal points. */

    for (j = 1; j < jmax-1; j++){

        i = 0;

        pnts[i][j].x_ksi = ( pnts[i+1][j].x - pnts[i][j].x ) / dksi;
        pnts[i][j].y_ksi = ( pnts[i+1][j].y - pnts[i][j].y ) / dksi;
        
        pnts[i][j].x_eta = 0.5 * ( (pnts[i][j+1].x - pnts[i][j-1].x) / deta );
        pnts[i][j].y_eta = 0.5 * ( (pnts[i][j+1].y - pnts[i][j-1].y) / deta );

    }

    /* Symmetry internal points. */

    for (i = 1; i < imax-1; i++){

        j = jmax-1;

        pnts[i][j].x_ksi = 0.5 * ( ( pnts[i+1][j].x - pnts[i-1][j].x ) / dksi );
        pnts[i][j].y_ksi = 0.5 * ( ( pnts[i+1][j].y - pnts[i-1][j].y ) / dksi );

        pnts[i][j].x_eta = ( pnts[i][j].x - pnts[i][j-1].x ) / dksi;
        pnts[i][j].y_eta = ( pnts[i][j].y - pnts[i][j-1].y ) / dksi;

    }

    /* Outlet internal points. */

    for (j = 1; j < jmax-1; j++){

        i = imax-1;

        pnts[i][j].x_ksi = ( pnts[i][j].x - pnts[i-1][j].x ) / dksi;
        pnts[i][j].y_ksi = ( pnts[i][j].y - pnts[i-1][j].y ) / dksi;

        pnts[i][j].x_eta = 0.5 * ( (pnts[i][j+1].x - pnts[i][j-1].x) / deta );
        pnts[i][j].y_eta = 0.5 * ( (pnts[i][j+1].y - pnts[i][j-1].y) / deta );

    }

    /* Wall internal points. */

    for (i = 1; i < imax-1; i++){

        j = 0;

        pnts[i][j].x_ksi = 0.5 * ( (pnts[i+1][j].x - pnts[i-1][j].x) / dksi );
        pnts[i][j].y_ksi = 0.5 * ( (pnts[i+1][j].y - pnts[i-1][j].y) / dksi );

        pnts[i][j].x_eta = ( pnts[i][j+1].x - pnts[i][j].x ) /deta;
        pnts[i][j].y_eta = ( pnts[i][j+1].y - pnts[i][j].y ) /deta;
    }

    /* Corner points. */

}

/* -------------------------------------------------------------------------- */
/*
 * This function is only for debug purposes. It tests the ranges of the loops
 * that will be used throughout the code in order to compute all the properties
 * of the field.  
 */

// void test_loop_ranges(T_DEFINE p_setup, T_POINTS ** pnts){

    // int i, j;

    // int imax = p_setup.imax;
    // int jmax = p_setup.jmax;

    /* Apply value through all nodes. */

    // for (i = 0; i < imax; i++){
        // for (j = 0; j < jmax; j++){
            // pnts[i][j].dummy = 1.0;
        // }
    // }

    /* Apply value only throgh the WALL. */

    // for (i = 0; i < imax; i++) pnts[i][0].dummy = 100.0;

    /* Apply values on the INLET. */

    // for (j = 0; j < jmax; j++) pnts[0][j].dummy = 100.0;

    /* Apply values on the upper SIDE. */

    // for (i = 0; i < imax; i++) pnts[i][jmax-1].dummy = 100.0;

    /* Apply values on the OUTLET. */

    // for (j = 0; j < jmax; j++) pnts[imax-1][j].dummy = 100.0;

// }
