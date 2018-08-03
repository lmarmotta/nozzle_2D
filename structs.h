#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/*
 * Problem definitions. This struct is not associated with each point, but with
 * each problem. 
 */

typedef struct T_DEFINE{

    int imax;
    int jmax;

    double T_t;
    double P_t;

} T_DEFINE;

/*
 * This struct holds the geometrical points read from a *.cgns file.
 */

typedef struct T_POINTS{

    /* Cartesian components. */

    double x;
    double y;

    /* Transformation jacobians. */

    double j_1;
    double jm1;

    double x_ksi;
    double y_ksi;

    double x_eta;
    double y_eta;

    /* Metric terms. */

    double ksi_x;
    double ksi_y;

    double eta_x;
    double eta_y;

} T_POINTS;
