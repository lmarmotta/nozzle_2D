#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/*
 * Problem definitions. This struct is not associated with each point, but with
 * each problem. 
 */

typedef struct T_DEFINE{

    int imax; /* Max number of points in psi */
    int jmax; /* Max number of points in eta */

    double T_t; /* Total temperature */
    double P_t; /* Total pressure */

} T_DEFINE;

/*
 * This struct holds the geometrical points read from a *.cgns file. */

typedef struct T_POINTS{

    double x;    // The x coordinate of the points. 
    double y;    // The y coordinate of the points.

    double dummy;

    double jm1;  // Metric Jacobian.
    double tau;  // Transformed time.

    double t_ksi; // 
    double x_ksi; // Metric terms in ksi direction.
    double y_ksi; //

    double t_eta; //
    double x_eta; // Metric terms in eta direction.
    double y_eta; //

    double rho;   // Fluid density.

} T_POINTS;
