#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/*
 * Problem definitions */

typedef struct T_DEFINE{

    int imax; /* Max number of points in psi */
    int jmax; /* Max number of points in eta */

    double T_t; /* Total temperature */
    double P_t; /* Total pressure */

} T_DEFINE;

/*
 * This struct holds the geometrical points read from a *.cgns file. */

typedef struct T_POINTS{

    double x; /* The x coordinate of the points. */
    double y; /* The y coordinate of the points. */

    double rho;

} T_POINTS;
