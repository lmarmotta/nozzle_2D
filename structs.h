#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/*
 * Problem definitions */

typedef struct p_define{

    int imax; /* Max number of points in psi */
    int jmax; /* Max number of points in eta */

    double T_t; /* Total temperature */

} p_define;

/*
 * This struct holds the geometrical points read from a *.cgns file. */

typedef struct t_points{

    double x; /* The x coordinate of the points. */
    double y; /* The y coordinate of the points. */

} t_points;
