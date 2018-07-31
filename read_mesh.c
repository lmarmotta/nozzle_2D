#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "setup.h"
#include "structs.h"

/* 
 * This function reads the size of the problem from a *.cgns file. */

void read_mesh_size(char * mesh_file_name, int imax, int jmax){
}

/*
 * This function reads the mesh from a *.cgns file and feeds the structured
 * data structures. */

void read_mesh_cgns(t_points * pnts){

    pnts[0].x = 1.0;

    printf("%lf", pnts[0].x);

}
