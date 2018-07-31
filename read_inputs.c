#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"
#include "cgnslib.h"

/*
 * This function reads the problem setup and feeds the proper struct */

void read_setup(char * setup_name, p_define * p_setup){

    char buf[30];

    FILE * fp_file = fopen(setup_name,"r");
    if (fp_file == NULL) {printf("\n ERROR: Unable to read setup file.\n"); exit(1);}

    /* Read total temperature */

    if (fscanf(fp_file,"%s %lf", buf, &p_setup->T_t) != 2)
    { printf("\nERROR: Non-conforming setup file.\n"); exit(1);}

    fclose(fp_file);

}

/* 
 * This function reads the size of the problem from a *.cgns file. */

void read_mesh_size(char * mesh_file_name, int * imax, int * jmax){

    /* Open cgns file and return the id */

    int cg_file;

    if ( cg_open(  mesh_file_name, CG_MODE_READ, &cg_file ) ) cg_error_exit();

    /* Close the cgns file index */

    cg_close(cg_file);

    /* Pass the number of maximun points to the function */

    *imax = -1; *jmax = -2;
}

/*
 * This function reads the mesh from a *.cgns file and feeds the structured
 * data structures. */

void read_mesh_cgns(t_points * pnts){

    pnts[0].x = 1.0;

}
