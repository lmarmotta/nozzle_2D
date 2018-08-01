#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"
#include "cgnslib.h"

/*
 * This function reads the problem setup and feeds the proper struct
 */

void read_setup(char * setup_name, T_DEFINE * p_setup){

    char buf[30];

    FILE * fp_file = fopen(setup_name,"r");
    if (fp_file == NULL) {printf("\n ERROR: Unable to read setup file.\n"); exit(1);}

    /* Read total temperature */

    if (fscanf(fp_file,"%s %lf", buf, &p_setup->T_t) != 2)
    { printf("\nERROR: Non-conforming setup file.\n"); exit(1);}

    /* Read total pressure. */

    if (fscanf(fp_file,"%s %lf", buf, &p_setup->P_t) != 2)
    { printf("\nERROR: Non-conforming setup file.\n"); exit(1);}

    /* Close the file */

    fclose(fp_file);
}

/* 
 * This function reads the size of the problem from a *.cgns file. 
 */

void read_mesh_size(char * mesh_file_name, int * imax, int * jmax){

    /* Open cgns file and return the id */

    int cg_file;

    /* Open cgns file and return the index */

    if (cg_open(mesh_file_name,CG_MODE_READ,&cg_file)) cg_error_exit();

    /* Set the indexes as our mesh is very simple. */

    int index_base = 1; int index_zone = 1;

    /* Read the zone information. */

    char zonename[33]; cgsize_t isize[3][3];

    cg_zone_read(cg_file,index_base,index_zone,zonename,isize[0]);

    /* Set the reference indexes */

    *imax = (int)isize[0][0];
    *jmax = (int)isize[0][1];

    /* Close the cgns file index */

    cg_close(cg_file);
}

/*
 * This function reads the mesh from a *.cgns file and feeds the structured
 * data structures. 
 */

void read_mesh_cgns(char * mesh_file_name, T_POINTS ** pnts){

    /* Open cgns file and return the id */

    int cg_file;

    /* Open cgns file and return the index */

    if (cg_open(mesh_file_name,CG_MODE_READ,&cg_file)) cg_error_exit();

    /* Set the indexes as our mesh is very simple. */

    int index_base = 1; int index_zone = 1;

    /* Read the zone information. */

    char zonename[33]; cgsize_t isize[2][2],irmin[2],irmax[2];;

    cg_zone_read(cg_file,index_base,index_zone,zonename,isize[0]);

    /* Set the reference indexes */

    irmin[0]=1;
    irmin[1]=1;

    irmax[0]=isize[0][0];
    irmax[1]=isize[0][1];

    /* The cgns procedures are not Ok with structs, set some aux arrays. */
    
    double x_coord[(int)irmax[0]][(int)irmax[1]];
    double y_coord[(int)irmax[0]][(int)irmax[1]];

    /* Read the coordinates. */

    cg_coord_read(cg_file,index_base,index_zone,"CoordinateX",CGNS_ENUMV(RealDouble),irmin,irmax,x_coord[0]);
    cg_coord_read(cg_file,index_base,index_zone,"CoordinateY",CGNS_ENUMV(RealDouble),irmin,irmax,y_coord[0]);

    /* Now feed the struct. */

    int i,j;

    for (i = 0; i<irmax[0]; i++){
        for (j = 0; j<irmax[1]; j++){
            pnts[i][j].x = x_coord[i][j];
            pnts[i][j].y = y_coord[i][j];
        }
    }

    /* Close the cgns file index */

    cg_close(cg_file);

}
