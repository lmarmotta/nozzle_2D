#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"
#include "cgnslib.h"

/*
 * This function reads the problem setup and feeds the proper struct
 */

void read_setup(char * setup_name, t_define * p_setup){

    char buf[30];

    FILE * fp_file = fopen(setup_name,"r");
    if (fp_file == NULL) {printf("\n ERROR: Unable to read setup file.\n"); exit(1);}

    /* Read the gamma constant. */

    if (fscanf(fp_file,"%s %lf", buf, &p_setup->gamma) != 2)
    { printf("\nERROR: Non-conforming setup file.\n"); exit(1);}

    /* Read the initial rho value. */

    if (fscanf(fp_file,"%s %lf", buf, &p_setup->i_rho) != 2)
    { printf("\nERROR: Non-conforming setup file.\n"); exit(1);}

    /* Read the initial rho*u value. */

    if (fscanf(fp_file,"%s %lf", buf, &p_setup->i_rhou) != 2)
    { printf("\nERROR: Non-conforming setup file.\n"); exit(1);}

    /* Read the initial rho*v value. */

    if (fscanf(fp_file,"%s %lf", buf, &p_setup->i_rhov) != 2)
    { printf("\nERROR: Non-conforming setup file.\n"); exit(1);}

    /* Read the initial e value. */

    if (fscanf(fp_file,"%s %lf", buf, &p_setup->i_e) != 2)
    { printf("\nERROR: Non-conforming setup file.\n"); exit(1);}

    /* Read the max number of iterations. */

    if (fscanf(fp_file,"%s %d", buf, &p_setup->n_max_iter) != 2)
    { printf("\nERROR: Non-conforming setup file.\n"); exit(1);}

    /* Read the DT. */

    if (fscanf(fp_file,"%s %lf", buf, &p_setup->dt) != 2)
    { printf("\nERROR: Non-conforming setup file.\n"); exit(1);}

    /* Read the BCIN_u. */

    if (fscanf(fp_file,"%s %lf", buf, &p_setup->BCIN_u) != 2)
    { printf("\nERROR: Non-conforming setup file.\n"); exit(1);}

    /* Read the BCIN_v. */

    if (fscanf(fp_file,"%s %lf", buf, &p_setup->BCIN_v) != 2)
    { printf("\nERROR: Non-conforming setup file.\n"); exit(1);}

    /* Read the BCIN_p. */

    if (fscanf(fp_file,"%s %lf", buf, &p_setup->BCIN_p) != 2)
    { printf("\nERROR: Non-conforming setup file.\n"); exit(1);}

    /* Read the BCIN_t. */

    if (fscanf(fp_file,"%s %lf", buf, &p_setup->BCIN_p) != 2)
    { printf("\nERROR: Non-conforming setup file.\n"); exit(1);}

    /* Close the file */

    fclose(fp_file);
}

void dump_setup(t_define p_setup){

    printf("\n--- Initial rho   : %lf\n",p_setup.i_rho);
    printf("\n--- Initial rho*u : %lf\n",p_setup.i_rhou);
    printf("\n--- Initial rho*v : %lf\n",p_setup.i_rhov);
    printf("\n--- Initial e     : %lf\n",p_setup.i_e);

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

void read_mesh_cgns(char * mesh_file_name, t_points ** pnts){

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

    /* The cgns procedure somehow do not like simple arrays allocated in
     * classical pointer to pointer fashion. I believe the misalignment of the
     * memory have something to do with it. Because of this limitation, a
     * different type of allocation is used here, where an array pointer is
     * used instead. It is only compatible with the -std=c99 as the concept of
     * variable-length arrays is used. */

    double(*x_coord)[(int)irmax[0]] = malloc(sizeof(double[(int)irmax[1]][(int)irmax[0]]));
    double(*y_coord)[(int)irmax[0]] = malloc(sizeof(double[(int)irmax[1]][(int)irmax[0]]));

    /* Read the coordinates. */

    cg_coord_read(cg_file,index_base,index_zone,"CoordinateX",CGNS_ENUMV(RealDouble),irmin,irmax,x_coord[0]);
    cg_coord_read(cg_file,index_base,index_zone,"CoordinateY",CGNS_ENUMV(RealDouble),irmin,irmax,y_coord[0]);

    /* Feed the structs. */

    int i,j;

    for (i = 0; i<(int)irmax[0]; i++){
        for (j = 0; j<(int)irmax[1]; j++){
            pnts[i][j].x = x_coord[j][i];
            pnts[i][j].y = y_coord[j][i];
        }
    }

    /* Close the cgns file index */

    cg_close(cg_file);

    /* We are not free from the free though... */

    free(x_coord);
    free(y_coord);

}
