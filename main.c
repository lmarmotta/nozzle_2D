#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "cgnslib.h"
#include "structs.h"

/*
 * Prototypes of the functions used.
 */

void read_mesh_cgns(char * mesh_file_name, T_POINTS ** pnts);
void read_mesh_size(char * mesh_file_name, int * imax, int * jmax);
void read_setup(char * setup_name, T_DEFINE * p_setup);
void calc_metric_relations(T_DEFINE p_setup, T_POINTS ** pnts);
void export_fields(T_POINTS ** pnts, T_DEFINE p_setup);
void test_loop_ranges(T_DEFINE p_setup, T_POINTS ** pnts);

/* 
 * Main function
 */

int main(int argc, char * argv[]){

    /* Get the input mesh and verify if it is ok ! */

    if (argc != 2){ 
        printf("ERROR: Problem with arguments.\n");
        printf("./nozzle <mesh_file>\n"); exit(1); 
    }

    /* Prompt the user about the procedure */

    printf("\n-Processing input file: %s\n",argv[1]);

    /* Read problem setup and feed the setup struct */

    T_DEFINE p_setup;

    char * setup_name = "input.in";

    read_setup(setup_name, &p_setup);

    printf("\n--Total Temperature : %lf\n", p_setup.T_t);
    printf("\n--Total Pressure    : %lf\n", p_setup.P_t);

    /* Prompt the user with respect to the mesh file */

    printf("\n-Processing mesh file: %s\n",argv[1]);

    /* Get the problem size */

    read_mesh_size(argv[1], &p_setup.imax, &p_setup.jmax);

    printf("\n--Mesh IMAX: %d Mesh JMAX: %d\n",p_setup.imax,p_setup.jmax);

    /* Allocate the main data-structure struct. */

    int i;

    T_POINTS ** pnts = (T_POINTS**)malloc(p_setup.imax*sizeof(T_POINTS*));

    if (pnts == NULL){ 
        printf("ERROR: Memory Allocation\n"); exit(1); 
    }

    for (i = 0; i<p_setup.imax; i++){
        pnts[i] = (T_POINTS*)malloc(p_setup.jmax*sizeof(T_POINTS)); 
        if (pnts[i] == NULL){
            printf("ERROR: Memory Allocation\n"); exit(1); 
        }
    }

    read_mesh_cgns(argv[1], pnts);

    printf("\n-Mesh read successfully.\n");

    /* Now compute the proper spatial transformations. */

    printf("\n-Computing the proper transformations.\n");

    calc_metric_relations(p_setup, pnts);

    /* Export post-processor file. */

    export_fields(pnts,p_setup);

    /* Free the main struct */

    for (i = 0; i<p_setup.imax; i++) free(pnts[i]); free(pnts);

    printf("\nSUCESS: Program finalized the run !\n"); return 0;

}
