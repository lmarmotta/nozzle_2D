#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "cgnslib.h"

#include "structs.h"

/*
 * Prototypes of the functions used.
 */

void read_mesh_cgns(T_POINTS * pnts);
void read_mesh_size(char * mesh_file_name, int * imax, int * jmax);
void read_setup(char * setup_name, T_DEFINE * p_setup);

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

    /* Set the pointer to the proper struct */

    T_POINTS * pnts = NULL;

    /* Get the problem size */

    read_mesh_size(argv[1], &p_setup.imax, &p_setup.jmax);

    printf("%d %d\n",p_setup.imax,p_setup.jmax);

    int size = 10;

    pnts = (T_POINTS*)malloc(size*sizeof(T_POINTS));

    read_mesh_cgns(pnts);

    /* Deallocate the structs */

    free(pnts);

    return 0;

}
