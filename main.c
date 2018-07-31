#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "cgnslib.h"

#include "structs.h"

/*
 * Prototypes of the functions used. */

void read_mesh_cgns(t_points * pnts);
void read_mesh_size(char * mesh_file_name, int * imax, int * jmax);
void read_setup(char * setup_name, p_define * p_setup);

int main(int argc, char * argv[]){

    /*
     * Get the input information and verify if it is ok ! */

    if (argc != 2){ 
        printf("ERROR: Problem with arguments.\n");
        printf("./nozzle <mesh_file>\n"); exit(1); 
    }

    /* Read problem setup and feed the setup struct */

    p_define p_setup;

    char * setup_name = "input.in";

    read_setup(setup_name, &p_setup);

    printf("%lf\n", p_setup.T_t);

    /* Prompt the user with respect to the mesh file */

    printf("\n - Processing mesh file: %s\n",argv[1]);

    /* Set the pointer to the proper struct */

    t_points * pnts = NULL;

    /* Get the problem size */

    read_mesh_size(argv[1], &p_setup.imax, &p_setup.jmax);

    printf("%d %d\n",p_setup.imax,p_setup.jmax);

    int size = 10;

    pnts = (t_points*)malloc(size*sizeof(t_points));

    read_mesh_cgns(pnts);

    /* 
     * Deallocate the structs */

    free(pnts);

    return 0;

}
