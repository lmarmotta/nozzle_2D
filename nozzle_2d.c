#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "cgnslib.h"

#include "setup.h"
#include "structs.h"

/*
 * Prototypes of the functions used. */

void read_mesh_cgns(t_points * pnts);
void read_mesh_size(char * mesh_file_name, int imax, int jmax);

int main(int argc, char * argv[]){

    /*
     * Get the input information and verify if it is ok ! */

    if (argc != 2){ 
        printf("ERROR: Problem with arguments.\n");
        printf("./nozzle <mesh_file>\n"); exit(1); 
    }

    /*
     * Gather the points for our mesh. */

    /* Prompt the user with respect to the mesh file */

    printf("\n - Processing mesh file: %s\n",argv[1]);

    /* Set the pointer to the proper struct */

    t_points * pnts = NULL;

    int size = 10;

    pnts = (t_points*)malloc(size*sizeof(t_points));

    read_mesh_cgns(pnts);

    printf("%lf", pnts[0].x);

    /* 
     * Deallocate the structs */

    free(pnts);

    return 0;

}
