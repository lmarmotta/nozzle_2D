#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "cgnslib.h"
#include "structs.h"

/* Allocate the main struct of the code. */

void alloc_struct_matrix(t_points *** pntss, int imax, int jmax){
    
    int i;

    /* Malloc the first row of the memory. */

    t_points ** pnts = (t_points**)calloc(imax,sizeof(t_points*));

    /* Check if the allocation of the first row is successfull. */

    if (pnts == NULL){ 
        printf("ERROR: Memory Allocation\n"); exit(1); 
    }

    /* Now, loop through pointers initializing the columns. */

    for (i = 0; i<imax; i++){

        pnts[i] = NULL;

        pnts[i] = (t_points*)calloc(jmax,sizeof(t_points)); 

        if (pnts[i] == NULL){
            printf("ERROR: Memory Allocation\n"); exit(1); 
        }
    }

    *pntss = pnts;
}

/* Free the main struct of the code. */

void free_struct_matrix(t_points ** pnts, int imax){

    int i;

    for (i = 0; i<imax; i++) free(pnts[i]);
        
    free(pnts);
}
