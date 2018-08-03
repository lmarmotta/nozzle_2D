#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "cgnslib.h"
#include "structs.h"

/* 
 * Safely allocate the main structs of the program. Needs to be done safely in
 * order to avoid valgrind errors.
 */

void alloc_struct_matrix(T_POINTS *** pntss, int imax, int jmax){
    
    int i;

    /* Malloc the first row of the memory. */

    T_POINTS ** pnts = (T_POINTS**)calloc(imax,sizeof(T_POINTS*));

    /* Check if the allocation of the first row is successfull. */

    if (pnts == NULL){ 
        printf("ERROR: Memory Allocation\n"); exit(1); 
    }

    /* Now, loop through pointers initializing the columns. */

    for (i = 0; i<imax; i++){

        pnts[i] = NULL;

        pnts[i] = (T_POINTS*)calloc(jmax,sizeof(T_POINTS)); 

        if (pnts[i] == NULL){
            printf("ERROR: Memory Allocation\n"); exit(1); 
        }
    }

    *pntss = pnts;

}
