#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"
#include "cgnslib.h"

void export_fields(T_POINTS ** pnts, T_DEFINE p_setup){

    /* Open the setup file */

    FILE * f_out = fopen("solution.dat","w");

    /* Print the header of the tecplot file. */

    fprintf(f_out,"TITLE = \" Projeto 01 \"\n");
    fprintf(f_out,"VARIABLES = \"X\",\"Y\"\n");
    fprintf(f_out,"ZONE T =\"Zone-one\", I= %d ,J= %d F=POINT\n",p_setup.imax,p_setup.jmax);
    
    int i, j;

    for (j = 0; j < p_setup.jmax; j++){
        for (i = 0; i < p_setup.imax; i++){
            fprintf(f_out,"%lf %lf\n",pnts[i][j].x,pnts[i][j].y);
        }
    }

    fclose(f_out);
}
