#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"
#include "cgnslib.h"

void export_fields(T_POINTS ** pnts, T_DEFINE p_setup){

    /* Open the setup file */

    FILE * f_out = NULL;
        
    f_out = fopen("solution.dat","w");

    /* Check the condition of the pointer. */

    if (f_out == NULL){
        printf("\nERROR: Output file cannot be opened !\n"); exit(1);
    }

    /* Print the header of the tecplot file. */

    fprintf(f_out,"TITLE = \" Projeto 01 \"\n");
    fprintf(f_out,"VARIABLES = \"X\" , \"Y\" , \"x_ksi\", \"y_ksi\", \"x_eta\", \"y_eta\" \n");
    fprintf(f_out,"ZONE T =\"Zone-one\", I= %d ,J= %d F=POINT\n",p_setup.imax,p_setup.jmax);
    
    int i, j;

    for (j = 0; j < p_setup.jmax; j++){
        for (i = 0; i < p_setup.imax; i++){
            fprintf(f_out,"%lf %lf %lf %lf %lf %lf\n",pnts[i][j].x, 
                                                      pnts[i][j].y, 
                                                      pnts[i][j].x_ksi, 
                                                      pnts[i][j].y_ksi, 
                                                      pnts[i][j].x_eta, 
                                                      pnts[i][j].y_eta);
        }
    }

    fclose(f_out);
}
