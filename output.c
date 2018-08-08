#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"
#include "cgnslib.h"

void export_fields(t_points ** pnts, t_define p_setup){

    /* Open the setup file */

    FILE * f_out = NULL;
        
    f_out = fopen("solution.dat","w");

    /* Check the condition of the pointer. */

    if (f_out == NULL){
        printf("\nERROR: Output file cannot be opened !\n"); exit(1);
    }

    /* Print the header of the tecplot file. */

    fprintf(f_out,"TITLE = \" Projeto 01 \"\n");
    fprintf(f_out,"VARIABLES = \"X\" , \"Y\" , \"Mach number\" , \"Pressure\" , \"Cov U\" , \"Cov V\"\n");
    fprintf(f_out,"ZONE T =\"Zone-one\", I= %d ,J= %d F=POINT\n",p_setup.imax,p_setup.jmax);
    
    int i, j;

    for (j = 0; j < p_setup.jmax; j++){
        for (i = 0; i < p_setup.imax; i++){
            fprintf(f_out,"%lf %lf %lf %lf %lf %lf\n",pnts[i][j].x, 
                                                      pnts[i][j].y,
                                                      pnts[i][j].mach,
                                                      pnts[i][j].pressure,
                                                      pnts[i][j].cov_u,
                                                      pnts[i][j].cov_v);
        }
    }

    fclose(f_out);
}

void comp_analysis(t_points ** pnts, t_define p_setup){

    int i, j;

    /* Separate field limits. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Bizu variables. */

    double rho, u, v, e;

    /* Loop and calc. */

    for (i = 0; i < imax; i++){
        for (j = 0; j < jmax; j++){

            /* Separate the needed properties. */

            rho = pnts[i][j].q[0];
            u   = pnts[i][j].q[1] / pnts[i][j].q[0];
            v   = pnts[i][j].q[2] / pnts[i][j].q[0];
            e   = pnts[i][j].q[3];

            /* Compute the static pressure. */

            pnts[i][j].pressure = (p_setup.gamma - 1.0)*(e - 0.5*rho*( pow(u,2.0) + pow(v,2.0)) );

            /* Compute the velocity of sound. */

            pnts[i][j].sound_speed = sqrt( (p_setup.gamma*pnts[i][j].pressure)/rho );

            /* Compute the Mach number. */

            double U = sqrt( pow(u,2.0) + pow(v,2.0) );

            pnts[i][j].mach = U/pnts[i][j].sound_speed;
        }
    }
}

