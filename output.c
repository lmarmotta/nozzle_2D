#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"
#include "cgnslib.h"
#include "externs.h"

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
    fprintf(f_out,"VARIABLES = \"X\" , \"Y\" , \"Mach number\" , \"Pressure\" , \"u\" , \"v\"\n");
    fprintf(f_out,"ZONE T =\"Zone-one\", I= %d ,J= %d F=POINT\n",p_setup.imax,p_setup.jmax);
    
    for (int j = 0; j < p_setup.jmax; j++){
        for (int i = 0; i < p_setup.imax; i++){

            /* Separate the properties. */

            double x = pnts[i][j].x;
            double y = pnts[i][j].y;

            double u = pnts[i][j].q[1]/pnts[i][j].q[0];
            double v = pnts[i][j].q[2]/pnts[i][j].q[0];

            double mach = pnts[i][j].mach;
            double p = pnts[i][j].pressure;

            /* Dump solution. */

            fprintf(f_out,"%lf %lf %lf %lf %lf %lf\n",x, y, mach, p, u, v);
        }
    }

    fclose(f_out);
}

void comp_analysis(t_points ** pnts, t_define p_setup){

    /* Separate field limits. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Bizu variables. */

    double rho, u, v, e;

    /* Loop and calc. */

    for (int i = 0; i < imax; i++){
        for (int j = 0; j < jmax; j++){

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

void dump_iteration(int iter){

    /* Dump residues to the screen. */

    printf(" -- iter: %10d | RHS(rho) = %lf | RHS(rhou) = %lf | RHS(rhov) = %lf | RHS(e) = %lf\n", iter, 
                                                                                                 log10(max_rhs_rho ), 
                                                                                                 log10(max_rhs_rhou), 
                                                                                                 log10(max_rhs_rhov), 
                                                                                                 log10(max_rhs_e));
}

void dump_residue_file(int iter, FILE ** res_output, t_define p_setup){

    fprintf(*res_output,"%10d %lf %lf %lf %lf\n",iter,log10(max_rhs_rho),log10(max_rhs_rhou),log10(max_rhs_rhov),log10(max_rhs_e));

    if (log10(max_rhs_rho)>p_setup.eps_blow){
        printf("\n BOOM: The code seens to be diverging. Make it better next time.\n");
        exit(1);
    }
}
