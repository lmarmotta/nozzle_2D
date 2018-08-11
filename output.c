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
    fprintf(f_out,"VARIABLES = \"X\" , \"Y\" , \"Mach number\" , \"Pressure\" , \"u\" , \"v\", \"T\" , \"RHS[0]\" , \"RHS[1]\" , \"RHS[2]\" , \"RHS[3]\"\n");
    fprintf(f_out,"ZONE T =\"Zone-one\", I= %d ,J= %d F=POINT\n",p_setup.imax,p_setup.jmax);
    
    for (int j = 0; j < p_setup.jmax; j++){
        for (int i = 0; i < p_setup.imax; i++){

            /* Separate the properties. */

            double x = pnts[i][j].x;
            double y = pnts[i][j].y;

            double u = pnts[i][j].q[1]/pnts[i][j].q[0];
            double v = pnts[i][j].q[2]/pnts[i][j].q[0];

            double mach = pnts[i][j].m;

            double p = pnts[i][j].p;

            double t = pnts[i][j].t;

            double rhs_0 = log10(pnts[i][j].RHS[0]);
            double rhs_1 = log10(pnts[i][j].RHS[1]);
            double rhs_2 = log10(pnts[i][j].RHS[2]);
            double rhs_3 = log10(pnts[i][j].RHS[3]);

            /* Dump solution. */

            fprintf(f_out,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                    x, y, mach, p, u, v, t,rhs_0,rhs_1,rhs_2,rhs_3);
        }
    }



    fclose(f_out);
}

void dump_iteration(int iter){

    /* Dump residues to the screen. */

    printf(" -- iter: %10d | RHS(rho) = %lf | RHS(rhou) = %lf | RHS(rhov) = %lf | RHS(e) = %lf\n", iter, 
                                                                                                 max_rhs_rho,
                                                                                                 max_rhs_rhou, 
                                                                                                 max_rhs_rhov, 
                                                                                                 max_rhs_e);
}

void dump_residue_file(int iter, FILE ** res_output, t_define p_setup){

    fprintf(*res_output,"%10d %lf %lf %lf %lf\n",iter,max_rhs_rho,max_rhs_rhou,max_rhs_rhov,max_rhs_e);

    if (log10(max_rhs_rho)>p_setup.eps_blow){
        printf("\n BOOM: The code seens to be diverging. Make it better next time.\n");
        exit(1);
    }
}
