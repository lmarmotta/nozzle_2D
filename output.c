#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"
#include "cgnslib.h"
#include "externs.h"

void export_fields(t_points ** pnts, t_define p_setup){

    /* Open solution file. */

    FILE * f_out = fopen("solution.dat", "w");

    /* Check the condition of the pointer. */

    if (f_out == NULL){
        printf("\nERROR: Output file cannot be opened !\n"); exit(1);
    }

    /* Print the header of the tecplot file. */

    fprintf(f_out,"TITLE = \" Projeto 01 \"\n");
    fprintf(f_out,"VARIABLES = \"X\" , \"Y\" , \"Z\" , \"Mach number\" , \"Pressure\" , \"u\" , \"v\" , \"T\" , \"jac\" , \"jac1\", \"dt\" , \"cov_u\" , \"cov_v\" , \"rhs_0\" , \"rhs_1\" , \"rhs_2\" , \"rhs_3\"\n");
    fprintf(f_out,"ZONE T =\"Zone-one\", I= %d ,J= %d , K=1, DATAPACKING=POINT\n",p_setup.imax,p_setup.jmax);
    
    for (int j = 0; j < p_setup.jmax; j++){
        for (int i = 0; i < p_setup.imax; i++){

            /* Separate the properties. */

            double x = pnts[i][j].x;
            double y = pnts[i][j].y;
            double z = 0.0;

            double u = pnts[i][j].q_hat[1]/pnts[i][j].q_hat[0];
            double v = pnts[i][j].q_hat[2]/pnts[i][j].q_hat[0];

            double mach = pnts[i][j].m;

            double p = pnts[i][j].p;

            double t = pnts[i][j].t;

            double jac  = pnts[i][j].J;
            double jac1 = pnts[i][j].J1;

            double dtt  = pnts[i][j].dt;

            double cu = pnts[i][j].cov_u;
            double cv = pnts[i][j].cov_v;

            double rhs_0 = pnts[i][j].RHS[0];
            double rhs_1 = pnts[i][j].RHS[1];
            double rhs_2 = pnts[i][j].RHS[2];
            double rhs_3 = pnts[i][j].RHS[3];

            /* Dump solution. */

            fprintf(f_out,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                    x, y, z, mach, p, u, v, t, jac, jac1, dtt, cu, cv, rhs_0, rhs_1, rhs_2, rhs_3);
        }
    }

    fclose(f_out);
}

void dump_iteration(int iter){

    /* Dump residues to the screen. */

    printf(" -- iter: %10d | RHS: %lf | %lf | %lf | %lf\n",iter, log10(max_rhs_rho), log10(max_rhs_rhou), log10(max_rhs_rhov), log10(max_rhs_e));

}

void dump_residue_file(int iter, FILE ** res_output, t_define p_setup){

    fprintf(*res_output,"%10d %lf %lf %lf %lf\n",iter,max_rhs_rho,max_rhs_rhou,max_rhs_rhov,max_rhs_e);

    /* Check rho residue. */

    if (log10(max_rhs_rho)>p_setup.eps_blow){
        printf("\n BOOM: The code seens to be diverging. Make it better next time.\n");
        exit(1);
    }

    /* Check rhou residue. */

    if (log10(max_rhs_rhou)>p_setup.eps_blow){
        printf("\n BOOM: The code seens to be diverging. Make it better next time.\n");
        exit(1);
    }

    /* Check rhov residue. */

    if (log10(max_rhs_rhov)>p_setup.eps_blow){
        printf("\n BOOM: The code seens to be diverging. Make it better next time.\n");
        exit(1);
    }

    /* Check e residue. */

    if (log10(max_rhs_rhov)>p_setup.eps_blow){
        printf("\n BOOM: The code seens to be diverging. Make it better next time.\n");
        exit(1);
    }
}

/* Converts an integer to string. */

char * to_string(int num){

    /* Get the size of the integer. */

    int num_length = floor(log10(abs(num))) + 1;

    /* Allocate the string size. */

    char * string = (char*)calloc(num_length,sizeof(char));

    /* Convert the string. */

    sprintf(string, "%d", num);

    return string;

    free(string);
}

/* Save debug files per iter. */

void save_for_gif(int num,t_points ** pnts, t_define p_setup){

    /* Get the length of the string. */
    
    int num_length = 0; num_length = floor(log10(abs(num))) + 1;

    /* Allocate a string. */

    char string[num_length];

    /* Convert iteration number string. */

    sprintf(string, "%d", num);

    /* Allocate space for the final file */

    int length_file_char = num_length + 10;

    char final_file[length_file_char];

    /* Initialize the char array. */

    memset(final_file, 0, sizeof final_file);

    strcat(final_file,"gif_file_");
    strcat(final_file,string);
    strcat(final_file,".dat");

    /* Open the file with the propwer name. */

    FILE * f_gif = fopen(final_file,"w");

    /* Check the condition of the pointer. */

    if (f_gif == NULL){
        printf("\nERROR: Output file cannot be opened !\n"); exit(1);
    }

    /* Print the header of the tecplot file. */

    fprintf(f_gif,"TITLE = \" Projeto 01 \"\n");
    fprintf(f_gif,"VARIABLES = \"X\" , \"Y\" , \"Z\" , \"Mach number\" , \"Pressure\" , \"u\" , \"v\" , \"T\" , \"jac\" , \"jac1\", \"dt\" , \"cov_u\" , \"cov_v\" , \"rhs_0\" , \"rhs_1\" , \"rhs_2\" , \"rhs_3\"\n");
    fprintf(f_gif,"ZONE T =\"Zone-one\", I= %d ,J= %d F=POINT\n",p_setup.imax,p_setup.jmax);
    
    for (int j = 0; j < p_setup.jmax; j++){
        for (int i = 0; i < p_setup.imax; i++){

            /* Separate the properties. */

            double x = pnts[i][j].x;
            double y = pnts[i][j].y;
            double z = 0.0;

            double u = pnts[i][j].q_hat[1]/pnts[i][j].q_hat[0];
            double v = pnts[i][j].q_hat[2]/pnts[i][j].q_hat[0];

            double mach = pnts[i][j].m;

            double p = pnts[i][j].p;

            double t = pnts[i][j].t;

            double jac  = pnts[i][j].J;
            double jac1 = pnts[i][j].J1;

            double dtt  = pnts[i][j].dt;

            double cu = pnts[i][j].cov_u;
            double cv = pnts[i][j].cov_v;

            double rhs_0 = pnts[i][j].RHS[0];
            double rhs_1 = pnts[i][j].RHS[1];
            double rhs_2 = pnts[i][j].RHS[2];
            double rhs_3 = pnts[i][j].RHS[3];

            /* Dump solution. */

            fprintf(f_gif,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                    x, y, z, mach, p, u, v, t, jac, jac1, dtt, cu, cv, rhs_0, rhs_1, rhs_2, rhs_3);

        }
    }

    /* Close the file. */

    fclose(f_gif);

}
