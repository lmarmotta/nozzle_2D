#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "cgnslib.h"
#include "structs.h"
#include "prototypes.h"

/* 
 * Main function
 */

int main(int argc, char * argv[]){

    /* Get the input mesh and verify if it is ok ! */

    if (argc != 3){ 
        printf("ERROR: Problem with arguments.\n");
        printf("./nozzle <input file> <mesh_file>\n"); exit(1); 
    }

    /* Read problem setup and feed the setup struct */

    t_define p_setup;

    /* Prompt the user about the procedure */

    printf("\n-Processing input file: %s\n",argv[1]);

    read_setup(argv[1], &p_setup);

    dump_setup(p_setup);

    /* Prompt the user with respect to the mesh file */

    printf("\n-Processing mesh file: %s\n",argv[2]);

    /* Get the problem size */

    read_mesh_size(argv[2], &p_setup.imax, &p_setup.jmax);

    printf("\n--Mesh IMAX: %d Mesh JMAX: %d\n",p_setup.imax,p_setup.jmax);

    /* Allocate the main data-structure struct. */

    t_points ** pnts = NULL;

    alloc_struct_matrix(&pnts,p_setup.imax,p_setup.jmax); 

    initialize_structs(p_setup, pnts);

    /* Read the whole mesh and feed the structs. */

    read_mesh_cgns(argv[2], p_setup, pnts);

    /* Apply the initial condition conditions. */

    printf("\n-Performing pre-processing computations....\n");

    /* Now compute the proper spatial transformations. */

    calc_metric_relations(p_setup, pnts);

    /* Apply initial condition. */

    apply_initial_condition(p_setup, pnts);

    /* Perform 0 iteration procedure. */

    printf("\n-Performing Zero iter procedures computations....\n");

    /* Build the fluxes. */

    build_fluxes(p_setup, pnts);

    /* Compute the RHS. */

    compute_rhs(p_setup, pnts);

    /* Apply the initial boundary conditions. */

    boundary_condition_euler(p_setup, pnts);

    /* Compute local time step for our mesh. */

    local_time(p_setup, pnts);

    /* Start the iterative procedure. */

    printf("\n-Starting iterations.\n\n");

    FILE * res_output = fopen("residues.dat", "w");

    int out_rate = p_setup.p_rate;
    int out_save = p_setup.n_save;

    for (int iter = 1; iter <= p_setup.n_max_iter; iter++){

        /* Now, build the fluxes. */

        build_fluxes(p_setup, pnts);

        /* Compute the RHS. */

        compute_rhs(p_setup, pnts);

        /* Apply dissipation. */

        art_dissip_2nd(p_setup, pnts);

        /* Marh in time. */

        explicitEuler(p_setup, pnts);

        /* Update the time. */

        local_time(p_setup, pnts);

        /* Update Boundary conditions. */

        boundary_condition_euler(p_setup, pnts);

        /* Dump a whole lotta of stuff. */

        dump_residue_file(iter, &res_output);

        if (iter == out_rate){
            dump_iteration(iter);
            out_rate += p_setup.p_rate;
        }

        if (iter == out_save){
            printf("\n Outputing solution. \n\n");
            if (p_setup.save_gif == 1) save_for_gif(iter, pnts, p_setup); 
            export_fields(pnts,p_setup);
            out_save += p_setup.n_save;
            
        }
    }

    
    /* Close the residue file. */

    fclose(res_output);

    /* Export post-processor file. */

    printf("\n-Output solution.\n");

    export_fields(pnts,p_setup);

    /* Free the main struct */

    printf("\n-Free memory.\n");

    free_struct_matrix(pnts, p_setup.imax);

    printf("\n\n +++ SUCESS: Program finalized the run ! +++\n\n"); return 0;

}
