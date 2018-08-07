#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "cgnslib.h"
#include "structs.h"
#include "prot.h"

/* 
 * Main function
 */

int main(int argc, char * argv[]){

    /* Get the input mesh and verify if it is ok ! */

    if (argc != 2){ 
        printf("ERROR: Problem with arguments.\n");
        printf("./nozzle <mesh_file>\n"); exit(1); 
    }

    /* Read problem setup and feed the setup struct */

    t_define p_setup;

    char * setup_name = "input.in";

    /* Prompt the user about the procedure */

    printf("\n-Processing input file: %s\n",setup_name);

    read_setup(setup_name, &p_setup);

    dump_setup(p_setup);

    /* Prompt the user with respect to the mesh file */

    printf("\n-Processing mesh file: %s\n",argv[1]);

    /* Get the problem size */

    read_mesh_size(argv[1], &p_setup.imax, &p_setup.jmax);

    printf("\n--Mesh IMAX: %d Mesh JMAX: %d\n",p_setup.imax,p_setup.jmax);

    /* Allocate the main data-structure struct. */

    t_points ** pnts = NULL;

    alloc_struct_matrix(&pnts,p_setup.imax,p_setup.jmax); 

    /* Read the whole mesh and feed the structs. */

    read_mesh_cgns(argv[1], pnts);

    /* Apply the initial condition. */

    printf("\n-Computing the initial condition.\n");

    apply_initial_condition(p_setup, pnts);

    /* Now compute the proper spatial transformations. */

    printf("\n-Computing the proper transformations.\n");

    calc_metric_relations(p_setup, pnts);


    printf("\n-Starting iterations.\n\n");

    int iter;

    /* Apply the initial boundary conditions. */

    boundary_condition(p_setup, pnts);

    /* Start the iterative procedure. */

    for (iter = 0; iter< p_setup.n_max_iter; iter++){

        /* Now, build the RHS_ij. */

        build_fluxes(p_setup, pnts);
        compute_rhs(p_setup, pnts);
        local_time(p_setup, pnts);
        boundary_condition(p_setup, pnts);

        /* Prompt the iteration number. */

        printf(" -- > iter: %d | RHS[1] = \n",iter);
    }
    

    /* Export post-processor file. */

    printf("\n-Output solution.\n");

    comp_analysis(pnts, p_setup);

    export_fields(pnts,p_setup);

    /* Free the main struct */

    printf("\n-Free memory.\n");

    free_struct_matrix(pnts, p_setup.imax);

    printf("\n\n +++ SUCESS: Program finalized the run ! +++\n\n"); return 0;

}
