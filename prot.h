#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/*
 * Prototypes of the functions used.
 */

void read_mesh_cgns(char * mesh_file_name, t_points ** pnts);
void read_mesh_size(char * mesh_file_name, int * imax, int * jmax);
void read_setup(char * setup_name, t_define * p_setup);
void calc_metric_relations(t_define p_setup, t_points ** pnts);
void export_fields(t_points ** pnts, t_define p_setup);
void alloc_struct_matrix(t_points *** pnts, int imax, int jmax);
void build_fluxes(t_define p_setup, t_points ** pnts);
void free_struct_matrix(t_points ** pnts, int imax);
void apply_initial_condition(t_define p_setup, t_points ** pnts);
void dump_setup(t_define p_setup);
void comp_analysis(t_points ** pnts, t_define p_setup);
void boundary_condition(t_define p_setup, t_points ** pnts);
void compute_rhs(t_define p_setup, t_points ** pnts);
void local_time(t_define p_setup, t_points ** pnts);
void rungeKuttaJST(t_define p_setup, t_points ** pnts);
void dump_iteration(int iter);
