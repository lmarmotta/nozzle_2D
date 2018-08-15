#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/*
 * Prototypes of the functions used.
 */

void export_fields(t_points ** pnts, t_define p_setup, FILE ** f_out);
void read_mesh_cgns(char * mesh_file_name, t_points ** pnts);
void read_mesh_size(char * mesh_file_name, int * imax, int * jmax);
void read_setup(char * setup_name, t_define * p_setup);
void calc_metric_relations(t_define p_setup, t_points ** pnts);
void alloc_struct_matrix(t_points *** pnts, int imax, int jmax);
void build_fluxes(t_define p_setup, t_points ** pnts);
void free_struct_matrix(t_points ** pnts, int imax);
void apply_initial_condition(t_define p_setup, t_points ** pnts);
void dump_setup(t_define p_setup);
void comp_properties(t_points ** pnts, t_define p_setup);
void boundary_condition_euler(t_define p_setup, t_points ** pnts);
void compute_rhs(t_define p_setup, t_points ** pnts);
void local_time(t_define p_setup, t_points ** pnts);
void explicitEuler(t_define p_setup, t_points ** pnts);
void dump_iteration(int iter);
void dump_residue_file(int iter, FILE ** res_output);
double ** alloc_double_matrix(int imax, int jmax);
void free_double_matrix(double ** matrix, int imax);
void explicitEuler(t_define p_setup, t_points ** pnts);
void art_dissip_2nd(t_define p_setup, t_points ** pnts);
