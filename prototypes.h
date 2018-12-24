#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/*
 * Prototypes of the functions used.
 */

void read_mesh_cgns(char * mesh_file_name, t_define p_setup, t_points ** pnts);
void read_mesh_size(char * mesh_file_name, int * imax, int * jmax);
void read_setup(char * setup_name, t_define * p_setup);
void calc_metric_relations(t_define p_setup, t_points ** pnts);
void alloc_struct_matrix(t_points *** pnts, int imax, int jmax);
void compute_fluxes(t_define p_setup, t_points ** pnts);
void free_struct_matrix(t_points ** pnts, int imax);
void apply_initial_condition(t_define p_setup, t_points ** pnts);
void dump_setup(t_define p_setup);
void comp_properties(t_points ** pnts, t_define p_setup);
void boundary_condition_euler(t_define p_setup, t_points ** pnts);
void compute_rhs(t_define p_setup, t_points ** pnts);
void local_time(t_define p_setup, t_points ** pnts);
void explicitEuler(t_define p_setup, t_points ** pnts);
void dump_iteration(int iter, double time);
void dump_residue_file(int iter, FILE ** res_output);
void free_double_matrix(double ** matrix, int imax);
void explicitEuler(t_define p_setup, t_points ** pnts);
void art_dissip(t_define p_setup, t_points ** pnts, int d_typ);
void export_fields(t_points ** pnts, t_define p_setup);
void save_for_gif(int num,t_points ** pnts, t_define p_setup);
void initialize_structs(t_define p_setup, t_points ** pnts);
void compute_jacobian(t_define p_setup, t_points ** pnts);
void inv(double ** matrix, int N);
void blk_tri(double *** main, double *** lower, double *** upper, int size_m, int num_m, double ** XB, double ** X);
void dmgss(double ** mA, double ** mB, double ** mC, int m, int n, int p);
void dmuls(double ** mA, double ** mB, double ** mC, int size);
void free_dcube(double *** cube, int imax, int jmax);
void free_dmatrix(double ** matrix, int imax);
void free_vector(double * vector);
double * alloc_dvector(int imax);
double ** alloc_dmatrix(int imax, int jmax);
double *** alloc_dcube(int imax, int jmax, int kmax);
void free_struct_matrix(t_points ** pnts, int imax);
void beam_warming(t_define p_setup, t_points ** pnts);
void sw_residue(t_define p_setup, t_points ** pnts);
void compute_split(t_define p_setup, t_points ** pnts);
