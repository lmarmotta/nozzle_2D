#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/*
 * ----------------------------------
 * Prototypes of the functions used.
 * ----------------------------------
 */


/* Pre-proc. */

void calc_metric_relations(t_define p_setup, t_points ** pnts);
void apply_initial_condition(t_define p_setup, t_points ** pnts);
void comp_properties(t_points ** pnts, t_define p_setup);
void boundary_condition_euler(t_define p_setup, t_points ** pnts);
void local_time(t_define p_setup, t_points ** pnts);
void initialize_structs(t_define p_setup, t_points ** pnts);

/* Linear algebra. */

void inv(double ** matrix, int N);
void dmgss(double ** mA, double ** mB, double ** mC, int m, int n, int p);
void dmuls(double ** mA, double ** mB, double ** mC, int size);
void blk_tri(double *** main, double *** lower, double *** upper, int size_m, 
        int num_m, double ** XB, double ** X);

/* Memory managment. */

void alloc_struct_matrix(t_points *** pnts, int imax, int jmax);
double * alloc_dvector(int imax);
double ** alloc_dmatrix(int imax, int jmax);
double *** alloc_dcube(int imax, int jmax, int kmax);
void free_dcube(double *** cube, int imax, int jmax);
void free_dmatrix(double ** matrix, int imax);
void free_vector(double * vector);
void free_struct_matrix(t_points ** pnts, int imax);

/* File and/or I/O */

void read_mesh_cgns(char * mesh_file_name, t_define p_setup, t_points ** pnts);
void read_mesh_size(char * mesh_file_name, int * imax, int * jmax);
void read_setup(char * setup_name, t_define * p_setup);
void dump_setup(t_define p_setup);
void dump_iteration(int iter, double time);
void save_for_gif(int num,t_points ** pnts, t_define p_setup);
void dump_residue_file(int iter, FILE ** res_output, t_define p_setup);
void export_fields(t_points ** pnts, t_define p_setup);

/* Centered specifics. */

void compute_fluxes(t_define p_setup, t_points ** pnts);
void compute_rhs(t_define p_setup, t_points ** pnts);
void explicitEuler(t_define p_setup, t_points ** pnts);
void art_dissip_2nd(t_define p_setup, t_points ** pnts);
void art_dissip_nlin(t_define p_setup, t_points ** pnts);

/* Beam Warming specifics. */

void compute_jacobian(t_define p_setup, t_points ** pnts);
void beam_warming(t_define p_setup, t_points ** pnts);
void compute_splited_jacobians(t_define p_setup, t_points ** pnts);

/* Steger Warming specifics. */

void compute_sw_fluxes(t_define p_setup, t_points ** pnts);
void compute_sw_residue_1sto(t_define p_setup, t_points ** pnts);
void compute_sw_residue_2ndo(t_define p_setup, t_points ** pnts);
void compute_sw_impicit_operator(t_define p_setup, t_points ** pnts);
