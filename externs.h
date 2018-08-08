#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

double max_rhs_rho;
double max_rhs_rhou;
double max_rhs_rhov;
double max_rhs_e;

extern void build_fluxes(t_define p_setup, t_points ** pnts);
extern void compute_rhs(t_define p_setup, t_points ** pnts);
