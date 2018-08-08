#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

extern void build_fluxes(t_define p_setup, t_points ** pnts);
extern void compute_rhs(t_define p_setup, t_points ** pnts);
