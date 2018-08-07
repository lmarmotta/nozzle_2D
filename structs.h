#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/*
 * Problem definitions. This struct is not associated with each point, but with
 * each problem. 
 */

typedef struct t_define{

    /* Reference dimensions of the mesh. */

    int imax, jmax;

    /* Properties which come from the input file. */

    int n_max_iter;

    double T_t, P_t, gamma, i_rho, i_rhou, i_rhov, i_e, CFL, F_Cv, F_R;

    /* Pressure inlet variables. */

    double BCIN_udir, BCIN_vdir, BCIN_pt, BCIN_tt, BCIN_p;

} t_define;

/*
 * This struct holds the properties associated with each point of the mesh.
 */

typedef struct t_points{

    /* Cartesian components. */

    double x, y;

    /* Transformation jacobians. */

    double j_1, jm1;

    double x_ksi, y_ksi;

    double x_eta, y_eta;

    /* Metric terms. */

    double ksi_x, ksi_y;

    double eta_x, eta_y;

    /* Covariant velocities. */

    double cov_u, cov_v;

    /* Non-transformed. */

    double q[4];

    /* Transformed fluxes. */

    double q_hat[4], e_hat[4], f_hat[4], RHS[4];

    /* Comprehensive variables. */

    double mach, sound_speed, pressure;

    /* Time step specifics. */

    double dt;

} t_points;
