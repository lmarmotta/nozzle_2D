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

    int n_max_iter, p_rate, n_save, save_gif;

    double T_t, P_t, gamma, i_rho, i_rhou, i_rhov, i_e, CFL, F_Cv, F_R, dissp2, dissp4;

    /* Pressure inlet variables. */

    double BCIN_udir, BCIN_vdir, BCIN_pt, BCIN_tt, BCIN_p;

    /* Maximun residue of the equation. */

    double eps_blow;

    /* Artificial dissipation selection. */

    /*
     * 1: Second difference dissipation.
     * 2: Fourth order dissipation.
     * 3: Non-Linear Pullian dissipation.
     */

    int d_typ;

    /* Type of numerical scheme:
     *
     * scheme = 1 - Explicit Euler with dissip.
     * scheme = 2 - Pullian implicit scheme.
     */

    int scheme;

} t_define;

/*
 * This struct holds the properties associated with each point of the mesh.
 */

typedef struct t_points{

    /* Cartesian components. */

    double x, y;

    /* Transformation jacobians. */

    double J, J1;

    double x_ksi, y_ksi;

    double x_eta, y_eta;

    /* Metric terms. */

    double ksi_x, ksi_y;

    double eta_x, eta_y;

    /* Covariant velocities. */

    double cov_u, cov_v;

    /* Transformed fluxes. */

    double q_hat[4], e_hat[4], f_hat[4], RHS[4];

    /* Comprehensive variables. */

    double a;  /* Speed of sound. */

    double m;  /* Mach number. */

    double p;  /* Static pressure. */

    double t;  /* Static temperature. */

    /* Time step specifics. */

    double dt;

    /* Inviscid Jacobians. */

    double A_hat[4][4], B_hat[4][4];

    /* Split matrices. */

    double T_ksi[4][4], T_eta[4][4];

    /* Inverse of the split matrices. */

    double T1_ksi[4][4], T1_eta[4][4];

    /* Artificial dissipation. */

    double diss_ksi[4], diss_eta[4];

    /* Fluxes separated for SW. */

    double f_plus[4], f_minu[4];
    double g_plus[4], g_minu[4];

} t_points;
