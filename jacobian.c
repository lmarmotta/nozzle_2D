#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"

/* Compute the Jacobian matrices. */

void compute_jacobian(t_define p_setup, t_points ** pnts){

    /* Separate the limits of the mesh. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Size of the jacobian matrices. */

    int nim = 4;

    /* Store the jacobians. */

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){

            /* Separate the primitives. */

            double rho = pnts[i][j].J * pnts[i][j].q_hat[0];
            double u   = pnts[i][j].q_hat[1] / pnts[i][j].q_hat[0];
            double v   = pnts[i][j].q_hat[2] / pnts[i][j].q_hat[0];
            double e   = pnts[i][j].J * pnts[i][j].q_hat[3];

            /* For A_hat, k is the derivatives in ksi. */

            double kx = pnts[i][j].ksi_x;
            double ky = pnts[i][j].ksi_y;

            /* Compute auxiliar variables. */

            double a1    = p_setup.gamma*(e/rho);
            double theta = kx*u + ky*v;
            double phi2  = 0.5*(p_setup.gamma - 1.0)*(pow(u,2.0) + pow(v,2.0));

            /* Store the A_hat, remember that [i][j] are i = lines and j = columns. */

            /* First line. */

            pnts[i][j].A_hat[0][0] = 0.0;
            pnts[i][j].A_hat[0][1] = kx;
            pnts[i][j].A_hat[0][2] = ky;
            pnts[i][j].A_hat[0][3] = 0.0;

            /* Second line. */

            pnts[i][j].A_hat[1][0] = -u*theta + kx*phi2;
            pnts[i][j].A_hat[1][1] =  theta - (p_setup.gamma - 2.0)*kx*u;
            pnts[i][j].A_hat[1][2] =  ky*u - (p_setup.gamma - 1.0)*kx*v;
            pnts[i][j].A_hat[1][3] =  (p_setup.gamma - 1.0)*kx;

            /* Third line. */

            pnts[i][j].A_hat[2][0] = -v*theta + ky*phi2;
            pnts[i][j].A_hat[2][1] = kx*v - (p_setup.gamma - 1.0)*ky*u;
            pnts[i][j].A_hat[2][2] = theta - (p_setup.gamma - 2.0)*ky*v;
            pnts[i][j].A_hat[2][3] = (p_setup.gamma - 1.0)*ky;

            /* Fourth line. */

            pnts[i][j].A_hat[3][0] = theta*(phi2 - a1);
            pnts[i][j].A_hat[3][1] = kx*a1 - (p_setup.gamma - 1.0)*u*theta;
            pnts[i][j].A_hat[3][2] = ky*a1 - (p_setup.gamma - 1.0)*v*theta;
            pnts[i][j].A_hat[3][3] = p_setup.gamma*theta;

            /* For B_hat, k is the derivatives in ksi. */

            kx = pnts[i][j].eta_x;
            ky = pnts[i][j].eta_y;

            /* Compute auxiliar variables. */

            a1    = p_setup.gamma*(e/rho);
            theta = kx*u + ky*v;
            phi2  = 0.5*(p_setup.gamma - 1.0)*(pow(u,2.0) + pow(v,2.0));

            /* Store the B_hat, remember that [i][j] are i = lines and j = columns. */

            /* First line. */

            pnts[i][j].B_hat[0][0] = 0.0;
            pnts[i][j].B_hat[0][1] = kx;
            pnts[i][j].B_hat[0][2] = ky;
            pnts[i][j].B_hat[0][3] = 0.0;

            /* Second line. */

            pnts[i][j].B_hat[1][0] = -u*theta + kx*phi2;
            pnts[i][j].B_hat[1][1] =  theta - (p_setup.gamma - 2.0)*kx*u;
            pnts[i][j].B_hat[1][2] =  ky*u - (p_setup.gamma - 1.0)*kx*v;
            pnts[i][j].B_hat[1][3] =  (p_setup.gamma - 1.0)*kx;

            /* Third line. */

            pnts[i][j].B_hat[2][0] = -v*theta + ky*phi2;
            pnts[i][j].B_hat[2][1] = kx*v - (p_setup.gamma - 1.0)*ky*u;
            pnts[i][j].B_hat[2][2] = theta - (p_setup.gamma - 2.0)*ky*v;
            pnts[i][j].B_hat[2][3] = (p_setup.gamma - 1.0)*ky;

            /* Fourth line. */

            pnts[i][j].B_hat[3][0] = theta*(phi2 - a1);
            pnts[i][j].B_hat[3][1] = kx*a1 - (p_setup.gamma - 1.0)*u*theta;
            pnts[i][j].B_hat[3][2] = ky*a1 - (p_setup.gamma - 1.0)*v*theta;
            pnts[i][j].B_hat[3][3] = p_setup.gamma*theta;

        }
    }

    /* Now, be consistent with the coordinate transformation. */

    for (int i = 0; i<imax; i++)
        for (int j = 0; j<jmax; j++)
            for (int ii = 0; ii<nim; ii++)
                for (int jj = 0; jj<nim; jj++){
                    pnts[i][j].A_hat[ii][jj] = pnts[i][j].J1 * pnts[i][j].A_hat[ii][jj]; 
                    pnts[i][j].B_hat[ii][jj] = pnts[i][j].J1 * pnts[i][j].B_hat[ii][jj]; 
                }
}


/* Compute the spliting matrices.. */

void compute_split(t_define p_setup, t_points ** pnts){

    /* Separate the limits of the mesh. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Size of the jacobian matrices. */

    int nim = 4;

    /* Store the jacobians. */

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){

            /* Separate the primitives. */

            double rho = pnts[i][j].J * pnts[i][j].q_hat[0];
            double u   = pnts[i][j].q_hat[1] / pnts[i][j].q_hat[0];
            double v   = pnts[i][j].q_hat[2] / pnts[i][j].q_hat[0];
            // double e   = pnts[i][j].J * pnts[i][j].q_hat[3];

            /* Compute the original split matrix. */

            double kx = pnts[i][j].ksi_x;
            double ky = pnts[i][j].ksi_y;

            /* Compute auxiliar variables. */

            double a     = pnts[i][j].a;
            double alpha = rho/(pow(2.0,0.5)*a);
            double beta  = 1.0/(pow(2.0,0.5)*rho*a);
            double ktx   = kx/(pow(kx*kx + ky*ky,0.5));
            double kty   = ky/(pow(kx*kx + ky*ky,0.5));
            double theta = ktx*u + kty*v;
            double phi2  = 0.5*(p_setup.gamma-1.0)*(pow(u,2.0) + pow(v,2.0));

            /* Compute the matrix ksi direction.. */

            pnts[i][j].T_ksi[0][0] = 1.0;
            pnts[i][j].T_ksi[0][1] = 0.0;
            pnts[i][j].T_ksi[0][2] = alpha;
            pnts[i][j].T_ksi[0][3] = alpha;

            pnts[i][j].T_ksi[1][0] = u;
            pnts[i][j].T_ksi[1][1] = kty*rho;
            pnts[i][j].T_ksi[1][2] = alpha*(u + ktx*a);
            pnts[i][j].T_ksi[1][3] = alpha*(u - ktx*a);

            pnts[i][j].T_ksi[2][0] = v;
            pnts[i][j].T_ksi[2][1] = -ktx*rho;
            pnts[i][j].T_ksi[2][2] = alpha*(v + kty*a);
            pnts[i][j].T_ksi[2][3] = alpha*(v - kty*a);

            pnts[i][j].T_ksi[3][0] = phi2/(p_setup.gamma - 1.0);
            pnts[i][j].T_ksi[3][1] = rho*(kty*u - ktx*v);
            pnts[i][j].T_ksi[3][2] = alpha*((phi2 + pow(a,2.0))/(p_setup.gamma - 1.0) + a*theta);
            pnts[i][j].T_ksi[3][3] = alpha*((phi2 + pow(a,2.0))/(p_setup.gamma - 1.0) - a*theta);

            /* Compute the original split matrix. */

            kx = pnts[i][j].eta_x;
            ky = pnts[i][j].eta_y;

            /* Compute auxiliar variables. */

            a     = pnts[i][j].a;
            alpha = rho/(pow(2.0,0.5)*a);
            beta  = 1.0/(pow(2.0,0.5)*rho*a);
            ktx   = kx/(pow(kx*kx + ky*ky,0.5));
            kty   = ky/(pow(kx*kx + ky*ky,0.5));
            theta = ktx*u + kty*v;
            phi2  = 0.5*(p_setup.gamma-1.0)*(pow(u,2.0) + pow(v,2.0));

            /* Compute the matrix eta direction.. */

            pnts[i][j].T_eta[0][0] = 1.0;
            pnts[i][j].T_eta[0][1] = 0.0;
            pnts[i][j].T_eta[0][2] = alpha;
            pnts[i][j].T_eta[0][3] = alpha;

            pnts[i][j].T_eta[1][0] = u;
            pnts[i][j].T_eta[1][1] = kty*rho;
            pnts[i][j].T_eta[1][2] = alpha*(u + ktx*a);
            pnts[i][j].T_eta[1][3] = alpha*(u - ktx*a);

            pnts[i][j].T_eta[2][0] = v;
            pnts[i][j].T_eta[2][1] = -ktx*rho;
            pnts[i][j].T_eta[2][2] = alpha*(v + kty*a);
            pnts[i][j].T_eta[2][3] = alpha*(v - kty*a);

            pnts[i][j].T_eta[3][0] = phi2/(p_setup.gamma - 1.0);
            pnts[i][j].T_eta[3][1] = rho*(kty*u - ktx*v);
            pnts[i][j].T_eta[3][2] = alpha*((phi2 + pow(a,2.0))/(p_setup.gamma - 1.0) + a*theta);
            pnts[i][j].T_eta[3][3] = alpha*((phi2 + pow(a,2.0))/(p_setup.gamma - 1.0) - a*theta);

            /* Now compute the inverse of the split matrices. */

            kx = pnts[i][j].ksi_x;
            ky = pnts[i][j].ksi_y;

            /* Compute auxiliar variables. */

            a     = pnts[i][j].a;
            alpha = rho/(pow(2.0,0.5)*a);
            beta  = 1.0/(pow(2.0,0.5)*rho*a);
            ktx   = kx/(pow(kx*kx + ky*ky,0.5));
            kty   = ky/(pow(kx*kx + ky*ky,0.5));
            theta = ktx*u + kty*v;
            phi2  = 0.5*(p_setup.gamma-1.0)*(pow(u,2.0) + pow(v,2.0));

            /* Compute the matrix ksi direction.. */

            pnts[i][j].T1_ksi[0][1] = (1.0 - phi2/pow(a,2.0));
            pnts[i][j].T1_ksi[0][2] = ((p_setup.gamma - 1.0)*u)/pow(a,2.0);
            pnts[i][j].T1_ksi[0][3] = ((p_setup.gamma - 1.0)*v)/pow(a,2.0);
            pnts[i][j].T1_ksi[0][4] = -(p_setup.gamma - 1.0)/pow(a,2.0);

            pnts[i][j].T1_ksi[1][1] = -(kty*u - ktx*v)/rho;
            pnts[i][j].T1_ksi[1][2] = kty/rho;
            pnts[i][j].T1_ksi[1][3] = -ktx/rho;
            pnts[i][j].T1_ksi[1][4] = 0.0;

            pnts[i][j].T1_ksi[2][1] = beta*(phi2 - a*theta);
            pnts[i][j].T1_ksi[2][2] = beta*(ktx*a - (p_setup.gamma - 1.0)*u);
            pnts[i][j].T1_ksi[2][3] = beta*(kty*a - (p_setup.gamma - 1.0)*v);
            pnts[i][j].T1_ksi[2][4] = beta*(p_setup.gamma - 1.0);

            pnts[i][j].T1_ksi[3][1] =   beta*(phi2 + a*theta);
            pnts[i][j].T1_ksi[3][2] = - beta*(ktx*a - (p_setup.gamma - 1.0)*u);
            pnts[i][j].T1_ksi[3][3] = - beta*(kty*a - (p_setup.gamma - 1.0)*v);
            pnts[i][j].T1_ksi[3][4] =   beta*(p_setup.gamma - 1.0);

            /* Now compute the inverse of the split matrices. */

            kx = pnts[i][j].eta_x;
            ky = pnts[i][j].eta_y;

            /* Compute auxiliar variables. */

            a     = pnts[i][j].a;
            alpha = rho/(pow(2.0,0.5)*a);
            beta  = 1.0/(pow(2.0,0.5)*rho*a);
            ktx   = kx/(pow(kx*kx + ky*ky,0.5));
            kty   = ky/(pow(kx*kx + ky*ky,0.5));
            theta = ktx*u + kty*v;
            phi2  = 0.5*(p_setup.gamma-1.0)*(pow(u,2.0) + pow(v,2.0));

            /* Compute the matrix ksi direction.. */

            pnts[i][j].T1_eta[0][1] = (1.0 - phi2/pow(a,2.0));
            pnts[i][j].T1_eta[0][2] = ((p_setup.gamma - 1.0)*u)/pow(a,2.0);
            pnts[i][j].T1_eta[0][3] = ((p_setup.gamma - 1.0)*v)/pow(a,2.0);
            pnts[i][j].T1_eta[0][4] = -(p_setup.gamma - 1.0)/pow(a,2.0);

            pnts[i][j].T1_eta[1][1] = -(kty*u - ktx*v)/rho;
            pnts[i][j].T1_eta[1][2] = kty/rho;
            pnts[i][j].T1_eta[1][3] = -ktx/rho;
            pnts[i][j].T1_eta[1][4] = 0.0;

            pnts[i][j].T1_eta[2][1] = beta*(phi2 - a*theta);
            pnts[i][j].T1_eta[2][2] = beta*(ktx*a - (p_setup.gamma - 1.0)*u);
            pnts[i][j].T1_eta[2][3] = beta*(kty*a - (p_setup.gamma - 1.0)*v);
            pnts[i][j].T1_eta[2][4] = beta*(p_setup.gamma - 1.0);

            pnts[i][j].T1_eta[3][1] =   beta*(phi2 + a*theta);
            pnts[i][j].T1_eta[3][2] = - beta*(ktx*a - (p_setup.gamma - 1.0)*u);
            pnts[i][j].T1_eta[3][3] = - beta*(kty*a - (p_setup.gamma - 1.0)*v);
            pnts[i][j].T1_eta[3][4] =   beta*(p_setup.gamma - 1.0);


        }
    }

    /* Now, be consistent with the coordinate transformation. */

    for (int i = 0; i<imax; i++)
        for (int j = 0; j<jmax; j++)
            for (int ii = 0; ii<nim; ii++)
                for (int jj = 0; jj<nim; jj++){

                    pnts[i][j].T_ksi[ii][jj] = pnts[i][j].J1 * pnts[i][j].T_ksi[ii][jj]; 
                    pnts[i][j].T_eta[ii][jj] = pnts[i][j].J1 * pnts[i][j].T_eta[ii][jj]; 

                    pnts[i][j].T1_ksi[ii][jj] = pnts[i][j].J1 * pnts[i][j].T1_ksi[ii][jj]; 
                    pnts[i][j].T1_eta[ii][jj] = pnts[i][j].J1 * pnts[i][j].T1_eta[ii][jj]; 

                }

}
