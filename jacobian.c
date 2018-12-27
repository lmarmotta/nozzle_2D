#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "structs.h"
#include "externs.h"
#include "prototypes.h"

/* Compute the Jacobian matrices. */

void compute_jacobian(t_define p_setup, t_points ** pnts){

    /* Separate the limits of the mesh. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

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

    int nim = 4;

    for (int i = 0; i<imax; i++)
        for (int j = 0; j<jmax; j++)
            for (int ii = 0; ii<nim; ii++)
                for (int jj = 0; jj<nim; jj++){
                    pnts[i][j].A_hat[ii][jj] = pnts[i][j].J1 * pnts[i][j].A_hat[ii][jj]; 
                    pnts[i][j].B_hat[ii][jj] = pnts[i][j].J1 * pnts[i][j].B_hat[ii][jj]; 
                }
}


/* Compute the spliting matrices.. */

void compute_splited_jacobians(t_define p_setup, t_points ** pnts){

    /* Separate the limits of the mesh. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Store Q matrices, or T matrices depending on your notation. */

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){

            /* Separate the primitives. */

            double rho = pnts[i][j].J * pnts[i][j].q_hat[0];
            double u   = pnts[i][j].q_hat[1] / pnts[i][j].q_hat[0];
            double v   = pnts[i][j].q_hat[2] / pnts[i][j].q_hat[0];

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

            beta  = 1.0/(pow(2.0,0.5)*rho*a);
            ktx   = kx/(pow(kx*kx + ky*ky,0.5));
            kty   = ky/(pow(kx*kx + ky*ky,0.5));
            theta = ktx*u + kty*v;
            phi2  = 0.5*(p_setup.gamma-1.0)*(pow(u,2.0) + pow(v,2.0));

            /* Compute the matrix ksi direction.. */

            pnts[i][j].T1_ksi[0][0] =  1.0 - (phi2/pow(a,2.0));
            pnts[i][j].T1_ksi[0][1] =  (p_setup.gamma - 1.0)*u/pow(a,2.0);
            pnts[i][j].T1_ksi[0][2] =  (p_setup.gamma - 1.0)*v/pow(a,2.0);
            pnts[i][j].T1_ksi[0][3] = -(p_setup.gamma - 1.0)/pow(a,2.0);

            pnts[i][j].T1_ksi[1][0] = -(kty*u - ktx*v)/rho;
            pnts[i][j].T1_ksi[1][1] =  kty/rho;
            pnts[i][j].T1_ksi[1][2] = -ktx/rho;
            pnts[i][j].T1_ksi[1][3] =  0.0;

            pnts[i][j].T1_ksi[2][0] =   beta*(phi2 - a*theta);
            pnts[i][j].T1_ksi[2][1] =   beta*(ktx*a - (p_setup.gamma - 1.0)*u);
            pnts[i][j].T1_ksi[2][2] =   beta*(kty*a - (p_setup.gamma - 1.0)*v);
            pnts[i][j].T1_ksi[2][3] =   beta*(p_setup.gamma - 1.0);

            pnts[i][j].T1_ksi[3][0] =   beta*(phi2 + a*theta);
            pnts[i][j].T1_ksi[3][1] = - beta*(ktx*a - (p_setup.gamma - 1.0)*u);
            pnts[i][j].T1_ksi[3][2] = - beta*(kty*a - (p_setup.gamma - 1.0)*v);
            pnts[i][j].T1_ksi[3][3] =   beta*(p_setup.gamma - 1.0);

            /* Now compute the inverse of the split matrices. */

            kx = pnts[i][j].eta_x;
            ky = pnts[i][j].eta_y;

            /* Compute auxiliar variables. */

            beta  = 1.0/(pow(2.0,0.5)*rho*a);
            ktx   = kx/(pow(kx*kx + ky*ky,0.5));
            kty   = ky/(pow(kx*kx + ky*ky,0.5));
            theta = ktx*u + kty*v;
            phi2  = 0.5*(p_setup.gamma-1.0)*(pow(u,2.0) + pow(v,2.0));

            /* Compute the matrix ksi direction.. */

            pnts[i][j].T1_eta[0][0] =  1.0 - phi2/pow(a,2.0);
            pnts[i][j].T1_eta[0][1] =  (p_setup.gamma - 1.0)*u/pow(a,2.0);
            pnts[i][j].T1_eta[0][2] =  (p_setup.gamma - 1.0)*v/pow(a,2.0);
            pnts[i][j].T1_eta[0][3] = -(p_setup.gamma - 1.0)/pow(a,2.0);

            pnts[i][j].T1_eta[1][0] = -(kty*u - ktx*v)/rho;
            pnts[i][j].T1_eta[1][1] =   kty/rho;
            pnts[i][j].T1_eta[1][2] = - ktx/rho;
            pnts[i][j].T1_eta[1][3] =   0.0;

            pnts[i][j].T1_eta[2][0] =   beta*(phi2 - a*theta);
            pnts[i][j].T1_eta[2][1] =   beta*(ktx*a - (p_setup.gamma - 1.0)*u);
            pnts[i][j].T1_eta[2][2] =   beta*(kty*a - (p_setup.gamma - 1.0)*v);
            pnts[i][j].T1_eta[2][3] =   beta*(p_setup.gamma - 1.0);

            pnts[i][j].T1_eta[3][0] =   beta*(phi2 + a*theta);
            pnts[i][j].T1_eta[3][1] = - beta*(ktx*a - (p_setup.gamma - 1.0)*u);
            pnts[i][j].T1_eta[3][2] = - beta*(kty*a - (p_setup.gamma - 1.0)*v);
            pnts[i][j].T1_eta[3][3] =   beta*(p_setup.gamma - 1.0);

        }
    }

    /* Now, do the multiplications and store the splited jacobians. */

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){

            /* Separate useful variables. */

            double a   = pnts[i][j].a;

            double k1 = pnts[i][j].ksi_x;
            double k2 = pnts[i][j].ksi_y;

            double kt1 = k1/(pow(k1*k1 + k2*k2,0.5));
            double kt2 = k2/(pow(k1*k1 + k2*k2,0.5));

            double eig[4], eig_p[4], eig_m[4];

            /* Compute the auxiliar and the needed variables to compute the fluxes. */

            eig[0] = kt1*pnts[i][j].cov_u;
            eig[1] = kt1*pnts[i][j].cov_u;
            eig[2] = kt1*pnts[i][j].cov_u + a*pow(k1*k1 + k2*k2,0.5);
            eig[3] = kt1*pnts[i][j].cov_u - a*pow(k1*k1 + k2*k2,0.5);

            /* Eq. 4.4 of SW original paper. */

            eig_p[0] = (eig[0] + fabs(eig[0])) / 2.0;
            eig_p[1] = (eig[1] + fabs(eig[1])) / 2.0;
            eig_p[2] = (eig[2] + fabs(eig[2])) / 2.0;
            eig_p[3] = (eig[3] + fabs(eig[3])) / 2.0;

            eig_m[0] = (eig[0] - fabs(eig[0])) / 2.0;
            eig_m[1] = (eig[1] - fabs(eig[1])) / 2.0;
            eig_m[2] = (eig[2] - fabs(eig[2])) / 2.0;
            eig_m[3] = (eig[3] - fabs(eig[3])) / 2.0;

            /* Check if the eigenvalues are consistent. */

            for (int ii = 0; ii<4; ii++){

                /* Test the positive eigenvalues. */

                if (fabs(eig[ii] - (eig_p[ii] + eig_m[ii])) > DBL_EPSILON){
                    printf("ERROR: Inconsistent Splited Eigenvalues\n"); 
                    printf("eig[%d] = %lf\n",ii,eig[ii]);
                    printf("eig_p[%d] = %lf\n",ii,eig_p[ii]);
                    printf("eig_m[%d] = %lf\n",ii,eig_m[ii]);
                    exit(1);
                }
            }

            /* Build a matrix with the eigenvalues. */

            double ** meig_p = alloc_dmatrix(4,4);
            double ** meig_m = alloc_dmatrix(4,4);

            /* Initialize. */

            for (int ii = 0; ii<4; ii++)
            for (int jj = 0; jj<4; jj++) meig_p[ii][jj] = 0.0;

            for (int ii = 0; ii<4; ii++)
            for (int jj = 0; jj<4; jj++) meig_m[ii][jj] = 0.0;

            /* Fill the matrices. */

            for (int ii = 0; ii<4; ii++) meig_p[ii][ii] = eig_p[ii];
            for (int ii = 0; ii<4; ii++) meig_m[ii][ii] = eig_m[ii];

            /* Build the first Jacobians in ksi direction. */

            double ** T    = alloc_dmatrix(4,4);
            double ** T1   = alloc_dmatrix(4,4);
            double ** aux1 = alloc_dmatrix(4,4);
            double ** aux2 = alloc_dmatrix(4,4);

            /* Copy the T matrices. */
            
            for (int ii = 0; ii<4; ii++){
                for (int jj = 0; jj<4; jj++){
                    T[ii][jj]  = pnts[i][j].T_ksi[ii][jj];
                    T1[ii][jj] = pnts[i][j].T1_ksi[ii][jj];
                } 
            }

            /* Multiply the first two parts (store in aux).*/

            dmgss(T,meig_p,aux1,4,4,4);

            /* Multiply the second two parts (store in the proper struct. */

            dmgss(aux1,T1,aux2,4,4,4);

            for (int ii = 0; ii<4; ii++)
            for (int jj = 0; jj<4; jj++) pnts[i][j].A_plus[ii][jj] = pnts[i][j].J1 * aux2[ii][jj];

            /* Multiply the first two parts (store in aux).*/

            dmgss(T,meig_m,aux1,4,4,4);

            /* Multiply the second two parts (store in the proper struct. */

            dmgss(aux1,T1,aux2,4,4,4);

            for (int ii = 0; ii<4; ii++)
            for (int jj = 0; jj<4; jj++) pnts[i][j].A_minu[ii][jj] = pnts[i][j].J1 * aux2[ii][jj];

            /* Compute the auxiliar and the needed variables to compute the fluxes. */

            k1 = pnts[i][j].eta_x;
            k2 = pnts[i][j].eta_y;

            kt1 = k1/(pow(k1*k1 + k2*k2,0.5));
            kt2 = k2/(pow(k1*k1 + k2*k2,0.5));

            /* Compute the auxiliar and the needed variables to compute the fluxes. */

            eig[0] = kt2*pnts[i][j].cov_v;
            eig[1] = kt2*pnts[i][j].cov_v;
            eig[2] = kt2*pnts[i][j].cov_v + a*pow(k1*k1 + k2*k2,0.5);
            eig[3] = kt2*pnts[i][j].cov_v - a*pow(k1*k1 + k2*k2,0.5);

            /* Eq. 4.4 of SW original paper. */

            eig_p[0] = (eig[0] + fabs(eig[0])) / 2.0;
            eig_p[1] = (eig[1] + fabs(eig[1])) / 2.0;
            eig_p[2] = (eig[2] + fabs(eig[2])) / 2.0;
            eig_p[3] = (eig[3] + fabs(eig[3])) / 2.0;

            eig_m[0] = (eig[0] - fabs(eig[0])) / 2.0;
            eig_m[1] = (eig[1] - fabs(eig[1])) / 2.0;
            eig_m[2] = (eig[2] - fabs(eig[2])) / 2.0;
            eig_m[3] = (eig[3] - fabs(eig[3])) / 2.0;

            /* Check if the eigenvalues are consistent. */

            for (int ii = 0; ii<4; ii++){

                /* Test the positive eigenvalues. */

                if (fabs(eig[ii] - (eig_p[ii] + eig_m[ii])) > DBL_EPSILON){
                    printf("ERROR: Inconsistent Splited Eigenvalues\n"); 
                    printf("eig[%d] = %lf\n",ii,eig[ii]);
                    printf("eig_p[%d] = %lf\n",ii,eig_p[ii]);
                    printf("eig_m[%d] = %lf\n",ii,eig_m[ii]);
                    exit(1);
                }
            }

            /* Copy the T matrices. */
            
            for (int ii = 0; ii<4; ii++){
                for (int jj = 0; jj<4; jj++){
                    T[ii][jj]  = pnts[i][j].T_eta[ii][jj];
                    T1[ii][jj] = pnts[i][j].T1_eta[ii][jj];
                } 
            }

            /* Multiply the first two parts (store in aux).*/

            dmgss(T,meig_p,aux1,4,4,4);

            /* Multiply the second two parts (store in the proper struct. */

            dmgss(aux1,T1,aux2,4,4,4);

            for (int ii = 0; ii<4; ii++)
            for (int jj = 0; jj<4; jj++) pnts[i][j].B_plus[ii][jj] = pnts[i][j].J1 * aux2[ii][jj];

            /* Multiply the first two parts (store in aux).*/

            dmgss(T,meig_m,aux1,4,4,4);

            /* Multiply the second two parts (store in the proper struct. */

            dmgss(aux1,T1,aux2,4,4,4);

            for (int ii = 0; ii<4; ii++)
            for (int jj = 0; jj<4; jj++) pnts[i][j].B_minu[ii][jj] = pnts[i][j].J1 * aux2[ii][jj];

            /* Check if the jacobians are consistent. */

//            for (int ii = 0; ii<4; ii++){
//                for (int jj = 0; jj<4; jj++){
//
//                    double elem1_A_hat = pnts[i][j].A_hat[ii][jj];
//                    double elem1_A_pos = pnts[i][j].A_plus[ii][jj];
//                    double elem1_A_neg = pnts[i][j].A_minu[ii][jj];
//
//                    double elem1_B_hat = pnts[i][j].B_hat[ii][jj];
//                    double elem1_B_pos = pnts[i][j].B_plus[ii][jj];
//                    double elem1_B_neg = pnts[i][j].B_minu[ii][jj];
//
//                    /* Test the jacobians of the ksi direction flux. */
//
//                    if (fabs(elem1_A_hat - (elem1_A_pos + elem1_A_neg)) > 2.0*DBL_EPSILON){
//                        printf("ERROR: Inconsistent Splited Jacobians\n"); 
//                        printf("pnts[%d][%d].A_hat[%d][%d] = %lf\n",i,j,ii,jj,pnts[i][j].A_hat[ii][jj]); 
//                        printf("pnts[%d][%d].A_plus[%d][%d] = %lf\n",i,j,ii,jj,pnts[i][j].A_plus[ii][jj]); 
//                        printf("pnts[%d][%d].A_minu[%d][%d] = %lf\n",i,j,ii,jj,pnts[i][j].A_minu[ii][jj]); 
//                        exit(1);
//                    }
//
//                    /* Test the jacobians of the eta direction flux. */
//
//                    if (fabs(elem1_B_hat - (elem1_B_pos + elem1_B_neg)) > 2.0*DBL_EPSILON){
//                        printf("ERROR: Inconsistent Splited Jacobians\n"); 
//                        printf("pnts[%d][%d].B_hat[%d][%d] = %lf\n",i,j,ii,jj,pnts[i][j].B_hat[ii][jj]); 
//                        printf("pnts[%d][%d].B_plus[%d][%d] = %lf\n",i,j,ii,jj,pnts[i][j].B_plus[ii][jj]); 
//                        printf("pnts[%d][%d].B_minu[%d][%d] = %lf\n",i,j,ii,jj,pnts[i][j].B_minu[ii][jj]); 
//                        exit(1);
//                    }
//                }
//            }

            /* Free the matrices. */

            free_dmatrix(T, 4);
            free_dmatrix(T1, 4);
            free_dmatrix(aux1, 4);
            free_dmatrix(aux2, 4);
            free_dmatrix(meig_p, 4);
            free_dmatrix(meig_m, 4);
        }
    }
}
