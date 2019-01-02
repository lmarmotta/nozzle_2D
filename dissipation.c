#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "structs.h"
#include "externs.h"
#include "prototypes.h"

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

/*
 * Second difference dissipation. 
 */

void art_dissip_2nd(t_define p_setup, t_points ** pnts){

    /* Separate bounds. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Compute the internal points. */

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            /* Ksi direction. */

            pnts[i][j].diss_ksi[0] = (pnts[i+1][j].J*pnts[i+1][j].q_hat[0]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[0]) + (pnts[i-1][j].J*pnts[i-1][j].q_hat[0]);
            pnts[i][j].diss_ksi[1] = (pnts[i+1][j].J*pnts[i+1][j].q_hat[1]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[1]) + (pnts[i-1][j].J*pnts[i-1][j].q_hat[1]);
            pnts[i][j].diss_ksi[2] = (pnts[i+1][j].J*pnts[i+1][j].q_hat[2]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[2]) + (pnts[i-1][j].J*pnts[i-1][j].q_hat[2]);
            pnts[i][j].diss_ksi[3] = (pnts[i+1][j].J*pnts[i+1][j].q_hat[3]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[3]) + (pnts[i-1][j].J*pnts[i-1][j].q_hat[3]);

            /* Eta direction. */

            pnts[i][j].diss_eta[0] = (pnts[i][j+1].J*pnts[i][j+1].q_hat[0]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[0]) + (pnts[i][j-1].J*pnts[i][j-1].q_hat[0]);
            pnts[i][j].diss_eta[1] = (pnts[i][j+1].J*pnts[i][j+1].q_hat[1]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[1]) + (pnts[i][j-1].J*pnts[i][j-1].q_hat[1]);
            pnts[i][j].diss_eta[2] = (pnts[i][j+1].J*pnts[i][j+1].q_hat[2]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[2]) + (pnts[i][j-1].J*pnts[i][j-1].q_hat[2]);
            pnts[i][j].diss_eta[3] = (pnts[i][j+1].J*pnts[i][j+1].q_hat[3]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[3]) + (pnts[i][j-1].J*pnts[i][j-1].q_hat[3]);

        }
    }

    /* Now, add the second difference dissipation to the residue. */

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            pnts[i][j].RHS[0] = pnts[i][j].RHS[0] - pnts[i][j].J1 * ( (p_setup.dissp2)*(pnts[i][j].diss_ksi[0] + pnts[i][j].diss_eta[0]) );
            pnts[i][j].RHS[1] = pnts[i][j].RHS[1] - pnts[i][j].J1 * ( (p_setup.dissp2)*(pnts[i][j].diss_ksi[1] + pnts[i][j].diss_eta[1]) );
            pnts[i][j].RHS[2] = pnts[i][j].RHS[2] - pnts[i][j].J1 * ( (p_setup.dissp2)*(pnts[i][j].diss_ksi[2] + pnts[i][j].diss_eta[2]) );
            pnts[i][j].RHS[3] = pnts[i][j].RHS[3] - pnts[i][j].J1 * ( (p_setup.dissp2)*(pnts[i][j].diss_ksi[3] + pnts[i][j].diss_eta[3]) );

        }
    }
}

/*
 * Fourth derivative dissipation.
 */

void art_dissip_4th(t_define p_setup, t_points ** pnts){

    /* Separate bounds. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* Temporalely takes away the coordinate transformation. */

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){

            pnts[i][j].q_hat[0] = pnts[i][j].J * pnts[i][j].q_hat[0];   
            pnts[i][j].q_hat[1] = pnts[i][j].J * pnts[i][j].q_hat[1];
            pnts[i][j].q_hat[2] = pnts[i][j].J * pnts[i][j].q_hat[2];
            pnts[i][j].q_hat[3] = pnts[i][j].J * pnts[i][j].q_hat[3];

        }
    }

    /* Compute the fourth difference dissipation for the internal points. */

    for (int i = 2; i<imax-2; i++){
        for (int j = 2; j<jmax-2; j++){

            /* Compute in ksi direction. */

            pnts[i][j].diss_ksi[0] = pnts[i+2][j].q_hat[0] - 4.0*pnts[i+1][j].q_hat[0] + 6.0*pnts[i][j].q_hat[0] - 4.0*pnts[i-1][j].q_hat[0] + pnts[i-2][j].q_hat[0];
            pnts[i][j].diss_ksi[1] = pnts[i+2][j].q_hat[1] - 4.0*pnts[i+1][j].q_hat[1] + 6.0*pnts[i][j].q_hat[1] - 4.0*pnts[i-1][j].q_hat[1] + pnts[i-2][j].q_hat[1];
            pnts[i][j].diss_ksi[2] = pnts[i+2][j].q_hat[2] - 4.0*pnts[i+1][j].q_hat[2] + 6.0*pnts[i][j].q_hat[2] - 4.0*pnts[i-1][j].q_hat[2] + pnts[i-2][j].q_hat[2];
            pnts[i][j].diss_ksi[3] = pnts[i+2][j].q_hat[3] - 4.0*pnts[i+1][j].q_hat[3] + 6.0*pnts[i][j].q_hat[3] - 4.0*pnts[i-1][j].q_hat[3] + pnts[i-2][j].q_hat[3];

            /* Compute the eta direction. */

            pnts[i][j].diss_eta[0] = pnts[i][j+2].q_hat[0] - 4.0*pnts[i][j+1].q_hat[0] + 6.0*pnts[i][j].q_hat[0] - 4.0*pnts[i][j-1].q_hat[0] + pnts[i][j-2].q_hat[0];
            pnts[i][j].diss_eta[1] = pnts[i][j+2].q_hat[1] - 4.0*pnts[i][j+1].q_hat[1] + 6.0*pnts[i][j].q_hat[1] - 4.0*pnts[i][j-1].q_hat[1] + pnts[i][j-2].q_hat[1];
            pnts[i][j].diss_eta[2] = pnts[i][j+2].q_hat[2] - 4.0*pnts[i][j+1].q_hat[2] + 6.0*pnts[i][j].q_hat[2] - 4.0*pnts[i][j-1].q_hat[2] + pnts[i][j-2].q_hat[2];
            pnts[i][j].diss_eta[3] = pnts[i][j+2].q_hat[3] - 4.0*pnts[i][j+1].q_hat[3] + 6.0*pnts[i][j].q_hat[3] - 4.0*pnts[i][j-1].q_hat[3] + pnts[i][j-2].q_hat[3];

        }
    }

    /* Compute the symmetry line. */

    for (int i = 2; i<imax-2; i++){

        int j = jmax - 2;

        /* Compute in ksi direction. */

        pnts[i][j].diss_ksi[0] = pnts[i+2][j].q_hat[0] - 4.0*pnts[i+1][j].q_hat[0] + 6.0*pnts[i][j].q_hat[0] - 4.0*pnts[i-1][j].q_hat[0] + pnts[i-2][j].q_hat[0];
        pnts[i][j].diss_ksi[1] = pnts[i+2][j].q_hat[1] - 4.0*pnts[i+1][j].q_hat[1] + 6.0*pnts[i][j].q_hat[1] - 4.0*pnts[i-1][j].q_hat[1] + pnts[i-2][j].q_hat[1];
        pnts[i][j].diss_ksi[2] = pnts[i+2][j].q_hat[2] - 4.0*pnts[i+1][j].q_hat[2] + 6.0*pnts[i][j].q_hat[2] - 4.0*pnts[i-1][j].q_hat[2] + pnts[i-2][j].q_hat[2];
        pnts[i][j].diss_ksi[3] = pnts[i+2][j].q_hat[3] - 4.0*pnts[i+1][j].q_hat[3] + 6.0*pnts[i][j].q_hat[3] - 4.0*pnts[i-1][j].q_hat[3] + pnts[i-2][j].q_hat[3];

        /* Compute the eta direction. */

        pnts[i][j].diss_eta[0] = pnts[i][j].q_hat[0] - 4.0*pnts[i][j-1].q_hat[0] + 6.0*pnts[i][j-2].q_hat[0] - 4.0*pnts[i][j-3].q_hat[0] + pnts[i][j-4].q_hat[0];
        pnts[i][j].diss_eta[1] = pnts[i][j].q_hat[1] - 4.0*pnts[i][j-1].q_hat[1] + 6.0*pnts[i][j-2].q_hat[1] - 4.0*pnts[i][j-3].q_hat[1] + pnts[i][j-4].q_hat[1];
        pnts[i][j].diss_eta[2] = pnts[i][j].q_hat[2] - 4.0*pnts[i][j-1].q_hat[2] + 6.0*pnts[i][j-2].q_hat[2] - 4.0*pnts[i][j-3].q_hat[2] + pnts[i][j-4].q_hat[2];
        pnts[i][j].diss_eta[3] = pnts[i][j].q_hat[3] - 4.0*pnts[i][j-1].q_hat[3] + 6.0*pnts[i][j-2].q_hat[3] - 4.0*pnts[i][j-3].q_hat[3] + pnts[i][j-4].q_hat[3];

    }


    /* Compute the lower surface. */

    for (int i = 2; i<imax-2; i++){

        int j = 1;

        /* Compute in ksi direction. */

        pnts[i][j].diss_ksi[0] = pnts[i+2][j].q_hat[0] - 4.0*pnts[i+1][j].q_hat[0] + 6.0*pnts[i][j].q_hat[0] - 4.0*pnts[i-1][j].q_hat[0] + pnts[i-2][j].q_hat[0];
        pnts[i][j].diss_ksi[1] = pnts[i+2][j].q_hat[1] - 4.0*pnts[i+1][j].q_hat[1] + 6.0*pnts[i][j].q_hat[1] - 4.0*pnts[i-1][j].q_hat[1] + pnts[i-2][j].q_hat[1];
        pnts[i][j].diss_ksi[2] = pnts[i+2][j].q_hat[2] - 4.0*pnts[i+1][j].q_hat[2] + 6.0*pnts[i][j].q_hat[2] - 4.0*pnts[i-1][j].q_hat[2] + pnts[i-2][j].q_hat[2];
        pnts[i][j].diss_ksi[3] = pnts[i+2][j].q_hat[3] - 4.0*pnts[i+1][j].q_hat[3] + 6.0*pnts[i][j].q_hat[3] - 4.0*pnts[i-1][j].q_hat[3] + pnts[i-2][j].q_hat[3];

        /* Compute the eta direction. */

        pnts[i][j].diss_eta[0] = pnts[i][j+4].q_hat[0] - 4.0*pnts[i][j+3].q_hat[0] + 6.0*pnts[i][j+2].q_hat[0] - 4.0*pnts[i][j+1].q_hat[0] + pnts[i][j].q_hat[0];
        pnts[i][j].diss_eta[1] = pnts[i][j+4].q_hat[1] - 4.0*pnts[i][j+3].q_hat[1] + 6.0*pnts[i][j+2].q_hat[1] - 4.0*pnts[i][j+1].q_hat[1] + pnts[i][j].q_hat[1];
        pnts[i][j].diss_eta[2] = pnts[i][j+4].q_hat[2] - 4.0*pnts[i][j+3].q_hat[2] + 6.0*pnts[i][j+2].q_hat[2] - 4.0*pnts[i][j+1].q_hat[2] + pnts[i][j].q_hat[2];
        pnts[i][j].diss_eta[3] = pnts[i][j+4].q_hat[3] - 4.0*pnts[i][j+3].q_hat[3] + 6.0*pnts[i][j+2].q_hat[3] - 4.0*pnts[i][j+1].q_hat[3] + pnts[i][j].q_hat[3];

    }

    /* Compute the inlet surface. */

    for (int j = 2; j<jmax-2; j++){

        int i = 1;

        /* Compute in ksi direction. */

        pnts[i][j].diss_ksi[0] = pnts[i+4][j].q_hat[0] - 4.0*pnts[i+3][j].q_hat[0] + 6.0*pnts[i+2][j].q_hat[0] - 4.0*pnts[i+1][j].q_hat[0] + pnts[i][j].q_hat[0];
        pnts[i][j].diss_ksi[1] = pnts[i+4][j].q_hat[1] - 4.0*pnts[i+3][j].q_hat[1] + 6.0*pnts[i+2][j].q_hat[1] - 4.0*pnts[i+1][j].q_hat[1] + pnts[i][j].q_hat[1];
        pnts[i][j].diss_ksi[2] = pnts[i+4][j].q_hat[2] - 4.0*pnts[i+3][j].q_hat[2] + 6.0*pnts[i+2][j].q_hat[2] - 4.0*pnts[i+1][j].q_hat[2] + pnts[i][j].q_hat[2];
        pnts[i][j].diss_ksi[3] = pnts[i+4][j].q_hat[3] - 4.0*pnts[i+3][j].q_hat[3] + 6.0*pnts[i+2][j].q_hat[3] - 4.0*pnts[i+1][j].q_hat[3] + pnts[i][j].q_hat[3];

        /* Compute the eta direction. */

        pnts[i][j].diss_eta[0] = pnts[i][j+2].q_hat[0] - 4.0*pnts[i][j+1].q_hat[0] + 6.0*pnts[i][j].q_hat[0] - 4.0*pnts[i][j-1].q_hat[0] + pnts[i][j-2].q_hat[0];
        pnts[i][j].diss_eta[1] = pnts[i][j+2].q_hat[1] - 4.0*pnts[i][j+1].q_hat[1] + 6.0*pnts[i][j].q_hat[1] - 4.0*pnts[i][j-1].q_hat[1] + pnts[i][j-2].q_hat[1];
        pnts[i][j].diss_eta[2] = pnts[i][j+2].q_hat[2] - 4.0*pnts[i][j+1].q_hat[2] + 6.0*pnts[i][j].q_hat[2] - 4.0*pnts[i][j-1].q_hat[2] + pnts[i][j-2].q_hat[2];
        pnts[i][j].diss_eta[3] = pnts[i][j+2].q_hat[3] - 4.0*pnts[i][j+1].q_hat[3] + 6.0*pnts[i][j].q_hat[3] - 4.0*pnts[i][j-1].q_hat[3] + pnts[i][j-2].q_hat[3];

    }

    /* Compute the outlet surface. */

    for (int j = 2; j<jmax-2; j++){

        int i = imax - 2;

        /* Compute in ksi direction. */

        pnts[i][j].diss_ksi[0] = pnts[i][j].q_hat[0] - 4.0*pnts[i-1][j].q_hat[0] + 6.0*pnts[i-2][j].q_hat[0] - 4.0*pnts[i-3][j].q_hat[0] + pnts[i-4][j].q_hat[0];
        pnts[i][j].diss_ksi[1] = pnts[i][j].q_hat[1] - 4.0*pnts[i-1][j].q_hat[1] + 6.0*pnts[i-2][j].q_hat[1] - 4.0*pnts[i-3][j].q_hat[1] + pnts[i-4][j].q_hat[1];
        pnts[i][j].diss_ksi[2] = pnts[i][j].q_hat[2] - 4.0*pnts[i-1][j].q_hat[2] + 6.0*pnts[i-2][j].q_hat[2] - 4.0*pnts[i-3][j].q_hat[2] + pnts[i-4][j].q_hat[2];
        pnts[i][j].diss_ksi[3] = pnts[i][j].q_hat[3] - 4.0*pnts[i-1][j].q_hat[3] + 6.0*pnts[i-2][j].q_hat[3] - 4.0*pnts[i-3][j].q_hat[3] + pnts[i-4][j].q_hat[3];

        /* Compute the eta direction. */

        pnts[i][j].diss_eta[0] = pnts[i][j+2].q_hat[0] - 4.0*pnts[i][j+1].q_hat[0] + 6.0*pnts[i][j].q_hat[0] - 4.0*pnts[i][j-1].q_hat[0] + pnts[i][j-2].q_hat[0];
        pnts[i][j].diss_eta[1] = pnts[i][j+2].q_hat[1] - 4.0*pnts[i][j+1].q_hat[1] + 6.0*pnts[i][j].q_hat[1] - 4.0*pnts[i][j-1].q_hat[1] + pnts[i][j-2].q_hat[1];
        pnts[i][j].diss_eta[2] = pnts[i][j+2].q_hat[2] - 4.0*pnts[i][j+1].q_hat[2] + 6.0*pnts[i][j].q_hat[2] - 4.0*pnts[i][j-1].q_hat[2] + pnts[i][j-2].q_hat[2];
        pnts[i][j].diss_eta[3] = pnts[i][j+2].q_hat[3] - 4.0*pnts[i][j+1].q_hat[3] + 6.0*pnts[i][j].q_hat[3] - 4.0*pnts[i][j-1].q_hat[3] + pnts[i][j-2].q_hat[3];

    }

    /* Now, lower left corner. */

    int i = 1; int j = 1;

    pnts[i][j].diss_ksi[0] = pnts[i+4][j].q_hat[0] - 4.0*pnts[i+3][j].q_hat[0] + 6.0*pnts[i+2][j].q_hat[0] - 4.0*pnts[i+1][j].q_hat[0] + pnts[i][j].q_hat[0];
    pnts[i][j].diss_ksi[1] = pnts[i+4][j].q_hat[1] - 4.0*pnts[i+3][j].q_hat[1] + 6.0*pnts[i+2][j].q_hat[1] - 4.0*pnts[i+1][j].q_hat[1] + pnts[i][j].q_hat[1];
    pnts[i][j].diss_ksi[2] = pnts[i+4][j].q_hat[2] - 4.0*pnts[i+3][j].q_hat[2] + 6.0*pnts[i+2][j].q_hat[2] - 4.0*pnts[i+1][j].q_hat[2] + pnts[i][j].q_hat[2];
    pnts[i][j].diss_ksi[3] = pnts[i+4][j].q_hat[3] - 4.0*pnts[i+3][j].q_hat[3] + 6.0*pnts[i+2][j].q_hat[3] - 4.0*pnts[i+1][j].q_hat[3] + pnts[i][j].q_hat[3];

    pnts[i][j].diss_eta[0] = pnts[i][j+4].q_hat[0] - 4.0*pnts[i][j+3].q_hat[0] + 6.0*pnts[i][j+2].q_hat[0] - 4.0*pnts[i][j+1].q_hat[0] + pnts[i][j].q_hat[0];
    pnts[i][j].diss_eta[1] = pnts[i][j+4].q_hat[1] - 4.0*pnts[i][j+3].q_hat[1] + 6.0*pnts[i][j+2].q_hat[1] - 4.0*pnts[i][j+1].q_hat[1] + pnts[i][j].q_hat[1];
    pnts[i][j].diss_eta[2] = pnts[i][j+4].q_hat[2] - 4.0*pnts[i][j+3].q_hat[2] + 6.0*pnts[i][j+2].q_hat[2] - 4.0*pnts[i][j+1].q_hat[2] + pnts[i][j].q_hat[2];
    pnts[i][j].diss_eta[3] = pnts[i][j+4].q_hat[3] - 4.0*pnts[i][j+3].q_hat[3] + 6.0*pnts[i][j+2].q_hat[3] - 4.0*pnts[i][j+1].q_hat[3] + pnts[i][j].q_hat[3];

    /* Upper left corner. */

    i = 1; j = jmax - 2;

    pnts[i][j].diss_ksi[0] = pnts[i+4][j].q_hat[0] - 4.0*pnts[i+3][j].q_hat[0] + 6.0*pnts[i+2][j].q_hat[0] - 4.0*pnts[i+1][j].q_hat[0] + pnts[i][j].q_hat[0];
    pnts[i][j].diss_ksi[1] = pnts[i+4][j].q_hat[1] - 4.0*pnts[i+3][j].q_hat[1] + 6.0*pnts[i+2][j].q_hat[1] - 4.0*pnts[i+1][j].q_hat[1] + pnts[i][j].q_hat[1];
    pnts[i][j].diss_ksi[2] = pnts[i+4][j].q_hat[2] - 4.0*pnts[i+3][j].q_hat[2] + 6.0*pnts[i+2][j].q_hat[2] - 4.0*pnts[i+1][j].q_hat[2] + pnts[i][j].q_hat[2];
    pnts[i][j].diss_ksi[3] = pnts[i+4][j].q_hat[3] - 4.0*pnts[i+3][j].q_hat[3] + 6.0*pnts[i+2][j].q_hat[3] - 4.0*pnts[i+1][j].q_hat[3] + pnts[i][j].q_hat[3];

    pnts[i][j].diss_eta[0] = pnts[i][j].q_hat[0] - 4.0*pnts[i][j-1].q_hat[0] + 6.0*pnts[i][j-2].q_hat[0] - 4.0*pnts[i][j-3].q_hat[0] + pnts[i][j-4].q_hat[0];
    pnts[i][j].diss_eta[1] = pnts[i][j].q_hat[1] - 4.0*pnts[i][j-1].q_hat[1] + 6.0*pnts[i][j-2].q_hat[1] - 4.0*pnts[i][j-3].q_hat[1] + pnts[i][j-4].q_hat[1];
    pnts[i][j].diss_eta[2] = pnts[i][j].q_hat[2] - 4.0*pnts[i][j-1].q_hat[2] + 6.0*pnts[i][j-2].q_hat[2] - 4.0*pnts[i][j-3].q_hat[2] + pnts[i][j-4].q_hat[2];
    pnts[i][j].diss_eta[3] = pnts[i][j].q_hat[3] - 4.0*pnts[i][j-1].q_hat[3] + 6.0*pnts[i][j-2].q_hat[3] - 4.0*pnts[i][j-3].q_hat[3] + pnts[i][j-4].q_hat[3];

    /* Lower right corner. */

    i = imax - 2; j = 1;

    pnts[i][j].diss_ksi[0] = pnts[i][j].q_hat[0] - 4.0*pnts[i-1][j].q_hat[0] + 6.0*pnts[i-2][j].q_hat[0] - 4.0*pnts[i-3][j].q_hat[0] + pnts[i-4][j].q_hat[0];
    pnts[i][j].diss_ksi[1] = pnts[i][j].q_hat[1] - 4.0*pnts[i-1][j].q_hat[1] + 6.0*pnts[i-2][j].q_hat[1] - 4.0*pnts[i-3][j].q_hat[1] + pnts[i-4][j].q_hat[1];
    pnts[i][j].diss_ksi[2] = pnts[i][j].q_hat[2] - 4.0*pnts[i-1][j].q_hat[2] + 6.0*pnts[i-2][j].q_hat[2] - 4.0*pnts[i-3][j].q_hat[2] + pnts[i-4][j].q_hat[2];
    pnts[i][j].diss_ksi[3] = pnts[i][j].q_hat[3] - 4.0*pnts[i-1][j].q_hat[3] + 6.0*pnts[i-2][j].q_hat[3] - 4.0*pnts[i-3][j].q_hat[3] + pnts[i-4][j].q_hat[3];

    pnts[i][j].diss_eta[0] = pnts[i][j+4].q_hat[0] - 4.0*pnts[i][j+3].q_hat[0] + 6.0*pnts[i][j+2].q_hat[0] - 4.0*pnts[i][j+1].q_hat[0] + pnts[i][j].q_hat[0];
    pnts[i][j].diss_eta[1] = pnts[i][j+4].q_hat[1] - 4.0*pnts[i][j+3].q_hat[1] + 6.0*pnts[i][j+2].q_hat[1] - 4.0*pnts[i][j+1].q_hat[1] + pnts[i][j].q_hat[1];
    pnts[i][j].diss_eta[2] = pnts[i][j+4].q_hat[2] - 4.0*pnts[i][j+3].q_hat[2] + 6.0*pnts[i][j+2].q_hat[2] - 4.0*pnts[i][j+1].q_hat[2] + pnts[i][j].q_hat[2];
    pnts[i][j].diss_eta[3] = pnts[i][j+4].q_hat[3] - 4.0*pnts[i][j+3].q_hat[3] + 6.0*pnts[i][j+2].q_hat[3] - 4.0*pnts[i][j+1].q_hat[3] + pnts[i][j].q_hat[3];

    /* Upper right corner. */

    i = imax - 2; j = jmax -2;

    pnts[i][j].diss_ksi[0] = pnts[i][j].q_hat[0] - 4.0*pnts[i-1][j].q_hat[0] + 6.0*pnts[i-2][j].q_hat[0] - 4.0*pnts[i-3][j].q_hat[0] + pnts[i-4][j].q_hat[0];
    pnts[i][j].diss_ksi[1] = pnts[i][j].q_hat[1] - 4.0*pnts[i-1][j].q_hat[1] + 6.0*pnts[i-2][j].q_hat[1] - 4.0*pnts[i-3][j].q_hat[1] + pnts[i-4][j].q_hat[1];
    pnts[i][j].diss_ksi[2] = pnts[i][j].q_hat[2] - 4.0*pnts[i-1][j].q_hat[2] + 6.0*pnts[i-2][j].q_hat[2] - 4.0*pnts[i-3][j].q_hat[2] + pnts[i-4][j].q_hat[2];
    pnts[i][j].diss_ksi[3] = pnts[i][j].q_hat[3] - 4.0*pnts[i-1][j].q_hat[3] + 6.0*pnts[i-2][j].q_hat[3] - 4.0*pnts[i-3][j].q_hat[3] + pnts[i-4][j].q_hat[3];

    pnts[i][j].diss_eta[0] = pnts[i][j].q_hat[0] - 4.0*pnts[i][j-1].q_hat[0] + 6.0*pnts[i][j-2].q_hat[0] - 4.0*pnts[i][j-3].q_hat[0] + pnts[i][j-4].q_hat[0];
    pnts[i][j].diss_eta[1] = pnts[i][j].q_hat[1] - 4.0*pnts[i][j-1].q_hat[1] + 6.0*pnts[i][j-2].q_hat[1] - 4.0*pnts[i][j-3].q_hat[1] + pnts[i][j-4].q_hat[1];
    pnts[i][j].diss_eta[2] = pnts[i][j].q_hat[2] - 4.0*pnts[i][j-1].q_hat[2] + 6.0*pnts[i][j-2].q_hat[2] - 4.0*pnts[i][j-3].q_hat[2] + pnts[i][j-4].q_hat[2];
    pnts[i][j].diss_eta[3] = pnts[i][j].q_hat[3] - 4.0*pnts[i][j-1].q_hat[3] + 6.0*pnts[i][j-2].q_hat[3] - 4.0*pnts[i][j-3].q_hat[3] + pnts[i][j-4].q_hat[3];

    /* Makes it coordinate transformed. */

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){

            pnts[i][j].q_hat[0] = pnts[i][j].J1 * pnts[i][j].q_hat[0];  
            pnts[i][j].q_hat[1] = pnts[i][j].J1 * pnts[i][j].q_hat[1];
            pnts[i][j].q_hat[2] = pnts[i][j].J1 * pnts[i][j].q_hat[2];
            pnts[i][j].q_hat[3] = pnts[i][j].J1 * pnts[i][j].q_hat[3];

        }
    }

    /* Now, add to the residue. */

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            pnts[i][j].RHS[0] = pnts[i][j].RHS[0] + pnts[i][j].J1 * p_setup.dissp2 * (pnts[i][j].diss_ksi[0] + pnts[i][j].diss_eta[0]);
            pnts[i][j].RHS[1] = pnts[i][j].RHS[1] + pnts[i][j].J1 * p_setup.dissp2 * (pnts[i][j].diss_ksi[1] + pnts[i][j].diss_eta[1]);
            pnts[i][j].RHS[2] = pnts[i][j].RHS[2] + pnts[i][j].J1 * p_setup.dissp2 * (pnts[i][j].diss_ksi[2] + pnts[i][j].diss_eta[2]);
            pnts[i][j].RHS[3] = pnts[i][j].RHS[3] + pnts[i][j].J1 * p_setup.dissp2 * (pnts[i][j].diss_ksi[3] + pnts[i][j].diss_eta[3]);

        }
    }
}

/*
 * Tau operator.
 */

double tau(double pf, double p, double pb){

    double tau_out = 0.0;

    tau_out = fabs(pf - 2.0*p + pb)/fabs(pf + 2.0*p + pb);

    return tau_out;
}

/*
 * Non linear dissipation eps2
 */

double e2(double tauf, double tau, double taub, double dt){

    double e2_out = 0.0;

    double k2 = 1.0/4.0;

    e2_out = k2*dt*max(tauf,max(tau,taub));

    return e2_out;
}

/*
 * Non linear dissipation eps4
 */

double e4(double e2, double dt){

    double e4_out = 0.0;

    double k4 = 1.0/100.0;

    e4_out = max(0.0, k4*dt - e2);

    return e4_out;
}

/*
 * Dissipation main operator.
 * q: is a vector with the forward and backward operators of q_hat.
 */

double nlim_op(double sigf, double sig, double jacf, double jac, double e2, double e4, double * q){

    double nlim_op_out = 0.0;

    /*
     * q[0] = q(i-2)
     * q[1] = q(i-1)
     * q[2] = q(i)
     * q[3] = q(i+1)
     * q[4] = q(i+2)
     */

    nlim_op_out = sigf*jacf + sig*jac;

    nlim_op_out = nlim_op_out*( e2*(q[3] - q[2]) - (e4*(q[4] - 3.0*q[3] + 3.0*q[2] - q[1])) );

    return nlim_op_out;
}

/*
 * Build the dissipation.
 */

void art_dissip_nli(t_define p_setup, t_points ** pnts){

    /* Separate bounds. */

    int imax = p_setup.imax;
    int jmax = p_setup.jmax;

    /* I will not compute the shifted operrators of the dissipation, so the
     * borders will be second derivatives. */

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            /* Ksi direction. */

            pnts[i][j].diss_ksi[0] = (pnts[i+1][j].J*pnts[i+1][j].q_hat[0]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[0]) + (pnts[i-1][j].J*pnts[i-1][j].q_hat[0]);
            pnts[i][j].diss_ksi[1] = (pnts[i+1][j].J*pnts[i+1][j].q_hat[1]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[1]) + (pnts[i-1][j].J*pnts[i-1][j].q_hat[1]);
            pnts[i][j].diss_ksi[2] = (pnts[i+1][j].J*pnts[i+1][j].q_hat[2]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[2]) + (pnts[i-1][j].J*pnts[i-1][j].q_hat[2]);
            pnts[i][j].diss_ksi[3] = (pnts[i+1][j].J*pnts[i+1][j].q_hat[3]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[3]) + (pnts[i-1][j].J*pnts[i-1][j].q_hat[3]);

            /* Eta direction. */

            pnts[i][j].diss_eta[0] = (pnts[i][j+1].J*pnts[i][j+1].q_hat[0]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[0]) + (pnts[i][j-1].J*pnts[i][j-1].q_hat[0]);
            pnts[i][j].diss_eta[1] = (pnts[i][j+1].J*pnts[i][j+1].q_hat[1]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[1]) + (pnts[i][j-1].J*pnts[i][j-1].q_hat[1]);
            pnts[i][j].diss_eta[2] = (pnts[i][j+1].J*pnts[i][j+1].q_hat[2]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[2]) + (pnts[i][j-1].J*pnts[i][j-1].q_hat[2]);
            pnts[i][j].diss_eta[3] = (pnts[i][j+1].J*pnts[i][j+1].q_hat[3]) - (2.0 * pnts[i][j].J*pnts[i][j].q_hat[3]) + (pnts[i][j-1].J*pnts[i][j-1].q_hat[3]);

        }
    }

    /* Compute the sigmas. */

    double ** sig = alloc_dmatrix(imax,jmax);

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){ 

            double a = pnts[i][j].a;

            double modU   = fabs(pnts[i][j].cov_u);
            double modV   = fabs(pnts[i][j].cov_v);

            double ksi_x2 = pnts[i][j].ksi_x*pnts[i][j].ksi_x;
            double ksi_y2 = pnts[i][j].ksi_y*pnts[i][j].ksi_y;

            double eta_x2 = pnts[i][j].eta_x*pnts[i][j].eta_x;
            double eta_y2 = pnts[i][j].eta_y*pnts[i][j].eta_y;

            sig[i][j] = modU + a*pow(ksi_x2 + ksi_y2, 0.5) + modV + a*pow(eta_x2 + eta_y2, 0.5);
        }
    }

    /* Compute the taus. */

    double ** tau_s = alloc_dmatrix(imax,jmax);

    /* Since the dissipation is operated, we have to store. */

    double *** nlim = alloc_dcube(imax, jmax, 4);

    /* Start storing the tau. */

    for (int i = 1; i<imax-1; i++)
    for (int j = 1; j<jmax-1; j++) tau_s[i][j] = tau(pnts[i+1][j].p, pnts[i][j].p, pnts[i-1][j].p); 

    /* Finalize the operator build in ksi. */

    for (int i = 2; i<imax-2; i++){
        for (int j = 2; j<jmax-2; j++){ 

            double sigf  = sig[i+1][j];
            double sig_v = sig[i][j];
            double jacf = pnts[i+1][j].J1;
            double jac  = pnts[i][j].J1;
            double e2_v   = e2(tau_s[i+1][j], tau_s[i][j], tau_s[i-1][j], pnts[i][j].dt);
            double e4_v   = e4(e2_v, pnts[i][j].dt); 

            /*
             * q[0] = q(i-2)
             * q[1] = q(i-1)
             * q[2] = q(i)
             * q[3] = q(i+1)
             * q[4] = q(i+2)
             */

            double * q = (double*)malloc(5*sizeof(double));

            for (int n = 0; n<4; n++){

                q[0] = pnts[i][j-2].J*pnts[i-2][j].q_hat[n];
                q[1] = pnts[i][j-1].J*pnts[i-1][j].q_hat[n];

                q[2] = pnts[i][j].J*pnts[i][j].q_hat[n];

                q[3] = pnts[i][j+1].J*pnts[i+1][j].q_hat[n];
                q[4] = pnts[i][j+2].J*pnts[i+2][j].q_hat[n];

                nlim[i][j][n] = nlim_op(sigf, sig_v, jacf, jac, e2_v, e4_v, q);
            }

            free(q);
        }
    }

    /* Finalize the operator build in ksi. */

    for (int i = 2; i<imax-2; i++)
        for (int j = 2; j<jmax-2; j++)
            for (int n = 0; n<4; n++) pnts[i][j].diss_ksi[n] = nlim[i][j][n] - nlim[i-1][j][n];

    /* Finalize the operator build in eta. */

    for (int i = 2; i<imax-2; i++){
        for (int j = 2; j<jmax-2; j++){ 

            double sigf = sig[i][j+1];
            double sig_v = sig[i][j];
            double jacf = pnts[i][j+1].J1;
            double jac  = pnts[i][j].J1;
            double e2_v = e2(tau_s[i][j+1], tau_s[i][j], tau_s[i][j-1], pnts[i][j].dt);
            double e4_v = e4(e2_v, pnts[i][j].dt); 

            /*
             * q[0] = q(i-2)
             * q[1] = q(i-1)
             * q[2] = q(i)
             * q[3] = q(i+1)
             * q[4] = q(i+2)
             */

            double * q = (double*)malloc(5*sizeof(double));

            for (int n = 0; n<4; n++){

                q[0] = pnts[i][j-2].J*pnts[i][j-2].q_hat[n];
                q[1] = pnts[i][j-1].J*pnts[i][j-1].q_hat[n];

                q[2] = pnts[i][j].J*pnts[i][j].q_hat[n];

                q[3] = pnts[i][j+1].J*pnts[i][j+1].q_hat[n];
                q[4] = pnts[i][j+2].J*pnts[i][j+2].q_hat[n];

                nlim[i][j][n] = nlim_op(sigf, sig_v, jacf, jac, e2_v, e4_v, q);
            }

            free(q);
        }
    }

    /* Finalize the operator build in eta. */

    for (int i = 2; i<imax-2; i++)
        for (int j = 2; j<jmax-2; j++)
            for (int n = 0; n<4; n++) pnts[i][j].diss_eta[n] = nlim[i][j][n] - nlim[i][j-1][n];

    /* Now, add to the residue. */

    for (int i = 1; i<imax-1; i++){
        for (int j = 1; j<jmax-1; j++){

            pnts[i][j].RHS[0] = pnts[i][j].RHS[0] - pnts[i][j].J1 * p_setup.dissp2 * (pnts[i][j].diss_ksi[0] + pnts[i][j].diss_eta[0]);
            pnts[i][j].RHS[1] = pnts[i][j].RHS[1] - pnts[i][j].J1 * p_setup.dissp2 * (pnts[i][j].diss_ksi[1] + pnts[i][j].diss_eta[1]);
            pnts[i][j].RHS[2] = pnts[i][j].RHS[2] - pnts[i][j].J1 * p_setup.dissp2 * (pnts[i][j].diss_ksi[2] + pnts[i][j].diss_eta[2]);
            pnts[i][j].RHS[3] = pnts[i][j].RHS[3] - pnts[i][j].J1 * p_setup.dissp2 * (pnts[i][j].diss_ksi[3] + pnts[i][j].diss_eta[3]);

        }
    }

    free_dmatrix(sig, imax);
    free_dmatrix(tau_s, imax);
    free_dcube(nlim, imax, jmax);
}
