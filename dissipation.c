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
 * Non linear dissipation eps2
 */

