/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2015  For the list of authors, see file AUTHORS.
 *
 * This file is part of SOFI2D.
 *
 * SOFI2D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 *
 * SOFI2D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SOFI2D See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/

/* $Id: update_s_elastic_PML.c 819 2015-04-17 11:07:06Z tmetz $*/
/*------------------------------------------------------------------------
 *   updating stress components at gridpoints of the CPML-frame (ABS=1 in the json file)
 *   by a staggered grid finite difference scheme of arbitrary (FDORDER) order accuracy in space
 *   and second order accuracy in time
 *   T. Bohlen
 *
 *   gx and gy are arrays with the locations of the boundary specified in subgrid_bounds.c
 *   for each subgrid
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void update_p_PML(int *gx, int *gy, float **rho, float **pi, float *hc,
                  float **vx, float **vy, float **sp, float **u,
                  float **psi_vxx, float **psi_vyy,
                  float *K_x, float *a_x, float *b_x, 
                  float *K_y, float *a_y, float *b_y, 
                  int nx2, int ny2, float flag){
/* flag = 1.0 : forward-propagation in clockwise direction;
   flag =-1.0 : backward-propagation in anticlockwise direction. */
    
    int i,j, h1,fdoh;
    float  vxx, vyy;
    extern int FDORDER, FW;
    extern float DT,DH;

    fdoh=FDORDER/2;

    /* left boundary */
    for ( j=gy[2]+1; j<=gy[3]; j++ ) {
        for ( i=gx[1]; i<=gx[2]; i++ ) {

            SSG_fd_v_ac(&vxx, &vyy, vx, vy, i, j, fdoh, hc);

            cpml_update_p_x ( i, j, &vxx, psi_vxx, K_x, a_x, b_x);

            wavefield_update_p(i, j, rho, pi, vxx, vyy, sp, u, DT, DH, flag);
        }
    }

    /* right boundary */
    for ( j=gy[2]+1; j<=gy[3]; j++ ) {
        for ( i=gx[3]+1; i<=gx[4]; i++ ) {

            SSG_fd_v_ac(&vxx, &vyy, vx, vy, i, j, fdoh, hc);
            
            h1 = ( i-nx2+2*FW );
            cpml_update_p_x ( h1, j, &vxx, psi_vxx, K_x, a_x, b_x);

            wavefield_update_p(i, j, rho, pi, vxx, vyy, sp, u, DT, DH, flag);
        }
    }

    /* top boundary */
    for ( j=gy[1]; j<=gy[2]; j++ ) {
        for ( i=gx[2]+1; i<=gx[3]; i++ ) {

            SSG_fd_v_ac(&vxx, &vyy, vx, vy, i, j, fdoh, hc);

            cpml_update_p_y ( i, j, &vyy, psi_vyy, K_y, a_y, b_y);

            wavefield_update_p(i, j, rho, pi, vxx, vyy, sp, u, DT, DH, flag);
        }
    }

    /* bottom boundary */
    for ( j=gy[3]+1; j<=gy[4]; j++ ) {
        for ( i=gx[2]+1; i<=gx[3]; i++ ) {

            SSG_fd_v_ac(&vxx, &vyy, vx, vy, i, j, fdoh, hc);

            h1 = ( j-ny2+2*FW );
            cpml_update_p_y ( i, h1, &vyy, psi_vyy, K_y, a_y, b_y);

            wavefield_update_p(i, j, rho, pi, vxx, vyy, sp, u, DT, DH, flag);
        }
    }

    /* corners */

    /*left-top*/
    for ( j=gy[1]; j<=gy[2]; j++ ) {
        for ( i=gx[1]; i<=gx[2]; i++ ) {

            SSG_fd_v_ac(&vxx, &vyy, vx, vy, i, j, fdoh, hc);

            cpml_update_p_x ( i, j, &vxx, psi_vxx, K_x, a_x, b_x);

            cpml_update_p_y ( i, j, &vyy, psi_vyy, K_y, a_y, b_y);

            wavefield_update_p(i, j, rho, pi, vxx, vyy, sp, u, DT, DH, flag);
        }
    }

    /*left-bottom*/
    for ( j=gy[3]+1; j<=gy[4]; j++ ) {
        for ( i=gx[1]; i<=gx[2]; i++ ) {

            SSG_fd_v_ac(&vxx, &vyy, vx, vy, i, j, fdoh, hc);

            cpml_update_p_x ( i, j, &vxx, psi_vxx, K_x, a_x, b_x);

            h1 = ( j-ny2+2*FW );
            cpml_update_p_y ( i, h1, &vyy, psi_vyy, K_y, a_y, b_y);

            wavefield_update_p(i, j, rho, pi, vxx, vyy, sp, u, DT, DH, flag);
        }
    }

    /* right-top */
    for ( j=gy[1]; j<=gy[2]; j++ ) {
        for ( i=gx[3]+1; i<=gx[4]; i++ ) {

            SSG_fd_v_ac(&vxx, &vyy, vx, vy, i, j, fdoh, hc);

            h1 = ( i-nx2+2*FW );
            cpml_update_p_x ( h1, j, &vxx, psi_vxx, K_x, a_x, b_x);

            cpml_update_p_y ( i, j, &vyy, psi_vyy, K_y, a_y, b_y);

            wavefield_update_p(i, j, rho, pi, vxx, vyy, sp, u, DT, DH, flag);
        }
    }

    /* right-bottom */
    for ( j=gy[3]+1; j<=gy[4]; j++ ) {
        for ( i=gx[3]+1; i<=gx[4]; i++ ) {

            SSG_fd_v_ac(&vxx, &vyy, vx, vy, i, j, fdoh, hc);

            h1 = ( i-nx2+2*FW );
            cpml_update_p_x ( h1, j, &vxx, psi_vxx, K_x, a_x, b_x);

            h1 = ( j-ny2+2*FW );
            cpml_update_p_y ( i, h1, &vyy, psi_vyy, K_y, a_y, b_y);

            wavefield_update_p(i, j, rho, pi, vxx, vyy, sp, u, DT, DH, flag);
        }
    }


}