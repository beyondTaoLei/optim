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

/* $Id: update_v_PML.c 819 2015-04-17 11:07:06Z tmetz $*/
/*------------------------------------------------------------------------
 *   updating particle velocities at gridpoints of the CPML-frame (ABS=1 in the json file)
 *   by a staggered grid finite difference scheme of FDORDER accuracy in space
 *   and second order accuracy in time
 *   T. Bohlen
 *
 *   gx and gy are arrays with the locations of the boundary specified in subgrid_bounds.c
 *   for each subgrid
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void update_v_PML(int *gx, int *gy, float  **rip, float **rjp, float *hc,
                  float **vx, float **vy, float **uvx, float **uvy, float **sp,
                  float **psi_sxx_x, float **psi_syy_y,
                  float *K_x_half, float *a_x_half, float *b_x_half, 
                  float *K_y_half, float *a_y_half, float *b_y_half,
                  int nx2, int ny2, float flag){
/* flag = 1.0 : forward-propagation in clockwise direction;
   flag =-1.0 : backward-propagation in anticlockwise direction. */

    int i, j, h1,fdoh;
    float sp_x, sp_y;

    extern int FDORDER; 
    extern int  FW;
    extern float DT,DH;

    fdoh=FDORDER/2;

    /* left boundary */
    for ( j=gy[2]+1; j<=gy[3]; j++ ) {
        for ( i=gx[1]; i<=gx[2]; i++ ) {

            SSG_fd_s_ac(&sp_x, &sp_y, sp, i, j, fdoh, hc);

            cpml_update_v_x ( i, j, &sp_x, psi_sxx_x, K_x_half, a_x_half, b_x_half);

            wavefield_update_v(i, j, rip, rjp, sp_x, sp_y, vx, vy, uvx, uvy, DT, DH, flag);

        }
    }

    /* right boundary */
    for ( j=gy[2]+1; j<=gy[3]; j++ ) {
        for ( i=gx[3]+1; i<=gx[4]; i++ ) {

            SSG_fd_s_ac(&sp_x, &sp_y, sp, i, j, fdoh, hc);

            h1 = ( i-nx2+2*FW );
            cpml_update_v_x ( h1, j, &sp_x, psi_sxx_x, K_x_half, a_x_half, b_x_half);
            
            wavefield_update_v(i, j, rip, rjp, sp_x, sp_y, vx, vy, uvx, uvy, DT, DH, flag);

        }
    }

    /* top boundary */
    for ( j=gy[1]; j<=gy[2]; j++ ) {
        for ( i=gx[2]+1; i<=gx[3]; i++ ) {

            SSG_fd_s_ac(&sp_x, &sp_y, sp, i, j, fdoh, hc);
            
            cpml_update_v_y ( i, j, &sp_y, psi_syy_y, K_y_half, a_y_half, b_y_half);
            
            wavefield_update_v(i, j, rip, rjp, sp_x, sp_y, vx, vy, uvx, uvy, DT, DH, flag);
        }
    }

    /* bottom boundary */
    for ( j=gy[3]+1; j<=gy[4]; j++ ) {
        for ( i=gx[2]+1; i<=gx[3]; i++ ) {

            SSG_fd_s_ac(&sp_x, &sp_y, sp, i, j, fdoh, hc);

            h1 = ( j-ny2+2*FW );
            cpml_update_v_y ( i, h1, &sp_y, psi_syy_y, K_y_half, a_y_half, b_y_half);

            wavefield_update_v(i, j, rip, rjp, sp_x, sp_y, vx, vy, uvx, uvy, DT, DH, flag);
        }
    }

    /* corners */

    /*left-top*/
    for ( j=gy[1]; j<=gy[2]; j++ ) {
        for ( i=gx[1]; i<=gx[2]; i++ ) {

            SSG_fd_s_ac(&sp_x, &sp_y, sp, i, j, fdoh, hc);

            cpml_update_v_x ( i, j, &sp_x, psi_sxx_x, K_x_half, a_x_half, b_x_half);
            
            cpml_update_v_y ( i, j, &sp_y, psi_syy_y, K_y_half, a_y_half, b_y_half);

            wavefield_update_v(i, j, rip, rjp, sp_x, sp_y, vx, vy, uvx, uvy, DT, DH, flag);
        }
    }

    /*left-bottom*/
    for ( j=gy[3]+1; j<=gy[4]; j++ ) {
        for ( i=gx[1]; i<=gx[2]; i++ ) {

            SSG_fd_s_ac(&sp_x, &sp_y, sp, i, j, fdoh, hc);

            cpml_update_v_x ( i, j, &sp_x, psi_sxx_x, K_x_half, a_x_half, b_x_half);
            
            h1 = ( j-ny2+2*FW );
            cpml_update_v_y ( i, h1, &sp_y, psi_syy_y, K_y_half, a_y_half, b_y_half);

            wavefield_update_v(i, j, rip, rjp, sp_x, sp_y, vx, vy, uvx, uvy, DT, DH, flag);
        }
    }

    /* right-top */
    for ( j=gy[1]; j<=gy[2]; j++ ) {
        for ( i=gx[3]+1; i<=gx[4]; i++ ) {

            SSG_fd_s_ac(&sp_x, &sp_y, sp, i, j, fdoh, hc);

            h1 = ( i-nx2+2*FW );
            cpml_update_v_x ( h1, j, &sp_x, psi_sxx_x, K_x_half, a_x_half, b_x_half);
            
            cpml_update_v_y ( i, j, &sp_y, psi_syy_y, K_y_half, a_y_half, b_y_half);

            wavefield_update_v(i, j, rip, rjp, sp_x, sp_y, vx, vy, uvx, uvy, DT, DH, flag);
        }
    }

    /* right-bottom */
    for ( j=gy[3]+1; j<=gy[4]; j++ ) {
        for ( i=gx[3]+1; i<=gx[4]; i++ ) {

            SSG_fd_s_ac(&sp_x, &sp_y, sp, i, j, fdoh, hc);

            h1 = ( i-nx2+2*FW );
            cpml_update_v_x ( h1, j, &sp_x, psi_sxx_x, K_x_half, a_x_half, b_x_half);

            h1 = ( j-ny2+2*FW );
            cpml_update_v_y ( i, h1, &sp_y, psi_syy_y, K_y_half, a_y_half, b_y_half);

            wavefield_update_v(i, j, rip, rjp, sp_x, sp_y, vx, vy, uvx, uvy, DT, DH, flag);
        }
    }

}