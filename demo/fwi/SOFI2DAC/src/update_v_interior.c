/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2016  For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS.
 * 
 * IFOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * IFOS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IFOS. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------
 *   updating particle velocities at gridpoints [nx1...nx2][ny1...ny2]
 *   by a staggered grid finite difference scheme of FDORDER accuracy in space
 *   and second order accuracy in time
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void update_v_interior(int *gx, int *gy, float  **rip, float **rjp, float *hc,
                       float **vx, float **vy, float **uvx, float **uvy, float **sp, float flag){

    int i, j,fdoh;
    float sp_x, sp_y;
    extern int FDORDER;
    extern float DT,DH;

    fdoh = FDORDER/2;

    /* ------------------------------------------------------------
    * Important!
    * rip and rjp are reciprocal values of averaged densities
    * ------------------------------------------------------------ */
    for (j=gy[2]+1; j<=gy[3]; j++){
    for (i=gx[2]+1; i<=gx[3]; i++){
        SSG_fd_s_ac(&sp_x, &sp_y, sp, i, j, fdoh, hc);
        
        wavefield_update_v(i, j, rip, rjp, sp_x, sp_y, vx, vy, uvx, uvy, DT, DH, flag);
    }
    }
}

// void update_v_interior(int *gx, int *gy, float  **rip, float **rjp, float *hc,
//                        float **vx, float **vy, float **uvx, float **uvy, float **sp){
// 
//     int i, j,fdoh;
//     float sp_x, sp_y;
//     extern int FDORDER;
// 
//     fdoh = FDORDER/2;
// 
//     /* ------------------------------------------------------------
//     * Important!
//     * rip and rjp are reciprocal values of averaged densities
//     * ------------------------------------------------------------ */
//     for (j=gy[2]+1; j<=gy[3]; j++){
//     for (i=gx[2]+1; i<=gx[3]; i++){
//         SSG_fd_s_ac(&sp_x, &sp_y, sp, i, j, fdoh, hc);
//         
//         wavefield_update_v(i, j, rip, rjp, sp_x, sp_y, vx, vy, uvx, uvy);
//     }
//     }
// }
// 
// void update_v_interior_back(int *gx, int *gy, float  **rip, float **rjp, float *hc,
//                        float **vx, float **vy, float **uvx, float **uvy, float **sp){
// 
//     int i, j,fdoh;
//     float sp_x, sp_y;
//     extern int FDORDER;
// 
//     fdoh = FDORDER/2;
// 
//     /* ------------------------------------------------------------
//     * Important!
//     * rip and rjp are reciprocal values of averaged densities
//     * ------------------------------------------------------------ */
//     for (j=gy[2]+1; j<=gy[3]; j++){
//     for (i=gx[2]+1; i<=gx[3]; i++){
//         SSG_fd_s_ac(&sp_x, &sp_y, sp, i, j, fdoh, hc);
//         
//         wavefield_update_v_b(i, j, rip, rjp, sp_x, sp_y, vx, vy, uvx, uvy);
//     }
//     }
// }
