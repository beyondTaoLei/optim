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
 *   updating stress components at gridpoints [nx1...nx2][ny1...ny2]
 *   by a staggered grid finite difference scheme of arbitrary (FDORDER) order accuracy in space
 *   and second order accuracy in time
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void update_p_interior(int *gx, int *gy, float **rho, float **pi, float *hc, 
                       float **vx, float **vy, float **sp, float **u, float flag){
	
	int i,j, fdoh;
	float  vxx, vyy;
	extern int FDORDER;
        extern float DT,DH;
	
	fdoh = FDORDER/2;
	
        for (j=gy[2]+1; j<=gy[3]; j++){
            for (i=gx[2]+1; i<=gx[3]; i++){
                SSG_fd_v_ac(&vxx, &vyy, vx, vy,i, j, fdoh, hc);
                
                wavefield_update_p(i, j, rho, pi, vxx, vyy, sp, u, DT, DH, flag);
            }
        }
}
/*
void update_p_interior(int *gx, int *gy, float **rho, float **pi, float *hc, 
                       float **vx, float **vy, float **sp, float **u){
	
	int i,j, fdoh;
	float  vxx, vyy;
	extern int FDORDER;
	
	fdoh = FDORDER/2;
	
        for (j=gy[2]+1; j<=gy[3]; j++){
            for (i=gx[2]+1; i<=gx[3]; i++){
                SSG_fd_v_ac(&vxx, &vyy, vx, vy,i, j, fdoh, hc);
                
                wavefield_update_p(i, j, rho, pi, vxx, vyy, sp, u);
            }
        }
}

void update_p_interior_back(int *gx, int *gy, float **rho, float **pi, float *hc, 
                       float **vx, float **vy, float **sp, float **u){
	
	int i,j, fdoh;
	float  vxx, vyy;
	extern int FDORDER;
	
	fdoh = FDORDER/2;
	
        for (j=gy[2]+1; j<=gy[3]; j++){
            for (i=gx[2]+1; i<=gx[3]; i++){
                SSG_fd_v_ac(&vxx, &vyy, vx, vy,i, j, fdoh, hc);
                
                wavefield_update_p_b(i, j, rho, pi, vxx, vyy, sp, u);
            }
        }
}
*/