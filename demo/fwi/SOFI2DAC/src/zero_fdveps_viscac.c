 /*-----------------------------------------------------------------------------------------
 * Copyright (C) 2013  For the list of authors, see file AUTHORS.
 *
 * This file is part of DENISE.
 * 
 * DENISE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * DENISE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with DENISE. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   zero wavefield
 *  
 *  
 *   last update 06/10/11, L. Groos
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void zero_fdveps_viscac(int ny1, int ny2, int nx1, int nx2, float ** vx, float ** vy, float ** sp, float ** vxp1, float ** vyp1,
			float ** psi_sxx_x, float ** psi_sxy_x, float ** psi_vxx, float ** psi_vyx, float ** psi_syy_y,
			float ** psi_sxy_y, float ** psi_vyy, float ** psi_vxy, float ** psi_vxxs, float ***pp){


	register int i, j, l;
	extern int FW, NX, NY, L;
	
	for (j=ny1;j<=ny2;j++){
		for (i=nx1;i<=nx2;i++){
			
			vx[j][i]=0.0;
			vy[j][i]=0.0;
			sp[j][i]=0.0;
			vxp1[j][i]=0.0;
			vyp1[j][i]=0.0;
			
		}
	}
	
	for (j=1;j<=NY;j++){
		for (i=1;i<=2*FW;i++){
			
			psi_sxx_x[j][i] = 0.0;
			psi_sxy_x[j][i] = 0.0;
			psi_vxx[j][i] = 0.0;
			psi_vxxs[j][i] = 0.0;
			psi_vyx[j][i] = 0.0;   
			
		}
	}
	
	for (j=1;j<=2*FW;j++){
		for (i=1;i<=NX;i++){
			
			psi_syy_y[j][i] = 0.0;
			psi_sxy_y[j][i] = 0.0;
			psi_vyy[j][i] = 0.0;
			psi_vxy[j][i] = 0.0;
			
		}
	}
	
	for (j=ny1;j<=ny2;j++){
		for (i=nx1;i<=nx2;i++){
			for (l=1;l<=L;l++){
				pp[j][i][l] = 0.0;
			}
		}
	}    
}
