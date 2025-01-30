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

/* $Id: update_s_visc_PML.c,v 1.1.1.1 2011/10/06 22:44:52 groos Exp $*/
/*------------------------------------------------------------------------
 *   updating stress components at gridpoints [nx1...nx2][ny1...ny2]
 *   by a staggered grid finite difference scheme of arbitrary (FDORDER) order accuracy in space
 *   and second order accuracy in time
 *   T. Bohlen
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void update_p_visc_PML(int nx1, int nx2, int ny1, int ny2, float ** vx, float ** vy, float ** sp, float ** pi, float **rho, float *hc, int infoout,
		       float ***p, float **g, float *bjm, float *cjm, float ***e,  
		       float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
		       float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
		       float ** psi_vxx, float ** psi_vyy, float ** psi_vxy, float ** psi_vyx){
	
	int i,j, m, fdoh, h, h1, l;
	float  vxx, vyy, vxy, vyx;
	float  dhi, dthalbe;	
	extern float DT, DH;
	extern int MYID, FDORDER, FW, L;
        extern int FREE_SURF, BOUNDARY;
	extern int NPROCX, NPROCY, POS[3];
	extern FILE *FP;
	double time1 = 0.0, time2;
	
	float sump=0.0;
	
	/*dhi = DT/DH;*/
	dhi=1.0/DH;
	fdoh = FDORDER/2;
	dthalbe = DT/2.0;
	
	if (infoout && (MYID==0)){
		time1=MPI_Wtime();
		fprintf(FP,"\n **Message from update_p_visc (printed by PE %d):\n",MYID);
		fprintf(FP," Updating stress components ...");
	}
	
	switch (FDORDER){
	
	case 2:
		for (j=ny1;j<=ny2;j++){
		for (i=nx1;i<=nx2;i++){
			vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1]))*dhi;
			
			vyx = (  hc[1]*(vy[j][i+1]-vy[j][i]))*dhi;
			
			vxy = (  hc[1]*(vx[j+1][i]-vx[j][i]))*dhi;
			
			vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i]))*dhi; 
			
			/* left boundary */                                         
			if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
				
				psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * vxx;
				vxx = vxx / K_x[i] + psi_vxx[j][i];
				
				psi_vyx[j][i] = b_x_half[i] * psi_vyx[j][i] + a_x_half[i] * vyx;
				vyx = vyx / K_x_half[i] + psi_vyx[j][i];                 
			}
			
			/* right boundary */                                         
			if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
				
				h1 = (i-nx2+2*FW);
				h = i;
				
				psi_vxx[j][h1] = b_x[h1] * psi_vxx[j][h1] + a_x[h1] * vxx;
				vxx = vxx / K_x[h1] + psi_vxx[j][h1];
				
				/*psi_vyx[j][h] = b_x_half[h] * psi_vyx[j][h] + a_x_half[h] * vyx;
				vyx = vyx / K_x_half[h] + psi_vyx[j][h];*/
				
				psi_vyx[j][h1] = b_x_half[h1] * psi_vyx[j][h1] + a_x_half[h1] * vyx;
				vyx = vyx / K_x_half[h1] + psi_vyx[j][h1];
			}
			
			/* top boundary */                                         
			if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
				
				psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * vyy;                                            
				psi_vxy[j][i] = b_y_half[j] * psi_vxy[j][i] + a_y_half[j] * vxy;
			
				vyy = vyy / K_y[j] + psi_vyy[j][i];
				vxy = vxy / K_y_half[j] + psi_vxy[j][i];
			}
			
			/* bottom boundary */                                         
			if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
				
				h1 = (j-ny2+2*FW);                                        
				h = j;
				
				psi_vyy[h1][i] = b_y[h1] * psi_vyy[h1][i] + a_y[h1] * vyy;                                            
				vyy = vyy / K_y[h1] + psi_vyy[h1][i];
				
				/*psi_vxy[j][i] = b_y_half[j] * psi_vxy[j][i] + a_y_half[j] * vxy;
				vxy = vxy / K_y_half[j] + psi_vxy[j][i];*/
				
				psi_vxy[h1][i] = b_y_half[h1] * psi_vxy[h1][i] + a_y_half[h1] * vxy;
				vxy = vxy / K_y_half[h1] + psi_vxy[h1][i];
			}
			
			/* computing sums of the old memory variables */
			sump=0.0;
			for (l=1;l<=L;l++){
				sump+=p[j][i][l];
			}
			
			/* updating components of the stress tensor, partially */
			sp[j][i] += (g[j][i]*(vxx+vyy))+(dthalbe*sump);
			
			/* now updating the memory-variables and sum them up*/
			sump=0.0;
			for (l=1;l<=L;l++){
				p[j][i][l] = bjm[l]*(p[j][i][l]*cjm[l]-(e[j][i][l]*(vxx+vyy)));
				sump += p[j][i][l];
			}
			
			
			/* and now the components of the stress tensor are
			   completely updated */
			sp[j][i]+=(dthalbe*sump);
		}
		}
	break;

	default:
		for (j=ny1;j<=ny2;j++){
		for (i=nx1;i<=nx2;i++){
			vxx = 0.0;
			vyy = 0.0;
			vyx = 0.0;
			vxy = 0.0;
			for (m=1; m<=fdoh; m++) {
				vxx += hc[m]*(vx[j][i+m-1] -vx[j][i-m]  );
				vyy += hc[m]*(vy[j+m-1][i] -vy[j-m][i]  );
				vyx += hc[m]*(vy[j][i+m]   -vy[j][i-m+1]);
				vxy += hc[m]*(vx[j+m][i]   -vx[j-m+1][i]);
			}
			vxx *= dhi;
			vyy *= dhi;
			vyx *= dhi;
			vxy *= dhi;
			
			sump=0.0;
			for (l=1;l<=L;l++){
				sump+=p[j][i][l];
			}
			
			sp[j][i] += (g[j][i]*(vxx+vyy))+(dthalbe*sump);

			sump=0.0;
			for (l=1;l<=L;l++){
				p[j][i][l] = bjm[l]*(p[j][i][l]*cjm[l]-(e[j][i][l]*(vxx+vyy)));
				sump += p[j][i][l];
			}
			
			sp[j][i]+=(dthalbe*sump);
		}
		}
	break;
	
	} /* end of switch(FDORDER) */
	
	if (infoout && (MYID==0)){
		time2=MPI_Wtime();
		fprintf(FP," finished (real time: %4.2f s).\n",time2-time1);
	}
}
