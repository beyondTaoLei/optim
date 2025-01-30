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


#include "fd.h"

void prepare_update_p(float *etajm, float *peta, float **ppi, float **prho, float **ptaup, float **g, float *bjm, float *cjm, float ***e){

	extern int NX, NY, L, MYID;//PARAMETERIZATION, MYID;
	extern float DT, *FL;
	int i, j, l;
	extern char  MFILE[STRING_SIZE];
// 	extern float F_REF;
		
	float sumpi, ws, *pts;
	float pi;
	
	
	/* vector for maxwellbodies */
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
	}
	
	/*ws=2.0*PI*FL[1];*/
// 	ws=2.0*PI*F_REF;
// 	if(MYID==0)printf("MYID %d: F_REF = %f\n",MYID,F_REF);
// 	
// 	sumpi=0.0;
// 	for (l=1;l<=L;l++){
// 		sumpi=sumpi+((ws*ws*pts[l]*pts[l])/(1.0+ws*ws*pts[l]*pts[l]));
// 	}
// 	
// 	for (l=1;l<=L;l++){
// 		etajm[l] = peta[l];
// 	}
// 	
// // 	if (PARAMETERIZATION==1){
// 		for (j=1;j<=NY;j++){
// 		for (i=1;i<=NX;i++){
// 			pi=(ppi[j][i]*ppi[j][i]*prho[j][i])/(1.0+sumpi*ptaup[j][i]);
// 			
// 			g[j][i] = pi*DT*(1.0+L*ptaup[j][i]);
// 			for (l=1;l<=L;l++){
// 				bjm[l] = 1.0/(1.0+(etajm[l]*0.5));
// 				cjm[l] = 1.0-(etajm[l]*0.5);
// 				e[j][i][l] = pi*etajm[l]*ptaup[j][i];
// 			}
// 		}
// 		}
// // 	}
// 	
// 	free_vector(pts,1,L);
}
