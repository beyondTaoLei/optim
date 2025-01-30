/*---------------------------------------------------------------------------------
 * Copyright (C) 2015 For the list of authors, see file AUTHORS.
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
 * along with SOFI2D. See file COPYING and/or
  * <http://www.gnu.org/licenses/gpl-2.0.html>.
---------------------------------------------------------------------------------*/

/* $Id: wavefield_update_s_el.c 819 2015-04-17 11:07:06Z tmetz $ */

/*Update Function of the stress-Wavefields in the elastic case*/

#include "ssg.h"

void wavefield_update_s_el ( int i, int j,float   vxx, float  vyx,float vxy,float  vyy, float **sxy,
                   float **sxx, float ** syy, float ** pi, float ** u, float ** uipjp )
{

	extern float DT;
	float fipjp, f, g;


	fipjp=uipjp[j][i]*DT;
	f = u[j][i]*DT;
	g = pi[j][i]*DT;

	/*Update  */
	sxy[j][i] += fipjp* ( vyx+vxy );
	sxx[j][i] += ( g* ( vxx+vyy ) )- ( 2.0*f*vyy );
	syy[j][i] += ( g* ( vxx+vyy ) )- ( 2.0*f*vxx );

}

void wavefield_update_s_el_new(int i, int j, float **pi, float **u, float **uipjp,
                               float vxx, float vyx, float vxy, float vyy, 
                               float **sxy, float **sxx, float **syy, 
                               float **uxy, float **ux, float ** uy,
                               float dt, float dh, float flag){
/* flag = 1.0 : forward-propagation in clockwise direction;
   flag =-1.0 : backward-propagation in anticlockwise direction. */

        /*Update  */
        uxy[j][i] = uipjp[j][i]*(vyx+vxy)/dh;
        ux[j][i] =  (pi[j][i]*(vxx+vyy)-2.0*u[j][i]*vyy)/dh;
        uy[j][i] =  (pi[j][i]*(vxx+vyy)-2.0*u[j][i]*vxx)/dh;
        
        sxy[j][i] += flag*uxy[j][i]*dt;
        sxx[j][i] += flag*ux[j][i]*dt;
        syy[j][i] += flag*uy[j][i]*dt;

}

// void wavefield_update_s_el_new2(int i, int j, float pi, float u, float uipjp,
//                                float vxx, float vyx, float vxy, float vyy, 
//                                float *sxy, float *sxx, float *syy, 
//                                float *uxy, float *ux, float * uy,
//                                float dt, float dh, float flag){
// /* flag = 1.0 : forward-propagation in clockwise direction;
//    flag =-1.0 : backward-propagation in anticlockwise direction. */
// 
//         /*Update  */
//         *uxy = uipjp*(vyx+vxy)/dh;
//         *ux =  (pi*(vxx+vyy)-2.0*u*vyy)/dh;
//         *uy =  (pi*(vxx+vyy)-2.0*u*vxx)/dh;
//         
//         *sxy += flag*(*uxy)*dt;
//         *sxx += flag*(*ux)*dt;
//         *syy += flag*(*uy)*dt;
// 
// }
