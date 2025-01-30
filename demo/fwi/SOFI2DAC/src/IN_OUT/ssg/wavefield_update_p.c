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

// void wavefield_update_p(int i, int j, float **rho, float **pi,
//                         float vxx, float vyy,  float **sp, float **u){
// 
// 	extern float DT,DH;
// 	float dtdh, g;
// //         dtdh = DT/DH;
//         
// 	g = pi[j][i];
// 	/*Update  */
//         u[j][i]   = g*(vxx+vyy)/DH;
//         sp[j][i] += u[j][i]*DT;
// 
// }
// 
// void wavefield_update_p_b(int i, int j, float **rho, float **pi,
//                         float vxx, float vyy,  float **sp, float **u){
// 
// 	extern float DT,DH;
// 	float dtdh, g;
// //         dtdh = DT/DH;
//         
// 	g = pi[j][i];
// 	/*Update  */
//         u[j][i]   = g*(vxx+vyy)/DH;
//         sp[j][i] -= u[j][i]*DT;
// 
// }

void wavefield_update_p(int i, int j, float **rho, float **pi,
                        float vxx, float vyy,  float **sp, float **u,
                        float dt, float dh, float flag){
    
/* flag = 1.0 : forward-propagation in clockwise direction;
   flag =-1.0 : backward-propagation in anticlockwise direction. 
   rho: density
   pi: bulk modulus kappa=vp*vp*rho
   */

	//extern float dt,dh;
	//float dtdh, g;
//         dtdh = dt/dh;
        
	//g = pi[j][i];
	/*Update  */
        u[j][i]   = pi[j][i]*(vxx+vyy)/dh;
        sp[j][i] += flag*u[j][i]*dt;

}
