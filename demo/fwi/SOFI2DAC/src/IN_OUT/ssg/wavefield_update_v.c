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

/* $Id: wavefield_update_v.c 819 2015-04-17 11:07:06Z tmetz $ */

/*Update Function of the particle velocity Wavefields */

#include "ssg.h"
/*
void wavefield_update_v(int i, int j, float **rip, float **rjp, 
                        float sp_x, float sp_y, float **vx, float **vy, 
                        float **uvx, float **uvy){
    float  dtdh;
    extern float DT,DH;
//  dtdh = DT/DH;
    uvx[j][i] =  sp_x *rip[j][i]/DH;
    uvy[j][i] =  sp_y *rjp[j][i]/DH;

    vx[j][i] +=  uvx[j][i]*DT;
    vy[j][i] +=  uvy[j][i]*DT;

}

void wavefield_update_v_b(int i, int j, float **rip, float **rjp, 
                        float sp_x, float sp_y, float **vx, float **vy, 
                        float **uvx, float **uvy){
    float  dtdh;
    extern float DT,DH;
//  dtdh = DT/DH;
    uvx[j][i] =  sp_x *rip[j][i]/DH;
    uvy[j][i] =  sp_y *rjp[j][i]/DH;

    vx[j][i] -=  uvx[j][i]*DT;
    vy[j][i] -=  uvy[j][i]*DT;

}
*/
void wavefield_update_v(int i, int j, float **rip, float **rjp, 
                        float sp_x, float sp_y, float **vx, float **vy, 
                        float **uvx, float **uvy, float dt, float dh, float flag){
/* flag = 1.0 : forward-propagation in clockwise direction;
   flag =-1.0 : backward-propagation in anticlockwise direction. */
    //float  dtdh;
    //extern float dt,dh;
//  dtdh = dt/dh;
    uvx[j][i] =  sp_x *rip[j][i]/dh;
    uvy[j][i] =  sp_y *rjp[j][i]/dh;

    vx[j][i] +=  flag*uvx[j][i]*dt;
    vy[j][i] +=  flag*uvy[j][i]*dt;

}