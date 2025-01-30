/*------------------------------------------------------------------------
 * Copyright (C) 2011 For the list of authors, see file AUTHORS.
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
--------------------------------------------------------------------------*/
/*------------------------------------------------------------------------
 *   spatial direvatives based on standard staggered grid(SSG) in 
 *   Cartesian coordinates as suggested by Virieux[1986] and Levander[1988]
 *
 *  ----------------------------------------------------------------------*/

#include "ssg.h"

/* regardless of efficiency, the following parts are useful when 
 * updating the wavefields around the boundary. */
void SSG_fd_v(float *vxx, float *vxy, float *vyy, float *vyx, 
           float **vx, float **vy, int i, int j, int fdoh, float *hc){
    /* when updating the stress components in elastic equation, one 
     * should calculate the particle velocity direvatives */
    
    *vxx = calc_vxx(vx, i, j, fdoh, hc);
    *vxy = calc_vxy(vx, i, j, fdoh, hc);
    *vyy = calc_vyy(vy, i, j, fdoh, hc);
    *vyx = calc_vyx(vy, i, j, fdoh, hc);
}

void SSG_fd_s(float *sxx_x, float *sxy_x, float *sxy_y, float *syy_y, 
           float ** sxx, float ** syy, float ** sxy, int i, int j, 
           int fdoh, float *hc){
    /* when updating the velocity components in elastic equation, one 
     * should calculate the stress direvatives */
    
    *sxx_x = calc_sxx_x(sxx, i, j, fdoh, hc);
    *sxy_x = calc_sxy_x(sxy, i, j, fdoh, hc);
    *sxy_y = calc_sxy_y(sxy, i, j, fdoh, hc);
    *syy_y = calc_syy_y(syy, i, j, fdoh, hc);
}

void SSG_fd_s_ac(float *sp_x, float *sp_y, float **sp, int i, int j, 
               int fdoh, float *hc){
    /* when updating the velocity components in acoustic equation, one 
     * should calculate the pressure direvatives */

    *sp_x = calc_sp_x(sp, i, j, fdoh, hc);
    *sp_y = calc_sp_y(sp, i, j, fdoh, hc);
}

void SSG_fd_v_ac(float *vxx, float *vyy, float **vx, float **vy, 
                 int i, int j, int fdoh, float *hc){
    /* when updating the stress components in elastic equation, one 
     * should calculate the particle velocity direvatives */
    
    *vxx = calc_vxx(vx, i, j, fdoh, hc);
    *vyy = calc_vyy(vy, i, j, fdoh, hc);
}

float calc_vxx(float **vx, int i, int j, int fdoh, float *hc){
    
    int m;
    float  sum1;
    
    sum1 = 0.0;
    for (m=1; m<=fdoh; m++) {
        sum1 += hc[m]*(vx[j][i+m-1] -vx[j][i-m]);
    }
        
    //sum1 *= dhi;
    return sum1;
    
}

float calc_vyy(float **vy, int i, int j, int fdoh, float *hc){
    
    int m;
    float  sum1;
    
    sum1 = 0.0;
    for (m=1; m<=fdoh; m++) {
        sum1 += hc[m]*(vy[j+m-1][i] -vy[j-m][i]);
    }
    //sum1 *= dhi;
    return sum1;
    
}

float calc_vxy(float **vx, int i, int j, int fdoh, float *hc){
    
    int m;
    float  sum1;
    
    sum1 = 0.0;
    for (m=1; m<=fdoh; m++) {
        sum1 += hc[m]*(vx[j+m][i] -vx[j-m+1][i]);
    }
        
    //sum1 *= dhi;
    return sum1;
    
}

float calc_vyx(float **vy, int i, int j, int fdoh, float *hc){
    
    int m;
    float  sum1;
    
    sum1 = 0.0;
    for (m=1; m<=fdoh; m++) {
        sum1 += hc[m]*(vy[j][i+m] -vy[j][i-m+1]);
    }
        
    //sum1 *= dhi;
    return sum1;
    
}

float calc_sxx_x(float **sxx, int i, int j, int fdoh, float *hc){
    
    int m;
    float  sum1;
    
    sum1 = 0.0;
    for (m=1; m<=fdoh; m++) {
        sum1 += hc[m]*(sxx[j][i+m] -sxx[j][i-m+1]);
    }
        
    //sum1 *= dhi;
    return sum1;
    
}

float calc_sxy_x(float **sxy, int i, int j, int fdoh, float *hc){
    
    int m;
    float  sum1;
    
    sum1 = 0.0;
    for (m=1; m<=fdoh; m++) {
        sum1 += hc[m]*(sxy[j][i+m-1] -sxy[j][i-m]);
    }
        
    //sum1 *= dhi;
    return sum1;
    
}

float calc_sxy_y(float **sxy, int i, int j, int fdoh, float *hc){
    
    int m;
    float  sum1;
    
    sum1 = 0.0;
    for (m=1; m<=fdoh; m++) {
        sum1 += hc[m]*(sxy[j+m-1][i] -sxy[j-m][i]);
    }
        
    //sum1 *= dhi;
    return sum1;
    
}

float calc_syy_y(float **syy, int i, int j, int fdoh, float *hc){
    
    int m;
    float  sum1;
    
    sum1 = 0.0;
    for (m=1; m<=fdoh; m++) {
        sum1 += hc[m]*(syy[j+m][i] -syy[j-m+1][i]);
    }
        
    //sum1 *= dhi;
    return sum1;
    
}

float calc_sp_x(float **sp, int i, int j, int fdoh, float *hc){
    
    int m;
    float  sum1;
    
    sum1 = 0.0;
    for (m=1; m<=fdoh; m++) {
        sum1 += hc[m]*(sp[j][i+m] -sp[j][i-m+1]);
    }
        
    //sum1 *= dhi;
    return sum1;
    
}

float calc_sp_y(float **sp, int i, int j, int fdoh, float *hc){
    
    int m;
    float  sum1;
    
    sum1 = 0.0;
    for (m=1; m<=fdoh; m++) {
        sum1 += hc[m]*(sp[j+m][i] -sp[j-m+1][i]);
    }
        
    //sum1 *= dhi;
    return sum1;
    
}