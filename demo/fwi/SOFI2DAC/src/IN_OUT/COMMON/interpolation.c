/*-------------------------------------------------------------------
 * Copyright (C) 2018  For the list of authors, see file AUTHORS.
 *
 * This file is part of SEISPLATFORM.
 * 
 * SEISPLATFORM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published 
 * by the Free Software Foundation, version 2.0 of the License only.
 * 
 * SEISPLATFORM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SEISPLATFORM. See file COPYING and/or 
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
-------------------------------------------------------------------*/

/* files to include */
#include "comm.h"

float linint(float x0, float y0, float x1, float y1, float x){
    /* linear interpolation */
    
    /* local varibales */
    float y;
    
    if(x0==x1){
        printf("error occurs in util_fun: linint\n");
        exit(1);
    }else{
        y = y0+(x-x0)*(y1-y0)/(x1-x0);
    }
    return y;
}

float linint_set(float *x, float *y, float xx, int n){
    /* linear interpolation on a data set */
    
    /* local varibales */
    int i, j;
    float yy=0.0;
    
    if(xx<x[0]){
        j=0;
        yy = linint(x[j], y[j], x[j+1], y[j+1], xx);
    }else if(xx>=x[n-1]){
        j=n-2;
        yy = linint(x[j], y[j], x[j+1], y[j+1], xx);
    }else{
        for(i=0;i<n;i++){
            if(xx>=x[i]&&xx<x[i+1]){
                j=i;
                yy = linint(x[j], y[j], x[j+1], y[j+1], xx);
                break;
            }
        }
    }
    
    return yy;
    
    
}
