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

/*----------------------------------------------------------------
 *   Calculate L2 adjoint source                                  
 *---------------------------------------------------------------*/

#include "adj.h"

float misfit_L2(float *obs, float *syn, float *res, float tstart, 
                float tend, int ns, float dt){

    /* declaration of variables */
    int j;
    float misfit=0.0;

    /* res will be set to zero */
    zero(&res[1],ns);
    
    /*************************************/
    /*         calculate residuals       */
    /*************************************/
    for(j=1;j<=ns;j++){
        /* calculate L2 residuals */
        res[j] = syn[j]-obs[j];
        misfit +=res[j]*res[j];
    }
    
    return misfit;
    
} /* end of function */
