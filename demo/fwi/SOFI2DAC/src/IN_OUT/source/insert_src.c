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

#include "source.h"
/* this function can be used in different cases, which is decided by flag.
 * Au=f
 * 
 * 1: as the normal source for forward modeling (FD) operator
 * 2: as the special source for FD operator
 * 3: as the boundary in receiver side for RTM or in the outer boundary 
 *    position for wavefield reconstruction used backward modeling (BD) operator 
 * //delete 4: as the adjoint source in receiver side for FWI used the backward operator
 * 0: no source in some cases
 * 
 */

float insert_src(float compnt, float sig, float b,int flag){
    
    /* local variables */
    float tmp;

    if(flag==1)
        /* adding source signals to any one component at source points */
        tmp = compnt+b*sig;
    else if(flag==2)
        /* replacing source signals to any one component at source points */
        tmp = b*sig;
    else if(flag==3)
        tmp = sig;
    else if(flag==4)
        /* adding source signals to any one component at source points */
        tmp = compnt-b*sig;
    else if(flag==0)
        /* do nothing */
        tmp = compnt;
    else{
        printf("I don't know which type of source you want to add.\n");
        exit(1);
    }
    
    return tmp;
}