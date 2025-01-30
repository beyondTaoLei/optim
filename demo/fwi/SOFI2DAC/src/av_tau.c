/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2016  For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS.
 * 
 * IFOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * IFOS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IFOS. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/

#include "fd.h"

void av_tau(float **taus, float **tausipjp){
    
    extern int NX, NY;
    int i, j;
    for (j=1;j<=NY;j++){
        for (i=1;i<=NX;i++){
            
            tausipjp[j][i] = 0.25*(taus[j][i]+taus[j][i+1]+taus[j+1][i]+taus[j+1][i+1]);
        }
    }
}
