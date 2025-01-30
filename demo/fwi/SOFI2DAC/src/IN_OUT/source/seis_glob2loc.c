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

#include <stdio.h>

void seis_glob2loc(int ntr,int ns, float **fulldate, float **section, int ** recpos_loc){
    
    int i, j, h=1;
    for(i=1;i<=ntr;i++){
        for(j=1;j<=ns;j++){
            section[h][j]=fulldate[recpos_loc[3][i]][j];
        }
        h++;
    }
}

void seis_glob2loc_vec(int ntr, float *fulldate, float *section, int ** recpos_loc){
    
    int i, h=1;
    for(i=1;i<=ntr;i++){
        section[h]=fulldate[recpos_loc[3][i]];
        h++;
    }
}