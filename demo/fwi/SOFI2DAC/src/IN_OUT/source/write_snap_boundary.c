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

void write_snap_boundary(FILE *fp1, float **sp, int **pos_loc, int num, int hin){
    /* saves the snapshots in the specified coordinates */

    int i,j,l;
    float *snap = NULL;
    
    //sprintf(snapfile_p,"%s.bin.%s.%i.%i",snap_file,ext,POS[1],POS[2]);
    
    snap = vector(1,num);
    
    for(l=1;l<=num;l++){
        i=(int)pos_loc[1][l];
        j=(int)pos_loc[2][l];
        snap[l] = sp[j][i];
    }
    write_wfd_file(fp1, &snap[1], num, hin);
    
    free_vector(snap, 1, num);

}


