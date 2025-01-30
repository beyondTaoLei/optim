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

int **splitpos(int iendx, int iendy, int *pos, int **recpos_g,int *ntr_loc, int ntr_glob){
    
    //extern int IENDX, IENDY, POS[4], REC1, REC2;
    int REC1, REC2;
    int a,b,i=0,j,k,found=0;
    int ** recpos_dummy, **recpos_local=NULL;
    recpos_dummy = imatrix(1,3,1,ntr_glob);
    
    REC1=0;
    REC2=0;
    
    for (j=1;j<=ntr_glob;j++) {
        
        a=(recpos_g[1][j]-1)/iendx;
        b=(recpos_g[2][j]-1)/iendy;
        
        if ((pos[1]==a)&&(pos[2]==b)) {
            
            found = 1;
            if(REC1==0){
                REC1=j;
            }
            
            i++;	/* Anzahl der Empfänger i jedes Prozesses ermitteln */
            recpos_dummy[1][i] = ((recpos_g[1][j]-1)%iendx)+1;
            recpos_dummy[2][i] = ((recpos_g[2][j]-1)%iendy)+1;
            recpos_dummy[3][i] = j;
            
        }
    }
    
    if (i>0) recpos_local = imatrix(1,3,1,i);
    REC2=REC1+(i-1);
    for (k=1;k<=i;k++){
        recpos_local[1][k] = recpos_dummy[1][k];
        recpos_local[2][k] = recpos_dummy[2][k];
        recpos_local[3][k] = recpos_dummy[3][k];
    }
    free_imatrix(recpos_dummy,1,3,1,ntr_glob);
    *ntr_loc=i;
    return recpos_local;
    
}
