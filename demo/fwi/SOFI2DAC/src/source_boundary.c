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

/*------------------------------------------------------------------------
 *   generate P-wave source at source nodes
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void source_boundary(float **sp, int **srcpos_loc, float **signals, int nt, int nsrc){
    /* local variables */
    int i, j, l;
    /* adding source wavelet to stress components 
    (explosive source) at source points */
    for (l=1;l<=nsrc;l++){
        i=(int)srcpos_loc[1][l];
        j=(int)srcpos_loc[2][l];
        sp[j][i]=signals[l][nt];
    }
    
}

void source_last_snap(int ishot, float **psp, int flag){
    
    extern int NX, NY, POS[3];
    
    /* local variables */
    int i,j;
    char snapfile[STRING_SIZE], ext[8];
    float val;
    FILE *fp=NULL;
    
    sprintf(ext,".bin");
    if(flag==1){
        sprintf(snapfile,"snapl%s.x.shot%i.%i.%i",ext,ishot,POS[1],POS[2]);
    }else if(flag==2){
        sprintf(snapfile,"snapl%s.y.shot%i.%i.%i",ext,ishot,POS[1],POS[2]);
    }else if(flag==3){
        sprintf(snapfile,"snapl%s.p.shot%i.%i.%i",ext,ishot,POS[1],POS[2]);
    }
    fp = fopen(snapfile, "r");
    for (i=1;i<=NX;i++){
    for (j=1;j<=NY;j++){
        if(feof(fp)){
            declare_error("last snap is to small. Check dimensions NX*NY of file.");
        }
        fread(&val, sizeof(float), 1, fp);
        psp[j][i]=val;
    }}
    fclose(fp);
}
