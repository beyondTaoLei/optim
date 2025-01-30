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

#include "fd.h"

int **receiver_boundary(FILE *fp, int *ntr){
    /* the sources are added for boundary saving strategy. 
     * One should read the signals from BIN-format file.
     * origin from file source/evaluate_wfd_coords.c
     */
    /* declaration of extern variables */
    extern int NXG, NYG, FREE_SURF, FW, FDORDER;
    
    /* local variables */
    int i,j,itr,num;
    int ix1,ix2,ix3,ix4;
    int iy1,iy2,iy3,iy4;
    int width = FDORDER/2;
    int **grids=NULL;
    
    /* signals from boundary in the case of boundary saving strategy */
    if(FREE_SURF){
        iy1=2;      iy2=iy1+width-1;
        iy4=NYG-FW; iy3=iy4-width+1;
    }else{
        iy1=FW+1;       iy2=iy1+width-1;
        iy4=NYG-FW;     iy3=iy4-width+1;
    }
    ix1=FW+1;       ix2=ix1+width-1;
    ix4=NXG-FW;     ix3=ix4-width+1;
    num=((ix4-ix1+1)+(iy3-iy2+1-2))*width*2;
    
    grids=imatrix(1,3,1,num);
    itr=1;
    for(i=ix1;i<=ix4;i++){
    for(j=iy1;j<=iy2;j++){
        grids[1][itr]=i;
        grids[2][itr]=j;
        itr++;
    }}
    for(i=ix1;i<=ix4;i++){
    for(j=iy3;j<=iy4;j++){
        grids[1][itr]=i;
        grids[2][itr]=j;
        itr++;
    }}
    for(i=ix1;i<=ix2;i++){
    for(j=iy2+1;j<=iy3-1;j++){
        grids[1][itr]=i;
        grids[2][itr]=j;
        itr++;
    }}
    for(i=ix3;i<=ix4;i++){
    for(j=iy2+1;j<=iy3-1;j++){
        grids[1][itr]=i;
        grids[2][itr]=j;
        itr++;
    }}
    if(itr-1!=num){
        fprintf(fp, "There are some mistakes in function evaluate_wfd_coords !\n");
        exit(0);
    }
    *ntr = num;
    
    return grids;

}