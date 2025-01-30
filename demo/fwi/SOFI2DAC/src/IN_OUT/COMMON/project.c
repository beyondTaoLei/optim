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

#include "comm.h"

// void create_win(int nx1, int nx2, int ny1, int ny2, int idx, int idy, 
//                 int *winx, int *winy){
    /* create the array to store grids, whose index is from 1 to 
     * 2*icount, as (ix1, iy1, ix2, iy2, ix3, iy3, ...) in order. you 
     * can create other types of windows based on your scheme. 
     */
    
    /* Local varibales */
/*    int i, j, ij, numx, numy, icount;
    
    numx = (nx2-nx1)/idx+1;
    numy = (ny2-ny1)/idy+1;
    icount = numx*numy;

    ij=1;
    for (j=ny1;j<=ny2;j=j+idy){
        for (i=nx1;i<=nx2;i=i+idx){
            winx[ij]=i;
            winy[ij]=j;
            ij = ij+1;
        }
    }
    if(ij-1 != icount){
        printf("ij=%d, icount=%d from function create_win\n", ij,icount);
        exit(1);
    }
}*/

void project2vec(float **m, float *v, int **win, int n1, int n2){
    /* the function is designed for 2D finite diffence method, in this 
     * case, one will project the 2D array [p] into one vector[v].
     */
    
    /* Local varibales */
    int ij;
    
    for (ij=n1;ij<=n2;ij++) v[ij] = m[win[2][ij]][win[1][ij]];
    
}



void project2vec2(float **p, float *v, int nx1, int nx2, int ny1,
                      int ny2, int idx, int idy){
    /* the function is designed for 2D finite diffence method, in this 
     * case, one will project the 2D array [p] into one vector[v].
     */
    
    /* Local varibales */
    int i, j, ij;
    
    ij=0;
    for (j=ny1;j<=ny2;j=j+idy){
        for (i=nx1;i<=nx2;i=i+idx){
            v[ij++] = p[j][i];
        }
    }
}

void project2arr(float **p, float *v, int nx1, int nx2, int ny1,
                      int ny2, int idx, int idy){
    /* the function is designed for 2D finite diffence method, in the 
     * case, one will project one vector[v] into the 2D array [p].
     */
    
    /* Local varibales */
    int i, j, ij;
    //numx=(nx2-nx1)/idx+1;
    //numy=(ny2-ny1)/idy+1;
    //num_of_v=numx*numy;
    
    ij=0;
    for (j=ny1;j<=ny2;j=j+idy){
        for (i=nx1;i<=nx2;i=i+idx){
            p[j][i] = v[ij++];
        }
    }
}

/*

gcc -c project.c  -lm -I ../include/
gcc -o project project.o ../COMMON/util.o
./project

 */
/*
int main(int argc, char **argv){
    
    int i,j,ic;
    int nx1, nx2, ny1, ny2, idx, idy;
    float **m, *v;
    int *win=NULL;
    
    nx1 =1;
    nx2 =10;
    ny1 =1;
    ny2 =20;
    idx =1;
    idy =1;
    m=matrix(ny1, ny2, nx1, nx2);
    for (j=ny1;j<=ny2;j++){
        for (i=nx1;i<=nx2;i++){
            m[j][i]= 1.0*(j*100+i);
    }}
    
    ic=count_elements(nx1, nx2, ny1, ny2, idx, idy);
    
    win= ivector(1, 2*ic);
    v= vector(1, ic);
    create_win(nx1, nx2, ny1, ny2, idx, idy, win);
    
    project2vec(m, v, win, ic);
    for (j=1;j<=ic;j++) printf("i,v=%d %f\n",j,v[j]);
}*/