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

void read_su_pos(char srcname[256], float dh, float *refrec, int nt_in,
                 float dt_in, int ***srcpos, int *nsrc_glob){        
    /*so far this software just supports the SU-format file.*/
    /* INPUT */
    /* 
     * char srcname[256]
     * float dh 
     * float *refrec
     * int nt_in
     * float dt_in
     */
    
    /* OUTPUT */
    /* 
     * int ***srcpos
     * int *nsrc
     */
    
    /* local variables */
    int nt, nsrc_tmp2;
    float dt;
    float **srcpos_tmp3d=NULL;
    int *repeated_tmp=NULL;
    
    srcpos_tmp3d = read_su_header(srcname, nsrc_glob, &nt, &dt);
    *srcpos=coord_grid(dh, refrec, srcpos_tmp3d, &repeated_tmp, *nsrc_glob, &nsrc_tmp2);
    if(*nsrc_glob!=nsrc_tmp2||nt!=nt_in||dt!=dt_in){
        printf("nsrc_su, nsrc_su2= %d %d\n",*nsrc_glob,nsrc_tmp2);
        printf("nt_in_su_file, nt_in_time_loop= %d %d\n",nt,nt_in);
        printf("dt_in_su_file, dt_in_time_loop= %f %f\n",dt,dt_in);
        printf("Please make sure the source positions are non-redundant, nt_in and dt_in are correct.\n");
        exit(0);
    }
    free_matrix(srcpos_tmp3d,1,3,1,*nsrc_glob);
    if(*nsrc_glob>nsrc_tmp2)  free_ivector(repeated_tmp,1,*nsrc_glob-nsrc_tmp2);
    
}

int **evaluate_wfd_coords(float dh, float *refrec, int fdorder, int free_surf,
                          int fw, int nxg, int nyg, float memory_depth,
                          int source_snap_type, int *nsrc_glob){
    /* If the sources are added as snapshots, which are useful for some 
     * algorithms, such as truncated Newton, boundary saving strategy. 
     * One should read the signals from BIN-format file.
     */
    
    /* local variables */
    int i,j,itr,num;
    int ix1,ix2,ix3,ix4;
    int iy1,iy2,iy3,iy4;
    int width = fdorder/2;
    int **grids=NULL;
    
    if(source_snap_type==1){
        /* signals from boundary in the case of boundary saving strategy */
        if(free_surf){
            iy1=2;      iy2=iy1+width-1;
            iy4=nyg-fw; iy3=iy4-width+1;
        }else{
            iy1=fw+1;       iy2=iy1+width-1;
            iy4=nyg-fw;     iy3=iy4-width+1;
        }
        ix1=fw+1;       ix2=ix1+width-1;
        ix4=nxg-fw;     ix3=ix4-width+1;
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
            printf("There are some mistakes in function evaluate_wfd_coords !\n");
            exit(0);
        }
    }else if(source_snap_type==2){
        /* signals from the specified depth */
        ix1=1;          ix2=nxg;
        iy1=(int)(memory_depth/dh); iy2=iy1 + width-1;
        num=nxg*width;
        grids=imatrix(1,3,1,num);
        itr=1;
        for(j=iy1;j<=iy2;j++){
        for(i=ix1;i<=ix2;i++){
            grids[1][itr]=i;
            grids[2][itr]=j;
            itr++;
        }}
        if(itr-1!=num){
            printf("There are some mistakes in function evaluate_wfd_coords !\n");
            exit(0);
        }
    }else if(source_snap_type==3){
        printf("There are some mistakes in function evaluate_wfd_coords !\n");
        exit(0);
       
    }else if(source_snap_type==4){
        /* signals from the whole model space */
        num=nxg*nyg;
        grids=imatrix(1,3,1,num);
        itr=1;
        for(i=1;i<=nxg;i++){
        for(j=1;j<=nyg;j++){
               grids[1][itr]=i;
               grids[2][itr]=j; 
               itr++;
        }}
    }
    
    //*srcpos_loc = splitpos(iendx, iendy, POS1, grids, nsrc_loc, num);
    //free_imatrix(grids,1,3,1,num);
    
    *nsrc_glob = num;
    
    return grids;

}