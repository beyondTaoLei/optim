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

void extract_snap_from_su(float **signal_sp, float *snap_sp, int nsrc_loc, int it){
    /* extracts current signals(one snapshot) from local source 
     * position(a time series) in current node */
    
    /* local variables */
    int i;
    for(i=1; i<=nsrc_loc; i++){
        snap_sp[i] = signal_sp[i][it];
        //printf("%d %f\n",it,snap_sp[i]);
    }
}

//FILE *open_wfd_file(char srcname[256], char flag[]){
    /* If the sources are added as snapshots, which are useful for some 
     * algorithms, such as truncated Newton, boundary saving strategy. 
     * One should read the signals from BIN-format file.
     */
/*    FILE *fp1;
    if(strcmp(flag,"r")==0){
        fp1 = fopen(srcname, "r");
    }else if(strcmp(flag,"w")==0){
        fp1 = fopen(srcname, "w");
    }
    
    return fp1;
    
}*/

//void close_wfd_file(FILE *fp1){
    /* If the sources are added as snapshots, which are useful for some 
     * algorithms, such as truncated Newton, boundary saving strategy. 
     * One should read the signals from BIN-format file.
     */
/*    
    fclose(fp1);
}*/

void read_wfd_file(FILE *fp1, float *wfd, int num, int hin){
    /* we also need to read wavefields from the disk. */
    
    fseek(fp1, hin*num*sizeof(float),SEEK_SET);
    fread(wfd,sizeof(float),num,fp1);
   
}

void write_wfd_file(FILE *fp1, float *wfd, int num, int hin){
    /* in order to save snapshots from the different time in one file, 
     * before that one should open the file, and after that one also 
     * should close the corresponding file.
     */
    
    fseek(fp1, hin*num*sizeof(float),SEEK_SET);
    fwrite(wfd,sizeof(float),num,fp1);
    
}