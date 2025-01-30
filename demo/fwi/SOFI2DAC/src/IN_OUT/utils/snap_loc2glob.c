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
/**
@brief rename the snapshots for sources and receivers
@details Because we name the wavefield snapshots with number, it's not 
convenient to assemble the kernels. Here we will rename the snapshots
with the global index, which means the source and receiver will use the 
same coordinate system. In future, we will name the source and snapshots
files by strings.

@file snap_loc2glob.c
@author Tao Lei (leit@sustech.edu.cn), Wei Zhang
@date 26 Sep 2019
@param[in]  verbose \b 1: print the details; \b 0: run in quiet pattern, default value;
@param[in]  fincf   the folder of snapshots from source sides <em>(optional & default: snap/inc)</em>;
@param[in]  fadjf   the folder of snapshots from receiver sides <em>(optional & default: snap/adj)</em>;
@param[in]  fsnap   the prefix of snapshots file name <em>(optional & default: test)</em>;
@param[in]  fspos  the relation between local index and global index 
            for sources <em>(optional & default: list/S_POS.dat)</em>;
@param[in]  frposl  the relation between local index and global index 
            for receivers, \b NOTING: the receivers exclude the ones in the 
            sources position <em>(optional & default: list/R_left_POS.dat)</em>;
@param[in,out] foutf the folder of output files <em>(optional & default: snap)</em>;
@code 
snap_loc2glob fincf=snap/inc fadjf=snap/adj fsnap=test
@endcode
*/
#include "comm.h"
#include "par.h"
#include<dirent.h>

/* Global variables */
int xargc; 
char **xargv;

int main(int argc, char *argv[]){

    /* IN and OUT */
    int verbose;
    char *fincf, *fadjf, *fsnap, *fspos, *frposl, *foutf;

    /* Local variables */
    int ns, ntr, ishot, ishot_g;
    float **pos_s, **pos_r;
    char oldname[256], newname[256], oldname2[256], newname2[256];
    DIR* dp;
    struct dirent* dirp;

    /* initialize getpar */
    initargs(argc,argv);
    
    /* 1. read formatted input list */
    verbose=0; getparint("verbose",&verbose);
    if(!getparstring("fincf", &fincf)) fincf="snap/inc";
    if(!getparstring("fadjf", &fadjf)) fadjf="snap/adj";
    if(!getparstring("fsnap", &fsnap)) fsnap="test";
    if(!getparstring("fspos", &fspos)) fspos="list/S_POS.dat";
    if(!getparstring("frposl", &frposl)) frposl="list/R_left_POS.dat";
    if(!getparstring("foutf", &foutf)) foutf="snap";
    
    if(verbose){
        printf("\n******\n");
        printf("%-18s| %s\n","fincf", fincf);
        printf("%-18s| %s\n","fadjf",fadjf);
        printf("%-18s| %s\n","fsnap",fsnap);
        printf("%-18s| %s\n","fspos",fspos);
        printf("%-18s| %s\n","frposl",frposl);
        printf("%-18s| %s\n","foutf",foutf);
        printf("******\n");
    }
    
    pos_s=read_text(2, 2, &ns,fspos);
    pos_r=read_text(2, 2, &ntr,frposl);
    /* source parts */
    if((dp=opendir(fincf))==NULL){
        printf("\ncannot open %s\n",fincf);
        exit(1);
    }
    while((dirp=readdir(dp))!=NULL){
        if((int)(strlen(dirp->d_name))>2){
            sprintf(oldname, "%s",dirp->d_name);
            ishot =extract_shot_num(oldname);
            ishot_g=(int)pos_s[ishot][2];
            update_shot_num2(oldname, ishot_g, newname);
            sprintf(oldname2,"%s/%s",fincf,oldname);
            sprintf(newname2,"%s/%s",foutf,newname);
            if(rename(oldname2, newname2) == 0){
                if(verbose) printf("%s--->%s\n",oldname2, newname2);
            }else{
                throw_error("Please the filename!\n");
            }

        }
        
    }
    closedir(dp);
    
    /* receiver parts */
    if((dp=opendir(fadjf))==NULL){
        printf("\ncannot open %s\n",fadjf);
        exit(1);
    }
    while((dirp=readdir(dp))!=NULL){
        if((int)(strlen(dirp->d_name))>2){
            sprintf(oldname, "%s",dirp->d_name);
            ishot =extract_shot_num(oldname);
            ishot_g=(int)pos_r[ishot][2];
            update_shot_num2(oldname, ishot_g, newname);
            sprintf(oldname2,"%s/%s",fadjf,oldname);
            sprintf(newname2,"%s/%s",foutf,newname);
            if(rename(oldname2, newname2) == 0){
                if(verbose) printf("%s--->%s\n",oldname2, newname2);
            }else{
                throw_error("Please the filename!\n");
            }
        }
    }
    closedir(dp);
    
    free_matrix(pos_s, 1, ns, 1, 2);
    free_matrix(pos_r, 1, ntr, 1, 2);
    
    return 0;
}
