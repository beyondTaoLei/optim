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
@brief Write the various wavelets into the file with the SU-format.
@file write_wavelet.c
@author Tao Lei (leit@sustech.edu.cn), Wei Zhang
@date 18 Sep 2019
@param[in] nt the time samples;
@param[in] delta the time interval;
@param[in] fsrc the source file in which you can set the corresponding 
parameters, <em>(optional & default: source.dat)</em>, each row is ordered as <em>[ x 0.0 y tshift fc amp ]</em>
@param[in] fpath the output path, <em>(optional & default: .)</em>
@param[in] shape the source shape, <em>(optional & default: 1)</em>
<em>

    ricker=1; fumue=2; from_SIGNAL_FILE=3; SIN**3=4; Gaussian_deriv=5; Spike=6; from_SIGNAL_FILE_in_su_format=7
</em>
@param[in] type  the source type, <em>(optional & default: 1)</em>
<em>

    explosive=1; force_in_x=2; force_in_y=3; rotated_force=4
</em>
@param[in] flag whether to save sources in the different files or not, <em>(optional & default: 1)</em>
<em>

    1: in several different files; 0: in one file.
</em>

@code 
write_wavelet nt=3000 delta=4.0e-04 fsrc=source.dat fpath=. shape=1 type=3 flag=1
@endcode
*/

#include "par.h"
#include "source.h"

/* Global variables */
int xargc; 
char **xargv;

int main(int argc, char **argv){
    
    int nshots, nt;
    int shape, type, flag;
    float  delta;
    
    /*local variables */
    int i,j,it,ishot;
    float **srcpos_tmp=NULL,**srcpos=NULL,**Sxy=NULL;
    float **signals_all=NULL, **signals=NULL;
    char file_out[256];
    
    char *fsrc, *fpath;
    
    /* 1. read formatted input list */
    /* initialize getpar */
    initargs(argc,argv);
    
    if(!getparstring("fsrc", &fsrc))  fsrc= "source.dat";
    if(!getparstring("fpath", &fpath)) fpath= ".";
    if(!getparint("nt",&nt)) throw_error("nt must be specified!\n");
    if(!getparfloat("delta",&delta)) throw_error("delta must be specified!\n");
    shape=1; getparint("shape",&shape);
    type=1; getparint("type",&type);
    flag=1; getparint("flag",&flag);
    printf("fsrc=%s\n",fsrc);
    printf("fpath=%s\n",fpath);
    printf("nt, delta=%d, %f\n",nt,delta);
    printf("shape, type, flag=%d, %d %d\n",shape,type,flag);

    srcpos_tmp = sources(fsrc, &nshots);
    srcpos =  matrix(1,3,1,nshots);
    Sxy =  matrix(1,3,1,1);
    signals=matrix(1,1,1,nt);
    for(i=1;i<=nshots;i++){
        for(j=1;j<=3;j++){
            srcpos[j][i]=srcpos_tmp[j][i];
    }}
    signals_all=wavelet(srcpos_tmp, nshots, shape, nt, delta);
    if(flag==1){
        for(ishot=1;ishot<=nshots;ishot++){
            sprintf(file_out,"%s/wavelet_shot%d.su",fpath,ishot);
            printf("%s\n",file_out);
            Sxy[1][1]=srcpos[1][ishot];
            Sxy[2][1]=srcpos[2][ishot];
            for(it=1;it<=nt;it++) signals[1][it]=signals_all[ishot][it];
            write_su(file_out, Sxy[1][1], Sxy[2][1], Sxy, signals, 1, nt, delta);
        }
    }else{
        sprintf(file_out,"%s/wavelet.su",fpath);
        write_su(file_out, 0.0, 0.0, srcpos, signals_all, nshots, nt, delta);
    }
    
    
    free_matrix(srcpos_tmp,1,8,1,nshots);
    
    free_matrix(srcpos,1,3,1,nshots);
    free_matrix(Sxy,1,3,1,1);
    free_matrix(signals,1,1,1,nt);
    
    return 0;
}





