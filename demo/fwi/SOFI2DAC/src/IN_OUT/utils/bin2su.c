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
@brief converts bin-format into SU-format.
@file bin2su.c
@author Tao Lei (leit@sustech.edu.cn), Wei Zhang
@date 09 July 2021
@param[in] sx source position x (in meter), <em>(optional & default: 0.0)</em>;
@param[in] sy source position y (in meter), <em>(optional & default: 0.0)</em>;
@param[in] delta the time interval (in second);
@param[in] fin the data seismogram in bin-format, which is in dimension of (ntr, nt);
@param[in] fout the output file, <em>(optional & default: out.su)</em>;
@param[in] frec the receiver file in ASCII format (nrec rows * 2 columns), as matrix [rx1, ry1; rx2, ry2; ...; rxn, ryn].
@code 
bin2su sx=10.0 sy=20.0 delta=1e-3 fin=obs_p_final.bin frec=receiver.dat
@endcode
*/

#include "par.h"
#include "source.h"

/* Global variables */
int xargc; 
char **xargv;

int main(int argc, char *argv[]){

    /*local variables */
    int itr, j, ntr, nt, size;
    float **recpos_tmp=NULL,**recpos=NULL, **data=NULL;
    char file_out[256];
    FILE *fp;
    
    float delta, sx, sy;
    char *fin, *fout, *frec;
    
    /* 1. read formatted input list */
    /* initialize getpar */
    initargs(argc,argv);
    
    if(!getparstring("fin", &fin))   throw_error("fin must be specified!\n");
    if(!getparstring("frec", &frec)) throw_error("frec must be specified!\n");
    if(!getparstring("fout", &fout)) fout= "out.su";
    if(!getparfloat("delta",&delta)) throw_error("delta must be specified!\n");
    if(!getparfloat("sx",&sx)) sx=0.0;
    if(!getparfloat("sy",&sy)) sy=0.0;
    printf("fin=%s\n",fin);
    printf("frec=%s\n",frec);
    printf("fout=%s\n",fout);
    printf("delta=%e s\n",delta);
    printf("sx=%e m\n",sx);
    printf("sy=%e m\n",sy);
    

    recpos_tmp =read_text(2, 2, &ntr,frec);
    recpos =matrix(1,2,1,ntr);
    for(itr=1;itr<=ntr;itr++){
        for(j=1;j<=2;j++){
            recpos[j][itr]=recpos_tmp[itr][j];
    }}
    fp = fopen(fin,"r");
    if (fp==NULL) {
        printf(" It Was not able to read BIN file: %s\n",fin);
        exit(1);
    }else{
        fseek(fp,0L,SEEK_END);
        size=ftell(fp);
        nt=size/4/ntr;
        data =matrix(1, ntr, 1, nt);
        V_read(fin,&data[1][1],ntr*nt);
    }
    write_su(fout, sx, sy, recpos, data, ntr, nt, delta);
    
    free_matrix(recpos_tmp,1,ntr,1,2);
    free_matrix(recpos,1,2,1,ntr);
    free_matrix(data,1, ntr, 1, nt);
    return 0;
}