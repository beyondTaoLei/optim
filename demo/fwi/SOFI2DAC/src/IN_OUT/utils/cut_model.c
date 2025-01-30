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
@brief Cut the model to extract one part of original model
@file cut_model.c
@author Tao Lei (leit@sustech.edu.cn), Wei Zhang
@date 22 Sep 2019
@param[in] verbose \b 1: print the details; \b 0(\e default): run in quiet pattern;
@param[in] n1   the number of grids in Y direction;
@param[in] n2   the number of grids in X direction <em>(optional)</em>;
@param[in] nc   the starting and end index in Y and X direciton 
                <em>(in order: ny1,ny2,nx1,nx2, delimited by commas)</em>;
@param[in] fin  the model file in binary file;
@param[in,out] fout the new model file in binary file <em>(optional & default: out.bin)</em>;
@code 
cut_model verbose=1 n1=350 nc=50,349,1,1700 fin=Marmousi2_10m.vp fout=new.bin
@endcode
*/

#include "par.h"

/* Global variables */
int xargc; 
char **xargv;

int main(int argc, char *argv[]) {
    /* IN and OUT */
    int verbose, n1, n2;
    char *fin, *fout;
    char* in_nyx[4]; //nc
    
    /* Local variables */
    int i, size1, nIN, n11, n12, n21, n22;
    long offset;
    float *media=NULL;
    
    FILE *fp, *fp_out;
    
    /* 1. read formatted input list */
    /* initialize getpar */
    initargs(argc,argv);

    verbose=0; getparint("verbose",&verbose);
    if(!getparstring("fin", &fin)) throw_error("fin must be specified!\n");
    if (!getparstring("fout", &fout))  fout= "out.bin";
    if(!getparint("n1",&n1)) throw_error("n1 must be specified!\n");
    if(!getparint("n2",&n2)){
        size1=get_file_size(fin);
        if(size1%(4*n1)==0){
            n2=size1/(4*n1);
        }else{
            throw_error("Make sure the size of model file is correct!\n");
        }
    }
    if(!(nIN=countparval( "nc" ))) 
        throw_error( "input nc(n11, n12, n21, n22) list not given" );
    if(nIN!=4) throw_error( "the number of input nc is not right" );
    for(i=0; i<nIN; i++) getnparstringarray(i+1, "nc", &(in_nyx[i]));
    n11 = atoi(in_nyx[0]);
    n12 = atoi(in_nyx[1]);
    n21 = atoi(in_nyx[2]);
    n22 = atoi(in_nyx[3]);
    
    if(verbose){
        printf("\n******\n");
        printf("n1, n2      \t| %d %d\n",n1, n2);
        printf("nc          \t| %d %d %d %d\n",n11, n12, n21, n22);
        printf("n1, n2(new) \t| %d %d\n",n12-n11+1, n22-n21+1);
        printf("fin         \t| %s\n",fin);
        printf("fout        \t| %s\n",fout);
        printf("******\n");
    }
    
    /* allocate memory */
    media = vector(1, n1);
    
    /* 2. read BIN-format file */
    fp=fopen(fin,"r");
    fp_out=fopen(fout,"w");
    for (i=n21;i<=n22;i++){
        offset=(long) (i-1)*n1*4;
        fseek(fp,offset,SEEK_SET);
        fread(&media[1],sizeof(float),n1,fp);
        fwrite(&media[n11],sizeof(float),n12-n11+1,fp_out);
    }
    fclose(fp);
    fclose(fp_out);
    
    /* deallocate memory */
    free_vector(media, 1, n1);
    
    return 0;
}

