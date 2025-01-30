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
@brief Add the seawater layer on the top of the model
@file add_seawater.c
@author Tao Lei (leit@sustech.edu.cn), Wei Zhang
@date 22 Sep 2019
@param[in] verbose \b 1: print the details; \b 0(\e default): run in quiet pattern;
@param[in] n1   the number of grids in Y direction;
@param[in] n2   the number of grids in X direction <em>(optional)</em>;
@param[in] nadd the number of water layer grids;
@param[in] val  the value of water layer<em>(optional & default: 1500.0)</em>;
@param[in] fin  the model file in binary file;
@param[in,out] fout the new model file in binary file <em>(optional & default: out.bin)</em>;
@code 
add_seawater verbose=1 n1=350 nadd=15 fin=Marmousi2_10m.vp fout=new.bin
@endcode
*/

#include "par.h"

/* Global variables */
int xargc; 
char **xargv;

int main(int argc, char *argv[]) {
    /* IN and OUT */
    int verbose, n1, n2, nadd;
    float val;
    char *fin, *fout;
    
    /* Local variables */
    int i, size1;
    float *water_layer=NULL, *media=NULL;
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
    if(!getparint("nadd",&nadd)) throw_error("nadd must be specified!\n");
    val=1500.0; getparfloat("val",&val);
    
    if(verbose){
        printf("\n******\n");
        printf("n2      \t| %d\n",n2);
        printf("n1_new  \t| %d\n",n1+nadd);
        printf("n1      \t| %d\n",n1);
        printf("nadd    \t| %d\n",nadd);
        printf("val     \t| %.3f\n",val);
        printf("fin     \t| %s\n",fin);
        printf("fout    \t| %s\n",fout);
        printf("******\n");
    }
    
    /* allocate memory */
    water_layer = vector(1, nadd);
    media = vector(1, n1);
    for (i=1;i<=nadd;i++){
        water_layer[i]=val;
    }
    
    /* 2. read BIN-format file */
    fp=fopen(fin,"r");
    fp_out=fopen(fout,"w");
    for (i=1;i<=n2;i++){
        fread(&media[1],sizeof(float),n1,fp);
        fwrite(&water_layer[1],sizeof(float),nadd,fp_out);
        fwrite(&media[1],sizeof(float),n1,fp_out);
    }
    fclose(fp);
    fclose(fp_out);
    
    /* deallocate memory */
    free_vector(water_layer, 1, nadd);
    free_vector(media, 1, n1);
    
    return 0;
}

