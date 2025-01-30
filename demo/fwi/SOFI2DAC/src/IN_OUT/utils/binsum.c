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

/*----------------------------------------------------------------
 *   sum up all the BIN-files                                  
 *---------------------------------------------------------------*/
#include "comm.h"
#include "par.h"

#define MAX_STACKS 1000

/* Global variables */
int xargc; 
char **xargv;

int main(int argc, char *argv[]) {
/*
    
    INPUT:
        fin: list of partial stack volume filenames for input (fields delimited by commas)
    OUTPUT:
        fout: the file of the summation.
*/

/* 
gcc -o binsum binsum.c -I../include -L../lib -lIN_OUT -lm
./binsum verbose=1 fin=data/obs1.bin,data/obs2.bin,data/obs3.bin fout=obs4.bin
*/
    /* IN and OUT */
    int verbose, opt;
    char *fin, *fout;
    
    /* Local variables */
    int nIN, size, num_bin, i, j;
    float *data=NULL, *sum1=NULL;
    char* inFilename[MAX_STACKS];
    FILE *fp;
    /* 1. read formatted input list */
    /* initialize getpar */
    initargs(argc,argv);
    
    verbose=0; getparint("verbose",&verbose);
    opt=1; getparint("opt",&opt);
    if(!(nIN=countparval( "fin" ))) err( "input file list not given" );
    if(!getparstring("fout", &fout)) fout= "out.bin";
    
    for(i=1; i<=nIN; i++) getnparstringarray(i, "fin", &(inFilename[i-1]));
    size=get_file_size(inFilename[0]);
    num_bin=size/4;

    for(i=2; i<=nIN; i++){
        if(size!=get_file_size(inFilename[i-1])) 
            err("the number of elements in input files is different!");
    }
    
    if(verbose){
        printf("\n******\n");
        printf("There are %d binary files in total.\n",nIN);
        printf("and, %d elements in every file.\n",num_bin);
        printf("fout: %s\n",fout);
        printf("******\n");
    }
    
    /* 2. read BIN-format file and sum up them*/
    data = vector(1, num_bin);
    sum1 = vector(1, num_bin);
    for(i=1; i<=nIN; i++){
        V_read(inFilename[i-1],&data[1],num_bin);
        for(j=1; j<=num_bin; j++) sum1[j] +=data[j];
    }
    
    /* 3. write adjoint source */
    V_write(fout, &sum1[1], num_bin);
    
    /* deallocate memory */
    free_vector(data, 1, num_bin);
    free_vector(sum1, 1, num_bin);

    return 0;
}

