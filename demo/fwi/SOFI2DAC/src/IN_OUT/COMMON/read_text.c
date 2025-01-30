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
/* be careful when you using!!! */
/*
 if(flag==1) the file is similar with source.dat
 if(flag==2) the file is similar with receiver.dat,
             and the last line should be empty.
 The output is matrix mat[nrow][ncol], which is [DIFFERENT] with the 
 one from source.c and receiver.c 
 */

float **read_text(int flag, int ncol, int *nrow, char filename[256]){

    int c;
    int nr, irow;
    float **mat = NULL, *vec = NULL;
    FILE *fp;

    if(ncol>6){
        printf("\n Reading values from file: \n\t%s\n",filename);
        printf(" So far we just read the file with no more than 6 columns!");
        exit(0);
    }
    /* read receiver positions from file */
    fp=fopen(filename,"r");
    if (fp==NULL){
        printf(" File %s could not be opened !",filename);
        exit(0);
    }
    
    nr=0;
    if(flag==1){/* the first line is the number of rows */
        /* read number of data */	
        fscanf(fp,"%i",&nr);
    }else{/* One should calculate the number of rows */
    /* The last line should be empty!!! */
        while ((c=fgetc(fp)) != EOF)
            if (c=='\n') ++(nr);
        rewind(fp);
    }

    vec=vector(1, ncol);
    mat=matrix(1, nr, 1, ncol);
    if(ncol==1){
        for(irow=1;irow<=nr;irow++){
            fscanf(fp,"%f",&vec[1]);
            mat[irow][1]=vec[1];
        }
    }else if(ncol==2){
        for(irow=1;irow<=nr;irow++){
            fscanf(fp,"%f %f",&vec[1],&vec[2]);
            mat[irow][1]=vec[1];
            mat[irow][2]=vec[2];
        }
    }else if(ncol==3){
        for(irow=1;irow<=nr;irow++){
            fscanf(fp,"%f %f %f",&vec[1],&vec[2],&vec[3]);
            mat[irow][1]=vec[1];
            mat[irow][2]=vec[2];
            mat[irow][3]=vec[3];
        }
    }else if(ncol==4){
        for(irow=1;irow<=nr;irow++){
            fscanf(fp,"%f %f %f %f",&vec[1],&vec[2],&vec[3],&vec[4]);
            mat[irow][1]=vec[1];
            mat[irow][2]=vec[2];
            mat[irow][3]=vec[3];
            mat[irow][4]=vec[4];
        }
    }else if(ncol==5){
        for(irow=1;irow<=nr;irow++){
            fscanf(fp,"%f %f %f %f %f",&vec[1],&vec[2],&vec[3],&vec[4],&vec[5]);
            mat[irow][1]=vec[1];
            mat[irow][2]=vec[2];
            mat[irow][3]=vec[3];
            mat[irow][4]=vec[4];
            mat[irow][5]=vec[5];
        }
    }else if(ncol==6){
        for(irow=1;irow<=nr;irow++){
            fscanf(fp,"%f %f %f %f %f %f",&vec[1],&vec[2],&vec[3],&vec[4],&vec[5],&vec[6]);
            mat[irow][1]=vec[1];
            mat[irow][2]=vec[2];
            mat[irow][3]=vec[3];
            mat[irow][4]=vec[4];
            mat[irow][5]=vec[5];
            mat[irow][6]=vec[6];
        }
    }

    fclose(fp);
    printf(" Number of data found: %i\n",nr);

    free_vector(vec, 1, ncol);
    *nrow= nr;
    
    return mat;
}