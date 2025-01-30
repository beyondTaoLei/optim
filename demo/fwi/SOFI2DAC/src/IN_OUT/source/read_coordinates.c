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

int **read_coordinates(int *ntr, float dh, float *refrec, char rec_file[256]){

    int **recpos1, **recpos = NULL;//, nxrec=0, nyrec=0, nzrec=0;
    int   itr=1, itr1=0, itr2=0, recflag=0, c;
    //int nxrec1, nxrec2, nyrec1, nyrec2, nzrec1=1, nzrec2=1;

    float xrec, yrec;
    FILE *fpr;

    /* read receiver positions from file */
    printf("\n Reading positions from file: \n\t%s\n",rec_file);
    fpr=fopen(rec_file,"r");
    if (fpr==NULL){
        printf(" Coordinates' file could not be opened !");
        exit(0);
    }
    *ntr=0;
    while ((c=fgetc(fpr)) != EOF)
            if (c=='\n') ++(*ntr);
    rewind(fpr);

    recpos1=imatrix(1,3,1,*ntr);
    for (itr=1;itr<=*ntr;itr++){
        fscanf(fpr,"%f%f\n",&xrec, &yrec);
        recpos1[1][itr]=iround((xrec+refrec[1])/dh);
        recpos1[2][itr]=iround((yrec+refrec[2])/dh);
        recpos1[3][itr]=iround((0.0+refrec[3])/dh);
    }
    fclose(fpr);
    printf(" Message from function read_coordinates_2d\n");/***/
    printf(" Number of receiver positions found: %i\n",*ntr);
    
    /* check if more than one receiver is located
                                at the same gridpoint */
    for (itr=1;itr<=(*ntr-1);itr++)
        for (itr1=itr+1;itr1<=*ntr;itr1++)
            if ((recpos1[1][itr]==recpos1[1][itr1])
                && (recpos1[2][itr]==recpos1[2][itr1])
                && (recpos1[3][itr]==recpos1[3][itr1]))
                    recpos1[1][itr1]=-(++recflag);

    recpos=imatrix(1,3,1,*ntr-recflag);
    for (itr=1;itr<=*ntr;itr++)
        if (recpos1[1][itr]>0){
            recpos[1][++itr2]=recpos1[1][itr];
            recpos[2][itr2]=recpos1[2][itr];
            recpos[3][itr2]=recpos1[3][itr];
        }

    *ntr=itr2;
    if (recflag>0){
        printf("\n\n");
        printf(" Warning:\n");
        printf(" Several receivers located at the same gridpoint !\n");
        printf(" Number of receivers reduced to %i\n", *ntr);
        printf("\n\n");
    }
    free_imatrix(recpos1,1,3,1,*ntr);

    return recpos;
}