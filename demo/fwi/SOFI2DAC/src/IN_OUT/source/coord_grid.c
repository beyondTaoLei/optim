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

/* run by all of the nodes in communicator MPI_COMM_WORLD or MPI_COMM_SHOT*/
int **coord_grid(float dh, float *refrec, float **coord, int **repeated_tr,
                 int ntr_prev, int *ntr){

    //refrec[4]
    
    int **recpos1, **recpos = NULL;//, nxrec=0, nyrec=0, nzrec=0;
    int   itr=1, itr1=0, itr2=0, recflag=0, i;
    //int nxrec1, nxrec2, nyrec1, nyrec2, nzrec1=1, nzrec2=1;
        
        recpos1=imatrix(1,3,1,ntr_prev);
        for (itr=1;itr<=ntr_prev;itr++){
                recpos1[1][itr]=iround((coord[1][itr]+refrec[1])/dh);
                recpos1[2][itr]=iround((coord[2][itr]+refrec[2])/dh);
                recpos1[3][itr]=iround((0.0+refrec[3])/dh);
        }
        /*if(myid==0){
            printf(" Message from function coord_grid (written by PE %d):\n",myid);
            printf(" Number of positions found: %i\n",ntr_prev);
        }*/

        /* recpos1[1][itr] X in m
        * recpos1[2][itr] Y in m
        *
        */

        /* check if more than one receiver is located
                                    at the same gridpoint */
        for (itr=1;itr<=(ntr_prev-1);itr++)
                for (itr1=itr+1;itr1<=ntr_prev;itr1++)
                        if ((recpos1[1][itr]==recpos1[1][itr1])
                            && (recpos1[2][itr]==recpos1[2][itr1])
                            && (recpos1[3][itr]==recpos1[3][itr1]))
                                recpos1[1][itr1]=-(++recflag);

        recpos=imatrix(1,3,1,ntr_prev-recflag);
        *repeated_tr=ivector(1,recflag);
        i=1;
        for (itr=1;itr<=ntr_prev;itr++)
                if (recpos1[1][itr]>0){
                        recpos[1][++itr2]=recpos1[1][itr];
                        recpos[2][itr2]=recpos1[2][itr];
                        recpos[3][itr2]=recpos1[3][itr];
                }else{
                    (*repeated_tr)[i++]=itr;
                }

        *ntr=itr2;
        if (recflag>0){
            printf("\n\n");
            printf(" Warning:\n");
            printf(" Several points located at the same gridpoint !\n");
            printf(" Number of receivers or sources reduced to %i\n", *ntr);
            printf("\n\n");
        }

        free_imatrix(recpos1,1,3,1,ntr_prev);
        
	return recpos;
}