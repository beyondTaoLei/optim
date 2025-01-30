/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2016  For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS.
 *
 * IFOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 *
 * IFOS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with IFOS. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
 -----------------------------------------------------------------------------------------*/

/*  Computation of local source coordinates
 *
 *
 */

#include "fd.h"
float **splitsrc(float **srcpos,int *nsrc_loc, int nsrc)
{
    
    extern int IENDX, IENDY, POS[3], MYID;
    extern float DH;
    
    int a,b,i=0,j,k,found=0;
    float ** srcpos_dummy, **srcpos_local=NULL;
    srcpos_dummy = matrix(1,8,1,nsrc);
    
    for (j=1;j<=nsrc;j++) {
        a=(iround(srcpos[1][j]/DH)-1)/IENDX;
        b=(iround(srcpos[2][j]/DH)-1)/IENDY;
        
        if ((POS[1]==a)&&(POS[2]==b)) {
            found = 1;
            i++;
            srcpos_dummy[1][i] = (float)(((iround(srcpos[1][j]/DH)-1)%IENDX)+1);
            srcpos_dummy[2][i] = (float)(((iround(srcpos[2][j]/DH)-1)%IENDY)+1);
            srcpos_dummy[3][i] = 0.0;
            srcpos_dummy[4][i] = srcpos[4][j];
            srcpos_dummy[5][i] = srcpos[5][j];
            srcpos_dummy[6][i] = srcpos[6][j];
            srcpos_dummy[7][i] = srcpos[7][j];
            srcpos_dummy[8][i] = srcpos[8][j];
        }
    }
    
    if (i>0) srcpos_local = matrix(1,8,1,i);
    
    for (k=1;k<=i;k++){
        srcpos_local[1][k] = srcpos_dummy[1][k]; // x-coordinate in gridpoints
        srcpos_local[2][k] = srcpos_dummy[2][k]; // y-coordinate in gridpoints
        srcpos_local[3][k] = srcpos_dummy[3][k]; // z-coordinate in gridpoints
        srcpos_local[4][k] = srcpos_dummy[4][k]; // tshift
        srcpos_local[5][k] = srcpos_dummy[5][k]; // fc
        srcpos_local[6][k] = srcpos_dummy[6][k]; // amp
        srcpos_local[7][k] = srcpos_dummy[7][k];
        srcpos_local[8][k] = srcpos_dummy[8][k];
    }
    free_matrix(srcpos_dummy,1,8,1,nsrc);
    
    //    if (i>0){
    //        printf("\n **Message from splitsrc:\n");
    //        printf(" Splitting of source positions from global to local grids finished.\n");
    //        printf(" MYID= %d \t \t no. of sources= %d\n",MYID,i);
    //        printf("\n **Message from splitsrc:\n");
    //        printf(" Table of local source positions (in gridpoints), time-shift, centre frequency and amplitude:\n");
    //        printf(" MYID\t  x\t  y\t  z\t  tshift  fc\t  amp\n");
    //        for (j=1;j<=i;j++)
    //            printf(" %3d\t%4.0f\t%4.0f\t%4.0f\t%6.2f\t%6.2f\t%6.2f\n",MYID,srcpos_local[1][j],srcpos_local[2][j],srcpos_local[3][j],
    //                   srcpos_local[4][j],srcpos_local[5][j],srcpos_local[6][j]);
    //    }
    
    *nsrc_loc=i;
    return srcpos_local;
    
}
