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
#include <sacio.h>
#include <sachd.h>
#include "source.h"


/* Define the maximum length of the data array */
float *read_sac(char filename[256], int *ns, float *dt){
    /* This routine gets the trace header information. */
    FILE *fp;
    long size=0;
    
    /* Define variables to be used in the call to rsac1() and getfhv() */
    int max, nerr;
    float *yarray=NULL , beg;
    
    fp = fopen(filename,"r");
    if (fp==NULL) {
        printf(" It was not able to read SAC file: %s\n",filename);
        exit(0);
    }
    fseek(fp,0L,SEEK_END);
    size=ftell(fp);
    max=(size-158*4)/4;
    fclose(fp);
    
    yarray=(float *)malloc((size_t) (max*sizeof(float)));
    /* Read in the SAC File */
    rsac1(filename, yarray, ns, &beg, dt, &max, & nerr, strlen( filename ) ) ;
    /* here *ns and max is same (max couldn't be less than ns )*/
    //printf("n=%d %d\n",*ns, max);
    
    return yarray;
}

//gcc -o sac_operators sac_operators.c -I/home/tao/seis_software/seisplatform/FD/SOFI2D/src/IN_OUT/include -L/home/tao/seis_software/SAC/lib -lsac -lsacio -lm
// int main(int argc, char *argv[]){
//     int ns;
//     float dt;
//     char filein[256];
//     sscanf(argv[1], "%s", filein);
//     
//     read_sac(filein, &ns, &dt);
//     
//     printf("%d %f\n",ns, dt);
// }
    
    