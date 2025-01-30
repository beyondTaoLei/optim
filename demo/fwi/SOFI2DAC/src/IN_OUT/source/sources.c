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

float **sources(char src_file[256], int *nsrc){
    
    float **srcpos = NULL;
    int   l, current_source=0, current_azimuth=0, nvarin=0;
    float xsrc, ysrc, zsrc, tshift, fc = 0.0;
    char  cline[256];
    FILE *fpsrc;
    
    int source_type=1;
    
    fpsrc=fopen(src_file,"r");
    if (fpsrc==NULL){
        printf(" Source file could not be opened !");
        exit(0);
    }
    *nsrc=0;


    /* read number of source positions */	
    fscanf(fpsrc,"%d",nsrc);
                            
    printf(" Number of source positions specified in %s : %d \n",src_file,*nsrc);

    srcpos=matrix(1,8,1,*nsrc);

    rewind(fpsrc);		
    fgets(cline,255,fpsrc);		/* Dummy fgets for ignoring first line */

    for (l=1;l<=*nsrc;l++){
            fgets(cline,255,fpsrc);

            nvarin=sscanf(cline,"%f%f%f%f%f%f%f%f",&xsrc, &zsrc, &ysrc, &tshift, &srcpos[5][l], &srcpos[6][l], &srcpos[7][l], &srcpos[8][l]);
            switch(nvarin){
                    case 0: xsrc=0.0;
                    case 1: zsrc=0.0;
                    case 2: ysrc=0.0;
                    case 3: printf(" No time shift defined for source %i in %s!\n",l, src_file);
                    case 4: printf(" No frequency defined for source %i in %s!\n",l, src_file);
                    case 5: printf(" No amplitude defined for source %i in %s!\n",l, src_file);
                    case 6: srcpos[7][l]=0.0;
                    case 7: srcpos[8][l]=source_type;

            }
            if ((srcpos[8][l]!=4) && (nvarin>6)) {
                current_source=(int)srcpos[8][l];
                printf(" source_type of source #%i is specified as %i, SOURCE_AZIMUTH is ignored.\n", l, current_source);
            }
            if ((srcpos[7][l]==0) && (srcpos[8][l]==4)) {
                current_azimuth=(int)srcpos[7][l];
                printf("\n WARNING: source_type of source #%i is specified as rotated force, but no SOURCE_AZIMUTH is specified in the source file! The SOURCE_AZIMUTH is set to %i degrees.\n",l,current_azimuth);
            }
            /* fscanf(fpsrc,"%f%f%f%f%f",&xsrc, &ysrc, &tshift, &fc, &amp); */ 
            srcpos[1][l]=xsrc;
            srcpos[2][l]=ysrc;
            srcpos[3][l]=0.0;
            srcpos[4][l]=tshift;
            fc=srcpos[5][l];
    }

    fclose(fpsrc);

    /* Compute maximum frequency */
    for (l=1;l<=*nsrc;l++)
            if (srcpos[5][l]>fc) fc=srcpos[5][l];
    printf(" Maximum frequency defined in %s: %6.2e Hz\n",src_file,fc);
    //TS=1.0/fc;

    /* outputs all sources per each subdomain / node*/


    printf(" number\t    x\t\t    y\t\t  tshift\t    fc\t\t   amp\t	source_azimuth\tsource_type\n");
    for (l=1;l<=*nsrc;l++)
            printf("    %i \t %6.2f \t %6.2f \t %6.2f \t %7.3f \t %6.2f   \t %6.2f  \t   %1.0f\n\n",
            l, srcpos[1][l],srcpos[2][l],srcpos[4][l],srcpos[5][l],srcpos[6][l],srcpos[7][l],srcpos[8][l]);

    return srcpos;
}