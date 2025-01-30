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

/*------------------------------------------------------------------------
 *   write local model to file
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void writemod(char modfile[STRING_SIZE], float ** array, int format){
    
    
    /* extern variables */
    extern int MYID, NX, NY, POS[3], IDX, IDY;//,VERBOSE;
    extern FILE *FP;
//     extern int VERBOSE;
    
    int i, j;
    FILE *fpmod;
    char file[STRING_SIZE];
    
//     if(VERBOSE) fprintf(FP,"\n\n PE %d is writing model to \n",MYID);
    sprintf(file,"%s.%i.%i",modfile,POS[1],POS[2]);
//     if(VERBOSE) fprintf(FP,"\t%s\n\n", file);
    fpmod=fopen(file,"w");
    for (i=1;i<=NX;i+=IDX)
        for (j=1;j<=NY;j+=IDY)
            writedsk(fpmod,array[j][i],format);
				
    fclose(fpmod);
    
    
}


