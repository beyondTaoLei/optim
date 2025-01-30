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
@brief Print the information of the SU-file
@file suinfo.c
@author Tao Lei (leit@sustech.edu.cn), Wei Zhang
@date 18 Sep 2019
@param[in]  fin the full path of su file;
@retval stdout \b ntr: the number of traces 
@retval stdout \b ns: the time samples, 
@retval stdout \b delta: time inteval (s).
@retval stdout \b time: the total time (s).
@retval stdout <b> Sx, Sy </b>: the source position.
@retval stdout <b> Rx, Ry </b>: the receivers' position.
@code 
./suinfo fin=data/data/toy_example_vy.su.shot1
@endcode
*/

#include "par.h"
#include "source.h"

/* Global variables */
int xargc; 
char **xargv;

int main(int argc, char *argv[]) {
    
    /* IN and OUT */
    char *fin;
    
    /* Local variables */
    int ns, ntr, itr;
    //int maxlag=NS/2;
    float delta, Sx, Sy;
    float **pos=NULL;
    
    /* 1. read formatted input list */
    /* initialize getpar */
    initargs(argc,argv);
    if(!getparstring("fin", &fin)) throw_error("fin must be specified!\n");

    /* 2. read SU-format file */
    pos=read_su_header_all(fin, &ntr, &ns, &delta, &Sx, &Sy);
    printf("\n******\n");
    printf("fin\t\t| %s\n",fin);
    printf("ntr, ns, delta\t| %d %d %.3e \n",ntr, ns, delta);
    printf("time\t\t| %.3e \n",ns*delta);
    printf("Sx, Sy\t\t| %.3e %.3e\n", Sx, Sy);
    printf("******\n\n");
    
    printf("*** Receivers ***\n");
    
    for(itr=1;itr<=ntr;itr++){
        printf("%.3e %.3e \n", pos[1][itr],pos[2][itr]);
    }
    
    /* deallocate memory */
    free_matrix(pos, 1, 3, 1, ntr);
    
    return 0;
}

