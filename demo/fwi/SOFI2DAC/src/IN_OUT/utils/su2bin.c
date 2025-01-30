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
@brief converts SU-format into bin-format.
@file su2bin.c
@author Tao Lei (leit@sustech.edu.cn), Wei Zhang
@date 12 July 2021
@param[in] fin the data seismogram in su-format, which is in dimension of (ntr, nt);
@param[in] fout the output file, <em>(optional & default: out.bin)</em>;
@code 
su2bin fin=obs_p_final.su
@endcode
*/

#include "par.h"
#include "source.h"

/* Global variables */
int xargc; 
char **xargv;

int main(int argc, char *argv[]){

    /*local variables */
    int ntr, nt;
    float **data=NULL;
    
    char *fin, *fout;
    
    /* 1. read formatted input list */
    /* initialize getpar */
    initargs(argc,argv);
    
    if(!getparstring("fin", &fin))   throw_error("fin must be specified!\n");
    if(!getparstring("fout", &fout)) fout= "out.bin";
    
    data = read_su(fin, &ntr, &nt);
    printf("fin=%s\n",fin);
    printf("fout=%s\n",fout);
    printf("ntr=%d\n",ntr);
    printf("nt=%d\n",nt);
    V_write(fout, &data[1][1], ntr*nt);

    free_matrix(data,1, ntr, 1, nt);
    return 0;
}