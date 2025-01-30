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
@brief Reverse the seismograms along the time axis (1:1:NT <--->NT:-1:1)
@file rtime.c
@author Tao Lei (leit@sustech.edu.cn), Wei Zhang
@date 18 Sep 2019
@param[in] verbose \b 1: print the details; \b 0(\e default): run in quiet pattern;
@param[in] fin the SU file;
@param[in] fout the reversed SU file <em>(optional & default: out.su)</em>;
@code 
rtime verbose=1 fin=data/syn/test_p_shot1.su fout=data/adj/test_p_shot1.su
@endcode
*/

#include "par.h"
#include "source.h"

/* Global variables */
int xargc; 
char **xargv;

int main(int argc, char *argv[]) {
    /* IN and OUT */
    int verbose;
    char *fin, *fout;
    
    /* Local variables */
    int ns, ntr, itr, it, invtime;
    //int maxlag=NS/2;
    float delta, Sx, Sy;
    float **pos=NULL, **data=NULL, *tr_data=NULL;
    
    /* 1. read formatted input list */
    /* initialize getpar */
    initargs(argc,argv);

    verbose=0; getparint("verbose",&verbose);
    if(!getparstring("fin", &fin))
        throw_error("fin must be specified!\n");
    if (!getparstring("fout", &fout))  fout= "out.su";
    if(verbose){
        printf("\n******\n");
        printf("fin\t\t| %s\n",fin);
        printf("fout\t\t| %s\n",fout);
        printf("******\n");
    }
    
    /* 2. read SU-format file */
    pos=read_su_header_all(fin, &ntr, &ns, &delta, &Sx, &Sy);
    data = read_su(fin, &ntr, &ns);
    tr_data = vector(1, ns);
    
    /* 3. reverse the seismogram */
    for(itr=1;itr<=ntr;itr++){
        for(it=1;it<=ns;it++){
            tr_data[it]=data[itr][it];
        }
        invtime=ns;
        for(it=1;it<=ns;it++){
            data[itr][invtime]=tr_data[it];
            invtime--;
        }
    }
    
    /* 4. write seismogram */
    write_su(fout, Sx, Sy, pos, data, ntr, ns, delta);
    
    /* deallocate memory */
    free_matrix(pos, 1, 3, 1, ntr);
    free_matrix(data, 1, ntr, 1, ns);
    free_vector(tr_data, 1, ns);
    
    if(verbose){
        printf("ntr, ns, delta: %d %d %.3e \n",ntr, ns, delta);
        printf("Sx, Sy: %.3e %.3e\n", Sx, Sy);
    }
    return 0;
}

