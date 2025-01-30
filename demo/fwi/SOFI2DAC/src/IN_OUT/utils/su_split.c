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
@brief Split up SU-format file into multiple files according to the number of traces. 
@details One can calculate the independent adjoint source at the 
        individual receiver position.
@file su_split.c
@author Tao Lei (leit@sustech.edu.cn), Wei Zhang
@date 23 Sep 2019
@param[in]  fin  the input path of original SU-format file;
@param[in]  fout the prefix of output path for each trace <em>(default: out)</em>;
@code 
su_split fin=data/toy_example_p.su.shot1 fout=out
@endcode
*/

#include "par.h"
#include "source.h"

/* Global variables */
int xargc; 
char **xargv;

int main(int argc, char *argv[]){
    
    /* IN and OUT */
    char *fin, *fout;
    
    /* Define variables */
    int ntr, ns, itr, it;
    float dt, Sx, Sy;
    float **pos=NULL, **data=NULL;
    float **pos_tr=NULL, **sig=NULL;
    char f_out[256];
    
    /* initialize getpar */
    initargs(argc,argv);
    if(!getparstring("fin", &fin)){
        printf("command as: ./su_split fin=data/toy_example_p.su.shot1 fout=out\n");
        throw_error("parameter [fin] must be specified!\n");}
    if (!getparstring("fout", &fout))  fout= "out";
    
    pos=read_su_header_all(fin, &ntr, &ns, &dt, &Sx, &Sy);
    data = read_su(fin, &ntr, &ns);
    pos_tr = matrix(1,3,1,1);
    sig = matrix(1, 1, 1, ns);
    
    for(itr=1;itr<=ntr;itr++){
        for(it=1;it<=ns;it++){
            sig[1][it]= data[itr][it];
        }
        pos_tr[1][1] = pos[1][itr];
        pos_tr[2][1] = pos[2][itr];
        pos_tr[3][1] = pos[3][itr];
        /* update the file name for each trace */
        sprintf(f_out, "%s_%d.su",fout, itr);
        /* write adjoint source */
        write_su(f_out, Sx, Sy, pos_tr, sig, 1, ns, dt);
    }
    
    free_matrix(pos, 1, 3, 1, ntr);
    free_matrix(data, 1, ntr, 1, ns);
    free_matrix(pos_tr, 1, 3, 1, 1);
    free_matrix(sig, 1, 1, 1, ns);

    return 0;
}