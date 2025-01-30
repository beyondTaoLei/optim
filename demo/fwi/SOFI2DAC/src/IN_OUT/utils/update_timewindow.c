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
@brief Update the time window according to the difference between the 
windowing data at previous iteration and raw data at current iteration.

@details It's not necessary to pick up the first break time at every
iteration, which is based on the fact that the seismograms are similar
during the adjacent iterations. It means we can move a bit the previous 
time window to obtain new one. we employe the cross correlation to find
this difference.
@file update_timewindow.c
@author Tao Lei (leit@sustech.edu.cn), Wei Zhang
@date 29 Sep 2019
@param[in]  verbose \b 1: print the details; \b 0(\e default): run in quiet pattern;
@param[in]  maxlag the maxinum shifted time [ seconds ] during calculating the correlation;
@param[in]  fref the full path of synthetic data at the reference iteration (such as 1st iteration);
@param[in]  fsyn the full path of raw synthetic data at current iteration;
@param[in]  frefw the input path of first-breaks picking time, which consists of two column, 
        1st column is the trace index, 2nd column is the picking time;
@param[out] fout the updated picking time window for current iteration;
@code 
update_timewindow fref=data/wavelet_shot1.su fsyn=data/wavelet_shot2.su frefw=picked_times/test_p_shot1.dat fout=picked_times/test_p_shot2.dat maxlag=0.4 verbose=1
@endcode
*/

#include "par.h"
#include "source.h"

/* Global variables */
int xargc; 
char **xargv;

int main(int argc, char *argv[]) {
/*
    the misfit function includes,
    L2 misfit
    wave-equation traveltime misfit
    
    INPUT:
        the file of the manipulated seismogram.
    OUTPUT:
        the file of the adjoint source and misfit.
*/
    /* IN and OUT */
    int verbose;
    float maxlag;
    char *fref, *fsyn, *frefw, *fout;
    
    /* Local variables */
    int itr, it, ntr, ns, nrow, maxlagi, offset, lag;
    float dt, **syn, **syn2, **twin, *corr, *tr_syn, *tr_syn2;
    FILE *fp=NULL;
    
    /* 1. read formatted input list */
    /* initialize getpar */
    initargs(argc,argv);

    verbose=0; getparint("verbose",&verbose);
    if(!getparfloat("maxlag",&maxlag)) throw_error("maxlag must be specified!\n");
    if(!getparstring("fref", &fref)) throw_error("fref must be specified!\n");
    if(!getparstring("fsyn", &fsyn)) throw_error("fsyn must be specified!\n");
    if(!getparstring("frefw", &frefw)) throw_error("frefw must be specified!\n");
    if(!getparstring("fout", &fout)) throw_error("fout must be specified!\n");
    
    if(verbose){
        printf("\n******\n");
        printf("maxlag      \t| %f (s)\n",maxlag);
        printf("fref        \t| %s\n",fref);
        printf("fsyn       \t| %s\n",fsyn);
        printf("frefw       \t| %s\n",frefw);
        printf("fout        \t| %s\n",fout);
        printf("******\n");
    }
    
    /* 2. read SU-format file */
    read_su_info(fref, &ntr, &ns, &dt);
    syn = read_su(fref, &ntr, &ns);
    syn2 = read_su(fsyn, &ntr, &ns);
    twin=read_text(2, 2, &nrow, frefw);
    if(nrow!=ntr) throw_error("Please check the picking-break time file again!\n");
    
    /* allocate memory */
    maxlagi=(int)(maxlag/dt);
    corr =vector(1, 2*maxlagi+1);
    tr_syn = vector(1, ns);
    tr_syn2 = vector(1, ns);
    
    /* 3. calculate the misfit and adjoint source */
    for(itr=1;itr<=ntr;itr++){
        /* a. initialize the array for extern subroutine */
        for(it=1;it<=ns;it++){
            tr_syn[it] = syn[itr][it];
            tr_syn2[it] = syn2[itr][it];
        }
        if(twin[itr][2]>0.0){
            /* 2. calculate traveltime difference */
            xcorr(&tr_syn[1], &tr_syn2[1], &corr[1], maxlagi, ns);
            offset = V_max_index(1, 2*maxlagi+1, corr)-1;
            lag=-maxlagi + offset;
            /*
            * obs: 1st input signal in xcorr function;
            * syn: 2nd input signal in xcorr function;
            * if lag >0 denotes obs is later than syn by lag samples, 
            * otherwise, obs is earlier than syn by -lag samples.
            */
            twin[itr][2] -= lag*dt;
        } 
    }
    
    /* write in the disk */
    fp = fopen(fout, "w");
    for(itr=1;itr<=ntr;itr++){
        fprintf(fp,"%8d %8.3e\n",itr, twin[itr][2]);
    }
    fclose(fp);
    
    /* deallocate memory */
    free_matrix(syn,  1, ntr, 1, ns);
    free_matrix(syn2, 1, ntr, 1, ns);
    free_matrix(twin, 1, ntr, 1, 2);
    free_vector(corr, 1, 2*maxlagi+1);
    free_vector(tr_syn, 1, ns);
    free_vector(tr_syn2, 1, ns);
    
    if(verbose) printf("finished...\n");
    return 0;
}

