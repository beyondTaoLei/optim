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
@file update_timewindow_one.c
@author Tao Lei (leit@sustech.edu.cn), Wei Zhang
@date 29 Sep 2019
@param[in]  verbose \b 1: print the details; \b 0(\e default): run in quiet pattern;
@param[in]  maxlag the maxinum shifted time [ seconds ] during calculating the correlation;
@param[in]  fref the full path of synthetic data at the reference iteration (such as 1st iteration);
@param[in]  fsyn the full path of raw synthetic data at current iteration;
@param[in]  bt the input path of first-breaks picking time, which consists of two column, 
        1st column is the trace index, 2nd column is the picking time;
@param[out] fout the updated picking time window for current iteration;
@code 
update_timewindow_one fref=data/wavelet_shot1.su fsyn=data/wavelet_shot2.su fout=picked_times/test_p_shot2.dat bt=0.2 maxlag=0.4 verbose=1
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
    float maxlag, bt;
    char *fref, *fsyn, *fout;
    
    /* Local variables */
    int it, ns, ntr, ntr2, maxlagi, offset, lag;
    float dt, **syn, **syn2, *corr, *tr_syn, *tr_syn2;
    FILE *fp=NULL;
    
    /* 1. read formatted input list */
    /* initialize getpar */
    initargs(argc,argv);

    verbose=0; getparint("verbose",&verbose);
    if(!getparfloat("bt",&bt)) throw_error("bt must be specified!\n");
    if(!getparfloat("maxlag",&maxlag)) throw_error("maxlag must be specified!\n");
    if(!getparstring("fref", &fref)) throw_error("fref must be specified!\n");
    if(!getparstring("fsyn", &fsyn)) throw_error("fsyn must be specified!\n");
    if(!getparstring("fout", &fout)) throw_error("fout must be specified!\n");
    
    if(verbose){
        printf("\n******\n");
        printf("bt          \t| %f (s)\n",bt);
        printf("maxlag      \t| %f (s)\n",maxlag);
        printf("fref        \t| %s\n",fref);
        printf("fsyn        \t| %s\n",fsyn);
        printf("fout        \t| %s\n",fout);
        printf("******\n");
    }
    
    /* 2. read SU-format file */
    read_su_info(fref, &ntr, &ns, &dt);
    syn = read_su(fref, &ntr, &ns);
    syn2 = read_su(fsyn, &ntr2, &ns);
    if(ntr!=ntr2||ntr!=1) throw_error("Please check the SU-file again!\n");
    
    /* allocate memory */
    maxlagi=(int)(maxlag/dt);
    corr =vector(1, 2*maxlagi+1);
    tr_syn = vector(1, ns);
    tr_syn2 = vector(1, ns);
    
    /* 3. calculate the misfit and adjoint source */
    /* a. initialize the array for extern subroutine */
    for(it=1;it<=ns;it++){
        tr_syn[it] = syn[1][it];
        tr_syn2[it] = syn2[1][it];
    }
    if(bt>0.0){
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
        bt -= lag*dt;
    } 
    
    /* write in the disk */
    fp = fopen(fout, "w");
    fprintf(fp,"%8.3e\n",bt);
    fclose(fp);
    
    /* deallocate memory */
    free_matrix(syn,  1, ntr, 1, ns);
    free_matrix(syn2, 1, ntr, 1, ns);
    free_vector(corr, 1, 2*maxlagi+1);
    free_vector(tr_syn, 1, ns);
    free_vector(tr_syn2, 1, ns);
    
    if(verbose) printf("finished...\n");
    return 0;
}

