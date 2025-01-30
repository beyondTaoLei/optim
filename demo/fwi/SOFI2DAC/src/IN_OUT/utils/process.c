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
@brief Preprocess the seismograms
@details This program includes some common manipulation operators 
applied to seismograms, such as, 

<em>

    filtering,
    
    integration with respect to time,
    
    time windowing, 
    
    normalization,
    
    trace killing.
</em>
@file process.c
@author Tao Lei (leit@sustech.edu.cn), Wei Zhang
@date 20 Jul 2019
@param[in]  verbose \b 1: print the details; \b 0: run in quiet pattern, default value;
@param[in]  fin  the input path of original seismograms;
@param[in]  fout the output path of manipulated seismograms;
//
@param[in]  filt prototype of filter(\e default: 0),
                \b 1: butterworth, 
                \b 2: Bessel, 
                \b 3: Chebyshev type I, 
                \b 4: Chebyshev type II, 
                \b 0: no filtering;
@param[in]  type type of filter (\e default: 1), 
                \b 1: bandpass, 
                \b 2: bandreject, 
                \b 3: lowpass, 
                \b 4: highpass;
@param[in]  order number of poles or order of the analog prototype (\e default: 6), 4-5 should be ample, it cannot exceed 10;
@param[in]  zero_phase  number of passes (\e default: 0), 
                \b 0: Forward filter only; 
                \b 1: Forward and reverse filtering;
@param[in]  fcmin low frequency cutoff [ hertz ], which is ignored on Sac_filter_lowpass;
@param[in]  fcmax high frequency cutoff [ hertz ], which is ignored on Sac_filter_highpass;
//
@param[in]  velocity integration of seismograms (\e default: 1); 
                \b 1: no integration; 
                \b 0: apply integration;
//
@param[in]  timewin  time window (\e default: 0), one can choose, 
                \b 0: no time window; 
                \b 1: apply time window with the same starting (\e t1) and end time (\e t2);
                \b 2: apply time window with first-breaks picking time;
@param[in]  t1 the starting time of time window [ seconds ], <em>for timewin=1</em>;
@param[in]  t2 the end time of  time window [ seconds ],  <em>for timewin=1</em>;
@param[in]  ltwin the length of time window [ seconds ], <em>for timewin=2</em>;
@param[in]  ftwin the input path of first-breaks picking time, which consists of two column, 
        1st column is the trace index, 2nd column is the picking time,  <em>for timewin=2</em>;

@param[in]  normalize normalization (\e default: 0); 
                \b 1: normalize all traces to their absolute maximum;
                \b 2: normalize all traces to the whole seismogram’s global maximum;
                \b 0: not normalize;
//
@param[in]  trkill apply trace killing (\e default: 0); 
                \b 1: kill traces with trace number between [ntr1, ntr2];
                \b 2: keep traces with trace number between [ntr1, ntr2];
                \b 3: kill traces with offsets bewteen [ntr1, ntr2];
                \b 4: keep traces with offsets bewteen [ntr1, ntr2];
@param[in]  ntr1 the starting trace number <em>(trkill: 1 or 2)</em> or offset <em>(trkill: 3 or 4)</em> of the killing window
@param[in]  ntr2 the end trace number <em>(trkill: 1 or 2)</em> or offset <em>(trkill: 3 or 4)</em> of the killing window
//
@param[in]  xtwin apply xtwin (\e default: 0)
@code
gcc -o process process.c -I../include -L../lib -L/home/tao/seis_software/SAC/lib -lIN_OUT -lsac -lsacio -lm -lfftw3

process verbose=0 fin=data/obs/toy_example_vy.su.shot1 fout=out.su filt=1 type=1 order=4 zero_phase=1 fcmin=15 fcmax=20 velocity=1 timewin=1 t1=50 t2=60 normalize=0 trkill=1 ntr1=10 ntr2=30 xtwin=0  

@endcode
*/

#include "par.h"
#include "source.h"
#include "process.h"

/* Global variables */
int xargc; 
char **xargv;

int main(int argc, char *argv[]) {
    
    /* IN and OUT */
    int filt, velocity, timewin, normalize, trkill, xtwin;
    int type, order, zero_phase, ntr1, ntr2, verbose;
    float t1, t2, fcmin, fcmax, ltwin;
    char *fxtwin, *fin, *fout, *ftwin;

    /* Local variables */
    int ns, ntr, itr, it, it1, it2, it0, it3, i, *ktable=NULL;
    int nrow, ipts=20, ltwin_i;
    float delta, Sx, Sy, offset, value, omega = PI/(2.*ipts);
    float **pos=NULL, **obs=NULL, *trace=NULL, **twin=NULL;
    
    /* initialize getpar */
    initargs(argc,argv);
    
    /* 1. read formatted input list */
    verbose=0; getparint("verbose",&verbose);
    if(!getparstring("fin", &fin)) throw_error("fin must be specified!\n");
    if(!getparstring("fout", &fout)) throw_error("fout must be specified!\n");
    /*filtering*/
    filt=0; getparint("filt",&filt);
    if(filt){
        type=1;         getparint("type",&type);
        order=6;        getparint("order",&order);
        zero_phase=0;   getparint("zero_phase",&zero_phase);
        if(!getparfloat("fcmin",&fcmin)) throw_error("fcmin must be specified!\n");
        if(!getparfloat("fcmax",&fcmax)) throw_error("fcmax must be specified!\n");
    }
    /*integration*/
    velocity=1;     getparint("velocity",&velocity);
    /*timewindow*/
    timewin=0;      getparint("timewin",&timewin);
    if(timewin==1){
        if(!getparfloat("t1",&t1)) throw_error("t1 must be specified!\n");
        if(!getparfloat("t2",&t2)) throw_error("t2 must be specified!\n");
    }else if(timewin==2){
        if(!getparfloat("ltwin",&ltwin)) throw_error("ltwin must be specified!\n");
        if(!getparstring("ftwin", &ftwin)) throw_error("ftwin must be specified!\n");
    }
    /*normalization*/
    normalize=0;     getparint("normalize",&normalize);
    /*tracekilling*/
    trkill=0;     getparint("trkill",&trkill);
    if(trkill){
        if(!getparint("ntr1",&ntr1)) throw_error("ntr1 must be specified!\n");
        if(!getparint("ntr2",&ntr2)) throw_error("ntr2 must be specified!\n");
    }
    /*xt_window*/
    xtwin=0;     getparint("xtwin",&xtwin);
    if(xtwin){
        if(!getparstring("fxtwin",&fxtwin)) throw_error("fxtwin must be specified!\n");
    }
    
    if(verbose){
        printf("\n******\n");
        printf("%-18s| %s\n","fin", fin);
        printf("%-18s| %s\n","fout",fout);
        printf("%-18s| %d\n","filt",filt);
        if(filt){
            printf("  %-16s| %d\n","type",type);
            printf("  %-16s| %d\n","order",order);
            printf("  %-16s| %d\n","zero_phase",zero_phase);
            printf("  %-16s| %f\n","fcmin",fcmin);
            printf("  %-16s| %f\n","fcmax",fcmax);
        }
        printf("%-18s| %d\n","velocity",velocity);
        
        printf("%-18s| %d\n","timewin",timewin);
        if(timewin==1){
            printf("  %-16s| %f\n","t1",t1);
            printf("  %-16s| %f\n","t2",t2);
        }else if(timewin==2){
            printf("  %-16s| %f\n","ltwin",ltwin);
            printf("  %-16s| %s\n","ftwin",ftwin);
        }
        printf("%-18s| %d\n","normalize",normalize);
        printf("%-18s| %d\n","trkill",trkill);
        if(trkill){
            printf("  %-16s| %d\n","ntr1",ntr1);
            printf("  %-16s| %d\n","ntr2",ntr2);
        }
        printf("%-18s| %d\n","xtwin",xtwin);
        if(xtwin){
            printf("  %-16s| %s\n","fxtwin",fxtwin);
        }
        printf("******\n");
    }
    
    /* 2. read SU-format file */
    pos=read_su_header_all(fin, &ntr, &ns, &delta, &Sx, &Sy);
    obs = read_su(fin, &ntr, &ns);
    trace = vector(1, ns);
    ktable = ivector(1, ntr);
    if(verbose) printf("ntr, ns, delta= %d %d %e\n", ntr, ns, delta);
    
    /* 3. preprocess the seismogram */
    //apply filtering
    if(filt){
        /* filt:: 1: butterworth, 2: bessel, 3: Chebyshev Type I, 4: Chebyshev Type II 
           type:: 1: BANDPASS, 2: BANDREJECT, 3: LOWPASS, 4: HIGHPASS */
        filtering(obs, ntr, ns, delta, fcmin, fcmax, order, zero_phase, filt, type);
    }
    
    //apply integration
    if(velocity!=1){
        for(itr=1;itr<=ntr;itr++){
            for(it=1;it<=ns;it++) trace[it] = obs[itr][it];
            Taper(&trace[1], 3, ns, 0.01);
            V_int(ns, 1., &trace[1], &trace[1]);
            for(it=1;it<=ns;it++) obs[itr][it] = trace[it];
        }
        
    }
    
    //apply timewindow
    if(timewin==1){// it keeps the seismogram from t1 to t2.
        it1 = (int)(t1/delta+1);
        it2 = (int)(t2/delta+1);
        //printf("%d %d\n",it1,it2);
        if(it1<1) it1=1;
        if(it2>ns-1) it2=ns-1;
        for(itr=1;itr<=ntr;itr++){
            for(it=1;it<it1;it++){
                obs[itr][it] = 0.0;
            }
            for(it=it2+1;it<=ns;it++){
                obs[itr][it] = 0.0;
            }
        }
    }else if(timewin==2){
        for(itr=1;itr<=ntr;itr++) ktable[itr]=1;
        twin=read_text(2, 2, &nrow, ftwin);
        ltwin_i=(int)(ltwin/delta);
        for(i=1;i<=nrow;i++){
            itr=(int) twin[i][1];
            if(twin[i][2]<0.0){/* kill this trace */
                ktable[itr]=1;
            }else{/* keep this trace */
                ktable[itr]=0;
                it1 = (int)(twin[i][2]/delta+1);
                it2 = it1+ltwin_i;
                if(it1<1) it1=1;
                if(it2>ns-1) it2=ns-1;
                it0=it1-ipts;
                it3=it2+ipts;
                if(it0<1) it0=1;
                if(it3>ns-1) it3=ns-1;
                for(it=1;it<it0;it++){
                    obs[itr][it] = 0.0;
                }
                for(it=it0;it<it1;it++){
                    value = 1.0-cos(omega*(it-it0));
                    obs[itr][it] *= value;
                }
                for(it=it2+1; it<it3+1;it++){
                    value = 1.0-cos(omega*(it3-it));
                    obs[itr][it] *= value;
                }
                for(it=it3+1;it<=ns;it++){
                    obs[itr][it] = 0.0;
                }
            }
        }
        for(itr=1;itr<=ntr;itr++){
            if(ktable[itr]==1){
                for(it=1;it<ns;it++){
                    obs[itr][it] = 0.0;
                }
            }
        }
        for(itr=1;itr<=ntr;itr++) ktable[itr]=0;
        free_matrix(twin, 1, nrow, 1, 2);
    }
    
    //apply normalization
    if(normalize==1){//normalize all traces to their absolute maximum
        for(itr=1;itr<=ntr;itr++){
            for(it=1;it<=ns;it++) trace[it] = obs[itr][it];
            V_normalize(1, ns, trace);
            for(it=1;it<=ns;it++) obs[itr][it] = trace[it];
        }
    }else if(normalize==2){//normalize all traces to the whole seismogram’s global maximum
        M_normalize(1, ns, 1, ntr, obs);
    }
    
    //apply tracekilling
    if(trkill){
        if(trkill==1){ //kill [ntr1, ntr2]
            if(ntr1<1) ntr1=1;
            if(ntr2>ntr) ntr2=ntr;
            for(itr=ntr1;itr<=ntr2;itr++) ktable[itr] = 1;
        }else if(trkill==2){//keep [ntr1, ntr2]
            for(itr=1;itr<ntr1;itr++) ktable[itr] = 1;
            for(itr=ntr2+1;itr<=ntr;itr++) ktable[itr] = 1;
        }else if(trkill==3||trkill==4){
            for(itr=1;itr<=ntr;itr++){//trkill==3 or 4
                offset = sqrt((pos[1][itr] - Sx)*(pos[1][itr] - Sx) + 
                        (pos[2][itr] - Sy)*(pos[2][itr] - Sy));
                //printf("%f %f %f %f %f\n",offset , pos[1][itr], pos[2][itr],Sx, Sy);
                if(offset>=ntr1&&offset<=ntr2) ktable[itr] = 1; //kill
            }
            if(trkill==4){
                for(itr=1;itr<=ntr;itr++){
                    ktable[itr] = 1-ktable[itr]; //keep
                }
            }
        }
        /* kill the specified traces */
        for(itr=1;itr<=ntr;itr++){
            if(ktable[itr]==1){
                for(it=1;it<=ns;it++) obs[itr][it] = 0.0;
            }
        }
    }
    
    //apply xtwin
    if(xtwin){
        
    }
    
    /* 4. write manipulated seismogram */
    write_su(fout, Sx, Sy, pos, obs, ntr, ns, delta);
    //////////////////////
    
    free_matrix(pos, 1, 3, 1, ntr);
    free_matrix(obs, 1, ntr, 1, ns);
    free_vector(trace,1, ns);
    free_ivector(ktable, 1, ntr);
    
    if(verbose) printf("finished...\n");
    return 0;
}
