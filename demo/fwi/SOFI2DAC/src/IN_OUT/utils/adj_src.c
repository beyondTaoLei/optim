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
@brief Calculate the adjoint sources and misfit measurement
@details This program is meant to be run to calculate the adjoint sources and misfit values, and one can choose the different measurement by changing the value of \e lnorm.
@file adj_src.c
@author Tao Lei (leit@sustech.edu.cn), Wei Zhang
@date 18 Sep 2019
@param[in]  verbose \b 1: print the details; \b 0(\e default): run in quiet pattern;
@param[in]  lnorm  \b 2(\e default): L2-norm misfit; \b 9: wave-equation traveltime misfit; 
                   \b 10: same with lnorm=9 but ignoring the time delay factor;
                   \b 11: multiple traveltime misfit adjoint sources;
                   \b 12: traveltime misfit based on Rytov approximation;
@param[in]  inv \b 1: output the adjoint source in the reverse time; \b 0(\e default): in the time sequence;
@param[in]  fsyn the full path of synthetic data;
@param[in]  fobs the full path of observed data;
@param[in]  order the SAC filtering parameter 1 for <em> lnorm=11, default= 4 </em>;
@param[in]  alpha the SAC filtering parameter 2 for <em> lnorm=11, default=0.2 </em>;
@param[in]  fmin the minimum frequency for \e lnorm=11;
@param[in]  fmax the maximum frequency for \e lnorm=11;
@param[in]  finc the frequency spacing for \e lnorm=11;
@param[out] fadj the full path of adjoint souce, <em> default: adj.su</em>;
@param[out] fmis the full path of misfit measurement, <em> default: misfit.bin </em>;
@param[out] ftime the full path of time delay, which is for the case of lnorm: 9 10 12, <em>default: timedelay.bin</em>;
@code 
adj_src verbose=1 lnorm=9 fsyn=data/syn/test_p_shot1.su fobs=data/obs/test_p_shot1.su fadj=data/adj/test_p_shot1.su

adj_src verbose=1 lnorm=11 fsyn=data/syn/test_p_shot1.su fobs=data/obs/test_p_shot1.su fadj=data/adj/test_p_shot1.su fmin=0.015 fmax=0.055 finc=0.005
@endcode
*/

#include "par.h"
#include "source.h"
#include "adj.h"

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
    int lnorm, inv, order, verbose;
    float alpha, fmin, fmax, finc;
    char *fobs, *fsyn, *fadj=NULL, *fmis, *ftime;
    
    /* Local variables */
    int ns, ntr, itr, it, invtime;
    int nfreq, ifreq, k_save_adj=1;
    //int maxlag=NS/2;
    float delta, Sx, Sy, Fc;
    float misfit=0.0, energy=0.0, stationary_lag;
    float **pos=NULL, **obs=NULL, **syn=NULL, **res=NULL, **freqs=NULL;
    float *tr_obs=NULL, *tr_syn=NULL, *tr_res=NULL, *misf_arr=NULL, *lag_arr=NULL;
    
    /* 1. read formatted input list */
    /* initialize getpar */
    initargs(argc,argv);
    
    lnorm=2; getparint("lnorm",&lnorm);
    inv=0; getparint("inv",&inv);
    verbose=0; getparint("verbose",&verbose);
    if(!getparstring("fsyn", &fsyn))
        throw_error("fsyn must be specified!\n");
    if(!getparstring("fobs", &fobs))
        throw_error("fobs must be specified!\n");
    if(!getparstring("fadj", &fadj)){
        fadj= "adj.su";
        k_save_adj=0;
    }
    if(!getparstring("fmis", &fmis))  fmis= "misfit.bin";
    if(lnorm==9||lnorm==10||lnorm==12||lnorm==13){
        if(!getparstring("ftime", &ftime))  ftime= "timedelay.bin";
    }else if(lnorm==11){
        order=4; getparint("order",&order);
        alpha=0.2; getparfloat("alpha",&alpha);
        if(!getparfloat("fmin",&fmin)) throw_error("fmin must be specified!\n");
        if(!getparfloat("fmax",&fmax)) throw_error("fmax must be specified!\n");
        if(!getparfloat("finc",&finc)) throw_error("finc must be specified!\n");
    }
    if(verbose){
        printf("\n******\n");
        printf("lnorm\t\t| %d\n",lnorm);
        printf("inv \t\t| %d\n",inv);
        printf("fsyn\t\t| %s\n",fsyn);
        printf("fobs\t\t| %s\n",fobs);
        if(k_save_adj==1) printf("fadj\t\t| %s\n",fadj);
        if(lnorm==9||lnorm==10||lnorm==12||lnorm==13){
            printf("ftime\t\t| %s\n",ftime);
        }else if(lnorm==11){
            printf("  %-16s| %d\n","order",order);
            printf("  %-16s| %d\n","nfreq",(int)((fmax -fmin)/finc+1));
            printf("  %-16s| %.3e\n","alpha",alpha);
            printf("  %-16s| %.3e\n","fmin",fmin);
            printf("  %-16s| %.3e\n","fmax",fmax);
            printf("  %-16s| %.3e\n","finc",finc);
        }
        printf("******\n");
    }
    
    /* special for lnorm=1 */
    if(lnorm==11){
        nfreq=(int)((fmax -fmin)/finc+1);
        freqs = matrix(1, nfreq, 1, 2);
        for(ifreq=1; ifreq<=nfreq; ifreq++){
            Fc=fmin+(ifreq-1.0)*finc;
            freqs[ifreq][1] = (1.0 -alpha)*Fc; 
            freqs[ifreq][2] = (1.0 +alpha)*Fc;
            printf("ifreq Tc Fc: %d %.2f %.5f\n", ifreq, 1.0/Fc, Fc);
        }
    }
    
    /* 2. read SU-format file */
    pos=read_su_header_all(fobs, &ntr, &ns, &delta, &Sx, &Sy);
    obs = read_su(fobs, &ntr, &ns);
    syn = read_su(fsyn, &ntr, &ns);
    res = matrix(1, ntr, 1, ns);
    tr_obs = vector(1, ns);
    tr_syn = vector(1, ns);
    tr_res = vector(1, ns);
    misf_arr = vector(1, ntr);
    lag_arr = vector(1, ns);
    
    /* 3. calculate the misfit and adjoint source */
    if(lnorm==1){
    }else if(lnorm==2){
        for(itr=1;itr<=ntr;itr++){
            /* a. initialize the array for extern subroutine */
            for(it=1;it<=ns;it++){
                tr_obs[it] = obs[itr][it];
                tr_syn[it] = syn[itr][it];
            }
            
            /* b. call extern subroutine */
            misf_arr[itr] = misfit_L2(tr_obs, tr_syn, tr_res, 0.0, ns*delta, ns, delta);
            
            /* c. obtain the residual */
            for(it=1;it<=ns;it++){
                res[itr][it]=tr_res[it];
            }
        }
    }else if(lnorm==3){
    }else if(lnorm==4){
    }else if(lnorm==5){
    }else if(lnorm==9){/* the adjoint source for calculating the misfit kernel */
        for(itr=1;itr<=ntr;itr++){
            /* a. initialize the array for extern subroutine */
            for(it=1;it<=ns;it++){
                tr_obs[it] = obs[itr][it];
                tr_syn[it] = syn[itr][it];
            }
            
            /* b. call extern subroutine */
            misf_arr[itr] = misfit_TT(tr_obs, tr_syn, tr_res, 0.0, ns*delta, ns, delta, 1, &stationary_lag);
            /* For debugging */
            lag_arr[itr]=stationary_lag;
            
            /* c. obtain the residual */
            for(it=1;it<=ns;it++){
                res[itr][it]=tr_res[it];
            }
        }
        V_write(ftime, &lag_arr[1], ntr);
        //V_write_txt(ftime, &lag_arr[1], ntr);
    }else if(lnorm==10){ /* the adjoint source for calculating the Banana-doughnut kernel */
        for(itr=1;itr<=ntr;itr++){
            /* a. initialize the array for extern subroutine */
            for(it=1;it<=ns;it++){
                tr_obs[it] = obs[itr][it];
                tr_syn[it] = syn[itr][it];
            }
            
            /* b. call extern subroutine */
            misf_arr[itr] = misfit_TT(tr_obs, tr_syn, tr_res, 0.0, ns*delta, ns, delta, 0, &stationary_lag);
            /* For debugging */
            lag_arr[itr]=stationary_lag;
            
            /* c. obtain the residual */
            for(it=1;it<=ns;it++){
                res[itr][it]=tr_res[it];
            }
        }
        V_write(ftime, &lag_arr[1], ntr);
        //V_write_txt(ftime, &lag_arr[1], ntr);
    }else if(lnorm==11){ /* the adjoint source generated by multiple traveltime misfit adjoint sources */
        for(itr=1;itr<=ntr;itr++){
            /* a. initialize the array for extern subroutine */
            for(it=1;it<=ns;it++){
                tr_obs[it] = obs[itr][it];
                tr_syn[it] = syn[itr][it];
            }
            
            /* b. call extern subroutine */
            misf_arr[itr] = misfit_TT_multi(tr_obs, tr_syn, tr_res, 0.0, ns*delta, ns, delta, freqs, nfreq, order);
            
            /* c. obtain the residual */
            for(it=1;it<=ns;it++){
                res[itr][it]=tr_res[it];
            }
        }
        
    }else if(lnorm==12){ /* the adjoint source generated based on Rytov approximation */
        for(itr=1;itr<=ntr;itr++){
            /* a. initialize the array for extern subroutine */
            for(it=1;it<=ns;it++){
                tr_obs[it] = obs[itr][it];
                tr_syn[it] = syn[itr][it];
            }
            
            /* b. call extern subroutine */
            misf_arr[itr] = misfit_Rytov(tr_obs, tr_syn, tr_res, 0.0, ns*delta, ns, delta, 0, &stationary_lag);
            /* For debugging */
            lag_arr[itr]=stationary_lag;
            
            /* c. obtain the residual */
            for(it=1;it<=ns;it++){
                res[itr][it]=tr_res[it];
            }
        }
        V_write(ftime, &lag_arr[1], ntr);
        
    }else if(lnorm==13){ /* the adjoint source for the measurement of a cross-correlation time shift */
        for(itr=1;itr<=ntr;itr++){
            /* a. initialize the array for extern subroutine */
            for(it=1;it<=ns;it++){
                tr_obs[it] = obs[itr][it];
                tr_syn[it] = syn[itr][it];
            }
            
            /* b. call extern subroutine */
            misf_arr[itr] = misfit_CC(tr_obs, tr_syn, tr_res, 0.0, ns*delta, ns, delta, &stationary_lag);
            /* For debugging */
            lag_arr[itr]=stationary_lag;
            printf("stationary_lag=%f\n",stationary_lag);
            
            /* c. obtain the residual */
            for(it=1;it<=ns;it++){
                res[itr][it]=tr_res[it];
            }
        }
        V_write(ftime, &lag_arr[1], ntr);
        //V_write_txt(ftime, &lag_arr[1], ntr);
    }
    
    /* d. sum up the misfit and energy from different traces */        
    for(itr=1;itr<=ntr;itr++){
        for(it=1;it<=ns;it++){
            energy += obs[itr][it]*obs[itr][it];
        }
        misfit += misf_arr[itr];
    }
    
    /* e. inverse the residual */
    if(inv){
        for(itr=1;itr<=ntr;itr++){
            for(it=1;it<=ns;it++){
                tr_res[it]=res[itr][it];
            }
            invtime=ns;
            for(it=1;it<=ns;it++){
                res[itr][invtime]=tr_res[it];
                invtime--;
            }
        }
    }
    
    /* 4. write adjoint source and misfit */
    if(k_save_adj==1) write_su(fadj, Sx, Sy, pos, res, ntr, ns, delta);
    V_write(fmis, &misfit, 1);
    
    /* deallocate memory */
    free_matrix(pos, 1, 3, 1, ntr);
    free_matrix(obs, 1, ntr, 1, ns);
    free_matrix(syn, 1, ntr, 1, ns);
    free_matrix(res, 1, ntr, 1, ns);
    free_vector(tr_obs, 1, ns);
    free_vector(tr_syn, 1, ns);
    free_vector(tr_res, 1, ns);
    free_vector(misf_arr, 1, ntr);
    free_vector(lag_arr, 1, ntr);
    if(lnorm==11) free_matrix(freqs, 1, nfreq, 1, 2);
    
    if(verbose) printf("misfit, energy = %e %e \n",misfit, energy);
    return 0;
}

