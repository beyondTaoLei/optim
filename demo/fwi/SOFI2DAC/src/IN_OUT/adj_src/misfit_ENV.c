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
/*------------------------------------------------------------------------
 *   Calculate Amplitude (in frequency-domain) adjoint source                                  
 *  ----------------------------------------------------------------------*/

/*----------------------------------------------------------------
 *   Calculate envelope adjoint source                                  
 *---------------------------------------------------------------*/
/*
@article{bozdag_misfit_2011,
	title = {Misfit functions for full waveform inversion based on instantaneous phase and envelope measurements},
	volume = {185},
	issn = {0956540X},
	url = {https://academic.oup.com/gji/article-lookup/doi/10.1111/j.1365-246X.2011.04970.x},
	doi = {10.1111/j.1365-246X.2011.04970.x},
	abstract = {},
	pages = {845--870},
	number = {2},
	journaltitle = {Geophysical Journal International},
	author = {Bozdağ, Ebru and Trampert, Jeannot and Tromp, Jeroen},
	urldate = {2019-01-04},
	date = {2011-05},
	langid = {english},
	keywords = {misfit, adjoint kernels},
	file = {}
}
*/

#include "comm.h"

float misfit_ENV(float *obs, float *syn, float *res, float tstart, 
                 float tend, int ns, float dt, float eps){
    /* 
     * the arrays are allocated by function vector
     * 
     * obs: observed data;
     * syn: synthetic data;
     * ns: the length of obs or syn;
     * tstart tend: the starting  and end time of measurement window;
     * dt: the time interval;
     * eps: stabilization factor, e.g. 0.1;
     * 
     * res: adjoint source.
     *
     */

    /* declaration of variables */
    int i, j, istart, iend, nlen, kof;
    float max_env, waterlevel, env2, coeff, misfit, tmp, syn_env_tmp;
    float *time_win, *syn_win, *obs_win;
    float *obs_env, *syn_hilb, *syn_env, *res_p1, *res_p2, *ln_env;
            
    istart = max(floor(tstart/dt), 1);
    iend   = min(floor(tend/dt), ns);
    printf("istart, iend= %d %d\n", istart,iend);
    nlen = iend-istart+1;
    
    time_win = vector(1, nlen);
    syn_win  = vector(1, nlen);
    obs_win  = vector(1, nlen);
    obs_env  = vector(1, nlen);
    syn_hilb = vector(1, nlen);
    syn_env  = vector(1, nlen);
    res_p1   = vector(1, nlen);
    res_p2   = vector(1, nlen);
    ln_env   = vector(1, nlen);
    
    /* 1. pre-processing */
    //time_window(&time_win[1], nlen);
    for (i=1; i<= nlen; i++){
        j=i+istart-1;
        syn_win[i] = syn[j];//*time_win[i];
        obs_win[i] = obs[j];//*time_win[i];
    }
    
    /*1. calculation of the Hilbert transform */
    envelope(&obs_win[1], &obs_env[1], nlen);
    envelope_hilbert(&syn_win[1], &syn_hilb[1], &syn_env[1], nlen);
    max_env = V_max_abs(1, nlen, syn_env);
    waterlevel=max_env*eps;
    printf("max_env=%e waterlevel=%e\n",max_env, waterlevel);
    V_write("./obs_env.bin",&obs_env[1],nlen);
    V_write("./syn_env.bin",&syn_env[1],nlen);
    V_write("./syn_hilb.bin",&syn_hilb[1], nlen);
    
    /*2. calculation of adjoint source and misfit */
    misfit=0.0;
    for (i=1;i<=nlen;i++){
        if(syn_env[i]>=waterlevel){
            ln_env[i] = log(obs_env[i]/syn_env[i]);
            env2=syn_env[i]*syn_env[i];
            res_p1[i] = ln_env[i]*syn_win[i]/env2;
            res_p2[i] = ln_env[i]*syn_hilb[i]/env2;
            kof=0;
        }else{
            coeff=pow(sin(0.5*PI*syn_env[i]/waterlevel),4);
            tmp = (1.0-coeff)*waterlevel;
            syn_env_tmp = syn_env[i]+tmp;
            ln_env[i] = coeff*log((obs_env[i]+tmp)/syn_env_tmp);
            
            env2=syn_env_tmp*syn_env_tmp;
            res_p1[i] = coeff*ln_env[i]*syn_win[i]/env2;
            res_p2[i] = coeff*ln_env[i]*syn_hilb[i]/env2;
            kof=1;
        }
        misfit += ln_env[i]*ln_env[i];
        
//         env2=syn_env[i]*syn_env[i]+waterlevel*waterlevel;
//         res_p1[i] = ln_env[i]*syn_win[i]/env2;
//         res_p2[i] = ln_env[i]*syn_hilb[i]/env2;
        printf("%d %d %e %e %e\n",i, kof, obs_env[i],syn_env[i],ln_env[i]);
    }
    misfit = 0.5*misfit;
    V_write("./ln_env.bin",&ln_env[1],nlen);
    V_write("./res_p1.bin",&res_p1[1],nlen);
    V_write("./res_p2.bin",&res_p2[1],nlen);
    
    /*3. calculation of the second part of adjoint source*/
    hilbert(&res_p2[1], &res_p2[1], nlen);
    V_write("./res_p2_2.bin",&res_p2[1],nlen);
    for (i=1;i<=nlen;i++){
        res_p2[i] = -res_p1[i] + res_p2[i];
    }
    for (i=1; i<= ns; i++){
        res[i] = 0.0;
    }
    for (i=1; i<= nlen; i++){
        j=i+istart-1;
        res[j] = res_p2[i];
    }
    
    free(time_win);
    free(syn_win);
    free(obs_win);
    free(obs_env);
    free(syn_hilb);
    free(syn_env);
    free(res_p1);
    free(res_p2);
    free(ln_env);
    
    return misfit;

}






