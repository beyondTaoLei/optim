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

/*----------------------------------------------------------------
 *   Calculate traveltime adjoint source                                  
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

#include "adj.h"

float misfit_CC(float *obs, float *syn, float *res, float tstart, 
                float tend, int ns, float dt, float *delay){
    /* 
     * the arrays are allocated by function vector
     * 
     * obs: observed data;
     * syn: synthetic data;
     * ns: the length of obs or syn;
     * tstart tend: the starting  and end time of measurement window;
     * dt: the time interval;
     * 
     * res: adjoint source.
     *
     */
    
    /* declaration of local variables */
    int offset, lag;
    float Nnorm, tmp;
    float *corr, *time_win, *syn_win, *obs_win, *syn_grad;//, *syn_int;
    //float *syn_env, *obs_env;
    
    int i, j, maxlag, istart, iend, nlen;
    
    float misfit;
    
    istart = max(floor(tstart/dt), 1);
    iend   = min(floor(tend/dt), ns);
    nlen = iend-istart+1;
    maxlag=nlen-1;
    /* maxlag: the maximum of time lag, limits the lag range from 
    –maxlag to maxlag. maxlag defaults to ns-1.
    One can set value from 0 to ns-1
    */
    
    corr    =vector(1, 2*maxlag+1);
    time_win=vector(1, nlen);
    syn_win =vector(1, nlen);
    obs_win =vector(1, nlen);
    syn_grad=vector(1, nlen);
    //syn_int =vector(1, nlen);
    //syn_env =vector(1, nlen);
    //obs_env =vector(1, nlen);
    
    /* 1. pre-processing */
    //time_window(&time_win[1], nlen);
    for (i=1; i<= nlen; i++){
        j=i+istart-1;
        syn_win[i] = syn[j];//*time_win[i];
        obs_win[i] = obs[j];//*time_win[i];
    }
    
    /* 2. calculate traveltime difference */
    xcorr(&syn_win[1], &obs_win[1], &corr[1], maxlag, nlen);/* change the input order */
    offset = V_max_index(1, 2*maxlag+1, corr)-1;
    lag=-maxlag + offset;
    /*
     * lag=Tsyn-Tobs, according to Fichtner's book
     */
    
    /* 3. calculate the normalization factor Nnorm */
    V_grad (1, nlen, syn_win, syn_grad);
    //V_integ(1, nlen, syn_win, syn_int);
    Nnorm = -dt*V_dotproduct(1, nlen, syn_win, syn_win);/* Nnorm<0 */

    /* 4. generate the traveltime adjoint source */
    zero(&res[1],ns);
    for (i=1; i<= nlen; i++){
        j=i+istart-1;
        res[j] = (syn_grad[i]/dt)/Nnorm;
    }
    
    tmp=lag*dt;
    *delay = tmp;
    misfit = tmp*tmp;
    
    /*
    if(lag > 0){
        printf("observed data delays by %d samples compared to synthetic data\n",lag);
        printf("Nnorm = %e, coeff = %e, misfit= %e\n", Nnorm, coeff, misfit);
    }else{
        printf("observed data arrives earlier by %d samples than synthetic data\n",-lag);
        printf("Nnorm = %e, coeff = %e, misfit= %e\n", Nnorm, coeff, misfit);
    }*/
    
    free_vector(corr, 1, 2*maxlag+1);
    free_vector(time_win,1, nlen);
    free_vector(syn_win, 1, nlen);
    free_vector(obs_win, 1, nlen);
    free_vector(syn_grad,1, nlen);
    //free_vector(syn_int, 1, nlen);
    //free_vector(syn_env, 1, nlen);
    //free_vector(obs_env, 1, nlen);
    
    return misfit;
}

