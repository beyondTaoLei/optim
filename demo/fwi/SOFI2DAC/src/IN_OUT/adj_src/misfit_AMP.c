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
 *   Calculate amplitude adjoint source                                  
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

float misfit_AMP(float *obs, float *syn, float *res, float tstart, 
                        float tend, int ns, float dt){
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
    float coeff, Mnorm, dlnAd, amp_syn2, amp_obs2;
    float *time_win, *syn_win, *obs_win, *syn_int, *obs_int;
    int i, j, istart, iend, nlen;
    float misfit;
    
    istart = max(floor(tstart/dt), 1);
    iend   = min(floor(tend/dt), ns);
    nlen = iend-istart+1;
    
    time_win=vector(1, nlen);
    syn_win =vector(1, nlen);
    obs_win =vector(1, nlen);
    syn_int =vector(1, nlen);
    obs_int =vector(1, nlen);
    
    /* 1. pre-processing */
    //time_window(&time_win[1], nlen);
    for (i=1; i<= nlen; i++){
        j=i+istart-1;
        syn_win[i] = syn[j];//*time_win[i];
        obs_win[i] = obs[j];//*time_win[i];
    }
    
    /* 2. calculate the normalization factor Mnorm */
    V_integ(1, nlen, syn_win, syn_int);
    V_integ(1, nlen, obs_win, obs_int);
    Mnorm = V_dotproduct(1, nlen, syn_int, syn_int);
    
    /* 3. calculate ln(Asyn/Aobs) */
    amp_syn2 = 0.0;
    amp_obs2 = 0.0;
    for (i=1; i<= nlen; i++){
        amp_syn2 += syn_int[i]*syn_int[i];
        amp_obs2 += obs_int[i]*obs_int[i];
    }
    dlnAd = 0.5*log(amp_syn2/amp_obs2);
    
    /* 4. generate the amplitude adjoint source */
    /* the adjoint source for calculating the Bonana-doughnut kernel */
    // coeff = 1.0/Mnorm;
    /* the adjoint source for calculating the misfit kernel */
    coeff = dlnAd/Mnorm;
    for (i=1; i<= ns; i++){
        res[i] = 0.0;
    }
    for (i=1; i<= nlen; i++){
        j=i+istart-1;
        res[j] = coeff*syn_win[i];
    }
    
    misfit = dlnAd*dlnAd;
    
//     printf("Mnorm = %e, coeff = %e, misfit= %e, dlnAd = %e\n", Mnorm, coeff, misfit, dlnAd);
//     printf("amp_syn2, amp_obs2, rate =%e, %e, %e\n",amp_syn2, amp_obs2, amp_syn2/amp_obs2);
    
    free_vector(time_win,1, nlen);
    free_vector(syn_win, 1, nlen);
    free_vector(obs_win, 1, nlen);
    free_vector(syn_int, 1, nlen);
    free_vector(obs_int,1, nlen);
    
    return misfit;
}

