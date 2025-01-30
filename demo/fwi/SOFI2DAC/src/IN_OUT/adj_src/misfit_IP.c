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

/* files to include */
#include "comm.h"

// #include <stdio.h>
// #include <math.h>
// #include <stdlib.h>
// #include <stddef.h>
// #include <string.h>
// #include <time.h>
// #include <fftw3.h>

//icc -O3 -c -o misfit_IP.o misfit_IP.c -I/home/tao/seis_software/seisplatform/FD/SOFI2D/src/IN_OUT/include/
float misfit_IP(float *obs, float *syn, float *res, int ns, float dt, float eps){
    /* INPUT */
    /*     
    * obs: the observed data  (velocity component)
    * syn: the synthetic data (velocity component)
    * ns: the length of the data
    * dt: the time sample (s)
    */
    /* OUTPUT */
    /*
     * res: adjoint source (velocity component)
     *    
     */

    /* declaration of variables */
    int i;
    float max_env, waterlevel, env2, misfit;
    float *obs_hilb, *obs_env, *obs_ph, *syn_hilb, *syn_env, *syn_ph;
    float *diff_ph, *res_p1, *res_p2;
            
    obs_hilb =(float *)malloc((size_t) (ns*sizeof(float)));
    obs_env  =(float *)malloc((size_t) (ns*sizeof(float)));
    obs_ph   =(float *)malloc((size_t) (ns*sizeof(float)));
    syn_hilb =(float *)malloc((size_t) (ns*sizeof(float)));
    syn_env  =(float *)malloc((size_t) (ns*sizeof(float)));
    syn_ph   =(float *)malloc((size_t) (ns*sizeof(float)));
    diff_ph  =(float *)malloc((size_t) (ns*sizeof(float)));
    
    res_p1  =(float *)malloc((size_t) (ns*sizeof(float)));
    res_p2  =(float *)malloc((size_t) (ns*sizeof(float)));
    
    /*1. run the Hilbert transform */
    envelope_hilbert(obs, obs_hilb, obs_env, ns);
    envelope_hilbert(syn, syn_hilb, syn_env, ns);
    max_env = V_max_abs(0, ns-1, syn_env);
    waterlevel=max_env*eps;
    
    /*2. calculate the misfit */
//     phase(syn, syn_hilb, syn_ph, ns);
//     V_write("/home/tao/seis_software/test/misfit/Bozdag_Fig5/fig/syn_hilb.bin",syn_hilb,ns);
//     V_write("/home/tao/seis_software/test/misfit/Bozdag_Fig5/fig/syn_env.bin",syn_env,ns);
//     phase(obs, obs_hilb, obs_ph, ns);
//     V_write("/home/tao/seis_software/test/misfit/Bozdag_Fig5/fig/obs_hilb.bin",obs_hilb,ns);
//     V_write("/home/tao/seis_software/test/misfit/Bozdag_Fig5/fig/obs_env.bin",obs_env,ns);
//     for (i=0;i<ns;i++){
//         diff_ph[i] = syn_ph[i]-obs_ph[i];
//         if(diff_ph[i] > 0.5*PI){
//             diff_ph[i]=diff_ph[i]-PI;
//         }else if(diff_ph[i] < -0.5*PI){
//             diff_ph[i]=diff_ph[i]+PI;
//         }
//         misfit += diff_ph[i]*diff_ph[i];
//     }
//     misfit = 0.5*misfit;
    
    
    phase_diff(obs, obs_hilb, syn, syn_hilb, diff_ph, ns);
    misfit=0.0;
    for (i=0;i<ns;i++){
        misfit += diff_ph[i]*diff_ph[i];
    }
    misfit = 0.5*misfit;
    
//     V_write("/home/tao/seis_software/test/misfit/IP_test2/ph_syn.bin",syn_ph,ns);
//     V_write("/home/tao/seis_software/test/misfit/IP_test2/ph_obs.bin",obs_ph,ns);
//     V_write("/home/tao/seis_software/test/misfit/IP_test2/ph_diff.bin",diff_ph,ns);
//     V_write("/home/tao/seis_software/test/misfit/Bozdag_Fig5/fig/ph_syn.bin",syn_ph,ns);
//     V_write("/home/tao/seis_software/test/misfit/Bozdag_Fig5/fig/ph_obs.bin",obs_ph,ns);
//     V_write("/home/tao/seis_software/test/misfit/Bozdag_Fig5/fig/ph_diff2.bin",diff_ph,ns);
    
    /*3. calculate of adjoint source */
    for (i=0;i<ns;i++){
        env2=syn_env[i]*syn_env[i]+waterlevel*waterlevel;
        res_p1[i] = diff_ph[i]*syn_hilb[i]/env2;
        res_p2[i] = diff_ph[i]*syn[i]/env2;
    }
    /* calculate the second part of adjoint source */
    hilbert(res_p2, res_p2, ns);
    for (i=0;i<ns;i++) res[i] = -res_p1[i] - res_p2[i];
    
    free(obs_hilb);
    free(obs_env);
    free(obs_ph);
    free(syn_hilb);
    free(syn_env);
    free(syn_ph);
    free(diff_ph);
    
    free(res_p1);
    free(res_p2);
    
    return misfit;

}






