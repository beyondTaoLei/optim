/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2016  For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS.
 * 
 * IFOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * IFOS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IFOS. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/

/* the routine is modified based on sac software */
   /*     Call xapiir ( Apply a IIR Filter ) 
     *        - yarray - Original Data 
     *        - nlen   - Number of points in yarray 
     *        - proto  - Prototype of Filter 
     *                 - SAC_FILTER_BUTTERWORTH       - Butterworth         "BU"
     *                 - SAC_FILTER_BESSEL            - Bessel              "BE"
     *                 - SAC_FILTER_CHEBYSHEV_TYPE_I  - Chebyshev Type I    "C1"
     *                 - SAC_FILTER_CHEBYSHEV_TYPE_II - Chebyshev Type II   "C2"
     *        - transition_bandwidth (Only for Chebyshev Filter) 
     *                 - Bandwidth as a fraction of the lowpass prototype 
     *                   cutoff frequency 
     *        - attenuation (Only for Chebyshev Filter) 
     *                 - Attenuation factor, equals amplitude reached at 
     *                   stopband egde 
     *        - order  - Number of poles or order of the analog prototype 
     *                   4 - 5 should be ample 
     *                   Cannot exceed 10 
     *        - type   - Type of Filter 
     *                 - SAC_FILTER_BANDPASS    "BP"
     *                 - SAC_FILTER_BANDREJECT  "BR"
     *                 - SAC_FILTER_LOWPASS     "LP"
     *                 - SAC_FILTER_HIGHPASS    "HP"
     *        - low    - Low Frequency Cutoff [ Hertz ] 
     *                   Ignored on SAC_FILTER_LOWPASS 
     *        - high   - High Frequency Cutoff [ Hertz ] 
     *                   Ignored on SAC_FILTER_HIGHPASS 
     *        - delta  - Sampling Interval [ seconds ] 
     *        - passes - Number of passes 
     *                 - 1 Forward filter only 
     *                 - 2 Forward and reverse (i.e. zero-phase) filtering 
     */                                 
#include "comm.h"
#include "sac.h"

//void Taper(float *data, int type, int nlen, float width);

void filtering(float ** data, int ntr, int ns, float dt, float fc1, 
               float fc2, int order, int zero_phase, int method, int type){

    /* 
    data    :   2-dimensional array containing seismograms (
    ntr     :   number of traces
    ns      :   number of samples
    method  :   definition of filter
                    1: "BU" + "BP"
                    2: "BE" + "BP"
                    3: 
    */

    /* declaration of local variables */
    int itr, j, passes;
    double delta_d, low, high;
    char prototype[3],passtype[3];
    float *seismogram;
    
    if(ntr &&method){
        if(method==1){ /* butterworth */
            sprintf(prototype,"BU");
        }else if(method==2){/* bessel */
            sprintf(prototype,"BE");
        }else if(method==3){/* Chebyshev Type I */
            sprintf(prototype,"C1");
        }else if(method==4){/* Chebyshev Type II */
            sprintf(prototype,"C2");
        }else{
            printf("So far method should be [1 4] in sac_filt.c! \n");
            exit(1);
        }
        
        if(type==1){ /* BANDPASS */
            sprintf(passtype,"BP");
        }else if(type==2){/* BANDREJECT */
            sprintf(passtype,"BR");
        }else if(type==3){/* LOWPASS */
            sprintf(passtype,"LP");
        }else if(type==4){/* HIGHPASS */
            sprintf(passtype,"HP");
        }else{
            printf("So far type should be [1 4] in sac_filt.c! \n");
            exit(1);
        }
        
        /* 1. initializing */
        delta_d = dt;
        low = fc1;
        high = fc2;
        if(zero_phase ==1){
            passes =2;
        }else{
            passes =1;
        }
        seismogram = vector(1, ns);
        
        /* 2. apply filtering */
        for (itr=1;itr<=ntr;itr++){
            for (j=1;j<=ns;j++){
                seismogram[j]=data[itr][j];
            }
            xapiir(&seismogram[1], ns, prototype, 0.0, 0.0, order, passtype,
                low, high, delta_d, passes);
            
            /* put zeros at each end*/
            Taper(&seismogram[1], 3, ns, 0.15);
            for (j=1;j<=ns;j++){
                data[itr][j]=seismogram[j];
            }
        }
        
        free_vector(seismogram,1, ns);
    
    }
    
} /* end of function */
