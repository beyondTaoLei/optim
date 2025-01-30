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
//modified from matlab codes
// function c = crosscorr(x,y,maxlag)


/* files to include */
// #include "comm.h"
// #include <stdio.h>
// #include <math.h>
// #include <stdlib.h>
// #include <stddef.h>
// #include <string.h>
// #include <time.h>
#include "comm.h"
#include <fftw3.h>
/*
 icc -c envelope.c  -I ../include/
*/
void hilbert(float *signals, float *hilb, int ns){
    /* signals and hilb can use same address if signals are not used again */
    /* declaration of variables */
    int i, npad;
    float  *h;
    double npadd;
            
    /* declaration of variables for FFTW3*/
    fftw_complex *in, *out;
    fftw_plan p1, p2;

    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    /* calculation of the Hilbert transform */
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */


    /* calculation of the FFT */
//     npad = (int)(pow(2.0, ceil(log((double)(ns))/log(2.0))+2.0) );  /* ns -> npad for usage in FFT*/
    npad=ns;
    npadd = (double)npad;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);

    /* calculation of the vector h */
    h =(float *)malloc((size_t) (npad*sizeof(float)));
    for(i=0;i<npad;i++) h[i]=0.0;
    
    if ((2*(npad/2))==npad) {
        /* npad is even */
        h[0]=1.0;
        h[npad/2]=1.0;
        for (i=1;i<=(npad/2-1);i++) h[i] = 2.0;
    } else {
        /* npad is odd */
        h[0]=1.0;
        for (i=1;i<=((npad+1)/2-1);i++) h[i]=2.0;
    }
                                    
    for (i=0;i<ns;i++){
        in[i][0]=signals[i];
        in[i][1]=0.0;
    }
    for (i=ns;i<npad;i++){
        in[i][0]=0.0;
        in[i][1]=0.0;
    }
    /* FFT */
    p1 = fftw_plan_dft_1d(npad, in, out, 1, FFTW_ESTIMATE);
    fftw_execute(p1);
    fftw_destroy_plan(p1);		    

    /* elementwise multiplication of FFT and h */
    for (i=0;i<npad;i++){
        out[i][0] *= h[i];
        out[i][1] *= h[i];
    }


    /* inverse FFT */
    p2 = fftw_plan_dft_1d(npad, out, in, -1, FFTW_ESTIMATE);
    fftw_execute(p2);
    fftw_destroy_plan(p2);

    /* %%%%%%%%%%%%%%%%%%%%%%% */
    /* calculation of envelope */
    /* %%%%%%%%%%%%%%%%%%%%%%% */
    for (i=0;i<ns;i++){
        hilb[i]=(1.0/npadd)*in[i][1]*(-1);
    }

    fftw_free(in);
    fftw_free(out);
    free(h);
	
}

void envelope(float *signals, float *env, int ns){
    /* signals and env can use same address if signals are not used again */
    int i;
    float *hilb;
    
    hilb    =(float *)malloc((size_t) (ns*sizeof(float)));
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    /* calculation of the Hilbert transform */
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    hilbert(signals, hilb, ns);
    
    /* %%%%%%%%%%%%%%%%%%%%%%% */
    /* calculation of envelope */
    /* %%%%%%%%%%%%%%%%%%%%%%% */
    for (i=0;i<ns;i++){
        env[i]=sqrt(signals[i]*signals[i]+hilb[i]*hilb[i]);
    }
    
    free(hilb);
}

void envelope_hilbert(float *signals, float *hilb, float *env, int ns){
    int i;
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    /* calculation of the Hilbert transform */
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    hilbert(signals, hilb, ns);
    
    /* %%%%%%%%%%%%%%%%%%%%%%% */
    /* calculation of envelope */
    /* %%%%%%%%%%%%%%%%%%%%%%% */
    for (i=0;i<ns;i++){
        env[i]=sqrt(signals[i]*signals[i]+hilb[i]*hilb[i]);
    }
    
}

void phase(float *signals, float *hilb, float *phase, int ns){
    
    int i;
    
    for (i=0;i<ns;i++){
        phase[i]=atan(hilb[i]/signals[i]);
        if(isnan(phase[i])) phase[i]=0.0;
    }
}

void phase2(float *signals, float *hilb, float *phase, int ns){
    
    int i;
    
    for (i=0;i<ns;i++){
        phase[i]=atan2(hilb[i],signals[i]);
    }
}

void phase_diff(float *signals_o, float *hilb_o, float *signals_s, 
                float *hilb_s, float *phase_diff, int ns){
    
    int i;
    float env_obs, env_syn, num, denom;
    
    for (i=0;i<ns;i++){
        env_obs=sqrt(signals_o[i]*signals_o[i]+hilb_o[i]*hilb_o[i]);
        env_syn=sqrt(signals_s[i]*signals_s[i]+hilb_s[i]*hilb_s[i]);
        
        num=signals_o[i]*hilb_s[i]-signals_s[i]*hilb_o[i];
        denom=env_obs*env_syn;
        phase_diff[i]=asin(num/denom);
        if(isnan(phase_diff[i])) phase_diff[i]=0.0;
    }
}

// void hilbert2(float *signals, float *hilb, int ns){
// 
//     /* declaration of variables */
//     int i, npad;
//     float  *h;
//     double npadd;
//             
//     /* declaration of variables for FFTW3*/
//     fftw_complex *in, *out;
//     fftw_plan p1, p2;
// 
//     /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
//     /* calculation of the Hilbert transform */
//     /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
// 
// 
//     /* calculation of the FFT */
// //     npad = (int)(pow(2.0, ceil(log((double)(ns))/log(2.0))+2.0) );  /* ns -> npad for usage in FFT*/
//     npad=ns;
//     npadd = (double)npad;
//     in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
//     out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
// 
//     /* calculation of the vector h */
//     h =(float *)malloc((size_t) (npad*sizeof(float)));
//     for(i=0;i<npad;i++) h[i]=0.0;
//     
//     if ((2*(npad/2))==npad) {
//         /* npad is even */
//         h[0]=1.0;
//         h[npad/2]=1.0;
//         for (i=1;i<=(npad/2-1);i++) h[i] = 2.0;
//     } else {
//         /* npad is odd */
//         h[0]=1.0;
//         for (i=1;i<=((npad+1)/2-1);i++) h[i]=2.0;
//     }
//                                     
//     for (i=0;i<ns;i++){
//         in[i][0]=signals[i];
//         in[i][1]=0.0;
//     }
//     for (i=ns;i<npad;i++){
//         in[i][0]=0.0;
//         in[i][1]=0.0;
//     }
//     /* FFT */
//     p1 = fftw_plan_dft_1d(npad, in, out, 1, FFTW_ESTIMATE);
//     fftw_execute(p1);
//     fftw_destroy_plan(p1);		    
// 
//     /* elementwise multiplication of FFT and h */
//     for (i=0;i<npad;i++){
//         out[i][0] *= h[i];
//         out[i][1] *= h[i];
//     }
// 
// 
//     /* inverse FFT */
//     p2 = fftw_plan_dft_1d(npad, out, in, -1, FFTW_ESTIMATE);
//     fftw_execute(p2);
//     fftw_destroy_plan(p2);
// 
//     /* %%%%%%%%%%%%%%%%%%%%%%% */
//     /* calculation of envelope */
//     /* %%%%%%%%%%%%%%%%%%%%%%% */
//     for (i=0;i<ns;i++){
//         hilb[i]=(1.0/npadd)*in[i][0];
//     }
// 
//     fftw_free(in);
//     fftw_free(out);
//     free(h);
// 	
// }