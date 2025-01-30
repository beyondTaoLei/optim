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
#include "comm.h"
#include <fftw3.h>
/*
 icc -c fft_use.c  -I ../include/
*/
void fft_amp_phase(float *signals, float *phase, int ns){
    
    /* the index of the vectors is between [0 ns-1] 
     * signals(INOUT): input signals in time series and return the amplitude in frequency series
     * phase(OUT):     the phase in frequency series
     */
    
    /* declaration of variables */
    int i, npad;
    //double npadd;
            
    /* declaration of variables for FFTW3*/
    fftw_complex *in, *out;
    fftw_plan p1;

    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    /* calculation of the Hilbert transform */
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */


    /* calculation of the FFT */
//     npad = (int)(pow(2.0, ceil(log((double)(ns))/log(2.0))+2.0) );  /* ns -> npad for usage in FFT*/
    npad=ns;
    //npadd = (double)npad;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
    
                                    
    for (i=0;i<ns;i++){
        in[i][0]=signals[i];
        in[i][1]=0.0;
    }
    for (i=ns;i<npad;i++){
        in[i][0]=0.0;
        in[i][1]=0.0;
    }
    /* FFT */
    p1 = fftw_plan_dft_1d(npad, in, out, -1, FFTW_ESTIMATE);
    fftw_execute(p1);
    fftw_destroy_plan(p1);		    

    /* elementwise multiplication of FFT and h */
    for (i=0;i<npad;i++){
        signals[i] = sqrt(out[i][0]*out[i][0]+out[i][1]*out[i][1])/(ns/2.0);
        /* the image part is negative of image part of matlab output,
         * and the later is right!!!  
         */
        phase[i] = atan2(out[i][1],out[i][0]);
    }
    signals[0] = signals[0]/2.0;
    
    fftw_free(in);
    fftw_free(out);
	
}

void fft_real_imag(float *signals, float *imag, int ns, int flag){
    
    /* the index of the vectors is between [0 ns-1] 
     * signals(INOUT): input signals in time series and 
     *                 return the real part in frequency series
     * imag(OUT):      the imaginary in frequency series
     flag: -1 forward transform; 
            1 backward DFT
     */
    
    /* declaration of variables */
    int i, npad;
    //double npadd;
            
    /* declaration of variables for FFTW3*/
    fftw_complex *in, *out;
    fftw_plan p1;

    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    /*          run fft transform           */
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */


    /* calculation of the FFT */
//     npad = (int)(pow(2.0, ceil(log((double)(ns))/log(2.0))+2.0) );  /* ns -> npad for usage in FFT*/
    npad=ns;
    //npadd = (double)npad;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
    
                                    
    for (i=0;i<ns;i++){
        in[i][0]=signals[i];
        in[i][1]=0.0;
    }
    for (i=ns;i<npad;i++){
        in[i][0]=0.0;
        in[i][1]=0.0;
    }
    /* FFT */
    p1 = fftw_plan_dft_1d(npad, in, out, flag, FFTW_ESTIMATE);
    fftw_execute(p1);
    fftw_destroy_plan(p1);		    

    /* elementwise multiplication of FFT and h */
    for (i=0;i<npad;i++){
        signals[i] = out[i][0];
        imag[i]    = out[i][1];/* same with matlab output */
    }
    
    fftw_free(in);
    fftw_free(out);
	
}


