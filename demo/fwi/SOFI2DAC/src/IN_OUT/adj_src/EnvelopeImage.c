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
 *   calculate envelope image, i.e., to obtain envelope at each T
 *   new code for group velocity analysis using frequency domain 
 *   Gaussian filter 
 *   copy from Yao Huajian                             
 *---------------------------------------------------------------*/
/*
@article{

}
*/

#include "comm.h"
#include <fftw3.h>

void EnvelopeImage(float *data, float *TPoint, int NumCtrT, int ns, 
                    float dt, float StaDist, float **env){
    /* 
     * the arrays are allocated by function vector
     * 
     * data     : one trace of seismogram, which contains one event;
     * TPoint   : the periods of interesting (s);
     * NumCtrT  : the number of TPoint;
     * ns       : the length of data;
     * dt       : the time interval (s);
     * StaDist  : the station offset from the source (kilometers); 
     * 
     * env      : envelope image.
     *
     */

    /* declaration of variables */
    int i, j, k, npad, Nfreq;
    float Fs, guassalfa=0.0, fc, CtrT;
    double npadd;
    float *Freqs, *win, *data2;
    float alfa[2][8] = {
                        {0, 100, 250, 500, 1000, 2000, 4000, 20000},
                        {5, 8,   12,  20,   25,  35,   50,   75}
                        };
            
    /* declaration of variables for FFTW3*/
    fftw_complex *in, *out, *out2;
    fftw_plan p1, p2;
    
    /* 1. preprocess work */
    npad = (int)(pow(2.0, ceil(log((double)(ns))/log(2.0))+1.0));
    npadd = (double)npad;
    Fs=1/dt;
    Nfreq = (int)ceil((npad+1)/2.0);
    
    /* 2. allocate memory */
    Freqs = vector(1, Nfreq);
    win   = vector(1, Nfreq);
    data2 = vector(1, ns);
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
    out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
    
    guassalfa = linint_set(alfa[0], alfa[1], StaDist, 8);
    printf("guassalfa=%f\n",guassalfa);
    if(guassalfa<=0.0){
        printf("Error occurs in function EnvelopeImage: guassalfa <= 0.0\n");
        exit(1);
    }
    for(i=1;i<=Nfreq;i++) Freqs[i]=(i-1)*Fs/npad;
    
    /* 3. initialize the work arrays */
    for (i=0;i<ns;i++){
        in[i][0]=data[i+1];
        in[i][1]=0.0;
    }
    for (i=ns;i<npad;i++){
        in[i][0]=0.0;
        in[i][1]=0.0;
    }
    /* 4. FFT */
    p1 = fftw_plan_dft_1d(npad, in, out, 1, FFTW_ESTIMATE);
    fftw_execute(p1);
    fftw_destroy_plan(p1);
    
    /* 5. apply filtering */
    for(j=1;j<=NumCtrT;j++){
        CtrT = TPoint[j];
        fc = 1/CtrT;
        /* firstly, we apply filter on the first half part,
         * which includes the zero frequency component and Nyquist component.
         */
        gauss_window(Nfreq, Freqs, win, guassalfa, fc);
        for(i=0; i<=Nfreq-1; i++){
//             h = exp(-guassalfa*pow((Freqs[i+1] - fc)/fc, 2));
            out2[i][0] = out[i][0]*win[i+1];
            out2[i][1] = out[i][1]*win[i+1];
        }
        
        /* secondly, we apply conjugate operator on the second half part. */
        k=1;
        for(i=npad-1; i>=Nfreq; i--){
            out2[i][0] =  out2[k][0];
            out2[i][1] = -out2[k][1];
            k++;
        }
        /* inverse FFT */
        p2 = fftw_plan_dft_1d(npad, out2, in, -1, FFTW_ESTIMATE);
        fftw_execute(p2);
        fftw_destroy_plan(p2);
        
        /* calculation of filtered data */
        for (i=0;i<ns;i++){
            data2[i+1]=(1.0/npadd)*in[i][0];
        }
        
        envelope(&data2[1], &env[j][1], ns);
    }
    
    
    free_vector(Freqs, 1, Nfreq);
    free_vector(data2, 1, ns);
    fftw_free(in);
    fftw_free(out);
    fftw_free(out2);
}





