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

void freqdomain_filt_vector(float *data, float *f_win, int nfreq, int ns){
    /* apply frequency domain filtering */
    /* 
     * the arrays are allocated by function vector
     * 
     * data (in/out)    : one trace of seismogram, which contains one event;
     * f_win            : a series of frequencies (Hz), PLEASE make sure f_win 
     *                    is generated corresponding to the appropriate rule;
     * nfreq           : the number of f_win;
     * ns               : the length of data;
     * [dt]             : the time interval (s);
     *
     */

    /* declaration of variables */
    int i, k, npad, Nfreq;
    //float Fs;
    double npadd;
    //float *Freqs;
            
    /* declaration of variables for FFTW3*/
    fftw_complex *in, *out, *out2;
    fftw_plan p1, p2;
    
    /* 1. preprocess work */
    /***************************************************************/
    /*********** generate the frequency series *********************/
    /***************************************************************/
    npad = (int)(pow(2.0, ceil(log((double)(ns))/log(2.0))+1.0));
    npadd = (double)npad;
    //Fs=1/dt;
    Nfreq = (int)ceil((npad+1)/2.0); //index is from 1
    if(Nfreq != nfreq){
        printf("Error occurs in function freqdomain_filt: Nfreq !=nfreq,\n");
        printf("and you should set fac to 1 in function generate_freqs.\n");
        exit(1);
    }
    //Freqs = vector(1, Nfreq);
    //for(i=1;i<=Nfreq;i++) Freqs[i]=(i-1)*Fs/npad;
    /***************************************************************/
    
    /* 2. allocate memory for fftw3 */
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
    out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
    
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
    /* firstly, we apply filter on the first half part,
    * which includes the zero frequency component and Nyquist component.
    */
    for(i=0; i<=Nfreq-1; i++){
        //h = exp(-guassalfa*pow((Freqs[i+1] - fc)/fc, 2));
        out2[i][0] = out[i][0]*f_win[i+1];
        out2[i][1] = out[i][1]*f_win[i+1];
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
        data[i+1]=(1.0/npadd)*in[i][0];
    }
    
    //free_vector(Freqs, 1, Nfreq);
    fftw_free(in);
    fftw_free(out);
    fftw_free(out2);
}





