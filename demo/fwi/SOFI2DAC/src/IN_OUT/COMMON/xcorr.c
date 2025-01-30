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
#include "comm.h"
#include <fftw3.h>
/*
 icc -c xcorr.c
*/
void xcorr(float * temp_TS, float * temp_TS1, float * temp_xcorr, int maxlag, int ns){
/* cross correlation */

    /* INPUT */
    /*     
    * temp_TS: the first time series
    * temp_TS1: the second time series
    * ns: the length of temp_TS or temp_TS1
    * maxlag: the maximum of time lag, limits the lag range from 
    *         –maxlag to maxlag. maxlag defaults to ns-1.
    *         One can set value from 0 to ns-1
    */
    
    /* OUTPUT */
    /*
     * temp_xcorr: cross correlation between temp_TS and temp_TS1, 
     *              whose length is 2*maxlag+1 
     * the correspecting time lag series is [-maxlag, ..., 0, ..., maxlag]
     *    
     */
    
    /* declaration of local variables */
    int i,j, npad;
    double npadd;
            
    /* declaration of variables for FFTW3*/
    fftw_complex *in1, *in2, *out1, *out2, *out3;
    fftw_plan p1, p2, p3;
    
    if(maxlag>ns-1){
        printf(" Variable maxlag, %d, should be less than ns, %d\n",maxlag,ns);
        exit(0);
    }
    
    npad = (int)(pow(2.0, ceil(log((double)(ns))/log(2.0))+2.0) );  /* ns -> npad for usage in FFT*/
    npadd = (double)npad;
//     (int)ceil(4.0)=4.0; (int)ceil(4.1)=5.0
    
    in1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
    out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
    
    in2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
    out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
    
    out3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
    
    for (j=0;j<ns;j++){
        in1[j][0]=temp_TS[j];
        in1[j][1]=0.0;
    }
    for (j=ns;j<npad;j++){
        in1[j][0]=0.0;
        in1[j][1]=0.0;
    }
    
    for (j=0;j<ns;j++){
        in2[j][0]=temp_TS1[j];
        in2[j][1]=0.0;
    }
    for (j=ns;j<npad;j++){
        in2[j][0]=0.0;
        in2[j][1]=0.0;
    }
    
    
    p1 = fftw_plan_dft_1d(npad, in1, out1, 1, FFTW_ESTIMATE);
    fftw_execute(p1);
    
    p2 = fftw_plan_dft_1d(npad, in2, out2, 1, FFTW_ESTIMATE);
    fftw_execute(p2);
            
    /* multiplication of complex vectors */
    for(i=0;i<npad;i++){
        out3[i][0]    = (out1[i][0] * out2[i][0]) + (out1[i][1] * out2[i][1]);
        out3[i][1] = -(out1[i][0] * out2[i][1]) + (out1[i][1] * out2[i][0]);
    }
    
    p3 = fftw_plan_dft_1d(npad, out3, in1, -1, FFTW_ESTIMATE);
    fftw_execute(p3);
            
    /* write output into data matrix */
    i=-1;
    for(j=npad-maxlag;j<npad;j++){
        i=i+1;
        temp_xcorr[i] = (1.0/npadd)*in1[j][0];
    }
    
    for(j=0;j<maxlag+1;j++){
        i=i+1;
        temp_xcorr[i] = (1.0/npadd)*in1[j][0];
    }
    
    fftw_free(in1);
    fftw_free(out1);
    
    fftw_free(in2);
    fftw_free(out2);
    
    fftw_free(out3);
    
    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_destroy_plan(p3);
}


// int max_index(float *temp_xcorr, int n){
// /* finds the index in which the element is maximum of whole array. */
//     int i, index=0;
//     float tmp, max1=0.0;
//     
//     for(i=0;i<n;i++){
//         tmp = temp_xcorr[i];
//         if(max1<=tmp){
//             max1=tmp;
//             index=i;
//         }
//     }
//     return index;
// }
// 
// int max_abs_index(float *temp_xcorr, int n){
// /* finds the index in which the element is maximum of whole array. */
//     int i, index=0;
//     float tmp, max1=0.0;
//     
//     for(i=0;i<n;i++){
//         tmp = fabs(temp_xcorr[i]);
//         if(max1<=tmp){
//             max1=tmp;
//             index=i;
//         }
//     }
//     return index;
// }

// gcc -c xcorr_test.c && gcc -o xcorr_test xcorr_test.o -lfftw3 -lm && ./xcorr_test
// int main(int argc, char **argv){
// 
//     /* declaration of local variables */
//     int i,j, npad,lag;
//     double npadd;
//     int ns=14;
//     int maxlag=3;//ns-1;
//     float temp_TS[]= {0,0,0,0,0,0,1,2,3,4,5,0,0,0};
//     float temp_TS1[]={0,0,0,1,2,3,4,5,0,0,0,0,0,0};
//     float temp_xcorr[2*maxlag+1];
//     
//     //lag: temp_TS lags temp_TS1 by lag samples
//     
//     xcorr(temp_TS, temp_TS1, temp_xcorr, maxlag, ns);
//     
//     printf("maxlag = %d \n",maxlag);
//     for(j=0;j<2*maxlag+1;j++){
//         printf("%d %.2f\n",j, temp_xcorr[j]);
//     }
//     i = max_index(temp_xcorr, 2*maxlag+1);
//     lag=-maxlag+i;
//     printf("lag = %d \n",lag);
//     
//     return 0;
// }

