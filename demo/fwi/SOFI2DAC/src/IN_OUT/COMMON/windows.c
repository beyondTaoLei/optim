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

/*
 icc -c envelope.c  -I ../include/
*/


// void time_window(float *time_win, int nlen){
//     /* extracted from wave2d_sub4.f90 in specfem3d software */
//     
//     int i, ipwr_t = 10;
//     float tmp;
//     
//     //sfac1 = (2./dble(nlen))**2;
//     for (i=0;i<nlen;i++){
//         //tmp = 1.                                   // boxcar window
//         //tmp = 1 - sfac1*(i - dble(nlen)/2.)**2     // welch window
//         tmp = 1. - pow(cos(PI*i/(nlen+1)), ipwr_t); // cosine window
// 
//         time_win[i] = tmp;
//     }
// }

float *generate_freqs(int ns, float dt, int fac, int *nfreq){
    
    /* declaration of variables */
    int i, npad, nfreq0;
    float Fs;
    float *freqs;
    //for frequency domain filtering, one SHOULD set fac to 1!!!
    if(fac >=0){
        npad = (int)(pow(2.0, ceil(log((double)(ns))/log(2.0))+fac*1.0));
    }else{
        npad = ns;
    }
    
    Fs=1/dt;
    nfreq0 = (int)ceil((npad+1)/2.0); //index is from 1
    //actually, it's not necessary to include the last point!!!
    
    freqs = vector(1, nfreq0);
    for(i=1;i<=nfreq0;i++) freqs[i]=(i-1)*Fs/npad;
    *nfreq = nfreq0;
    return freqs;
}

void gauss_window(int nfreq, float *freqs, float *win, float alfa, float fc){
    
    /* declaration of variables */
    int i;
    
    if(alfa<=0.0){
        printf("error occurs in util_fun: gauss_window\n");
        exit(1);
    }
    for(i=1; i<=nfreq; i++){
        win[i] = exp(-alfa*pow((freqs[i] - fc)/fc, 2));
    }
}

void cosine_window(int nfreq, float *freqs, float *win, float *fc){
    /* the filter is designed by Renat Shigapov, in page 26 */
    
    /* declaration of variables */
    int i;
    
    for(i=1; i<=nfreq; i++){
        if(freqs[i]>=fc[0] && freqs[i]<=fc[1]){
            win[i] = (1.0-cos(PI*(freqs[i]-fc[0])/(fc[1]-fc[0])))/2.0;
        }else if(freqs[i]>=fc[2] && freqs[i]<=fc[3]){
            win[i] = (1.0+cos(PI*(freqs[i]-fc[2])/(fc[3]-fc[2])))/2.0;
        }else if(freqs[i]>fc[1] && freqs[i]<fc[2]){
            win[i] = 1.0;
        }else{
            win[i] = 0.0;
        }
    }
}

void cosine_window2(int nfreq, float *freqs, float *win, float fc){
    /* the filter is designed by Renat Shigapov, in page 26 */
    
    /* declaration of variables */
    int i;
    float fc0=1.0, fc1=fc, fc2=fc, fc3=2.0*fc;
    
    for(i=1; i<=nfreq; i++){
        if(freqs[i]>=fc0 && freqs[i]<=fc1){
            win[i] = (1.0-cos(PI*(freqs[i]-fc0)/(fc1-fc0)))/2.0;
        }else if(freqs[i]>=fc2 && freqs[i]<=fc3){
            win[i] = (1.0+cos(PI*(freqs[i]-fc2)/(fc3-fc2)))/2.0;
        }else if(freqs[i]>fc1 && freqs[i]<fc2){
            win[i] = 1.0;
        }else{
            win[i] = 0.0;
        }
    }
}

void Taper(float *data, int type, int nlen, float width){
    /* apply taper action at each end of an array */
    /* modified from SAC software */
    /* declaration of variables */
    /* 
     * TYPE 
     *      1 SINE:    Apply a sine taper. 
     *      2 COSINE:  Apply a cosine taper.
     *      3 HANNING: Apply a Hanning taper.
            4 HAMMING: Apply a Hamming taper.
     */
    
    int i, ipts;
    float f0, f1, tpts, omega, value;
    
    /* the width should be between 0 and 0.5 */
    if(width > 0.5){
        printf("Errors occurs in func: %s \n", __func__);
        exit(1);
    }
    
    /* -- Determine number of points for taper. */
    tpts = width*(float)( nlen + 1 );
    ipts = max( 2, (int)( tpts ) );
                
    /* -- Taper each end. */
    if(type==1){/* Sine taper. */
        omega = PI/(2.*(float)( ipts ));
        /* ndxl = ndx1 + *nlen - 1; */
        for( i = 0; i <= (ipts - 1); i++ ){
            value = sin( omega*(float)( i ) );
            data[i] *= value;
            data[nlen-1-i] *= value;
        }
    }else if(type==2){ /* Cosine taper. */
        omega = PI/(2.*(float)( ipts ));
        /* ndxl = ndx1 + *nlen - 1; */
        for( i = 0; i <= (ipts - 1); i++ ){
            value = 1.0 - cos( omega*(float)( i ) );
            data[i] *= value;
            data[nlen-1-i] *= value;
        }
    }else if(type==3 ||type==4){
    /* --- Hanning and Hamming tapers differ only in coefficients. */    
        if(type==3){/* Hanning taper */
            f0 = 0.50;
            f1 = 0.50; 
        }else{/* Hamming taper */
            f0 = 0.54;
            f1 = 0.46;
        }
        omega = PI/(float)( ipts );
        /* ndxl = ndx1 + *nlen - 1; */
        for( i = 0; i <= (ipts - 1); i++ ){
            value = f0 - f1*cos( omega*(float)( i ) );
            data[i] *= value;
            data[nlen-1-i] *= value;
        }
    }
    
}

