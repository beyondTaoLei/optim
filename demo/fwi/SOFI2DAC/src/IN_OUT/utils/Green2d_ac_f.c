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
/**
@brief Calculate the gradient in the frequency domain for 2D acoustic case
@details This program is meant to be run to calculate the gradient for 2D 
acoustic case. One can apply discrete Fourier transfomation(DFT) on the 
fly when time domain forward modelling is employed.
@file Green2d_ac_f.c
@author Tao Lei (leit@sustech.edu.cn), Wei Zhang
@date 18 Sep 2019
@param[in]  verbose \b 1: print the details; \b 0: run in quiet pattern, default value;
@param[in]  omega1 the omega1-th index of the frequency sample;
@param[in]  omega_inc the increment of frequency sample, omega=omega1+omega_inc*(i-1);
@param[in]  fwav_s the full path of wavelet <em>(clockwise)</em> which is added at the receiver position;
@param[in]  fadj_s the actual adjoint source <em>(clockwise)</em>;
@param[in]  fmod    the full path of vp model, the path of density model 
                    is derived from vp model, the following [fin, fadj and fout_vp] are similar;
@param[in]  fin    the full path of event name (from the source side);
@param[in]  fadj    the full path of station name (from the receiver side);
@param[in,out] fout_vp the full path of vp gradient;
@code 
mpicc -o Green2d_ac_f Green2d_ac_f.c -I ../include/ -L../lib \
-L/home/tao/seis_software/SAC/lib -lIN_OUT -lsac -lsacio -lm -lfftw3

mpirun -np 4 Green2d_ac_f verbose=1 omega1=2 omega_inc=1 \
fwav_s=source/glob/135.sac fadj_s=data/syn_pre/40/135.p.sac fmod=./model/inv/inv.vp \
fin=snap/inc/40 fadj=snap/adj/135 fout_vp=gradient/40.135.vp  

@endcode
*/

#include "par.h"
#include "source.h"
#define STRING_SIZE 256

/* Global variables */
int xargc; 
char **xargv;

int main(int argc, char *argv[]){
    
    /* IN and OUT */
    int verbose, omega1, omega_inc,num, nfreq;
    char *fwav_s, *fin, *fout;
    
    /* Local variables */
    int nlen, it,i,j;
    float dt, corr_r, corr_i, max1, df,omega;
    int *irecord=NULL;
    float *vec, *wavelet_r, *wavelet_i, *wavelet_abs;
    float *sp_inc, *sp_i_inc;
    char f_sp_inc[256], f_sp_i_inc[256]; 
    
    /* 1. read formatted input list */
    /* initialize getpar */
    initargs(argc,argv);
    
    /* these lines are useful for SEISPLATFORM software */
    verbose=0; getparint("verbose",&verbose);
    if(!getparint("omega1", &omega1)) throw_error("omega1 must be specified!\n");
    if(!getparint("omega_inc", &omega_inc)) throw_error("omega_inc must be specified!\n");
    if(!getparint("num", &num)) throw_error("num must be specified!\n");
    if(!getparint("nfreq", &nfreq)) throw_error("nfreq must be specified!\n");
    /* the file or path */
    if(!getparstring("fwav_s", &fwav_s)) throw_error("fwav_s must be specified!\n");
    if(!getparstring("fin", &fin)) throw_error("fin must be specified!\n");
    if(!getparstring("fout", &fout)) throw_error("fout must be specified!\n");

    /* allocate the memory */
    /* adjoint source */
    irecord = ivector(1,nfreq);
    /* wavefield space */
    sp_inc = vector(1, num*nfreq);
    sp_i_inc = vector(1, num*nfreq);
    
    /* read the wavefields */
    sprintf(f_sp_inc,"%s.sp.bin",fin);
    sprintf(f_sp_i_inc,"%s.sp.i.bin",fin);
    V_read(f_sp_inc, &sp_inc[1], num*nfreq);
    V_read(f_sp_i_inc, &sp_i_inc[1], num*nfreq);
    
    /*3. calculate the correlation */
    /* wavelet in sac */
    vec = read_sac(fwav_s, &nlen, &dt);
    dt=1.0;
    wavelet_r = vector(1, nlen);
    wavelet_i = vector(1, nlen);
    wavelet_abs = vector(1, nlen);
    for(it=1; it<nlen-1; it++) wavelet_r[it+1]=(vec[it+1]-vec[it-1])/(2.0*dt);
    //for(it=2; it<=nlen; it++) wavelet_r[it]=wavelet_r[it-1]+vec[it-1];
    free(vec);
    
    /*3. DFT */
    fft_real_imag(&wavelet_r[1], &wavelet_i[1], nlen, -1);
    max1=-1.0;
    for(it=2; it<=nlen; it++){/*ignore the zero component */
        wavelet_abs[it]=sqrt(wavelet_r[it]*wavelet_r[it]+wavelet_i[it]*wavelet_i[it]);
        if(wavelet_abs[it]>max1) max1=wavelet_abs[it];
    }
    //if(verbose) printf("max(abs(wavelet))=%e\n",max1);
    if(max1==0.0) throw_error("Please check wavelet file!\n");

    df=2.0*PI/(nlen*dt);
    for(it=1; it<=nfreq; it++){
    //for(it=41; it<=41; it++){
        irecord[it]=0;
        i=omega1+omega_inc*(it-1);
        printf("it, if=%d %d\n",it, i);
        omega=(i-1)*df;
        if(wavelet_abs[i]>=0.001*max1){
            irecord[it]=1;
        }
    }

    /* remove the wavelet factor */ 
    for(it=1;it<=nfreq;it++){
        /* remove the compoments with small amplitude */
        if(irecord[it]==0){
            printf("%d\n", it);
            continue;
        }
        
        /* calculate the zero-lag cross-correlation */
        for(i=1;i<=num;i++){
            j=(it-1)*num+i;
            corr_r =(sp_inc[j]*wavelet_r[it] + sp_i_inc[j]*wavelet_i[it])/wavelet_abs[it];
            corr_i =(-sp_inc[j]*wavelet_i[it]+ sp_i_inc[j]*wavelet_r[it])/wavelet_abs[it];
            sp_inc[j]=corr_r;
            sp_i_inc[j]=corr_i;
        }
    }
    sprintf(f_sp_inc,"%s.sp.bin",fout);
    sprintf(f_sp_i_inc,"%s.sp.i.bin",fout);
    V_write(f_sp_inc, &sp_inc[1], num*nfreq);
    V_write(f_sp_i_inc, &sp_i_inc[1], num*nfreq);
    
    /* deallocate memory */
    free_vector(wavelet_r, 1, nlen);
    free_vector(wavelet_i, 1, nlen);
    free_vector(wavelet_abs, 1, nlen);
    free_ivector(irecord ,1,nfreq);
    free_vector(sp_inc ,1, num*nfreq);
    free_vector(sp_i_inc ,1, num*nfreq);
    
    if(verbose){
        printf("\n******\n");
        printf("dt        \t| %.3e\n",dt);
        printf("num       \t| %d\n",num);
        printf("nlen      \t| %d\n",nlen);
        printf("nfreq     \t| %d\n",nfreq);
        printf("omega1    \t| %d\n",omega1);
        printf("omega_inc \t| %d\n",omega_inc);
        printf("fin      \t| %s\n",fin);
        printf("fwav_s    \t| %s\n",fwav_s);
        printf("******\n");
    }
    
    return 0;
}

