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
/* files to include */
#include "comm.h"
#include <fftw3.h>

void V_der(int ns, float dt, float *x, float *r){
    /* the gradient of an N-dimensional array 
       where index range of x is 0: ns-1
     */
    int i, nfreq0, npad;
    float coeff, tmp, tmpi, *freqs=NULL;
    
    /* declaration of variables for FFTW3*/
    fftw_complex *in, *out;
    fftw_plan p1, p2;
    
    npad=ns;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
    
    coeff=2.0*PI/(npad*dt);
    freqs =(float *)malloc((size_t) (npad*sizeof(float)));
    for(i=0;i<npad;i++) freqs[i]=0.0;
    //L=ns*dt;
    //For *even* N, k = (2*pi/L) * [0:N/2-1, 0, -N/2+1:-1];
    //For *odd* N, use k = (2*pi/L) * [0:(N-1)/2, -(N-1)/2:-1];
    nfreq0 = (int)ceil(npad/2.0);
    for(i=1;i<nfreq0;i++){
        freqs[i]=coeff*i;
        freqs[npad-i]=-freqs[i];
    }
    
    for (i=0;i<ns;i++){
        in[i][0]=x[i];
        in[i][1]=0.0;
        //printf("%d %f\n",i, freqs[i]/coeff);
    }
    
    /* FFT */
    p1 = fftw_plan_dft_1d(npad, in, out, -1, FFTW_ESTIMATE);
    fftw_execute(p1);
    fftw_destroy_plan(p1);		    

    /* elementwise multiplication of FFT and (i*omega) */
    for (i=0;i<npad;i++){
        tmp=out[i][0];
        tmpi=out[i][1];
        out[i][0] = -tmpi*freqs[i];
        out[i][1] =  tmp *freqs[i];
    }

    /* inverse FFT */
    p2 = fftw_plan_dft_1d(npad, out, in, 1, FFTW_ESTIMATE);
    fftw_execute(p2);
    fftw_destroy_plan(p2);

    /* derivative */
    for (i=0;i<ns;i++){
        r[i]=(1.0/npad)*in[i][0];
    }

    fftw_free(in);
    fftw_free(out);
    free(freqs);
}

void V_int(int ns, float dt, float *x, float *r){
    /* the integration of an N-dimensional array 
       where index range of x is 0: ns-1
       noting, the 0-th frequency component is removed!
     */
    int i, nfreq0, npad;
    float coeff, tmp, tmpi, *freqs=NULL;
    
    /* declaration of variables for FFTW3*/
    fftw_complex *in, *out;
    fftw_plan p1, p2;
    
    npad=ns;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
    
    coeff=2.0*PI/(npad*dt);
    freqs =(float *)malloc((size_t) (npad*sizeof(float)));
    for(i=0;i<npad;i++) freqs[i]=0.0;
    //L=ns*dt;
    //For *even* N, k = (2*pi/L) * [0:N/2-1, 0, -N/2+1:-1];
    //For *odd* N, use k = (2*pi/L) * [0:(N-1)/2, -(N-1)/2:-1];
    nfreq0 = (int)ceil(npad/2.0);
    for(i=1;i<nfreq0;i++){
        freqs[i]=coeff*i;
        freqs[npad-i]=-freqs[i];
    }
    
    for (i=0;i<ns;i++){
        in[i][0]=x[i];
        in[i][1]=0.0;
        //printf("%d %f\n",i, freqs[i]/coeff);
    }
    
    /* FFT */
    p1 = fftw_plan_dft_1d(npad, in, out, -1, FFTW_ESTIMATE);
    fftw_execute(p1);
    fftw_destroy_plan(p1);		    

    /* elementwise multiplication of FFT and (i*omega) */
    for (i=0;i<npad;i++){
        tmp=out[i][0];
        tmpi=out[i][1];
        out[i][0] = tmpi/freqs[i];
        out[i][1] = -tmp/freqs[i];
    }
    out[0][0]=0.0; out[0][1]=0.0;
    if(npad%2==0){
        i=(int)npad/2;
        out[i][0]=0.0; out[i][1]=0.0;
    }

    /* inverse FFT */
    p2 = fftw_plan_dft_1d(npad, out, in, 1, FFTW_ESTIMATE);
    fftw_execute(p2);
    fftw_destroy_plan(p2);

    /* derivative */
    for (i=0;i<ns;i++){
        r[i]=(1.0/npad)*in[i][0];
    }

    fftw_free(in);
    fftw_free(out);
    free(freqs);
}

/*
int main(int argc, char **argv){
*/    
    /* Local varibales */
/*    int i, j,n1, n2;
    float a[10],b[10],c[10],d[10],e[10],f[10],tmp;
    n1=0;
    n2=9;
    
    
    for (i=n1;i<=n2;i++){
        a[i]=i*1.0;
        b[i]=i*2.0;
        c[i]=i*3.0;
        d[i]=i*4.0;
        e[i]=i*5.0;
        f[i]=i*6.0;
    }*/
    //a[1]=-1;a[4]=-20;
    
    //V_copy(n1, n2, a, c);
    //V_add(n1, n2, a, b, c);
    //V_subtract(n1, n2, a, b, c);
    //V_multiply(n1, n2, a, b, c);
    //V_divide(n1, n2, a, b, a);
    //V_square_sum(n1, n2, a, b, c);
    //V_square_sum_a(n1, n2, a, b, c);
//     V_AXPY(n1, n2, 5.0, a, b, c);
    //V_AXPb(n1, n2, 10.0, a,1, c);
    
    //V_scaling(n1, n2, a,10.0);
    //V_sqrt(n1, n2, a,c);
    //V_multiply_sqrt(n1, n2, a, b, c);
    //V_multiply_a(n1, n2, a, b, c);
    /*V_multiply_a2(n1, n2, a, b, c,d,e,f);*/
/*    printf("\n");
    for (i=n1;i<=n2;i++) printf("%5.1f ",c[i]);
    printf("\n");
    
    //V_write("1.bin",e,10);
    
    //V_read("/home/tao/seis_software/proj01/precond/approx_hessian_it1_shot1.bin",a,1700*350);
    //tmp=global_percentile(a, n1, n2, 0.03);
    //printf("%e \n",tmp);
    //printf("%f \n",V_normL2(n1, n2, a));
    //i=find_min_positive_num1(n1,n2,a);
    //printf("%d \n",i);
    //for (i=n1;i<=n2;i++) printf("%5.1f ",f[i]);
    
}*/

