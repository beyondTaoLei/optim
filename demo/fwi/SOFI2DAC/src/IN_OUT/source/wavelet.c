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

#include "source.h"

float ** wavelet(float ** srcpos_loc, int nsrc, int src_shape, int nt, float dt){
    
    /*local variables */
    int nts, it, k;
    float *psource=NULL, tshift, amp=0.0, a, fc, tau, t, ts, ag;
    float ** signals;

        
    signals=matrix(1,nsrc,1,nt);
    
    for (it=1;it<=nt;it++){
        t=(float)it*dt;
        
        for (k=1;k<=nsrc;k++) {
            tshift=srcpos_loc[4][k];
            fc=srcpos_loc[5][k];
            a=srcpos_loc[6][k];
            ts=1.0/fc;
            
            switch (src_shape){
                case 1 :
                    /* New Ricker Wavelet, equal to SOFI2D */
                    tau=PI*(t-1.5*ts-tshift)/(ts);
                    amp=(((1.0-2.0*tau*tau)*exp(-tau*tau)));
                    break;
                case 2 :
                    if ((t<tshift) || (t>(tshift+ts))) amp=0.0;
                    else amp=((sin(2.0*PI*(t-tshift)*fc)
                                -0.5*sin(4.0*PI*(t-tshift)*fc)));
                    break;
                case 3 :
                    /* source wavelet from file SOURCE_FILE */
                    /*if (it<=nts) amp=psource[it];
                    else amp=0.0;*/
                    break;
                case 4 :
                    /* sinus raised to the power of three */
                    if ((t<tshift) || (t>(tshift+ts))) amp=0.0;
                    else amp=pow(sin(PI*(t+tshift)/ts),3.0);
                    break;
                    
                    break;
                case 5 :
                    /* first derivative of a Gaussian */
                    ts=1.2/fc;
                    ag  = PI*PI*fc*fc;
                    amp = - 2.0 * ag * (t-ts) * exp(-ag*(t-ts)*(t-ts));
                    break;
                case 6 :
                    /* Bandlimited Spike */
                    amp=0.0;
                    if(it==1+iround(tshift/dt)){
                        amp = 1.0;}
                    break;
                case 7 :
                    /* source wavelet from file SOURCE_FILE */
                    //amp=psource[it];
                    break;
                case 8 :
                    /* integral of sinus raised to the power of three */
                    if (t<tshift) {
                        amp=0.0;}
                    if ((t>=tshift) && (t<=(tshift+ts))){
                        amp=(ts/(0.75*PI))*(0.5-0.75*cos(PI*(t-tshift)/ts)+0.25*pow(cos(PI*(t-tshift)/ts),3.0));}
                    if (t>(tshift+ts))
                    {amp=ts/(0.75*PI);}
                    break;
                default :
                    printf("Which source-wavelet ? ");
            }
            
            
            signals[k][it]=amp*a;
        }
    }
    
    return signals;	
    
}