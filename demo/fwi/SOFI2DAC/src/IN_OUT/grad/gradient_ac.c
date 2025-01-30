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

#include "grad.h"

/* The function is called in the time loop */
/*void eprecond_time(int n1, int n2, float *vx, float *vy, float *W){
    
    V_square_sum_a(n1, n2, vx, vy, W);
    
}*/

/* The function is called in the time loop */
void calc_correlation_ac_time(int n1,int n2,int inv_vp,int inv_rho,
        float *fvx, float *fvy, float *fsp, float *bvx, float *bvy, 
        float *bsp, float *cvx, float *cvy, float *csp){
    /* before calculating the gradient, one has to calculate zero-lag
     * cross correlation of incident wavefield and adjoint wavefield.
     * 
     * For time-domain simulation, one will extract the snapshots of 
     * specified components for FWI gradient at any time step, this 
     * function is used to set the pointer that can point to these 
     * kinds of components.
     */
    /* Local varibales */
    
    
    /* INPUT */
    /* int n1, n2: the start and end of index,  
     * 
     * st_wfd *wfdi: includes all of the required incident wavefields
     * st_wfd *wfda: includes all of the required adjoint wavefields
     * 
     */
    /* OUTPUT */
    /* st_grad *grad: includes all of the gradients(vp/vs/rho...)
     * 
     */
    
    if(inv_vp||inv_rho)
        V_multiply_a(n1, n2,fsp, bsp,csp);
    if(inv_rho){
        V_multiply_a(n1, n2,fvx, bvx,cvx);
        V_multiply_a(n1, n2,fvy, bvy,cvy);
    }
}

/* The following is called after the time loop */
void gradient_ac(int n1, int n2, float dt, int inv_vp, int inv_rho, 
                 float *Mrho, float *Mvp, float *cvx,float *cvy,float *csp,
                 float *grad_lam, float *grad_rho_s,float *grad_vp, float *grad_rho){
    /* converts zero-lag cross correlation into gradient */
    
    /* Local varibales */
    int i;
    float lam, lam2;
    
    /* INPUT */
    /*float *cvx;
    float *cvy;
    float *csp;*/
    /* OUTPUT */
    /*
    float *grad_lam;
    float *grad_rho_s;
    float *grad_vp;
    float *grad_rho;*/
    
    /* ----------------------------------------- */
    /* calculate gradient direction lamda rho    */
    /* ----------------------------------------- */
    if(inv_vp||inv_rho){
        for(i=n1;i<=n2;i++){
            lam = Mrho[i]*(Mvp[i]*Mvp[i]);
            lam2=lam*lam;
            /* calculate gradient direction lamda */
            grad_lam[i] =-dt*(0.25*csp[i]/lam2);
        }
    }
    
    if(inv_rho){
        /* calculate gradient direction rho */
        for (i=n1;i<=n2;i++) grad_rho_s[i]=-dt*(cvx[i]+cvy[i]);
    }
    
    /* ----------------------------------------- */
    /* calculate complete gradient direction     */
    /* ----------------------------------------- */
    if(inv_vp){
        for (i=n1;i<=n2;i++) grad_vp[i]=2.0*Mvp[i]*Mrho[i]*grad_lam[i];
    }
    
    if(inv_rho){
        for (i=n1;i<=n2;i++)
            grad_rho[i] =Mvp[i]*Mvp[i]*grad_lam[i] + grad_rho_s[i];
    }
    
}


void save_grad_ac(int n1, int n2, int inv_vp, int inv_rho, 
               float *grad_vp, float *grad_rho, int ishot){
    
    /* Local varibales */
    char str_vp[225], str_vs[225], str_rho[225];
    
    sprintf(str_vp,"./gradient/grad_vp_shot%d.bin",ishot);
    sprintf(str_rho,"./gradient/grad_rho_shot%d.bin",ishot);
    if(inv_vp)  V_write(str_vp,&grad_vp[n1], n2-n1+1);
    if(inv_rho) V_write(str_rho,&grad_rho[n1], n2-n1+1);
    
    
}
