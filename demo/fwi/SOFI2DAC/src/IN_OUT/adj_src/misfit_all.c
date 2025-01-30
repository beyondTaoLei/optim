/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2016  For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS.
 *
 * IFOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 *
 * IFOS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with IFOS. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
 -----------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   Calculate Data Residuals and Misfit
 *  ----------------------------------------------------------------------*/
#include "fd.h"
#include "adj.h"

void calc_misfit_res_energy_kernel(float **section, float **sectiondata, int *kill_vector, float **sectiondiff, float *misfit, float *energy){
    
//     #include "globvar1_extern.h"
//     #include "globvar2_extern.h"
    extern int NS, NTR, LNORM, TRKILL;
    extern float DT;
    
    extern int POW_LNORM8;
    extern float WATERLEVEL_LNORM8;
    extern float FC_source;

    /* declaration of variables */
    int ntr, ns;
    int i,j,invtime, umax=0;
    float rms, signL1;
    float l2=0.0, energy2=0.0;
    float abs_section, abs_sectiondata, sectiondata_mult_section;
    float *abs_section_arr=NULL, *abs_sectiondata_arr=NULL, *sectiondata_mult_section_arr=NULL;
    float **intseis_sectiondata_envelope=NULL, **intseis_section_envelope=NULL, **intseis_sectiondata_hilbert=NULL, **intseis_section_hilbert=NULL;
    float **env_res=NULL, **dummy_1=NULL, **dummy_2=NULL;
    char seis_save[STRING_SIZE];
    float stationary_lag;
    float *tr_obs=NULL, *tr_syn=NULL, *tr_res=NULL, *lag_arr=NULL, *misf_arr=NULL;
    int maxlag=NS/2;
    
    ntr=NTR;
    ns=NS;

    /* sectiondiff will be set to zero */
    umax=ntr*ns;
    zero(&sectiondiff[1][1],umax);
    
    /*rms=sqrt(rms/Lcount);*/
    rms=1.0;
    
                    /*******************************************/
                    /*              data manipulation          */
                    /*******************************************/
    /* Integration of measured and synthetic data  */
    if(VELOCITY==0){
        matrix_integration(NTR, NS, DT, section, section);
        matrix_integration(NTR, NS, DT, sectiondata, sectiondata);
    }
    
    /* TIME WINDOWING */
    if(TIMEWIN){
        matrix_operate(1, NS, 1, NTR, 1, 1, section, timewin, section, 3);
        matrix_operate(1, NS, 1, NTR, 1, 1, sectiondata, timewin, sectiondata, 3);
    }
    /* NORMALIZE TRACES */
    if(NORMALIZE==1){
        normalize_data(section,NTR, NS);
        normalize_data(sectiondata,NTR, NS);
    }
    
                    /*******************************************/
                    /*         prepare for special LNORM       */
                    /*******************************************/
    /*calculate energy for (LNORM==5)||(LNORM==7)*/
    if((LNORM==5)||(LNORM==7)){
        abs_section_arr=vector(1,ntr);
        abs_sectiondata_arr=vector(1,ntr);
        sectiondata_mult_section_arr=vector(1,ntr);
    
        for(i=1;i<=ntr;i++){
            if((TRKILL!=0)&&(kill_vector[i]==1))    continue;
            abs_sectiondata=0.0;
            abs_section=0.0;
            sectiondata_mult_section=0.0;
            for(j=1;j<=ns;j++){
                abs_sectiondata+=sectiondata[i][j]*sectiondata[i][j];
                abs_section+=section[i][j]*section[i][j];
                sectiondata_mult_section+=sectiondata[i][j]*section[i][j]; /* calculation of dot product for measured (section) and synthetic (sectiondata) data*/
            }
            /* keep in memory */
            abs_sectiondata_arr[i]=sqrt(abs_sectiondata);
            abs_section_arr[i]=sqrt(abs_section);
            sectiondata_mult_section_arr[i]=sectiondata_mult_section;
        }
    }
    /* calculate envelope and hilbert transform for LNORM==8 */
    if(LNORM==8){
        intseis_section_envelope = matrix(1,ntr,1,ns);
        intseis_sectiondata_envelope = matrix(1,ntr,1,ns);
        intseis_section_hilbert = matrix(1,ntr,1,ns);
        intseis_sectiondata_hilbert = matrix(1,ntr,1,ns);
        env_res = matrix(1,ntr,1,ns);
        dummy_1 = matrix(1,ntr,1,ns);
        dummy_2 = matrix(1,ntr,1,ns);
        
        /*calc_envelope(sectiondata,intseis_sectiondata_envelope,ns,ntr);
        calc_envelope(section,intseis_section_envelope,ns,ntr);
        calc_hilbert(section,intseis_section_hilbert,ns,ntr);*/
        calc_hilbert_su(sectiondata, intseis_sectiondata_hilbert, ns, ntr);
        matrix_energy(1, ns, 1, ntr, 1, 1,sectiondata, intseis_sectiondata_hilbert, intseis_sectiondata_envelope);
        calc_hilbert_su(section, intseis_section_hilbert, ns, ntr);
        matrix_energy(1, ns, 1, ntr, 1, 1,section, intseis_section_hilbert, intseis_section_envelope);
        ////////////////
        //sprintf(seis_save,"%s_L8.su",SEIS_FILE);
        //saveseis_local(1, intseis_sectiondata_hilbert, seis_save);
        ////////////////
        if(POW_LNORM8==1){
            for(i=1;i<=ntr;i++){
            for(j=1;j<=ns;j++){
                env_res[i][j] = intseis_section_envelope[i][j] - intseis_sectiondata_envelope[i][j];
            }}
        }else if(POW_LNORM8==2){
            for(i=1;i<=ntr;i++){
            for(j=1;j<=ns;j++){
                env_res[i][j] = intseis_section_envelope[i][j]*intseis_section_envelope[i][j] - intseis_sectiondata_envelope[i][j]*intseis_sectiondata_envelope[i][j];
            }}
        }
        if(POW_LNORM8==1){
            for(i=1;i<=ntr;i++){
            for(j=1;j<=ns;j++){
                dummy_1[i][j] = env_res[i][j] / (intseis_section_envelope[i][j]+WATERLEVEL_LNORM8) * intseis_section_hilbert[i][j];
            }}
        }else if(POW_LNORM8==2){
            for(i=1;i<=ntr;i++){
            for(j=1;j<=ns;j++){
                dummy_1[i][j] = env_res[i][j] * intseis_section_hilbert[i][j];
            }}
        }
        //calc_hilbert(dummy_1,dummy_2,ns,ntr);
        calc_hilbert_su(dummy_1, dummy_2, ns, ntr);
    }
    /* end of calculate envelope and hilbert transform for LNORM==8 */
    
    /* allocate memory for WTI calculation if LNORM==9 */
    if(LNORM==9){
        tr_obs=vector(1, ns);
        tr_syn=vector(1, ns);
        tr_res=vector(1, ns);
//         lag_arr=vector(1, ntr);
        misf_arr=vector(1, ntr);
    }
    
                    /*************************************/
                    /*         calculate residuals       */
                    /*************************************/
    for(i=1;i<=ntr;i++){
        if((TRKILL!=0)&&(kill_vector[i]==1))    continue;
        
        /*reverse time direction */
        invtime=ns;
        for(j=1;j<=ns;j++){
                    
            /* calculate L1 residuals */
            if(LNORM==1){
                if(((sectiondata[i][j]-section[i][j])/rms)>0){signL1=1.0;}
                if(((sectiondata[i][j]-section[i][j])/rms)<0){signL1=-1.0;}
                if(((sectiondata[i][j]-section[i][j])/rms)==0){signL1=0.0;}
                sectiondiff[i][invtime]=signL1/rms;
            }
            
            /* calculate L2 residuals */
            if(LNORM==2){
                sectiondiff[i][invtime] = section[i][j]-sectiondata[i][j];
            }
            
            /* calculate Cauchy residuals */  /* NOT UP TO DATE */
            if(LNORM==3){
                sectiondiff[i][invtime]=((sectiondata[i][j]-section[i][j])/rms)/(1+(((sectiondata[i][j]-section[i][j])/rms)*((sectiondata[i][j]-section[i][j])/rms)))/rms;
            }
            
            /* calculate sech residuals */   /* NOT UP TO DATE */
            if(LNORM==4){
                sectiondiff[i][invtime]=(tanh((sectiondata[i][j]-section[i][j])/rms))/rms;
            }
            
            /* calculate LNORM 5 or LNORM 7 residuals */
            if((LNORM==5) || (LNORM==7)){
                sectiondiff[i][invtime]=((section[i][j]*sectiondata_mult_section_arr[i])/(abs_section_arr[i]*abs_section_arr[i]*abs_section_arr[i]*abs_sectiondata_arr[i])) - (sectiondata[i][j]/(abs_section_arr[i]*abs_sectiondata_arr[i]));
            }
            
            /* calculate LNORM 8 residuals */
            if (LNORM==8){
                if(POW_LNORM8==1){
                    sectiondiff[i][invtime]=section[i][j]/(intseis_section_envelope[i][j]+WATERLEVEL_LNORM8)*env_res[i][j]-dummy_2[i][j];
                }else if(POW_LNORM8==2){
                    sectiondiff[i][invtime]=section[i][j]*env_res[i][j]-dummy_2[i][j];
                }
            }
            
            invtime--;          /* reverse time direction */
        }
        
        /* calculate WTI residuals */
        if(LNORM==9){
//             sectiondiff[i][invtime] = section[i][j]-sectiondata[i][j];
            for(j=1;j<=ns;j++){
                tr_obs[j]=sectiondata[i][j];
                tr_syn[j]=section[i][j];
            }
//             WTI_misfit(tr_obs, tr_syn, tr_res, maxlag, ns, DT, &stationary_lag);
//             WTI_misfit_syn(tr_obs, tr_syn, tr_res, maxlag, ns, DT, &stationary_lag);
//             lag_arr[i]=stationary_lag;
            
            misf_arr[i] = misfit_TT(tr_obs, tr_syn, tr_res, 0.0, ns*DT, ns, DT, &stationary_lag);
            //if(MYID==0) printf("%d %f\n",i, stationary_lag);
            
            invtime=ns;
            for(j=1;j<=ns;j++){
                sectiondiff[i][invtime]=tr_res[j];
                invtime--;
            }
        }
    }
                    /*******************************************/
                    /*         special work for LNORM = 8      */
                    /*******************************************/
    /* because of noise of high frequency */
    if(LNORM==8)    timedomain_filt(sectiondiff,FC_source,2,ntr,ns,1);

                    /************************************************************/
                    /* calculate misfit between measured data and modeling data */
                    /************************************************************/
    for(i=1;i<=ntr;i++){
        if((TRKILL!=0)&&(kill_vector[i]==1))    continue;
        for(j=1;j<=ns;j++){
            /* calculate norm for different LNORM */
            if(LNORM==2){
                l2+=sectiondiff[i][j]*sectiondiff[i][j];
            }
            
            if((LNORM==5)||(LNORM==7)){
                l2-=(sectiondata[i][j]*section[i][j])/(abs_sectiondata_arr[i]*abs_section_arr[i]);
            }
            
            if(LNORM==8){
                //l2+=(intseis_section_envelope[i][j]-intseis_sectiondata_envelope[i][j]) * (intseis_section_envelope[i][j]-intseis_sectiondata_envelope[i][j]);
                l2+=env_res[i][j] *env_res[i][j];
            }
            
            /*if(sws==2){
                l2+=fabs(sectiondiff[i][j])*fabs(sectiondiffold[i][j]);
            }*/
            
            /*l2+=sectiondiff[i][j];*/
            
        }
        if(LNORM==9){
//             l2+=lag_arr[i]*lag_arr[i];
            l2 += misf_arr[i];
        }
    }                  
    
                    /*************************************/
                    /* calculate energy of measured data */
                    /*************************************/
    for(i=1;i<=ntr;i++){
        if((TRKILL!=0)&&(kill_vector[i]==1))    continue;
        for(j=1;j<=ns;j++){
            energy2 += sectiondata[i][j]*sectiondata[i][j];
        }
    }
    *misfit=l2;
    *energy=energy2;
    
    if((LNORM==5)||(LNORM==7)){
        free_vector(abs_section_arr,1,ntr);
        free_vector(abs_sectiondata_arr,1,ntr);
        free_vector(sectiondata_mult_section_arr,1,ntr);
    }
    if(LNORM==8){
        free_matrix(intseis_sectiondata_envelope,1,ntr,1,ns);
        free_matrix(intseis_section_envelope,1,ntr,1,ns);
        free_matrix(intseis_section_hilbert,1,ntr,1,ns);
        free_matrix(intseis_sectiondata_hilbert,1,ntr,1,ns);
        free_matrix(env_res,1,ntr,1,ns);
        free_matrix(dummy_1,1,ntr,1,ns);
        free_matrix(dummy_2,1,ntr,1,ns);
    }
    if(LNORM==9){
        free_vector(tr_obs, 1, ns);
        free_vector(tr_syn, 1, ns);
        free_vector(tr_res, 1, ns);
//         free_vector(lag_arr, 1, ntr);
        free_vector(misf_arr, 1, ntr);
    }
    
} /* end of function */
