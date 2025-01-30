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
@file tomo_el2d_fast.c
@author Tao Lei (leit@sustech.edu.cn), Wei Zhang
@date 18 Sep 2019
@param[in]  verbose \b 1: print the details; \b 0: run in quiet pattern, default value;
@param[in]  omega1 the omega1-th index of the frequency sample;
@param[in]  omega_inc the increment of frequency sample, omega=omega1+omega_inc*(i-1);
@param[in]  fmod the full path of vp model, the path of density and 
                 shear velocity model are derived from vp model;
@param[in]  fevents the file of event name list;
@param[in]  fstations the file of station name list;
@param[in]  fincf the folder of snapshots from the source side;
@param[in]  fadjf the folder of snapshots from the receiver side;
@param[in]  fwavf the folder of wavelet <em>(clockwise)</em> which is added at the receiver position;
@param[in]  fadj_sf the folder of actual adjoint source or seismograms <em>(clockwise)</em>, which supports more than one measurement;
@param[in,out] cmpt which component one will use, such as vx, vy;
@param[in,out] foutf the folder of vp gradient;
@code 
mpicc -o tomo_el2d_fast tomo_el2d_fast.c -I ../include/ -L../lib \
-L/home/tao/seis_software/SAC/lib -lIN_OUT -lsac -lsacio -lm

mpirun -np 40 ../seisplatform/bin/tomo_el2d_fast omega1=2 omega_inc=1 \
mem=20 fmod=model/inv/inv.vp fevents=list/eventname.list \
fstations=list/stationname.list fincf=snap/inc fadjf=snap/adj \
fwavf=source/glob fadj_sf=data/adj,data/adj_2 cmpt=vy foutf=gradient

@endcode
*/

#include "par.h"
#include "source.h"
#include "grad.h"
#include"mpi.h"
#define BLOCK_LOW(rank,size,n)    ((rank)*(n)/(size))
#define BLOCK_HIGH(rank,size,n)    (BLOCK_LOW((rank)+1,size,n)-1)
#define BLOCK_SIZE(rank,size,n)    (BLOCK_HIGH(rank,size,n)-BLOCK_LOW(rank,size,n)+1)
#define STRING_SIZE 256
#define SHIFT_IND 1
#define FREE_ARGUMENT char *
#define MAX_STACKS 10

/* Global variables */
int xargc; 
char **xargv;

void check_file_exist(char **event_list, char **station_list, char *fincf, 
                      char *fadjf, char *f_vp, char *f_vs, char *f_rho, 
                      char **phase_list, int nshots, int nrec, int nphase);
void read_wfd(char *finc, float *sxx_inc, float *sxx_i_inc, float *sxy_inc, 
              float *sxy_i_inc, float *syy_inc, float *syy_i_inc,
              int nfreq, int num, int mynum, int view_offset);
void read_wfd_one(char *f_sp_inc, float *sp_inc, int nfreq, int num, 
                  int mynum, int view_offset);
void read_mod(char *filename, float *mod, int view_offset, int mynum);
void write_mod3(MPI_File mpi_fh, float *mod, float *mod2, long long view_offset_g,
                long long view_offset_g2, int mynum);
void stress2strain2D(float *sxx_inc, float *sxx_i_inc, float *sxy_inc, 
                   float *sxy_i_inc, float *syy_inc, float *syy_i_inc,
                   float *Mvp, float *Mvs, float *Mrho, int mynum, int nfreq);
void cal_kernel(int mynum, int nfreq, int nlen, int nphase,
                float **fac_r_all, float **fac_i_all, float *sxx_inc, 
                float *sxx_i_inc, float *sxy_inc, float *sxy_i_inc, 
                float *syy_inc, float *syy_i_inc, float *sxx_adj, 
                float *sxx_i_adj, float *sxy_adj, float *sxy_i_adj, 
                float *syy_adj, float *syy_i_adj, float *Mvp, float *Mvs, 
                float *Mrho, float *grad_vp, float *grad_vs);

int main(int argc, char *argv[]){

    /* IN and OUT */
    int verbose, inv_vs, omega1, omega_inc;
    float mem;
    char *fmod, *fevents, *fstations, *fincf, *fadjf, *fwavf, *cmpt, *foutf;//, *flags , *fadj_sf
    
    /* Local variables */
    int size_wfd, size_wfd2, size_mod, view_offset, view_offset2;
    int num, nfreq, nlen, nphase, nshots, nrec, ntask, n_max;//, nphase2
    int mynum, mynum_loc, idx_g1, idx_l1;
    int it, ii, jj, iph, offset_s, offset_r, offset_p, itask;//im, 
    long long ic, jc, ip, idx, ibase, view_offset_g, view_offset_g2;
    float dt, mem_total;
    
    float *vec;
    float *Mvp, *Mvs, *Mrho, *grad_vp, *grad_vs;
    float *fac_all, **fac_r_all, **fac_i_all;
    float *sxx_inc, *sxy_inc, *syy_inc;
    float *sxx_i_inc, *sxy_i_inc, *syy_i_inc;
    float *sxx_adj, *sxy_adj, *syy_adj;
    float *sxx_i_adj, *sxy_i_adj, *syy_i_adj;
    char **event_list = NULL, **station_list =NULL;
    char finc[STRING_SIZE], fadj[STRING_SIZE], str1[STRING_SIZE];
    char f_vp[STRING_SIZE], f_vs[STRING_SIZE], f_rho[STRING_SIZE];
    /*MPI parameters */
    int NP,MYID;
    MPI_File *fp_kernel;
    /* multiple separated phases*/
    //int flag_list[MAX_STACKS];
    char *phase_list[MAX_STACKS];//*strf[MAX_STACKS], 
    
    /* Initialize MPI environment */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &MYID);
    MPI_Comm_size(MPI_COMM_WORLD, &NP);
    
    /* 1. read formatted input list */
    /* initialize getpar */
    initargs(argc,argv);
    
    /* these lines are useful for SEISPLATFORM software */
    verbose=0; getparint("verbose",&verbose);
    //inv_vs=0; getparint("inv_vs",&inv_vs);
    if(!getparint("omega1", &omega1)) throw_error("omega1 must be specified!\n");
    if(!getparint("omega_inc", &omega_inc)) throw_error("omega_inc must be specified!\n");
    if(!getparfloat("mem",&mem)) throw_error("mem must be specified!\n");
    /* the file or path */
    if(!getparstring("fmod", &fmod)) throw_error("fmod must be specified!\n");
    if(!getparstring("fevents", &fevents)) throw_error("fevents must be specified!\n");
    if(!getparstring("fstations", &fstations)) throw_error("fstations must be specified!\n");
    if(!getparstring("fincf", &fincf)) throw_error("fincf must be specified!\n");
    if(!getparstring("fadjf", &fadjf)) throw_error("fadjf must be specified!\n");
    if(!getparstring("fwavf", &fwavf)) throw_error("fwavf must be specified!\n");
    //if(!getparstring("fadj_sf", &fadj_sf)) throw_error("fadj_sf must be specified!\n");
    if(!(nphase=countparval("fadj_sf" ))) throw_error( "fadj_sf list must be specified!" );
    //if(!(nphase2=countparval("flags" ))) throw_error( "flags list must be specified!" );
    //if(nphase!=nphase2) throw_error( "fadj_sf and flags should have the same number elements!" );
    inv_vs=0;
    for(iph=1; iph<=nphase; iph++){
        getnparstringarray(iph, "fadj_sf", &(phase_list[iph]));
        //getnparstringarray(iph, "flags", &(strf[iph]));
        //flag_list[iph] = atoi(strf[iph]);
        //if(flag_list[iph]<1 ||flag_list[iph]>3){/* 1(vp+vs), 2(vp) or 3(vs) */
        //    throw_error( "the element of flag_list must be 1, 2 or 3!");
        //}
        //if(flag_list[iph]!=2) inv_vs=1;
    }
    if(!getparstring("cmpt", &cmpt)) throw_error("cmpt must be specified!\n");
    if(!getparstring("foutf", &foutf)) throw_error("foutf must be specified!\n");
    MPI_Barrier(MPI_COMM_WORLD);
    
    /* check the file size */
    event_list = read_string(&nshots, fevents);
    station_list = read_string(&nrec, fstations);
    if(MYID==0){
        sprintf(str1,"%s/%s.sxx.bin",fincf, event_list[0]);
        size_wfd=get_file_size(str1);
        sprintf(str1,"%s/%s.sxx.bin",fadjf, station_list[0]);
        size_wfd2=get_file_size(str1);
        size_mod=get_file_size(fmod);
        if((size_wfd%4!=0)||(size_mod%4!=0)||(size_wfd2!=size_wfd)
            ||(size_wfd%size_mod!=0)){
            printf("size_wfd, size_wfd2, size_mod=%d %d %d\n",
                size_wfd, size_wfd2, size_mod);
            throw_error("Please check the size of the files!\n");
        }else{
            num=size_mod/4;
            nfreq=size_wfd/size_mod;
        }
        sprintf(str1,"%s/%s.sac",fwavf,station_list[0]);
        vec = read_sac(str1, &nlen, &dt);
        free(vec);
    }
    /* bcast integers */
    MPI_Bcast(&num,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&nfreq,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&nlen,1,MPI_INT,0,MPI_COMM_WORLD);
    /* bcast floats */
    MPI_Bcast(&dt,1,MPI_FLOAT,0,MPI_COMM_WORLD);
    dt=1.0;
    /* bcast string */
    MPI_Barrier(MPI_COMM_WORLD);
    
    mynum = BLOCK_SIZE(MYID,NP,num);
    idx_g1=BLOCK_LOW(MYID,NP,num);
    mem_total=num/1024.*(nshots+nrec)/1024.*nfreq*6*4/1024.;
    ntask=(int) ceil(mem_total/mem);
    
    /* the maximum size for every reading in the model space */
    n_max=0;
    for(itask=0; itask<ntask; itask++){
        mynum_loc = BLOCK_SIZE(itask,ntask,mynum);
        if(n_max < mynum_loc) n_max = mynum_loc;
    }
    if(verbose&&(MYID==1||NP==1)){
        printf("\n******\n");
        printf("mem_total \t| %.2f G\n",mem_total);
        printf("mem       \t| %.2f G\n",mem);
        printf("ntask     \t| %d\n",ntask);
        printf("dt        \t| %.3e\n",dt);
        printf("inv_vs    \t| %d\n",inv_vs);
        printf("num       \t| %d\n",num);
        printf("nlen      \t| %d\n",nlen);
        printf("nfreq     \t| %d\n",nfreq);
        printf("omega1    \t| %d\n",omega1);
        printf("omega_inc \t| %d\n",omega_inc);
        printf("nshots    \t| %d\n",nshots);
        printf("nrec      \t| %d\n",nrec);
        printf("fmod      \t| %s\n",fmod);
        printf("fevents   \t| %s\n",fevents);
        printf("fstations \t| %s\n",fstations);
        printf("fincf     \t| %s\n",fincf);
        printf("fadjf     \t| %s\n",fadjf);
        printf("fwavf     \t| %s\n",fwavf);
        printf("nphase    \t| %d\n",nphase);
        //printf("fadj_sf   \t| %s\n",fadj_sf);
        //printf("flags     \t| %s\n",flags);
        printf("cmpt      \t| %s\n",cmpt);
        printf("foutf     \t| %s\n",foutf);
        printf("******\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    sprintf(f_vp,"%s",fmod);
    strrpc2(f_vp,"vp","vs",f_vs);
    strrpc2(f_vp,"vp","rho",f_rho);
    if(MYID==0) check_file_exist(event_list, station_list, fincf, fadjf, f_vp, 
                     f_vs, f_rho, phase_list, nshots, nrec, nphase);
    
    /* allocate the memory */
    /* model space */
    Mvp     = vector(1, n_max);
    Mvs     = vector(1, n_max);
    Mrho    = vector(1, n_max);
    grad_vp = vector(1, n_max*nphase);
    grad_vs = vector(1, n_max*nphase);
    /* adjoint source */
    fac_all     = vector(1, nfreq*2*nrec*nshots*nphase);
    fac_r_all   = matrix(1,nphase,1, nfreq);
    fac_i_all   = matrix(1,nphase,1, nfreq);
    /* wavefield space */
    sxx_inc     = vector(1, n_max*nfreq*nshots);
    sxy_inc     = vector(1, n_max*nfreq*nshots);
    syy_inc     = vector(1, n_max*nfreq*nshots);
    sxx_i_inc   = vector(1, n_max*nfreq*nshots);
    sxy_i_inc   = vector(1, n_max*nfreq*nshots);
    syy_i_inc   = vector(1, n_max*nfreq*nshots);
    sxx_adj     = vector(1, n_max*nfreq*nrec);
    sxy_adj     = vector(1, n_max*nfreq*nrec);
    syy_adj     = vector(1, n_max*nfreq*nrec);
    sxx_i_adj   = vector(1, n_max*nfreq*nrec);
    sxy_i_adj   = vector(1, n_max*nfreq*nrec);
    syy_i_adj   = vector(1, n_max*nfreq*nrec);
    /* file handles */
    fp_kernel =(MPI_File *)malloc((size_t) (nphase*sizeof(MPI_File)));
    
    /* open file for writing kernels of all phases */
    for(iph=1;iph<=nphase;iph++){
        sprintf(str1,"%s/%s.bin",foutf, &phase_list[iph][5]); /* data/adj */
        MPI_File_open(MPI_COMM_WORLD, str1, MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fp_kernel[iph-1]);
        MPI_File_set_size(fp_kernel[iph-1],(long long)num*2*nrec*nshots*sizeof(float));
    }
    
    /**************** read the correction factors *******************/
    /* extract the factors for all frequencies from all source-receiver pairs */
    /*  index order: nfreq --> 2 --> nrec --> nshots --> nphase, means,
        storing the real (nfreq) and imaginary part (nfreq) for 1st receiver 
        of 1st source, then for 2nd receiver of 1st source, until the 
        last receiver of 1st source, repeatedly for the 2nd source until 
        the last source. Here we consider the multi. phases are picked. */
    if(MYID==0){
        for(iph=1;iph<=nphase;iph++){
            sprintf(str1,"%s/factor.bin",phase_list[iph]);
            V_read(str1, &fac_all[1+(iph-1)*nfreq*2*nrec*nshots],nfreq*2*nrec*nshots);
        }
    }
    MPI_Bcast(&fac_all[1],nphase*nfreq*2*nrec*nshots,MPI_FLOAT,0,MPI_COMM_WORLD);
    
    for(itask=0; itask<ntask; itask++){
        mynum_loc = BLOCK_SIZE(itask,ntask,mynum);
        idx_l1=idx_g1+BLOCK_LOW(itask,ntask,mynum);
        view_offset=idx_l1*sizeof(float);/* global offset in the file */
        view_offset2=(idx_l1+num)*sizeof(float);
        /* models */
        read_mod(f_vp,  Mvp,  view_offset, mynum_loc);
        read_mod(f_vs,  Mvs,  view_offset, mynum_loc);
        read_mod(f_rho, Mrho, view_offset, mynum_loc);
        
        /******************* Prepare the wavefields *********************/
        /* read wavefields from source side */
        /* the wavefield file is in the order of (num, nfreq), the 
         * wavefield array is in the order of (mynum_loc, nfreq, nshots)*/
        for(ii=0;ii<nshots;ii++){
            sprintf(finc,"%s/%s",fincf,event_list[ii]);
            offset_s=ii*mynum_loc*nfreq;/* offset in the array */
            read_wfd(finc, &sxx_inc[offset_s], &sxx_i_inc[offset_s], 
                     &sxy_inc[offset_s], &sxy_i_inc[offset_s], 
                     &syy_inc[offset_s], &syy_i_inc[offset_s], 
                     nfreq, num, mynum_loc, view_offset);
            /* convert stress to strain tensor */
            stress2strain2D(&sxx_inc[offset_s], &sxx_i_inc[offset_s], 
                            &sxy_inc[offset_s], &sxy_i_inc[offset_s], 
                            &syy_inc[offset_s], &syy_i_inc[offset_s], 
                            Mvp, Mvs, Mrho, mynum_loc, nfreq);
        }
        
        /* read wavefields from receiver side */
        /* the wavefield file is in the order of (num, nfreq), the 
         * wavefield array is in the order of (mynum_loc, nfreq, nrec)*/
        for(jj=0;jj<nrec;jj++){
            sprintf(fadj,"%s/%s",fadjf,station_list[jj]);
            offset_r=jj*mynum_loc*nfreq;/* offset in the array */
            read_wfd(fadj, &sxx_adj[offset_r], &sxx_i_adj[offset_r], 
                     &sxy_adj[offset_r], &sxy_i_adj[offset_r], 
                     &syy_adj[offset_r], &syy_i_adj[offset_r], 
                     nfreq, num, mynum_loc, view_offset);
            /* convert stress to strain tensor */
            stress2strain2D(&sxx_adj[offset_r], &sxx_i_adj[offset_r], 
                            &sxy_adj[offset_r], &sxy_i_adj[offset_r], 
                            &syy_adj[offset_r], &syy_i_adj[offset_r], 
                            Mvp, Mvs, Mrho, mynum_loc, nfreq);
        }
        if(verbose&&MYID==0){
            printf("Finished loading the wavefield in itask: %d\n", itask);
        }
        /******************* Kernel calculation *********************/
        /* assemble the kernel for the source-receiver pair one by one */
        for(ii=0;ii<nshots;ii++){
            ic=(long long)ii*nfreq*2*nrec;
            for(jj=0;jj<nrec;jj++){
                jc=(long long)jj*nfreq*2;
                /*extract the factor for current source-receiver pair */
                for(iph=1;iph<=nphase;iph++){
                    ip=(long long)(iph-1)*nfreq*2*nrec*nshots;
                    for(it=1;it<=nfreq;it++){
                        /* OLD fac_all: (nfreq, 2, nrec, nshots, nphase) */
                        /*     index    (it,   1/2, jj,   ii,     iph) */
                        /* NEW fac_r_all:  (nfreq, nphase)--->(it, iph)*/
                        idx=it+jc+ic+ip;/* only for real part in the file */
                        fac_r_all[iph][it]=fac_all[idx];
                        fac_i_all[iph][it]=fac_all[idx+nfreq];
                    }
                }
                /*3. calculate the kernels */
                /* the wavefield array is in the order as (mynum_loc, nfreq, nshots)*/
                offset_s=ii*mynum_loc*nfreq;
                offset_r=jj*mynum_loc*nfreq;
                cal_kernel(mynum_loc, nfreq, nlen, nphase, fac_r_all, fac_i_all,
                           &sxx_inc[offset_s], &sxx_i_inc[offset_s], 
                           &sxy_inc[offset_s], &sxy_i_inc[offset_s], 
                           &syy_inc[offset_s], &syy_i_inc[offset_s], 
                           &sxx_adj[offset_r], &sxx_i_adj[offset_r], 
                           &sxy_adj[offset_r], &sxy_i_adj[offset_r], 
                           &syy_adj[offset_r], &syy_i_adj[offset_r], 
                           Mvp, Mvs, Mrho, grad_vp, grad_vs);
                
                /*4. write gradient vp */
                for(iph=1;iph<=nphase;iph++){
                    offset_p=(iph-1)*mynum_loc;
//                     if(flag_list[iph]==2){/* only gradient vp */
//                         for(im=1;im<=mynum_loc;im++){
//                             grad_vs[offset_p+im]=0.0;
//                         }
//                     }else if(flag_list[iph]==3){/* only gradient vs */
//                         for(im=1;im<=mynum_loc;im++){
//                             grad_vp[offset_p+im]=0.0;
//                         }
//                     }
                    ///////////////
                    /* the saved file is as (num, 2, nrec, nshots), 
                       and the index is as  (im, 1/2, jj,  ii) */
                    ibase=((long long)jj*num*2+(long long)ii*num*2*nrec)*sizeof(float);
                    view_offset_g=view_offset+ibase; /* only for vp */
                    view_offset_g2=view_offset2+ibase;
                    write_mod3(fp_kernel[iph-1], &(grad_vp[offset_p]), &(grad_vs[offset_p]), 
                               view_offset_g, view_offset_g2, mynum_loc);
                    //MPI_Barrier(MPI_COMM_WORLD);
                }
            }/* nrec */
            if(verbose&&MYID==0){ 
                printf("Finished writing in shot %d th of %d\n", ii, nshots);
            }
        }/* nshots */
        if(verbose&&MYID==0){
            printf("Finished kernel calculation (itask): %d\n", itask);
        }
    }/* ntask */
    
    /* close file for writing kernels of all phases */
    for(iph=1;iph<=nphase;iph++){
        MPI_File_close(&fp_kernel[iph-1]);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    /* deallocate the memory */
    /* model space */
    free_vector(Mvp     ,1, n_max);
    free_vector(Mvs     ,1, n_max);
    free_vector(Mrho    ,1, n_max);
    free_vector(grad_vp ,1, n_max*nphase);
    free_vector(grad_vs ,1, n_max*nphase);
    /* adjoint source */
    free_vector(fac_all     ,1, nfreq*2*nrec*nshots*nphase);
    free_matrix(fac_r_all   ,1,nphase,1, nfreq);
    free_matrix(fac_i_all   ,1,nphase,1, nfreq);
    /* wavefield space */
    free_vector(sxx_inc     ,1, n_max*nfreq*nshots);
    free_vector(sxy_inc     ,1, n_max*nfreq*nshots);
    free_vector(syy_inc     ,1, n_max*nfreq*nshots);
    free_vector(sxx_i_inc   ,1, n_max*nfreq*nshots);
    free_vector(sxy_i_inc   ,1, n_max*nfreq*nshots);
    free_vector(syy_i_inc   ,1, n_max*nfreq*nshots);
    free_vector(sxx_adj     ,1, n_max*nfreq*nrec);
    free_vector(sxy_adj     ,1, n_max*nfreq*nrec);
    free_vector(syy_adj     ,1, n_max*nfreq*nrec);
    free_vector(sxx_i_adj   ,1, n_max*nfreq*nrec);
    free_vector(sxy_i_adj   ,1, n_max*nfreq*nrec);
    free_vector(syy_i_adj   ,1, n_max*nfreq*nrec);
    /* file handles */
    free(fp_kernel);
    
    MPI_Finalize();
    return 0;
}

/* subroutines */
void check_file_exist(char **event_list, char **station_list, char *fincf, 
                      char *fadjf, char *f_vp, char *f_vs, char *f_rho, 
                      char **phase_list, int nshots, int nrec, int nphase){
    
    /* Local variables */
    int ii, jj, iph;
    char finc[STRING_SIZE], fadj[STRING_SIZE], str1[STRING_SIZE];
    FILE *fp= NULL;
    /* vp */
    fp = fopen(f_vp, "r");
    if(fp == NULL) {
        printf("%s:", f_vp); 
        throw_error("Unable to open file!\n");
    }else 
        fclose(fp);
    /* vs */
    fp = fopen(f_vs, "r");
    if(fp == NULL) {
        printf("%s:", f_vs); 
        throw_error("Unable to open file!\n");
    }else 
        fclose(fp);
    /* rho */
    fp = fopen(f_rho, "r");
    if(fp == NULL) {
        printf("%s:", f_rho); 
        throw_error("Unable to open file!\n");
    }else 
        fclose(fp);
    /* factor */
    for(iph=1;iph<=nphase;iph++){
        sprintf(str1,"%s/factor.bin",phase_list[iph]);
        fp = fopen(str1, "r");
        if(fp == NULL) {
            printf("%s:", str1); 
            throw_error("Unable to open file!\n");
        }else 
            fclose(fp);
    }
    /* source wavefield */
    for(ii=0;ii<nshots;ii++){
        sprintf(finc,"%s/%s",fincf,event_list[ii]);
        sprintf(str1,"%s.sxx.bin",finc);
        fp = fopen(str1, "r");
        if(fp == NULL) {
            printf("%s:", str1); 
            throw_error("Unable to open file!\n");
        }else 
            fclose(fp);
    }
    /* receiver wavefield */
    for(jj=0;jj<nrec;jj++){
        sprintf(fadj,"%s/%s",fadjf,station_list[jj]);
        sprintf(str1,"%s.sxx.bin",fadj);
        fp = fopen(str1, "r");
        if(fp == NULL) {
            printf("%s:", str1); 
            throw_error("Unable to open file!\n");
        }else 
            fclose(fp);
    }                     
}

void read_wfd(char *finc, float *sxx_inc, float *sxx_i_inc, float *sxy_inc, 
              float *sxy_i_inc, float *syy_inc, float *syy_i_inc,
              int nfreq, int num, int mynum, int view_offset){
    
    /* Local variables */
    char f_sxx_inc[STRING_SIZE], f_sxy_inc[STRING_SIZE], f_syy_inc[STRING_SIZE];
    char f_sxx_i_inc[STRING_SIZE], f_sxy_i_inc[STRING_SIZE], f_syy_i_inc[STRING_SIZE];

    sprintf(f_sxx_inc,"%s.sxx.bin",finc);
    sprintf(f_sxy_inc,"%s.sxy.bin",finc);
    sprintf(f_syy_inc,"%s.syy.bin",finc);
    sprintf(f_sxx_i_inc,"%s.sxx.i.bin",finc);
    sprintf(f_sxy_i_inc,"%s.sxy.i.bin",finc);
    sprintf(f_syy_i_inc,"%s.syy.i.bin",finc);
    read_wfd_one(f_sxx_inc,       sxx_inc,   nfreq, num, mynum, view_offset);
    read_wfd_one(f_sxx_i_inc,     sxx_i_inc, nfreq, num, mynum, view_offset);
    read_wfd_one(f_syy_inc,       syy_inc,   nfreq, num, mynum, view_offset);
    read_wfd_one(f_syy_i_inc,     syy_i_inc, nfreq, num, mynum, view_offset);
    read_wfd_one(f_sxy_inc,       sxy_inc,   nfreq, num, mynum, view_offset);
    read_wfd_one(f_sxy_i_inc,     sxy_i_inc, nfreq, num, mynum, view_offset);
}

void read_wfd_one(char *f_sp_inc, float *sp_inc, int nfreq, int num, 
                  int mynum, int view_offset){
    
    /* Local variables */
    int view_offset2, it, save_offset;
    MPI_File mpi_fh;
    MPI_Status status;
    
    /* real part */
    MPI_File_open(MPI_COMM_WORLD,f_sp_inc, MPI_MODE_RDONLY,MPI_INFO_NULL,&mpi_fh);
    for(it=1;it<=nfreq;it++){
        view_offset2=(it-1)*num*sizeof(float)+view_offset;
        save_offset=(it-1)*mynum+1;
        /* real part */
        MPI_File_read_at_all(mpi_fh, view_offset2, &sp_inc[save_offset], mynum, MPI_FLOAT, &status);
    }
    /* close the files */
    MPI_File_close(&mpi_fh);
}

void read_mod(char *filename, float *mod, int view_offset, int mynum){
    
    /* Local variables */
    //int rc;
    MPI_File mpi_fh;
    MPI_Status status;
    
    MPI_File_open(MPI_COMM_WORLD,filename, MPI_MODE_RDONLY,MPI_INFO_NULL,&mpi_fh);
    MPI_File_read_at_all(mpi_fh, view_offset, &mod[1], mynum, MPI_FLOAT, &status);
    MPI_File_close(&mpi_fh);
    //MPI_Barrier(MPI_COMM_WORLD);
}

// void write_mod(char *filename, float *mod, int view_offset, int mynum){
//     
//     /* Local variables */
//     MPI_Status status;
//     MPI_File mpi_fh;
//     
//     MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&mpi_fh);
//     MPI_File_set_view(mpi_fh, view_offset, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL); 
//     MPI_File_write_all(mpi_fh, &mod[1], mynum, MPI_FLOAT, &status);
//     MPI_File_close(&mpi_fh);
//     //MPI_Barrier(MPI_COMM_WORLD);
//     
// }

// void write_mod2(char *filename, float *mod, float *mod2, int view_offset, 
//                 int view_offset2, int mynum){
//     
//     /* Local variables */
//     MPI_Status status;
//     MPI_File mpi_fh;
//     
//     MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&mpi_fh);
//     /* vp */
//     MPI_File_set_view(mpi_fh, view_offset, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL); 
//     MPI_File_write_all(mpi_fh, &mod[1], mynum, MPI_FLOAT, &status);
//     /* vs */
//     MPI_File_set_view(mpi_fh, view_offset2, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL); 
//     MPI_File_write_all(mpi_fh, &mod2[1], mynum, MPI_FLOAT, &status);
//     MPI_File_close(&mpi_fh);
//     //MPI_Barrier(MPI_COMM_WORLD);
//     
// }

void write_mod2(char *filename, float *mod, float *mod2, int view_offset, 
                int view_offset2, int mynum){
    
    /* Local variables */
    MPI_Status status;
    MPI_File mpi_fh;
    
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&mpi_fh);
    /* vp */
    //MPI_File_set_view(mpi_fh, view_offset, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL); MPI_STATUS_IGNORE
    MPI_File_write_at_all(mpi_fh, view_offset, &mod[1], mynum, MPI_FLOAT, &status);
    /* vs */
    //MPI_File_set_view(mpi_fh, view_offset2, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL); 
    MPI_File_write_at_all(mpi_fh, view_offset2, &mod2[1], mynum, MPI_FLOAT, &status);
    MPI_File_close(&mpi_fh);
    //MPI_Barrier(MPI_COMM_WORLD);
    
}

void write_mod3(MPI_File mpi_fh, float *mod, float *mod2, long long view_offset_g,
                long long view_offset_g2, int mynum){
    
    /* Local variables */
    
    /* vp */
    MPI_File_write_at_all(mpi_fh, view_offset_g, &mod[1], mynum, MPI_FLOAT, MPI_STATUS_IGNORE);
    /* vs */
    MPI_File_write_at_all(mpi_fh, view_offset_g2, &mod2[1], mynum, MPI_FLOAT, MPI_STATUS_IGNORE);
    //MPI_Barrier(MPI_COMM_WORLD);
}

void stress2strain2D(float *sxx_inc, float *sxx_i_inc, float *sxy_inc, 
                   float *sxy_i_inc, float *syy_inc, float *syy_i_inc,
                   float *Mvp, float *Mvs, float *Mrho, int mynum, int nfreq){
    /* modified from Zhang Wei's codes FD3DspherEw */
    
    /* Local variables */
    int i, j, it;
    float lam, mu, E0, E2, E3;//, E1
    
    for(i=1;i<=mynum;i++){
        if(Mvs[i]>0.1){/* means vs is not zero */
            mu =Mrho[i]*Mvs[i]*Mvs[i];
            lam =Mrho[i]*Mvp[i]*Mvp[i]-2.0*mu;
            E2=-lam/(4.0*mu*(lam+mu));
            E3=0.5/mu;
        }else{/* means vs is close to zero */
            lam =Mrho[i]*Mvp[i]*Mvp[i];
            E2=0.25/lam; /* here set the trace(strain) is same */
            E3=0.0;
        }
        for(it=1;it<=nfreq;it++){
            j=(it-1)*mynum+i;
            /* sxx, syy, sxy */
            E0=E2*(sxx_inc[j]+syy_inc[j]);
            sxx_inc[j] =E0+E3*sxx_inc[j];
            syy_inc[j] =E0+E3*syy_inc[j];
            sxy_inc[j] =E3*sxy_inc[j];
            /* sxx_i, syy_i, sxy_i */
            E0=E2*(sxx_i_inc[j]+syy_i_inc[j]);
            sxx_i_inc[j] =E0+E3*sxx_i_inc[j];
            syy_i_inc[j] =E0+E3*syy_i_inc[j];
            sxy_i_inc[j] =E3*sxy_i_inc[j];
        }
    }
}

void cal_kernel(int mynum, int nfreq, int nlen, int nphase,
                float **fac_r_all, float **fac_i_all, float *sxx_inc, 
                float *sxx_i_inc, float *sxy_inc, float *sxy_i_inc, 
                float *syy_inc, float *syy_i_inc, float *sxx_adj, 
                float *sxx_i_adj, float *sxy_adj, float *sxy_i_adj, 
                float *syy_adj, float *syy_i_adj, float *Mvp, float *Mvs, 
                float *Mrho, float *grad_vp, float *grad_vs){
    
    /* Local variables */
    int i, j, it, iph, ii;
    float corr_r, corr_i, corr_r2, corr_i2, grad1, grad2;
    
    /* initialize the memory */
    for(ii=1;ii<=mynum*nphase;ii++){
        grad_vp[ii]=0.0;
        grad_vs[ii]=0.0;
    }
    
    /*3. calculate the correlation */
    for(it=1;it<=nfreq;it++){
        /* remove the compoments with small amplitude */
        
        /* calculate the zero-lag cross-correlation */
        for(i=1;i<=mynum;i++){
            j=(it-1)*mynum+i;
            corr_r =(sxx_inc[j]+syy_inc[j])*(sxx_adj[j]+syy_adj[j])
                   -(sxx_i_inc[j]+syy_i_inc[j])*(sxx_i_adj[j]+syy_i_adj[j]);
            corr_i =(sxx_inc[j]+syy_inc[j])*(sxx_i_adj[j]+syy_i_adj[j])
                   +(sxx_i_inc[j]+syy_i_inc[j])*(sxx_adj[j]+syy_adj[j]);
            corr_r2 =2.0*(sxx_inc[j]*sxx_adj[j]+syy_inc[j]*syy_adj[j])
                -2.0*(sxx_i_inc[j]*sxx_i_adj[j]+syy_i_inc[j]*syy_i_adj[j])
                +4.0*(sxy_inc[j]*sxy_adj[j])
                -4.0*(sxy_i_inc[j]*sxy_i_adj[j]);
            corr_i2 =2.0*(sxx_inc[j]*sxx_i_adj[j]+syy_inc[j]*syy_i_adj[j])
                +2.0*(sxx_i_inc[j]*sxx_adj[j]+syy_i_inc[j]*syy_adj[j])
                +4.0*(sxy_inc[j]*sxy_i_adj[j])
                +4.0*(sxy_i_inc[j]*sxy_adj[j]);
            for(iph=1;iph<=nphase;iph++){
                ii=(iph-1)*mynum+i;
                grad_vp[ii] += corr_r*fac_r_all[iph][it] - corr_i*fac_i_all[iph][it];
                grad_vs[ii] += corr_r2*fac_r_all[iph][it] - corr_i2*fac_i_all[iph][it];
            }
            
        }
    }
    /* the relative perturbation kernels */
    for(iph=1;iph<=nphase;iph++){
        for(i=1;i<=mynum;i++){
            ii=(iph-1)*mynum+i;
            grad1= 2.0*Mrho[i]*Mvp[i]*Mvp[i]*grad_vp[ii];
            grad2= 2.0*Mrho[i]*Mvs[i]*Mvs[i]*(grad_vs[ii]-2.0*grad_vp[ii]);
            grad_vp[ii]=grad1/nlen*2.0;
            grad_vs[ii]=grad2/nlen*2.0;
        }
    }
}
