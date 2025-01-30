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
@file tomo_el2d.c
@author Tao Lei (leit@sustech.edu.cn), Wei Zhang
@date 18 Sep 2019
@param[in]  verbose \b 1: print the details; \b 0: run in quiet pattern, default value;
@param[in]  inv_vs  \b 1: inverse parameter vs; \b 0: not inverse it;
@param[in]  omega1 the omega1-th index of the frequency sample;
@param[in]  omega_inc the increment of frequency sample, omega=omega1+omega_inc*(i-1);
@param[in]  fwav_s the full path of wavelet <em>(clockwise)</em> which is added at the receiver position;
@param[in]  fadj_s the actual adjoint source <em>(clockwise)</em>;
@param[in]  fmod    the full path of vp model, the path of density model 
                    is derived from vp model, the following [finc, fadj and fout_vp] are similar;
@param[in]  finc    the full path of event name (from the source side);
@param[in]  fadj    the full path of station name (from the receiver side);
@param[in,out] fout_vp the full path of vp gradient;
@code 
mpicc -o tomo_el2d tomo_el2d.c -I ../include/ -L../lib \
-L/home/tao/seis_software/SAC/lib -lIN_OUT -lsac -lsacio -lm -lfftw3

mpirun -np 4 tomo_el2d verbose=1 inv_vs=1 omega1=2 omega_inc=1 \
fwav_s=source/glob/135.sac fadj_s=data/syn_pre/40/135.vy.sac fmod=./model/inv/inv.vp \
finc=snap/inc/40 fadj=snap/adj/135 fout_vp=gradient/40.135.vp

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

/* Global variables */
int xargc; 
char **xargv;
void read_wfd(char *f_sxx_inc, char *f_sxx_i_inc, 
              char *f_sxy_inc, char *f_sxy_i_inc, 
              char *f_syy_inc, char *f_syy_i_inc, 
              float *sxx_inc, float *sxx_i_inc,
              float *sxy_inc, float *sxy_i_inc,
              float *syy_inc, float *syy_i_inc,
              int nfreq, int num, int myid, int np, int inv_vs);
void read_wfd_one(char *f_sp_inc, float *sp_inc, int nfreq, int num, 
                  int myid, int np);
void read_mod(char *filename, float *mod, int num, int myid, int np);
void write_mod(char *filename, float *mod, int num, int myid, int np);
void generate_factor(char *fwav_s, char *fadj_s, int omega1, int omega_inc, 
                     int nfreq, int myid, int *irecord, float *fac_r, float *fac_i);
void stress2strain2D(float *sxx_inc, float *sxx_i_inc, float *sxy_inc, 
                   float *sxy_i_inc, float *syy_inc, float *syy_i_inc,
                   float *Mvp, float *Mvs, float *Mrho, int mynum, int nfreq);
void cal_kernel(int mynum, int nfreq, int nlen, int *irecord, 
                float *fac_r, float *fac_i, float *sxx_inc, 
                float *sxx_i_inc, float *sxy_inc, float *sxy_i_inc, 
                float *syy_inc, float *syy_i_inc, float *sxx_adj, 
                float *sxx_i_adj, float *sxy_adj, float *sxy_i_adj, 
                float *syy_adj, float *syy_i_adj, float *Mvp, float *Mvs, 
                float *Mrho, float *grad_vp, float *grad_vs, int inv_vs);

int main(int argc, char *argv[]){
    
    /* IN and OUT */
    int verbose, inv_vs, omega1, omega_inc;
    char *fwav_s, *fadj_s, *fmod, *finc, *fadj, *fout_vp;
    
    /* Local variables */
    int size_wfd, size_wfd2, size_mod, num;
    int mynum, nfreq, nlen;
    float dt;
    int *irecord=NULL;
    float *vec, *fac_r, *fac_i, *Mvp, *Mvs, *Mrho, *grad_vp, *grad_vs;
    float *sxx_inc, *sxy_inc, *syy_inc;
    float *sxx_i_inc, *sxy_i_inc, *syy_i_inc;
    float *sxx_adj, *sxy_adj, *syy_adj;
    float *sxx_i_adj, *sxy_i_adj, *syy_i_adj;
    char f_vp[STRING_SIZE], f_vs[STRING_SIZE], f_rho[STRING_SIZE], str1[STRING_SIZE];
    char f_sxx_inc[STRING_SIZE], f_sxy_inc[STRING_SIZE], f_syy_inc[STRING_SIZE];
    char f_sxx_i_inc[STRING_SIZE], f_sxy_i_inc[STRING_SIZE], f_syy_i_inc[STRING_SIZE];
    char f_sxx_adj[STRING_SIZE], f_sxy_adj[STRING_SIZE], f_syy_adj[STRING_SIZE];
    char f_sxx_i_adj[STRING_SIZE], f_sxy_i_adj[STRING_SIZE], f_syy_i_adj[STRING_SIZE];
    /*MPI parameters */
    int NP, MYID;
    
    /* Initialize MPI environment */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &MYID);
    MPI_Comm_size(MPI_COMM_WORLD, &NP);
    
    /* 1. read formatted input list */
    /* initialize getpar */
    initargs(argc,argv);
    
    /* these lines are useful for SEISPLATFORM software */
    verbose=0; getparint("verbose",&verbose);
    inv_vs=0; getparint("inv_vs",&inv_vs);
    if(!getparint("omega1", &omega1)) throw_error("omega1 must be specified!\n");
    if(!getparint("omega_inc", &omega_inc)) throw_error("omega_inc must be specified!\n");
    /* the file or path */
    if(!getparstring("fmod", &fmod)) throw_error("fmod must be specified!\n");
    if(!getparstring("finc", &finc)) throw_error("finc must be specified!\n");
    if(!getparstring("fadj", &fadj)) throw_error("fadj must be specified!\n");
    if(!getparstring("fwav_s", &fwav_s)) throw_error("fwav_s must be specified!\n");
    if(!getparstring("fadj_s", &fadj_s)) throw_error("fadj_s must be specified!\n");
    if(!getparstring("fout_vp", &fout_vp)) throw_error("fout_vp must be specified!\n");
    MPI_Barrier(MPI_COMM_WORLD);
    
    /* check the file size */
    if(MYID==0){
        sprintf(str1,"%s.sxx.bin",finc);
        size_wfd=get_file_size(str1);
        sprintf(str1,"%s.sxx.bin",fadj);
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
        vec = read_sac(fwav_s, &nlen, &dt);
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
    
    if(verbose&&(MYID==1||NP==1)){
        printf("\n******\n");
        printf("inv_vs    \t| %d\n",inv_vs);
        printf("dt        \t| %.3e\n",dt);
        printf("num       \t| %d\n",num);
        printf("nlen      \t| %d\n",nlen);
        printf("nfreq     \t| %d\n",nfreq);
        printf("omega1    \t| %d\n",omega1);
        printf("omega_inc \t| %d\n",omega_inc);
        printf("fmod      \t| %s\n",fmod);
        printf("finc      \t| %s\n",finc);
        printf("fadj      \t| %s\n",fadj);
        printf("fwav_s    \t| %s\n",fwav_s);
        printf("fadj_s    \t| %s\n",fadj_s);
        printf("fout_vp   \t| %s\n",fout_vp);
        printf("******\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    /* allocate the memory */
    /* model space */
    Mvp     = vector(1, mynum);
    Mvs     = vector(1, mynum);
    Mrho    = vector(1, mynum);
    grad_vp = vector(1, mynum);
    grad_vs = vector(1, mynum);
    /* adjoint source */
    irecord = ivector(1,nfreq);
    fac_r   = vector(1, nfreq);
    fac_i   = vector(1, nfreq);
    /* wavefield space */
    sxx_inc     = vector(1, mynum*nfreq);
    sxy_inc     = vector(1, mynum*nfreq);
    syy_inc     = vector(1, mynum*nfreq);
    sxx_i_inc   = vector(1, mynum*nfreq);
    sxy_i_inc   = vector(1, mynum*nfreq);
    syy_i_inc   = vector(1, mynum*nfreq);
    sxx_adj     = vector(1, mynum*nfreq);
    sxy_adj     = vector(1, mynum*nfreq);
    syy_adj     = vector(1, mynum*nfreq);
    sxx_i_adj   = vector(1, mynum*nfreq);
    sxy_i_adj   = vector(1, mynum*nfreq);
    syy_i_adj   = vector(1, mynum*nfreq);
   
    /* models */
    sprintf(f_vp,"%s",fmod);
    strrpc2(f_vp,"vp","vs",f_vs);
    strrpc2(f_vp,"vp","rho",f_rho);
    read_mod(f_vp,  Mvp, num, MYID,NP);
    read_mod(f_vs,  Mvs, num, MYID,NP);
    read_mod(f_rho, Mrho, num, MYID,NP);
    
    /* calculate the factor */
    generate_factor(fwav_s, fadj_s, omega1, omega_inc, nfreq, MYID,
                    irecord, fac_r, fac_i);
    
    /* read the wavefields */
    sprintf(f_sxx_inc,"%s.sxx.bin",finc);
    sprintf(f_sxy_inc,"%s.sxy.bin",finc);
    sprintf(f_syy_inc,"%s.syy.bin",finc);
    sprintf(f_sxx_i_inc,"%s.sxx.i.bin",finc);
    sprintf(f_sxy_i_inc,"%s.sxy.i.bin",finc);
    sprintf(f_syy_i_inc,"%s.syy.i.bin",finc);
    sprintf(f_sxx_adj,"%s.sxx.bin",fadj);
    sprintf(f_sxy_adj,"%s.sxy.bin",fadj);
    sprintf(f_syy_adj,"%s.syy.bin",fadj);
    sprintf(f_sxx_i_adj,"%s.sxx.i.bin",fadj);
    sprintf(f_sxy_i_adj,"%s.sxy.i.bin",fadj);
    sprintf(f_syy_i_adj,"%s.syy.i.bin",fadj);
    read_wfd(f_sxx_inc, f_sxx_i_inc, f_sxy_inc, f_sxy_i_inc, f_syy_inc, 
             f_syy_i_inc, sxx_inc, sxx_i_inc, sxy_inc, sxy_i_inc, syy_inc, 
             syy_i_inc, nfreq, num, MYID, NP, inv_vs);
    read_wfd(f_sxx_adj, f_sxx_i_adj, f_sxy_adj, f_sxy_i_adj, f_syy_adj, 
             f_syy_i_adj, sxx_adj, sxx_i_adj, sxy_adj, sxy_i_adj, syy_adj, 
             syy_i_adj, nfreq, num, MYID, NP, inv_vs);
    
    /* convert stress to strain tensor */
    stress2strain2D(sxx_inc, sxx_i_inc, sxy_inc, sxy_i_inc, syy_inc, 
                    syy_i_inc, Mvp, Mvs, Mrho, mynum, nfreq);
    stress2strain2D(sxx_adj, sxx_i_adj, sxy_adj, sxy_i_adj, syy_adj, 
                    syy_i_adj, Mvp, Mvs, Mrho, mynum, nfreq);
    /*3. calculate the kernels */
    cal_kernel(mynum, nfreq, nlen, irecord, fac_r, fac_i, 
               sxx_inc, sxx_i_inc, sxy_inc, sxy_i_inc, syy_inc, syy_i_inc, 
               sxx_adj, sxx_i_adj, sxy_adj, sxy_i_adj, syy_adj, syy_i_adj,
               Mvp, Mvs, Mrho, grad_vp, grad_vs, inv_vs);
    
    /* write gradient vp */
    write_mod(fout_vp, grad_vp, num, MYID, NP);
    if(inv_vs){
        strrpc2(fout_vp,"vp","vs",str1);
        write_mod(str1, grad_vs, num, MYID, NP);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    /* deallocate the memory */
    /* model space */
    free_vector(Mvp     ,1, mynum);
    free_vector(Mvs     ,1, mynum);
    free_vector(Mrho    ,1, mynum);
    free_vector(grad_vp ,1, mynum);
    free_vector(grad_vs ,1, mynum);
    /* adjoint source */
    free_ivector(irecord ,1,nfreq);
    free_vector(fac_r   ,1, nfreq);
    free_vector(fac_i   ,1, nfreq);
    /* wavefield space */
    free_vector(sxx_inc     ,1, mynum*nfreq);
    free_vector(sxy_inc     ,1, mynum*nfreq);
    free_vector(syy_inc     ,1, mynum*nfreq);
    free_vector(sxx_i_inc   ,1, mynum*nfreq);
    free_vector(sxy_i_inc   ,1, mynum*nfreq);
    free_vector(syy_i_inc   ,1, mynum*nfreq);
    free_vector(sxx_adj     ,1, mynum*nfreq);
    free_vector(sxy_adj     ,1, mynum*nfreq);
    free_vector(syy_adj     ,1, mynum*nfreq);
    free_vector(sxx_i_adj   ,1, mynum*nfreq);
    free_vector(sxy_i_adj   ,1, mynum*nfreq);
    free_vector(syy_i_adj   ,1, mynum*nfreq);
    
    MPI_Finalize();
    return 0;
}

/* subroutines */
void read_wfd(char *f_sxx_inc, char *f_sxx_i_inc, 
              char *f_sxy_inc, char *f_sxy_i_inc, 
              char *f_syy_inc, char *f_syy_i_inc, 
              float *sxx_inc, float *sxx_i_inc,
              float *sxy_inc, float *sxy_i_inc,
              float *syy_inc, float *syy_i_inc,
              int nfreq, int num, int myid, int np, int inv_vs){
    //////////?????????????????????[1]
    read_wfd_one(f_sxx_inc,       sxx_inc,   nfreq, num, myid, np);
    read_wfd_one(f_sxx_i_inc,     sxx_i_inc, nfreq, num, myid, np);
    read_wfd_one(f_syy_inc,       syy_inc,   nfreq, num, myid, np);
    read_wfd_one(f_syy_i_inc,     syy_i_inc, nfreq, num, myid, np);
    if(inv_vs){
        read_wfd_one(f_sxy_inc,       sxy_inc,   nfreq, num, myid, np);
        read_wfd_one(f_sxy_i_inc,     sxy_i_inc, nfreq, num, myid, np);
    }
}

void read_wfd_one(char *f_sp_inc, float *sp_inc, int nfreq, int num, 
                  int myid, int np){
    
    /* Local variables */
    int rc, mynum, view_offset, it, save_offset;
    MPI_File fp_sp_inc;
    MPI_Status status;
    
    /* real part */
    rc = MPI_File_open(MPI_COMM_WORLD,f_sp_inc, MPI_MODE_RDONLY,MPI_INFO_NULL,&fp_sp_inc);
    if(rc){
        printf("%s:", f_sp_inc);
        throw_error("Unable to open file!\n");
    }
    mynum = BLOCK_SIZE(myid,np,num);
    for(it=1;it<=nfreq;it++){
        /* remove the compoments with small amplitude */
        //if(irecord[it]==0) continue;
        
        view_offset=((it-1)*num+BLOCK_LOW(myid,np,num))*sizeof(float);
        save_offset=(it-1)*mynum+1;
        /* real part */
        MPI_File_read_at_all(fp_sp_inc, view_offset, &sp_inc[save_offset], mynum, MPI_FLOAT, &status);
    }
    /* close the files */
    MPI_File_close(&fp_sp_inc);
}

void read_mod(char *filename, float *mod, int num, int myid, int np){
    
    /* Local variables */
    int rc, mynum, view_offset;
    MPI_File fp_vp;
    MPI_Status status;
    
    /* models */
    rc = MPI_File_open(MPI_COMM_WORLD,filename, MPI_MODE_RDONLY,MPI_INFO_NULL,&fp_vp);
    if(rc){
        printf("%s:", filename);
        throw_error("Unable to open file!\n");
    }
    /*2. read models */
    mynum = BLOCK_SIZE(myid,np,num);
    view_offset=BLOCK_LOW(myid,np,num)*sizeof(float);
    MPI_File_read_at_all(fp_vp, view_offset, &mod[1], mynum, MPI_FLOAT, &status);
    MPI_File_close(&fp_vp);
    MPI_Barrier(MPI_COMM_WORLD);
}

void write_mod(char *filename, float *mod, int num, int myid, int np){
    
    /* Local variables */
    int mynum, view_offset;
    MPI_Status status;
    MPI_File mpi_fh;
    
    mynum = BLOCK_SIZE(myid,np,num);
    view_offset=BLOCK_LOW(myid,np,num)*sizeof(float);
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&mpi_fh);
    MPI_File_set_view(mpi_fh, view_offset, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL); 
    MPI_File_write_all(mpi_fh, &mod[1], mynum, MPI_FLOAT, &status);
    MPI_File_close(&mpi_fh);
    //MPI_Barrier(MPI_COMM_WORLD);
    
}

void generate_factor(char *fwav_s, char *fadj_s, int omega1, int omega_inc, 
                     int nfreq, int myid, int *irecord, float *fac_r, float *fac_i){
    /* INPUT */
    /* 
     * fwav_s: the path of wavelet
     * fadj_s: the path of particle velocity, instead of the one of displacement
     * 
     */
    
    /* Local variables */
    int i, it, nlen, nlen2;
    float max1, max2, dt, tmp, p1=0.0;//, df, omega
    //double fac_r1, fac_i1;
    float *vec;//, *weight;
    float *wavelet_r, *wavelet_i, *wavelet_abs; 
    float *seismo_r, *seismo_i, *seismo_abs;
    
    
    /*1. read the wavelet and seismogram */
    /*a. calculate the factor for adjoint source 
    this removes the footprint of the source function and inserts 
    the ajoint source on the receiver sides */
    if(myid==0){
        /* wavelet in sac */
        vec = read_sac(fwav_s, &nlen, &dt);
        wavelet_r = vector(1, nlen);
        //for(it=2; it<=nlen; it++) wavelet_r[it]=(vec[it]-vec[it-1])/dt;
        //for(it=2; it<=nlen; it++) wavelet_r[it]=wavelet_r[it-1]+vec[it-1];
        for(it=1; it<=nlen; it++) wavelet_r[it]=vec[it-1];
        free(vec);
        /* adjoint source in sac */
        vec = read_sac(fadj_s, &nlen2, &dt);
        seismo_r = vector(1, nlen);
        for(it=1; it<=nlen; it++) seismo_r[it]=vec[it-1];
        for(it=1; it<=nlen; it++) p1 += seismo_r[it]*seismo_r[it];
        p1=p1*dt;
        free(vec);
        if(nlen!=nlen2||(omega1+omega_inc*(nfreq-1)>nlen)){
            throw_error("Please check the size of the *.su file!\n");
        }
        
        /*2. allocate memory */
        wavelet_i   = vector(1, nlen);
        wavelet_abs = vector(1, nlen);
        seismo_i   = vector(1, nlen);
        seismo_abs = vector(1, nlen);
        //weight      = vector(1, nlen);
        
        /*3. DFT */
        fft_real_imag(&wavelet_r[1], &wavelet_i[1], nlen, -1);
        fft_real_imag(&seismo_r[1], &seismo_i[1], nlen, -1);
        max1=-1.0;
        for(it=2; it<=nlen; it++){/*ignore the zero component */
            wavelet_abs[it]=sqrt(wavelet_r[it]*wavelet_r[it]+wavelet_i[it]*wavelet_i[it]);
            if(wavelet_abs[it]>max1) max1=wavelet_abs[it];
        }
        //if(verbose) printf("max(abs(wavelet))=%e\n",max1);
        if(max1==0.0) throw_error("Please check wavelet file!\n");
        
        max2=-1.0;
        for(it=2; it<=nlen; it++){/*ignore the zero component */
            seismo_abs[it]=sqrt(seismo_r[it]*seismo_r[it]+seismo_i[it]*seismo_i[it]);
            if(seismo_abs[it]>max2) max2=seismo_abs[it];
        }
        //if(verbose) printf("max(abs(wavelet))=%e\n",max2);
        if(max2==0.0) throw_error("Please check seismogram file!\n");
        
//         /* calculate the weight */
//         df=2.0*PI/(nlen*dt);
//         tmp=0.0;
//         weight[1]=0.0;
//         for(it=2; it<=nlen/2; it++){/*ignore the zero component */
//             omega=(it-1)*df;
//             weight[it]=omega*omega*seismo_abs[it]*seismo_abs[it];
//             tmp+=weight[it];
//         }
//         for(it=2; it<=nlen/2; it++){
//             weight[it]/=tmp;
//             //printf("it: %d %.3e\n", it, weight[it]);
//         }
        
        /* record the useful frequency points */
        for(it=1; it<=nfreq; it++){
        //for(it=41; it<=41; it++){
            irecord[it]=0;
            i=omega1+omega_inc*(it-1);
            //omega=(i-1)*df;
            if(wavelet_abs[i]>=0.001*max1&&seismo_abs[i]>=0.001*max2){
                irecord[it]=1;
                //printf("2irecord[it]=%d %d\n",irecord[it],i);
                tmp=-p1*wavelet_abs[i]*wavelet_abs[i];
                //tmp=weight[i]/(seismo_abs[i]*seismo_abs[i]*wavelet_abs[i]*wavelet_abs[i]);
                fac_r[it]=(wavelet_r[i]*seismo_r[i]-wavelet_i[i]*seismo_i[i])/tmp;
                fac_i[it]=-(wavelet_r[i]*seismo_i[i]+wavelet_i[i]*seismo_r[i])/tmp;
                /* same with matlab result */
            }
        }
        /* deallocate memory */
        free_vector(wavelet_r, 1, nlen);
        free_vector(wavelet_i, 1, nlen);
        free_vector(wavelet_abs, 1, nlen);
        free_vector(seismo_r, 1, nlen);
        free_vector(seismo_i, 1, nlen);
        free_vector(seismo_abs, 1, nlen);
        //free_vector(weight, 1, nlen);
    }
    /* bcast integers */
    MPI_Bcast(&irecord[1],nfreq,MPI_INT,0,MPI_COMM_WORLD);
    /* bcast floats */
    MPI_Bcast(&fac_r[1],nfreq,MPI_FLOAT,0,MPI_COMM_WORLD);
    MPI_Bcast(&fac_i[1],nfreq,MPI_FLOAT,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
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

void cal_kernel(int mynum, int nfreq, int nlen, int *irecord, 
                float *fac_r, float *fac_i, float *sxx_inc, 
                float *sxx_i_inc, float *sxy_inc, float *sxy_i_inc, 
                float *syy_inc, float *syy_i_inc, float *sxx_adj, 
                float *sxx_i_adj, float *sxy_adj, float *sxy_i_adj, 
                float *syy_adj, float *syy_i_adj, float *Mvp, float *Mvs, 
                float *Mrho, float *grad_vp, float *grad_vs, int inv_vs){
    
    /* Local variables */
    int i, j, it;
    float corr_r, corr_i, grad1, grad2;
    
    /* initialize the memory */
    for(i=1;i<=mynum;i++){
        grad_vp[i]=0.0;
        grad_vs[i]=0.0;
    }
    
    /*3. calculate the correlation */
    for(it=1;it<=nfreq;it++){
        /* remove the compoments with small amplitude */
        if(irecord[it]==0){
            continue;
        }
        
        /* calculate the zero-lag cross-correlation */
        for(i=1;i<=mynum;i++){
            j=(it-1)*mynum+i;
            corr_r =(sxx_inc[j]+syy_inc[j])*(sxx_adj[j]+syy_adj[j])
                   -(sxx_i_inc[j]+syy_i_inc[j])*(sxx_i_adj[j]+syy_i_adj[j]);
            corr_i =(sxx_inc[j]+syy_inc[j])*(sxx_i_adj[j]+syy_i_adj[j])
                   +(sxx_i_inc[j]+syy_i_inc[j])*(sxx_adj[j]+syy_adj[j]);
            grad_vp[i] += corr_r*fac_r[it] - corr_i*fac_i[it];
        }
        if(inv_vs){
          for(i=1;i<=mynum;i++){
            j=(it-1)*mynum+i;
            corr_r =2.0*(sxx_inc[j]*sxx_adj[j]+syy_inc[j]*syy_adj[j])
                -2.0*(sxx_i_inc[j]*sxx_i_adj[j]+syy_i_inc[j]*syy_i_adj[j])
                +4.0*(sxy_inc[j]*sxy_adj[j])
                -4.0*(sxy_i_inc[j]*sxy_i_adj[j]);
            corr_i =2.0*(sxx_inc[j]*sxx_i_adj[j]+syy_inc[j]*syy_i_adj[j])
                +2.0*(sxx_i_inc[j]*sxx_adj[j]+syy_i_inc[j]*syy_adj[j])
                +4.0*(sxy_inc[j]*sxy_i_adj[j])
                +4.0*(sxy_i_inc[j]*sxy_adj[j]);
            grad_vs[i] += corr_r*fac_r[it] - corr_i*fac_i[it];
          }
        }
    }
    /* the relative perturbation kernels */
    for(i=1;i<=mynum;i++){
        grad1= 2.0*Mrho[i]*Mvp[i]*Mvp[i]*grad_vp[i];
        grad2= 2.0*Mrho[i]*Mvs[i]*Mvs[i]*(grad_vs[i]-2.0*grad_vp[i]);
        grad_vp[i]=grad1/nlen*2.0;
        grad_vs[i]=grad2/nlen*2.0;
    }
}