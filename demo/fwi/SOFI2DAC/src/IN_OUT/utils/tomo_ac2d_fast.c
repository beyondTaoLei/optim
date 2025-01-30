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
@file tomo_ac2d_fast.c
@author Tao Lei (leit@sustech.edu.cn), Wei Zhang
@date 18 Sep 2019
@param[in]  verbose \b 1: print the details; \b 0: run in quiet pattern, default value;
@param[in]  omega1 the omega1-th index of the frequency sample;
@param[in]  omega_inc the increment of frequency sample, omega=omega1+omega_inc*(i-1);
@param[in]  nincs the shot number in one group;
@param[in]  nincr the receiver number in one group;
@param[in]  fmod    the full path of vp model, the path of density model 
                    is derived from vp model;
@param[in]  fevents the file of event name list;
@param[in]  fstations the file of station name list;
@param[in]  fincf the folder of snapshots from the source side;
@param[in]  fadjf the folder of snapshots from the receiver side;
@param[in]  fwavf the folder of wavelet <em>(clockwise)</em> which is added at the receiver position;
@param[in]  fadj_sf the folder of actual adjoint source <em>(clockwise)</em>;
@param[in,out] foutf the folder of vp gradient;
@code 
mpicc -o tomo_ac2d_fast tomo_ac2d_fast.c -I ../include/ -L../lib \
-L/home/tao/seis_software/SAC/lib -lIN_OUT -lsac -lsacio -lm -lfftw3

mpirun -np 60 tomo_ac2d_fast omega1=2 omega_inc=1 \
nincs=95 nincr=95 fmod=./model/inv/inv.vp fevents=list/eventname.list \
fstations=list/stationname.list fincf=snap/inc fadjf=snap/adj \
fwavf=source/glob fadj_sf=data/adj foutf=grad verbose=1

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

/* Global variables */
int xargc; 
char **xargv;

void read_wfd(char *f_sp_inc, char *f_sp_i_inc, float *sp_inc, float *sp_i_inc, 
              int nfreq, int num, int myid, int np);
void read_wfd_one(char *f_sp_inc, float *sp_inc, int nfreq, int num, 
                  int myid, int np);
void read_mod(char *filename, float *mod, int num, int myid, int np);
void write_mod(char *filename, float *mod, int num, int myid, int np);
void generate_factor(char *fwav_s, char *fadj_s, int omega1, int omega_inc, 
                     int nfreq, int myid, int *irecord, float *fac_r, float *fac_i);
void cal_kernel(int mynum, int nfreq, int nlen, int *irecord, float *fac_r,
                float *fac_i, float *sp_inc, float *sp_i_inc, float *sp_adj, 
                float *sp_i_adj, float *Mvp, float *csp);

int main(int argc, char *argv[]){

    /* IN and OUT */
    int verbose, omega1, omega_inc, nincs, nincr;
    char *fmod, *fevents, *fstations, *fincf, *fadjf, *fwavf, *fadj_sf, *foutf;
    
    /* Local variables */
    int size_wfd, size_wfd2, size_mod, num, it;
    int mynum, nfreq, nlen;
    int nshots, nrec, ngroup_s, ngroup_r;
    int ii, jj, kk, igs, igr, istart_s, istart_r, num_s, num_r, offset_s, offset_r;
    float dt;
    int *irecord, **irecord_all;
    float *vec, *Mvp, *csp;
    float *fac_r, *fac_i, **fac_r_all, **fac_i_all;
    float *sp_inc, *sp_adj, *sp_i_inc, *sp_i_adj;
    char **event_list = NULL, **station_list =NULL;
    char **finc_list, **fadj_list, **event_list_cur, **station_list_cur;
    char **wavelet_list_cur, **adj_src_list_cur, **kername_list_cur;
    char f_sp_inc[STRING_SIZE], f_sp_i_inc[STRING_SIZE];
    char f_sp_adj[STRING_SIZE], f_sp_i_adj[STRING_SIZE];
    char str1[STRING_SIZE];
    FILE *fp=NULL;
    /*MPI parameters */
    int NP,MYID;
    
    /* Initialize MPI environment */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &MYID);
    MPI_Comm_size(MPI_COMM_WORLD, &NP);
    
    /* 1. read formatted input list */
    /* initialize getpar */
    initargs(argc,argv);
    
    /* these lines are useful for SEISPLATFORM software */
    verbose=0; getparint("verbose",&verbose);
    if(!getparint("omega1", &omega1)) throw_error("omega1 must be specified!\n");
    if(!getparint("omega_inc", &omega_inc)) throw_error("omega_inc must be specified!\n");
    if(!getparint("nincs", &nincs)) throw_error("nincs must be specified!\n");
    if(!getparint("nincr", &nincr)) throw_error("nincr must be specified!\n");
    /* the file or path */
    if(!getparstring("fmod", &fmod)) throw_error("fmod must be specified!\n");
    if(!getparstring("fevents", &fevents)) throw_error("fevents must be specified!\n");
    if(!getparstring("fstations", &fstations)) throw_error("fstations must be specified!\n");
    if(!getparstring("fincf", &fincf)) throw_error("fincf must be specified!\n");
    if(!getparstring("fadjf", &fadjf)) throw_error("fadjf must be specified!\n");
    if(!getparstring("fwavf", &fwavf)) throw_error("fwavf must be specified!\n");
    if(!getparstring("fadj_sf", &fadj_sf)) throw_error("fadj_sf must be specified!\n");
    if(!getparstring("foutf", &foutf)) throw_error("foutf must be specified!\n");
    MPI_Barrier(MPI_COMM_WORLD);
    
    /* check the file size */
    event_list = read_string(&nshots, fevents);
    station_list = read_string(&nrec, fstations);
    if(MYID==0){
        sprintf(str1,"%s/%s.sp.bin",fincf, event_list[0]);
        size_wfd=get_file_size(str1);
        sprintf(str1,"%s/%s.sp.bin",fadjf, station_list[0]);
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
    ngroup_s=(int)ceil(1.0*nshots/nincs);
    ngroup_r=(int)ceil(1.0*nrec/nincr);
    
    if(verbose&&(MYID==1||NP==1)){
        printf("\n******\n");
        printf("dt        \t| %.3e\n",dt);
        printf("num       \t| %d\n",num);
        printf("nlen      \t| %d\n",nlen);
        printf("nfreq     \t| %d\n",nfreq);
        printf("omega1    \t| %d\n",omega1);
        printf("omega_inc \t| %d\n",omega_inc);
        printf("nshots    \t| %d\n",nshots);
        printf("nrec      \t| %d\n",nrec);
        printf("nincs     \t| %d\n",nincs);
        printf("nincr     \t| %d\n",nincr);
        printf("ngroup_s  \t| %d\n",ngroup_s);
        printf("ngroup_r  \t| %d\n",ngroup_r);
        printf("fmod      \t| %s\n",fmod);
        printf("fevents   \t| %s\n",fevents);
        printf("fstations \t| %s\n",fstations);
        printf("fincf     \t| %s\n",fincf);
        printf("fadjf     \t| %s\n",fadjf);
        printf("fwavf     \t| %s\n",fwavf);
        printf("fadj_sf   \t| %s\n",fadj_sf);
        printf("foutf     \t| %s\n",foutf);
        printf("******\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    /* allocate the memory */
    /* model space */
    Mvp    = vector(1, mynum);
    csp    = vector(1, mynum);
    /* adjoint source */
    irecord = ivector(1,nfreq);
    fac_r   = vector(1, nfreq);
    fac_i   = vector(1, nfreq);
    irecord_all = imatrix(1,nincs*nincr,1,nfreq);
    fac_r_all   = matrix(1,nincs*nincr,1, nfreq);
    fac_i_all   = matrix(1,nincs*nincr,1, nfreq);
    /* wavefield space */
    sp_inc      = vector(1, mynum*nfreq*nincs);
    sp_i_inc    = vector(1, mynum*nfreq*nincs);
    sp_adj      = vector(1, mynum*nfreq*nincr);
    sp_i_adj    = vector(1, mynum*nfreq*nincr);
    /* kernel list space */
    event_list_cur      = chmatrix(0, nincs-1, 0, 255);
    finc_list           = chmatrix(0, nincs-1, 0, 255);
    station_list_cur    = chmatrix(0, nincr-1, 0, 255);
    fadj_list           = chmatrix(0, nincr-1, 0, 255);
    wavelet_list_cur    = chmatrix(0, nincr-1, 0, 255);
    adj_src_list_cur    = chmatrix(0, nincs*nincr-1, 0, 255);
    kername_list_cur    = chmatrix(0, nincs*nincr-1, 0, 255);
    
    /* models */
    read_mod(fmod,  Mvp, num, MYID,NP);
    
    /* source group */
    for(igs=0;igs<ngroup_s;igs++){
        istart_s=BLOCK_LOW(igs,ngroup_s,nshots); //iend_s  =BLOCK_HIGH(igs,ngroup_s,nshots);
        num_s   =BLOCK_SIZE(igs,ngroup_s,nshots);//num_s not more than nincs
        /* initialize current source list */
        for(ii=0;ii<num_s;ii++){
            sprintf(event_list_cur[ii],"%s",event_list[istart_s+ii]);
            sprintf(finc_list[ii],"%s/%s",fincf,event_list_cur[ii]);
            //if(MYID==0)printf("%s %s\n",event_list_cur[ii],finc_list[ii]);
        }
        
        /* read wavefields from source side */
        for(ii=0;ii<num_s;ii++){
            /*wavefields: the seeds are f_sp_inc and f_sp_adj*/
            sprintf(f_sp_inc,"%s.sp.bin",finc_list[ii]);
            sprintf(f_sp_i_inc,"%s.sp.i.bin",finc_list[ii]);
            offset_s=ii*mynum*nfreq;
            read_wfd(f_sp_inc, f_sp_i_inc, &sp_inc[offset_s], &sp_i_inc[offset_s], nfreq, num, MYID, NP);
        }
            
        /* receiver group */
        for(igr=0;igr<ngroup_r;igr++){
            istart_r=BLOCK_LOW(igr,ngroup_r,nrec);//iend_r  =BLOCK_HIGH(igr,ngroup_r,nrec);
            num_r   =BLOCK_SIZE(igr,ngroup_r,nrec); //num_r not more than nincr
            if(MYID==0)printf("%d %d %d %d %d %d\n",igs, istart_s, num_s, igr, istart_r, num_r);
            /* initialize current receiver list */
            for(jj=0;jj<num_r;jj++){
                sprintf(station_list_cur[jj],"%s",station_list[istart_r+jj]);
                sprintf(fadj_list[jj],"%s/%s",fadjf,station_list_cur[jj]);
                /* wavelet */
                sprintf(wavelet_list_cur[jj],"%s/%s.sac",fwavf,station_list[istart_r+jj]);
                //if(MYID==0)printf("%s %s %s\n",station_list_cur[jj], fadj_list[jj],wavelet_list_cur[jj]);
                /* adjoint source and kernel name*/
                for(ii=0;ii<num_s;ii++){
                    kk=ii*num_r+jj;
                    sprintf(adj_src_list_cur[kk],"%s/%s/%s.p.sac",fadj_sf,event_list_cur[ii],station_list_cur[jj]);
                    sprintf(kername_list_cur[kk],"%s/%s.%s",foutf,event_list_cur[ii],station_list_cur[jj]);
                    //if(MYID==0)printf("%d %s %s\n",kk,adj_src_list_cur[kk], kername_list_cur[kk]);
                }
                //if(MYID==0)printf("\n\n");
            }
            
            /* calculate the factor conj(F/conj(S)) in the kernel formula */
            for(ii=0;ii<num_s;ii++){
            for(jj=0;jj<num_r;jj++){
                kk=ii*num_r+jj;
                fp=fopen(adj_src_list_cur[kk],"r");
                if (fp==NULL){
                    if(verbose && MYID==0) printf("No adjoint source: %s!\n",adj_src_list_cur[kk]);
                    for(it=1; it<=nfreq; it++) irecord_all[kk+1][it]=0;
                    continue;
                }else{
                    fclose(fp);
                    /* calculate the factor */
                     /*1. prepare the files for acoustic case */
                    /*a. calculate the factor for adjoint source 
                    this removes the footprint of the source function and inserts 
                    the ajoint source on the receiver sides */
                    generate_factor(wavelet_list_cur[jj], adj_src_list_cur[kk], 
                                    omega1, omega_inc, nfreq, MYID, irecord, fac_r, fac_i);
                    for(it=1; it<=nfreq; it++){
                        irecord_all[kk+1][it] = irecord[it];
                        fac_r_all[kk+1][it]   = fac_r[it];
                        fac_i_all[kk+1][it]   = fac_i[it];
                    }
                }
            }}

            /* read wavefields from receiver side */
            for(jj=0;jj<num_r;jj++){
                sprintf(f_sp_adj,"%s.sp.bin",fadj_list[jj]);
                sprintf(f_sp_i_adj,"%s.sp.i.bin",fadj_list[jj]);
                offset_r=jj*mynum*nfreq;
                read_wfd(f_sp_adj, f_sp_i_adj, &sp_adj[offset_r], &sp_i_adj[offset_r], nfreq, num, MYID, NP);
            }
            
            /* start to assemble the kernel from source side (num_s) 
             * and receiver side (num_r) one by one */
            for(ii=0;ii<num_s;ii++){
            for(jj=0;jj<num_r;jj++){
                kk=ii*num_r+jj;
                offset_s=ii*mynum*nfreq;
                offset_r=jj*mynum*nfreq;
                /////////////////////////
                for(it=1; it<=nfreq; it++){
                    irecord[it] = irecord_all[kk+1][it];
                    fac_r[it]   = fac_r_all[kk+1][it];
                    fac_i[it]   = fac_i_all[kk+1][it];
                }
                cal_kernel(mynum, nfreq, nlen, irecord, fac_r, fac_i, 
                           &sp_inc[offset_s], &sp_i_inc[offset_s], 
                           &sp_adj[offset_r], &sp_i_adj[offset_r], 
                           Mvp, csp);
                /* write gradient vp */
                sprintf(str1,"%s.vp.bin",kername_list_cur[kk]);
                write_mod(str1, csp, num, MYID, NP);
                //MPI_Barrier(MPI_COMM_WORLD);
            }}
            MPI_Barrier(MPI_COMM_WORLD);
            //if(MYID==0)printf("\n\n");
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    /* deallocate the memory */
    /* model space */
    free_vector(Mvp    ,1, mynum);
    free_vector(csp    ,1, mynum);
    /* adjoint source */
    free_ivector(irecord ,1,nfreq);
    free_vector(fac_r   ,1, nfreq);
    free_vector(fac_i   ,1, nfreq);
    free_imatrix(irecord_all ,1,nincs*nincr,1,nfreq);
    free_matrix(fac_r_all   ,1,nincs*nincr,1, nfreq);
    free_matrix(fac_i_all   ,1,nincs*nincr,1, nfreq);
    /*  wavefield space */
    free_vector(sp_inc      ,1, mynum*nfreq*nincs);
    free_vector(sp_i_inc    ,1, mynum*nfreq*nincs);
    free_vector(sp_adj      ,1, mynum*nfreq*nincr);
    free_vector(sp_i_adj    ,1, mynum*nfreq*nincr);
    /* kernel list space */
    free_chmatrix(event_list_cur      ,0, nincs-1, 0, 255);
    free_chmatrix(finc_list           ,0, nincs-1, 0, 255);
    free_chmatrix(station_list_cur    ,0, nincr-1, 0, 255);
    free_chmatrix(fadj_list           ,0, nincr-1, 0, 255);
    free_chmatrix(wavelet_list_cur    ,0, nincr-1, 0, 255);
    free_chmatrix(adj_src_list_cur    ,0, nincs*nincr-1, 0, 255);
    free_chmatrix(kername_list_cur    ,0, nincs*nincr-1, 0, 255);
    
    MPI_Finalize();
    return 0;
}

/* subroutines */
void read_wfd(char *f_sp_inc, char *f_sp_i_inc, float *sp_inc, float *sp_i_inc, 
              int nfreq, int num, int myid, int np){
    //////////?????????????????????[1]
    read_wfd_one(f_sp_inc,       sp_inc,   nfreq, num, myid, np);
    read_wfd_one(f_sp_i_inc,     sp_i_inc, nfreq, num, myid, np);
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
    /* Local variables */
    int i, it, nlen, nlen2;
    float max1, max2, dt, df, omega, tmp;
    //double fac_r1, fac_i1;
    float *vec, *weight;
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
        for(it=2; it<=nlen; it++) wavelet_r[it]=wavelet_r[it-1]+vec[it-1];
        free(vec);
        /* adjoint source in sac */
        vec = read_sac(fadj_s, &nlen2, &dt);
        seismo_r = vector(1, nlen);
        for(it=1; it<=nlen; it++) seismo_r[it]=vec[it-1];
        free(vec);
        if(nlen!=nlen2||(omega1+omega_inc*(nfreq-1)>nlen)){
            throw_error("Please check the size of the *.su file!\n");
        }
        
        /*2. allocate memory */
        wavelet_i   = vector(1, nlen);
        wavelet_abs = vector(1, nlen);
        seismo_i   = vector(1, nlen);
        seismo_abs = vector(1, nlen);
        weight      = vector(1, nlen);
        
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
        
        /* calculate the weight */
        df=2.0*PI/(nlen*dt);
        tmp=0.0;
        weight[1]=0.0;
        for(it=2; it<=nlen/2; it++){/*ignore the zero component */
            omega=(it-1)*df;
            weight[it]=omega*omega*seismo_abs[it]*seismo_abs[it];
            tmp+=weight[it];
        }
        for(it=2; it<=nlen/2; it++){
            weight[it]/=tmp;
            //printf("it: %d %.3e\n", it, weight[it]);
        }
        
        /* record the useful frequency points */
        for(it=1; it<=nfreq; it++){
        //for(it=41; it<=41; it++){
            irecord[it]=0;
            i=omega1+omega_inc*(it-1);
            omega=(i-1)*df;
            if(wavelet_abs[i]>=0.001*max1&&seismo_abs[i]>=0.001*max2){
                irecord[it]=1;
                //printf("2irecord[it]=%d %d\n",irecord[it],i);
                tmp=weight[i]*omega/(seismo_abs[i]*seismo_abs[i]*wavelet_abs[i]*wavelet_abs[i]);
                //tmp=weight[i]/(seismo_abs[i]*seismo_abs[i]*wavelet_abs[i]*wavelet_abs[i]);
                fac_r[it]=(wavelet_r[i]*seismo_r[i]-wavelet_i[i]*seismo_i[i])*tmp;
                fac_i[it]=-(wavelet_r[i]*seismo_i[i]+wavelet_i[i]*seismo_r[i])*tmp;
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
        free_vector(weight, 1, nlen);
    }
    /* bcast integers */
    MPI_Bcast(&irecord[1],nfreq,MPI_INT,0,MPI_COMM_WORLD);
    /* bcast floats */
    MPI_Bcast(&fac_r[1],nfreq,MPI_FLOAT,0,MPI_COMM_WORLD);
    MPI_Bcast(&fac_i[1],nfreq,MPI_FLOAT,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
}

void cal_kernel(int mynum, int nfreq, int nlen, int *irecord, float *fac_r,
                float *fac_i, float *sp_inc, float *sp_i_inc, float *sp_adj, 
                float *sp_i_adj, float *Mvp, float *csp){
    
    /* Local variables */
    int i, j, k, it;
    float corr_r, corr_i;
    
    /* initialize the memory */
    for(i=1;i<=mynum;i++){
        csp[i]=0.0;
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
            k=(it-1)*mynum+i;
            corr_r =sp_inc[j]*sp_adj[k] - sp_i_inc[j]*sp_i_adj[k];
            corr_i =sp_inc[j]*sp_i_adj[k]+sp_i_inc[j]*sp_adj[k];
            //Kr_im[i] += corr_r*fac_r[it] - corr_i*fac_i[it];
            csp[i] += corr_r*fac_i[it] + corr_i*fac_r[it];
        }
    }
    for(i=1;i<=mynum;i++){
        csp[i]= csp[i]/nlen*2.0/(Mvp[i]*Mvp[i]);
    }
}
