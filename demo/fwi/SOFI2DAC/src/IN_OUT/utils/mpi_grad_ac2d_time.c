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
@brief Calculate the gradient in the time domain for 2D acoustic case
@details This program is meant to be run to calculate the gradient for 2D 
acoustic case. One can employ time domain forward modelling.
@file mpi_grad_ac2d_time.c
@author Tao Lei (leit@sustech.edu.cn), Wei Zhang
@date 18 Sep 2019
@param[in]  verbose 1: print the details; 0: run in quiet pattern, default value;
@param[in]  inv_vp  1: inverse parameter vp; 0: not inverse it;
@param[in]  inv_rho 1: inverse parameter rho; 0: not inverse it;
@param[in]  inv_adj 1: reverse the snapshot when one calculates the adjoint wavefields
                    with clockwise TDFD operator; 0: the adjoint wavefields are 
                    calculated by time-reversal TDFD operator;
@param[in]  fmod    the full path of vp model, the path of density model 
                    is derived from vp model, the following [finc, fadj and fout_vp] are similar;
@param[in]  finc    the full path of vx component <em>(clockwise)</em> excited from the source side;
@param[in]  fadj    the full path of vx component excited from the receiver side, please check the parameter <em>inv_adj</em>;
@param[in,out] fout_vp the full path of vp gradient;
@code 
mpicc -o mpi_grad_ac2d_time mpi_grad_ac2d_time.c -I../include/ -L../lib -lIN_OUT -lm

mpirun -np 4 mpi_grad_ac2d_time verbose=1 inv_vp=1 inv_rho=1 fmod=./model/inv/inv.vp finc=./snap/inc/test.vx.shot1.bin \
fadj=./snap/adj/test.vx.shot1.bin fout_vp=./gradient/grad_vp_shot1.bin inv_adj=1
@endcode
*/

#include "par.h"
#include "comm.h"
#include "grad.h"
#include"mpi.h"
#define BLOCK_LOW(rank,size,n)    ((rank)*(n)/(size))
#define BLOCK_HIGH(rank,size,n)    (BLOCK_LOW((rank)+1,size,n)-1)
#define BLOCK_SIZE(rank,size,n)    (BLOCK_HIGH(rank,size,n)-BLOCK_LOW(rank,size,n)+1)
#define STRING_SIZE 256

/* Global variables */
int xargc; 
char **xargv;
/* 
 rank: the number of the current processor
 size: the total number of processors
 n   : the total number of jobs or elements
 job0, job1, job2, job3, job4, job5, job6, job7, job8
 here n=9, if size=3, then there are 3 blocks, as,
 BLOCK INDEX                    | BLOCK_LOW BLOCK_HIGH BLOCK_SIZE
 block0: job0, job1, job2       | 0 2 3
 block1: job3, job4, job5       | 3 5 3
 block2: job6, job7, job8       | 6 8 3
 it doesn't matter if n can't be devided by size.
 if n=8 and size=3, then blocks are as: 0 1   |  2 3 4  |  5 6 7
*/
int main(int argc, char *argv[]){
    int NP,MYID;
    int rc;
    
    MPI_Status status;
    
    /* IN and OUT */
    int verbose, inv_vp, inv_rho, inv_adj, eprecond;
    float epsilon, dt;
    char *fmod, *finc, *fadj, *fout_vp;
    
    /* Local variables */
    int size_wfd, size_wfd2, size_mod, num, nt;
    int mynum, it, view_offset, view_offset2;
    float *Mvp, *Mrho;
    float *sp_inc, *vx_inc, *vy_inc, *sp_adj, *vx_adj, *vy_adj;
    float *cvx, *cvy, *csp, *grad_lam, *grad_rho_s, *grad_vp, *grad_rho;
    char f_mod[256], f_vp[256], f_rho[256]; 
    char f_vx_inc[256], f_vy_inc[256], f_sp_inc[256];
    char f_vx_adj[256], f_vy_adj[256], f_sp_adj[256];
    char f_out_vp[256], f_out_rho[256];
    MPI_File fp_vp, fp_rho;
    MPI_File fp_sp_inc, fp_vx_inc, fp_vy_inc, fp_sp_adj, fp_vx_adj, fp_vy_adj;
    MPI_File mpi_fh;
    
    /* Initialize MPI environment */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &MYID);
    MPI_Comm_size(MPI_COMM_WORLD, &NP);
    
    /* 1. read formatted input list */
    /* initialize getpar */
    initargs(argc,argv);
    
    if(MYID==0){
        verbose=0; getparint("verbose",&verbose);
        inv_vp=0; getparint("inv_vp",&inv_vp);
        inv_rho=0; getparint("inv_rho",&inv_rho);
        if(inv_vp==0 && inv_rho==0){
            throw_error("at least one class of parameters should be inverted!\n");
        }
        inv_adj=0; getparint("inv_adj",&inv_adj);
        /* preconditioning */
        eprecond=0; getparint("eprecond",&eprecond);
        epsilon=0.03; getparfloat("epsilon",&epsilon);
        dt=1.0; getparfloat("dt",&dt);
        
        /* the files of incident and adjoint wavefields */
        if(!getparstring("fmod", &fmod)) throw_error("fmod must be specified!\n");
        if(!getparstring("finc", &finc)) throw_error("finc must be specified!\n");
        if(!getparstring("fadj", &fadj)) throw_error("fadj must be specified!\n");
        if(!getparstring("fout_vp", &fout_vp)) throw_error("fout_vp must be specified!\n");
        size_wfd=get_file_size(finc);
        size_wfd2=get_file_size(fadj);
        size_mod=get_file_size(fmod);
        if((size_wfd%4!=0)||(size_mod%4!=0)||(size_wfd2!=size_wfd)
            ||(size_wfd%size_mod!=0)){
            printf("size_wfd, size_wfd2, size_mod=%d %d %d\n",
                size_wfd, size_wfd2, size_mod);
            throw_error("Please check the size of the files!\n");
        }else{
            num=size_mod/4;
            nt=size_wfd/size_mod;
            strcpy(f_mod, fmod);
            strcpy(f_vx_inc, finc);
            strcpy(f_vx_adj, fadj);
            strcpy(f_out_vp, fout_vp);
        } 
    }
    
    /* bcast integers */
    MPI_Bcast(&verbose,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&inv_vp,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&inv_rho,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&inv_adj,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&eprecond,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&num,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&nt,1,MPI_INT,0,MPI_COMM_WORLD);
    /* bcast floats */
    MPI_Bcast(&epsilon,1,MPI_FLOAT,0,MPI_COMM_WORLD);
    MPI_Bcast(&dt,1,MPI_FLOAT,0,MPI_COMM_WORLD);
    /* bcast strings */
    MPI_Bcast(&f_mod,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(&f_vx_inc,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(&f_vx_adj,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(&f_out_vp,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    sprintf(f_vp,"%s",f_mod);
    strrpc2(f_vp,"vp","rho",f_rho); //model
    strrpc2(f_out_vp,"vp","rho",f_out_rho); //gradient
    strrpc2(f_vx_inc,".vx.",".vy.",f_vy_inc); //wavefields
    strrpc2(f_vx_inc,".vx.",".sp.",f_sp_inc);
    strrpc2(f_vx_adj,".vx.",".vy.",f_vy_adj);
    strrpc2(f_vx_adj,".vx.",".sp.",f_sp_adj);
    
    if(MYID==1||NP==1){
        //printf("\nMYID=%d, NP=%d\n",MYID, NP);
        if(verbose){
            printf("\n******\n");
            printf("inv_vp    \t| %d\n",inv_vp);
            printf("inv_rho   \t| %d\n",inv_rho);
            printf("inv_adj   \t| %d\n",inv_adj);
            printf("eprecond  \t| %d\n",eprecond);
            printf("epsilon   \t| %.2e\n",epsilon);
            printf("dt        \t| %.3e\n",dt);
            printf("num       \t| %d\n",num);
            printf("nt        \t| %d\n",nt);
            printf("f_vp      \t| %s\n",f_vp);
            printf("f_rho     \t| %s\n",f_rho);
            printf("f_vx_inc  \t| %s\n",f_vx_inc);
            printf("f_vy_inc  \t| %s\n",f_vy_inc);
            printf("f_sp_inc  \t| %s\n",f_sp_inc);
            printf("f_vx_adj  \t| %s\n",f_vx_adj);
            printf("f_vy_adj  \t| %s\n",f_vy_adj);
            printf("f_sp_adj  \t| %s\n",f_sp_adj);
            printf("f_out_vp  \t| %s\n",f_out_vp);
            printf("f_out_rho  \t| %s\n",f_out_rho);
            printf("******\n");
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    mynum = BLOCK_SIZE(MYID,NP,num);
    //printf("MYID,BLOCK_LOW,BLOCK_HIGH,mynum=%d %d %d %d\n", MYID, BLOCK_LOW(MYID,NP,num), BLOCK_HIGH(MYID,NP,num), mynum);
    
    /* allocate the memory */
    Mvp    = vector(1, mynum);
    Mrho   = vector(1, mynum);
    //fullmod= vector(1, mynum);
    sp_inc = vector(1, mynum);
    vx_inc = vector(1, mynum);
    vy_inc = vector(1, mynum);
    sp_adj = vector(1, mynum);
    vx_adj = vector(1, mynum);
    vy_adj = vector(1, mynum);
    cvx    = vector(1, mynum);
    cvy    = vector(1, mynum);
    csp    = vector(1, mynum);
    grad_lam    = vector(1, mynum);
    grad_rho_s  = vector(1, mynum);
    grad_vp     = vector(1, mynum);
    grad_rho    = vector(1, mynum);
    
    /*1. prepare the files for acoustic case */
    /* models */
    rc = MPI_File_open(MPI_COMM_WORLD,f_vp, MPI_MODE_RDONLY,MPI_INFO_NULL,&fp_vp);
    if(rc) throw_error("Unable to open file: f_vp!\n");
    rc = MPI_File_open(MPI_COMM_WORLD,f_rho, MPI_MODE_RDONLY,MPI_INFO_NULL,&fp_rho);
    if(rc) throw_error("Unable to open file: f_rho!\n");
    /* incident wavefields */
    rc = MPI_File_open(MPI_COMM_WORLD,f_sp_inc, MPI_MODE_RDONLY,MPI_INFO_NULL,&fp_sp_inc);
    if(rc) throw_error("Unable to open file: f_sp_inc!\n");
    rc = MPI_File_open(MPI_COMM_WORLD,f_vx_inc, MPI_MODE_RDONLY,MPI_INFO_NULL,&fp_vx_inc);
    if(rc) throw_error("Unable to open file: f_vx_inc!\n");
    rc = MPI_File_open(MPI_COMM_WORLD,f_vy_inc, MPI_MODE_RDONLY,MPI_INFO_NULL,&fp_vy_inc);
    if(rc) throw_error("Unable to open file: f_vy_inc!\n");
    /* adjoint wavefields */
    rc = MPI_File_open(MPI_COMM_WORLD,f_sp_adj, MPI_MODE_RDONLY,MPI_INFO_NULL,&fp_sp_adj);
    if(rc) throw_error("Unable to open file: f_sp_adj!\n");
    rc = MPI_File_open(MPI_COMM_WORLD,f_vx_adj, MPI_MODE_RDONLY,MPI_INFO_NULL,&fp_vx_adj);
    if(rc) throw_error("Unable to open file: f_vx_adj!\n");
    rc = MPI_File_open(MPI_COMM_WORLD,f_vy_adj, MPI_MODE_RDONLY,MPI_INFO_NULL,&fp_vy_adj);
    if(rc) throw_error("Unable to open file: f_vy_adj!\n");
    
    /*2. read models */
    view_offset=BLOCK_LOW(MYID,NP,num)*sizeof(float);
    MPI_File_read_at_all(fp_vp, view_offset, &Mvp[1], mynum, MPI_FLOAT, &status);
    MPI_File_read_at_all(fp_rho, view_offset, &Mrho[1], mynum, MPI_FLOAT, &status);
    MPI_File_close(&fp_vp);
    MPI_File_close(&fp_rho);
    MPI_Barrier(MPI_COMM_WORLD);
    
    /*3. calculate the correlation */
    for(it=1;it<=nt;it++){
        view_offset=((it-1)*num+BLOCK_LOW(MYID,NP,num))*sizeof(float);
        if(inv_adj){
            view_offset2=((nt-it)*num+BLOCK_LOW(MYID,NP,num))*sizeof(float);
        }else{
            view_offset2=view_offset;
        }
        /* read the wavefields */
        /* incident parts */
        MPI_File_read_at_all(fp_sp_inc, view_offset, &sp_inc[1], mynum, MPI_FLOAT, &status);
        MPI_File_read_at_all(fp_vx_inc, view_offset, &vx_inc[1], mynum, MPI_FLOAT, &status);
        MPI_File_read_at_all(fp_vy_inc, view_offset, &vy_inc[1], mynum, MPI_FLOAT, &status);
        /* adjoint parts */
        MPI_File_read_at_all(fp_sp_adj, view_offset2, &sp_adj[1], mynum, MPI_FLOAT, &status);
        MPI_File_read_at_all(fp_vx_adj, view_offset2, &vx_adj[1], mynum, MPI_FLOAT, &status);
        MPI_File_read_at_all(fp_vy_adj, view_offset2, &vy_adj[1], mynum, MPI_FLOAT, &status);
        
        /* calculate the zero-lag cross-correlation */
        calc_correlation_ac_time(1, mynum, inv_vp, inv_rho, vx_inc, vy_inc, 
                                 sp_inc, vx_adj, vy_adj, sp_adj, cvx, cvy, csp);
        
    }
    /* close the files */
    MPI_File_close(&fp_sp_inc);
    MPI_File_close(&fp_vx_inc);
    MPI_File_close(&fp_vy_inc);
    MPI_File_close(&fp_sp_adj);
    MPI_File_close(&fp_vx_adj);
    MPI_File_close(&fp_vy_adj);
    MPI_Barrier(MPI_COMM_WORLD);
    
    /*4. calculate gradident */
    gradient_ac(1, mynum, dt, inv_vp, inv_rho, Mrho, Mvp, cvx, cvy, csp, 
                grad_lam, grad_rho_s, grad_vp, grad_rho);
    
    view_offset=BLOCK_LOW(MYID,NP,num)*sizeof(float);
    /* write gradient vp */
    //printf("%s\n%s\n",f_out_vp,f_out_rho);
    MPI_File_open(MPI_COMM_WORLD, f_out_vp, MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&mpi_fh);
    MPI_File_set_view(mpi_fh, view_offset, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL); 
    MPI_File_write_all(mpi_fh, &grad_vp[1], mynum, MPI_FLOAT, &status);
    MPI_File_close(&mpi_fh);
    MPI_Barrier(MPI_COMM_WORLD);
    /* write gradient rho */
    MPI_File_open(MPI_COMM_WORLD, f_out_rho, MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&mpi_fh);
    MPI_File_set_view(mpi_fh, view_offset, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL); 
    MPI_File_write_all(mpi_fh, &grad_rho[1], mynum, MPI_FLOAT, &status);
    MPI_File_close(&mpi_fh);

    /* deallocate the memory */
    free_vector(Mvp, 1, mynum);
    free_vector(Mrho, 1, mynum);
    //free_vector(fullmod, 1, mynum);
    free_vector(sp_inc, 1, mynum);
    free_vector(vx_inc, 1, mynum);
    free_vector(vy_inc, 1, mynum);
    free_vector(sp_adj, 1, mynum);
    free_vector(vx_adj, 1, mynum);
    free_vector(vy_adj, 1, mynum);
    free_vector(cvx, 1, mynum);
    free_vector(cvy, 1, mynum);
    free_vector(csp, 1, mynum);
    free_vector(grad_lam, 1, mynum);
    free_vector(grad_rho_s, 1, mynum);
    free_vector(grad_vp, 1, mynum);
    free_vector(grad_rho, 1, mynum);
    
    MPI_Finalize();
    return 0;
}
