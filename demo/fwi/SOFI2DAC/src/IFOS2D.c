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

/* ----------------------------------------------------------------------
 * This is program IFOS Version 2.0.3
 * Inversion of Full Observerd Seismograms
 *
 *  ----------------------------------------------------------------------*/


#include "fd.h"           /* general include file for viscoelastic FD programs */

#include "globvar.h"      /* definition of global variables  */
#include "cseife.h"

#include "stfinv/stfinv.h" /* libstfinv - inversion for source time function */

int main(int argc, char **argv){
    /* variables in main */
    int ns, nseismograms=0, nt, nd, fdo3, j, i, iter, h, infoout, SHOTINC,  hin, hin1, do_stf=0;
    int NTDTINV, nxny, nxnyi, imat, imat1, imat2, IDXI, IDYI, hi, NTST, NTSTI;
    int lsnap, nsnap=0, lsamp=0, buffsize,  swstestshot, snapseis, snapseis1;
    int ntr=0, ntr_loc=0, ntr_glob=0, nsrc=0, nsrc_loc=0, nsrc_glob=0, ishot, irec, nshots=0, nshots1, Lcount, itest, itestshot;
    
    float muss, lamss;
    float memdyn, memmodel, memseismograms, membuffer, memtotal, eps_scale;
    float fac1, fac2;
    float opteps_vp, opteps_vs, opteps_rho, Vp_avg, C_vp, Vs_avg, C_vs, rho_avg, C_rho;
    float memfwt, memfwt1, memfwtdata;
    char *buff_addr, ext[10], *fileinp;
    char jac[225];
    
    double time1, time2, time3, time4, time5, time6, time7, time8,
    time_av_v_update=0.0, time_av_s_update=0.0, time_av_v_exchange=0.0,
    time_av_s_exchange=0.0, time_av_timestep=0.0;
    
    float L2, L2sum, L2_all_shots, L2sum_all_shots, *L2t, alphanom, alphadenom;
    int sum_killed_traces=0, sum_killed_traces_testshots=0, killed_traces=0, killed_traces_testshots=0;
    int *ptr_killed_traces=&killed_traces, *ptr_killed_traces_testshots=&killed_traces_testshots;
    
    float energy, energy_sum, energy_all_shots, energy_sum_all_shots = 0.0;
    float energy_SH, energy_sum_SH, energy_all_shots_SH, energy_sum_all_shots_SH;
    float L2_SH, L2sum_SH, L2_all_shots_SH, L2sum_all_shots_SH;
    
    // Pointer for dynamic wavefields:
    float  **  psxx, **  psxy, **  psyy, **  psxz, **  psyz, **psp, ** ux, ** uy, ** uxy, ** uyx, ** u, ** Vp0, ** uttx, ** utty, ** Vs0, ** Rho0;
    float  **  pvx, **  pvy, **  pvz, **waveconv, **waveconv_lam, **waveconv_mu, **waveconv_rho, **waveconv_rho_s, **waveconv_u, **waveconvtmp, **wcpart, **wavejac,**waveconv_rho_s_z,**waveconv_u_z,**waveconv_rho_z;
    float **waveconv_shot, **waveconv_u_shot, **waveconv_rho_shot, **waveconv_u_shot_z, **waveconv_rho_shot_z;
    float  **  pvxp1, **  pvyp1, **  pvzp1, **  pvxm1, **  pvym1, **  pvzm1;
    float ** gradg, ** gradp,** gradg_rho, ** gradp_rho, ** gradg_u, ** gradp_u, ** gradp_u_z,** gradp_rho_z;
    float  **  prho,**  prhonp1, **prip=NULL, **prjp=NULL, **pripnp1=NULL, **prjpnp1=NULL, **  ppi, **  pu, **  punp1, **  puipjp, **  ppinp1;
    float  **  vpmat, ***forward_prop_x, ***forward_prop_y, ***forward_prop_rho_x, ***forward_prop_u, ***forward_prop_rho_y, ***forward_prop_p;
    
    float ***forward_prop_z_xz,***forward_prop_z_yz,***forward_prop_rho_z,**waveconv_mu_z;
    float ** uxz, ** uyz;
    
    float  ** sectionvx=NULL, ** sectionvy=NULL, ** sectionvz=NULL, ** sectionp=NULL, ** sectionpnp1=NULL,
    ** sectioncurl=NULL, ** sectiondiv=NULL, ** sectionvxdata=NULL, ** sectionvydata=NULL, ** sectionvzdata=NULL, ** sectionvxdiff=NULL, ** sectionvzdiff=NULL, ** sectionvxdiffold=NULL, ** sectionvydiffold=NULL, ** sectionvzdiffold=NULL,** sectionpdata=NULL, ** sectionpdiff=NULL, ** sectionpdiffold=NULL,
    ** sectionvydiff=NULL, ** sectionpn=NULL, ** sectionread=NULL, ** sectionvy_conv=NULL, ** sectionvy_obs=NULL, ** sectionvx_conv=NULL,** sectionvx_obs=NULL, ** sectionvz_conv=NULL,** sectionvz_obs=NULL,
    ** sectionp_conv=NULL,** sectionp_obs=NULL, * source_time_function=NULL;
    float  **  absorb_coeff, ** taper_coeff, * epst1, * epst2,  * epst3, * picked_times;
    float  ** srcpos=NULL, **srcpos_loc=NULL, ** srcpos1=NULL, **srcpos_loc_back=NULL, ** signals=NULL,** signals_SH=NULL,  *hc=NULL;
    int   ** recpos=NULL, ** recpos_loc=NULL;
    /*int   ** tracekill=NULL, TRKILL, DTRKILL;*/
    int * DTINV_help;
    
    float ** bufferlef_to_rig,  ** bufferrig_to_lef, ** buffertop_to_bot, ** bufferbot_to_top;
    
    /* PML variables */
    float * d_x, * K_x, * alpha_prime_x, * a_x, * b_x, * d_x_half, * K_x_half, * alpha_prime_x_half, * a_x_half, * b_x_half, * d_y, * K_y, * alpha_prime_y, * a_y, * b_y, * d_y_half, * K_y_half, * alpha_prime_y_half, * a_y_half, * b_y_half;
    float ** psi_sxx_x, ** psi_syy_y, ** psi_sxy_y, ** psi_sxy_x, ** psi_vxx, ** psi_vyy, ** psi_vxy, ** psi_vyx, ** psi_vxxs;
    float ** psi_sxz_x, ** psi_syz_y, ** psi_vzx, ** psi_vzy;
    
    /* Variables for viscoelastic modeling */
    float **ptaus=NULL, **ptaup=NULL, *etaip=NULL, *etajm=NULL, *peta=NULL, **ptausipjp=NULL, **fipjp=NULL, ***dip=NULL, *bip=NULL, *bjm=NULL;
    float *cip=NULL, *cjm=NULL, ***d=NULL, ***e=NULL, ***pr=NULL, ***pp=NULL, ***pq=NULL, **f=NULL, **g=NULL;
    float ***pt=NULL, ***po=NULL; // SH Simulation
    
    /* Variables for step length calculation */
    int step1, step2, step3=0, itests, iteste, stepmax, countstep;
    float scalefac;
    
    int RECINC, ntr1;
    int SOURCE_SHAPE_OLD=0;
    int SOURCE_SHAPE_OLD_SH=0;
    
    /* Variables for conjungate gradient */
    int PCG_iter_start=1;
    
    /* Variables for L-BFGS */
    int LBFGS_NPAR=3;
    int LBFGS_iter_start=1;
    float **s_LBFGS,**y_LBFGS, *rho_LBFGS;
    int l=0;
    int m=0;
    
    /* Check wolfe */
    int steplength_search=0;
    int FWI_run=1;
    int gradient_optimization=1;
    float alpha_SL_min=0, alpha_SL_max=0, alpha_SL=1.0;
    float alpha_SL_old;
    float ** waveconv_old,** waveconv_u_old,** waveconv_rho_old;
    float ** waveconv_up,** waveconv_u_up,** waveconv_rho_up;
    float L2_SL_old=0, L2_SL_new=0;
    float c1_SL=1e-4, c2_SL=0.9;
    int wolfe_status;
    int wolfe_sum_FWI=0;
    int wolfe_found_lower_L2=0;
    float alpha_SL_FS;
    float L2_SL_FS;
    int use_wolfe_failsafe=0;
    int wolfe_SLS_failed=0;
    
    /* Variables for energy weighted gradient */
    float ** Ws, **Wr, **We;
    float ** Ws_SH, **Wr_SH, **We_SH;
    float ** We_sum,** We_sum_SH;
    float We_sum_max1;
    float We_max_SH,We_max;
    
    int * recswitch=NULL;
    float ** fulldata=NULL, ** fulldata_vx=NULL, ** fulldata_vy=NULL, ** fulldata_vz=NULL, ** fulldata_p=NULL, ** fulldata_curl=NULL, ** fulldata_div=NULL;
    
    /*vector for abort criterion*/
    float * L2_hist=NULL;
    
    /* help variable for MIN_ITER */
    int min_iter_help=0;
    
    float ** workflow=NULL;
    int workflow_lines;
    char workflow_header[STRING_SIZE];
    int change_wavetype_iter=-10; /* Have to be inialized negative */
    int wavetype_start; /* We need this due to MPI Comm */
    int buf1=0, buf2=0;
    WORKFLOW_STAGE=1;
    
    /* variable for time domain filtering */
    float F_LOW_PASS;
    float *F_LOW_PASS_EXT=NULL;
    int nfrq=0;
    int FREQ_NR=1;

    float JOINT_EQUAL_PSV=0.0, JOINT_EQUAL_SH=0.0;
    float JOINT_EQUAL_PSV_all=0.0, JOINT_EQUAL_SH_all=0.0;
    int JOINT_EQUAL_new_max=1;
    
    FILE *fprec, *FPL2;
    
    FILE *FPL2_JOINT;
    char L2_joint_log[STRING_SIZE];
    
    /* General parameters */
    int nt_out;
    
    MPI_Request *req_send, *req_rec;
    MPI_Status  *send_statuses, *rec_statuses;
    
    /* Initialize MPI environment */
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&NP);
    MPI_Comm_rank(MPI_COMM_WORLD,&MYID);
    
    setvbuf(stdout, NULL, _IONBF, 0);
    
    if (MYID == 0){
        time1=MPI_Wtime();
        clock();
    }
    
    /* print program name, version etc to stdout*/
    if (MYID == 0) info(stdout);
    
    /* read parameters from parameter-file (stdin) */
    fileinp=argv[1];
    FP=fopen(fileinp,"r");
    if(FP==NULL) {
        if (MYID == 0){
            printf("\n==================================================================\n");
            printf(" Cannot open IFOS input file %s \n",fileinp);
            printf("\n==================================================================\n\n");
            declare_error(" --- ");
        }
    }
    
    /* read json formatted input file */
    read_par_json(stdout,fileinp);
    
    exchange_par();
    
    wavetype_start=WAVETYPE;
    if (MYID == 0) note(stdout);
    
    
    /* open log-file (each PE is using different file) */
    /*	fp=stdout; */
    sprintf(ext,".%i",MYID);
    strcat(LOG_FILE,ext);
    
    /* If Verbose==0, no PE will write a log file */
    if(!VERBOSE) sprintf(LOG_FILE,"/dev/null");
    
    if ((MYID==0)) FP=stdout;
    else {
        FP=fopen(LOG_FILE,"w");
    }
    fprintf(FP," This is the log-file generated by PE %d \n\n",MYID);
    
    /* domain decomposition */
    initproc();
    
    NT=iround(TIME/DT);  	  /* number of timesteps */
    /*ns=iround(NT/NDT);*/           /* number of samples per trace */
    ns=NT;	/* in a FWI one has to keep all samples of the forward modeled data
             at the receiver positions to calculate the adjoint sources and to do
             the backpropagation; look at function saveseis_glob.c to see that every
             NDT sample for the forward modeled wavefield is written to su files*/
    
    /* output of parameters to log-file or stdout */
    if (MYID==0) write_par(FP);
    
    
    /* NXG, NYG denote size of the entire (global) grid */
    NXG=NX;
    NYG=NY;
    
    /* In the following, NX and NY denote size of the local grid ! */
    NX = IENDX;
    NY = IENDY;
        
    /* Reading source positions from SOURCE_FILE */
    srcpos=sources(&nsrc);
    nsrc_glob=nsrc;
    ishot=0;    
    
    if (SEISMO){
        recpos=receiver(&ntr, srcpos, ishot);
        recswitch = ivector(1,ntr);
        recpos_loc = splitrec(recpos,&ntr_loc, ntr, recswitch);
        ntr_glob=ntr;
        ntr=ntr_loc;
    }
    
    /* memory allocation for abort criterion*/
    L2_hist = vector(1,1000);
    
    if(INV_STF) fulldata = matrix(1,ntr_glob,1,NT);
    
    /* estimate memory requirement of the variables in megabytes*/
    
    switch (SEISMO){
        case 1 : /* particle velocities only */
            nseismograms=2;
            break;
        case 2 : /* pressure only */
            nseismograms=1;
            break;
        case 3 : /* curl and div only */
            nseismograms=2;
            break;
        case 4 : /* everything */
            nseismograms=5;
            break;
        case 5 : /* everything except curl and div */
            nseismograms=3;
            break;
    }
    
    /* use only every DTINV time sample for the inversion */
    /*DTINV=15;*/
    DTINV_help=ivector(1,NT);
    NTDTINV=ceil((float)NT/(float)DTINV);		/* round towards next higher integer value */
    
    /* save every IDXI and IDYI spatial point during the forward modelling */
    IDXI=1;
    IDYI=1;
    
    /*allocate memory for dynamic, static and buffer arrays */
    fac1=(NX+FDORDER)*(NY+FDORDER);
    fac2=sizeof(float)*pow(2.0,-20.0);
    
    nd = FDORDER/2 + 1;
    
    // decide how much space for exchange is needed
    switch (WAVETYPE) {
        case 1:
            fdo3 = 2*nd;
            break;
        case 2:
            fdo3 = 1*nd;
            break;
        case 3:
            fdo3 = 3*nd;
            break;
        default:
            fdo3 = 2*nd;
            break;
    }
    
    
    if (L){
        memdyn=(5.0+3.0*(float)L)*fac1*fac2;
        memmodel=(12.0+3.0*(float)L)*fac1*fac2;
        
    } else {
        memdyn=5.0*fac1*fac2;
        memmodel=6.0*fac1*fac2;
    }
    memseismograms=nseismograms*ntr*ns*fac2;
    
    memfwt=5.0*((NX/IDXI)+FDORDER)*((NY/IDYI)+FDORDER)*NTDTINV*fac2;
    memfwt1=20.0*NX*NY*fac2;
    memfwtdata=6.0*ntr*ns*fac2;
    
    membuffer=2.0*fdo3*(NY+NX)*fac2;
    buffsize=2.0*2.0*fdo3*(NX+NY)*sizeof(MPI_FLOAT);
    memtotal=memdyn+memmodel+memseismograms+memfwt+memfwt1+memfwtdata+membuffer+(buffsize*pow(2.0,-20.0));
    
    
    if (MYID==0 && WAVETYPE == 1){
        fprintf(FP,"\n **Message from main (printed by PE %d):\n",MYID);
        fprintf(FP," Size of local grids: NX=%d \t NY=%d\n",NX,NY);
        fprintf(FP," Each process is now trying to allocate memory for:\n");
        fprintf(FP," Dynamic variables: \t\t %6.2f MB\n", memdyn);
        fprintf(FP," Static variables: \t\t %6.2f MB\n", memmodel);
        fprintf(FP," Seismograms: \t\t\t %6.2f MB\n", memseismograms);
        fprintf(FP," Buffer arrays for grid exchange:%6.2f MB\n", membuffer);
        fprintf(FP," Network Buffer for MPI_Bsend: \t %6.2f MB\n", buffsize*pow(2.0,-20.0));
        fprintf(FP," ------------------------------------------------ \n");
        fprintf(FP," Total memory required: \t %6.2f MB.\n\n", memtotal);
    }
    
    
    /* allocate buffer for buffering messages */
    buff_addr=malloc(buffsize);
    if (!buff_addr) declare_error("allocation failure for buffer for MPI_Bsend !");
    MPI_Buffer_attach(buff_addr,buffsize);
    
    /* allocation for request and status arrays */
    req_send=(MPI_Request *)malloc(REQUEST_COUNT*sizeof(MPI_Request));
    req_rec=(MPI_Request *)malloc(REQUEST_COUNT*sizeof(MPI_Request));
    send_statuses=(MPI_Status *)malloc(REQUEST_COUNT*sizeof(MPI_Status));
    rec_statuses=(MPI_Status *)malloc(REQUEST_COUNT*sizeof(MPI_Status));
    
    
    /* memory allocation for dynamic (wavefield) arrays */
    if(!ACOUSTIC){
        switch (WAVETYPE) {
            case 1: // P and SV Waves
                psxx =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
                psxy =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
                psyy =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
                break;
                
            case 2: // SH Waves
                psxz =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
                psyz =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
                break;
                
            case 3: // P, SH and SV Waves
                psxx =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
                psxy =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
                psyy =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
                psxz =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
                psyz =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
                break;
        }
    }else{
        psp  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    }
    
    if(GRAD_METHOD==2) {
        /* Allocate memory for L-BFGS */
        
        if(WAVETYPE==2) LBFGS_NPAR=2;
        
        s_LBFGS=fmatrix(1,N_LBFGS,1,LBFGS_NPAR*NX*NY);
        
        y_LBFGS=fmatrix(1,N_LBFGS,1,LBFGS_NPAR*NX*NY);
        
        rho_LBFGS=vector(1,N_LBFGS);
        
        for(l=1;l<=N_LBFGS;l++){
            for(m=1;m<=LBFGS_NPAR*NX*NY;m++){
                s_LBFGS[l][m]=0.0;
                y_LBFGS[l][m]=0.0;
            }
            rho_LBFGS[l]=0.0;
        }
    }
    
    if(!ACOUSTIC){
        if(WAVETYPE==1||WAVETYPE==3){
            ux   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            uy   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            uxy  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            uyx  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            uttx   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            utty   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        }
        if(WAVETYPE==2||WAVETYPE==3){
            uxz   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            uyz   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        }
    }else{
        u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    }
    
    switch (WAVETYPE) {
        case 1: // P and SV Waves
            pvx  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            pvy  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            pvxp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            pvyp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            pvxm1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            pvym1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            break;
            
        case 2: // SH Waves
            pvz  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            pvzp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            pvzm1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            break;
            
        case 3: // P and SV Waves
            pvx  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            pvy  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            pvxp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            pvyp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            pvxm1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            pvym1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            pvz  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            pvzp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            pvzm1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            break;
    }
    
    Vp0  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    if(!ACOUSTIC)
        Vs0  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    Rho0  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    
    /* memory allocation for static (model) arrays */
    prho =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    prhonp1 =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    prip =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    prjp =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    pripnp1 =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    prjpnp1 =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    ppi  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    ppinp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    if(!ACOUSTIC){
        pu   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        punp1   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        puipjp   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    }
    vpmat   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    
    
    if((EPRECOND==1)||(EPRECOND==3)){
        if(WAVETYPE==1 || WAVETYPE==3) {
            We_sum = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            Ws = matrix(-nd+1,NY+nd,-nd+1,NX+nd); /* total energy of the source wavefield */
            Wr = matrix(-nd+1,NY+nd,-nd+1,NX+nd); /* total energy of the receiver wavefield */
            We = matrix(-nd+1,NY+nd,-nd+1,NX+nd); /* total energy of source and receiver wavefield */
        }
        if(WAVETYPE==2 || WAVETYPE==3) {
            We_sum_SH = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            Ws_SH = matrix(-nd+1,NY+nd,-nd+1,NX+nd); /* total energy of the source wavefield */
            Wr_SH = matrix(-nd+1,NY+nd,-nd+1,NX+nd); /* total energy of the receiver wavefield */
            We_SH = matrix(-nd+1,NY+nd,-nd+1,NX+nd); /* total energy of source and receiver wavefield */
        }
    }
    
    if (L) {
        /* dynamic (wavefield) arrays for viscoelastic modeling */
        pr = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        pp = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        pq = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        /* memory allocation for static arrays for viscoelastic modeling */
        dip = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        d =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        e =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        ptaus =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        ptausipjp =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        if(WAVETYPE==2 || WAVETYPE==3) {
            pt = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
            po = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        }
        ptaup =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        fipjp =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        f =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        g =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        peta =  vector(1,L);
        etaip =  vector(1,L);
        etajm =  vector(1,L);
        bip =  vector(1,L);
        bjm =  vector(1,L);
        cip =  vector(1,L);
        cjm =  vector(1,L);
    }
    
    NTST=20;
    nxnyi=(NX/IDXI)*(NY/IDYI);
    
    /* Parameters for step length calculations */
    stepmax = STEPMAX; /* number of maximum misfit calculations/steplength 2/3*/
    scalefac = SCALEFAC; /* scale factor for the step length */
    
    if(FORWARD_ONLY==0){
        waveconv = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        waveconv_lam = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        waveconv_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        
        waveconvtmp = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        wcpart = matrix(1,3,1,3);
        wavejac = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        
        if(!ACOUSTIC){
            forward_prop_x =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
            forward_prop_y =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
        }else{
            forward_prop_p =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
        }
        gradg = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        gradp = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        
        if(WAVETYPE==1 || WAVETYPE==3){
            forward_prop_rho_x =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
            forward_prop_rho_y =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
        }
        if(WAVETYPE==2 || WAVETYPE==3){
            forward_prop_rho_z =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
            forward_prop_z_xz =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
            forward_prop_z_yz =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
            waveconv_rho_shot_z = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            waveconv_u_shot_z = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            waveconv_mu_z = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            waveconv_rho_s_z = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            waveconv_u_z = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            waveconv_rho_z = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            gradp_u_z = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            gradp_rho_z = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        }
        
        gradg_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        gradp_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        waveconv_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        waveconv_rho_s = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        waveconv_rho_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        
        if(WOLFE_CONDITION){
            
            c1_SL=WOLFE_C1_SL;
            c2_SL=WOLFE_C2_SL;
            
            waveconv_old= matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            if(!ACOUSTIC) waveconv_u_old= matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            waveconv_rho_old= matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            
            waveconv_up= matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            if(!ACOUSTIC) waveconv_u_up= matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            waveconv_rho_up= matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        }
        
        if(!ACOUSTIC){
            forward_prop_u =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
            gradg_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            gradp_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            waveconv_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            waveconv_mu = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            waveconv_u_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        }
        
    }
    
    /* Allocate memory for boundary */
    if(FW>0){
        d_x = vector(1,2*FW);
        K_x = vector(1,2*FW);
        alpha_prime_x = vector(1,2*FW);
        a_x = vector(1,2*FW);
        b_x = vector(1,2*FW);
        
        d_x_half = vector(1,2*FW);
        K_x_half = vector(1,2*FW);
        alpha_prime_x_half = vector(1,2*FW);
        a_x_half = vector(1,2*FW);
        b_x_half = vector(1,2*FW);
        
        d_y = vector(1,2*FW);
        K_y = vector(1,2*FW);
        alpha_prime_y = vector(1,2*FW);
        a_y = vector(1,2*FW);
        b_y = vector(1,2*FW);
        
        d_y_half = vector(1,2*FW);
        K_y_half = vector(1,2*FW);
        alpha_prime_y_half = vector(1,2*FW);
        a_y_half = vector(1,2*FW);
        b_y_half = vector(1,2*FW);
        
        if (WAVETYPE==1||WAVETYPE==3){
            psi_sxx_x =  matrix(1,NY,1,2*FW);
            psi_syy_y =  matrix(1,2*FW,1,NX);
            psi_sxy_y =  matrix(1,2*FW,1,NX);
            psi_sxy_x =  matrix(1,NY,1,2*FW);
            psi_vxx   =  matrix(1,NY,1,2*FW);
            psi_vxxs  =  matrix(1,NY,1,2*FW);
            psi_vyy   =  matrix(1,2*FW,1,NX);
            psi_vxy   =  matrix(1,2*FW,1,NX);
            psi_vyx   =  matrix(1,NY,1,2*FW);
        }
        if(WAVETYPE==2||WAVETYPE == 3 ){
            psi_sxz_x =  matrix(1,NY,1,2*FW);
            psi_syz_y =  matrix(1,2*FW,1,NX);
            psi_vzx   =  matrix(1,NY,1,2*FW);
            psi_vzy   =  matrix(1,2*FW,1,NX);
        }
    }
    
    taper_coeff=  matrix(1,NY,1,NX);
    
    
    /* memory allocation for buffer arrays in which the wavefield
     information which is exchanged between neighbouring PEs is stored */
    bufferlef_to_rig = matrix(1,NY,1,fdo3);
    bufferrig_to_lef = matrix(1,NY,1,fdo3);
    buffertop_to_bot = matrix(1,NX,1,fdo3);
    bufferbot_to_top = matrix(1,NX,1,fdo3);
    
    /* Allocate memory to save full seismograms */
    switch (SEISMO){
        case 1 : /* particle velocities only */
            switch (WAVETYPE) {
                case 1:
                    fulldata_vx = matrix(1,ntr_glob,1,NT);
                    fulldata_vy = matrix(1,ntr_glob,1,NT);
                    break;
                    
                case 2:
                    fulldata_vz = matrix(1,ntr_glob,1,NT);
                    break;
                    
                case 3:
                    fulldata_vx = matrix(1,ntr_glob,1,NT);
                    fulldata_vy = matrix(1,ntr_glob,1,NT);
                    fulldata_vz = matrix(1,ntr_glob,1,NT);
                    break;
            }
            break;
        case 2 : /* pressure only */
            fulldata_p = matrix(1,ntr_glob,1,NT);
            break;
        case 3 : /* curl and div only */
            fulldata_div = matrix(1,ntr_glob,1,NT);
            fulldata_curl = matrix(1,ntr_glob,1,NT);
            break;
        case 4 : /* everything */
            switch (WAVETYPE) {
                case 1:
                    fulldata_vx = matrix(1,ntr_glob,1,NT);
                    fulldata_vy = matrix(1,ntr_glob,1,NT);
                    break;
                    
                case 2:
                    fulldata_vz = matrix(1,ntr_glob,1,NT);
                    break;
                    
                case 3:
                    fulldata_vx = matrix(1,ntr_glob,1,NT);
                    fulldata_vy = matrix(1,ntr_glob,1,NT);
                    fulldata_vz = matrix(1,ntr_glob,1,NT);
                    break;
            }
            fulldata_p = matrix(1,ntr_glob,1,NT);
            fulldata_div = matrix(1,ntr_glob,1,NT);
            fulldata_curl = matrix(1,ntr_glob,1,NT);
            break;
        case 5 : /* everything except curl and div*/
            switch (WAVETYPE) {
                case 1:
                    fulldata_vx = matrix(1,ntr_glob,1,NT);
                    fulldata_vy = matrix(1,ntr_glob,1,NT);
                    break;
                    
                case 2:
                    fulldata_vz = matrix(1,ntr_glob,1,NT);
                    break;
                    
                case 3:
                    fulldata_vx = matrix(1,ntr_glob,1,NT);
                    fulldata_vy = matrix(1,ntr_glob,1,NT);
                    fulldata_vz = matrix(1,ntr_glob,1,NT);
                    break;
            }
            fulldata_p = matrix(1,ntr_glob,1,NT);
            break;
            
    }
    if (ntr>0){
      alloc_sections(ntr,ns,&sectionvx,&sectionvy,&sectionvz,&sectionp,&sectionpnp1,&sectionpn,&sectioncurl,&sectiondiv,
	&sectionpdata,&sectionpdiff,&sectionpdiffold,&sectionvxdata,&sectionvxdiff,&sectionvxdiffold,&sectionvydata,
	&sectionvydiff,&sectionvydiffold,&sectionvzdata,&sectionvzdiff,&sectionvzdiffold);
    }
    
    /* Memory for seismic data */
    sectionread=matrix(1,ntr_glob,1,ns);
        
    /* Memory for inversion for source time function */
    if((INV_STF==1)||(TIME_FILT==1) || (TIME_FILT==2)){
        sectionp_conv=matrix(1,ntr_glob,1,NT);
        sectionp_obs=matrix(1,ntr_glob,1,NT);
        source_time_function = vector(1,NT);
        switch (WAVETYPE) {
            case 1:
                sectionvy_conv=matrix(1,ntr_glob,1,NT);
                sectionvy_obs=matrix(1,ntr_glob,1,NT);
                sectionvx_conv=matrix(1,ntr_glob,1,NT);
                sectionvx_obs=matrix(1,ntr_glob,1,NT);
                break;
                
            case 2:
                sectionvz_conv=matrix(1,ntr_glob,1,NT);
                sectionvz_obs=matrix(1,ntr_glob,1,NT);
                break;
                
            case 3:
                sectionvy_conv=matrix(1,ntr_glob,1,NT);
                sectionvy_obs=matrix(1,ntr_glob,1,NT);
                sectionvx_conv=matrix(1,ntr_glob,1,NT);
                sectionvx_obs=matrix(1,ntr_glob,1,NT);
                sectionvz_conv=matrix(1,ntr_glob,1,NT);
                sectionvz_obs=matrix(1,ntr_glob,1,NT);
                break;
        }
    }
    
    /* memory for source position definition */
    srcpos1=fmatrix(1,8,1,1);
    
    /* memory of L2 norm */
    L2t = vector(1,4);
    epst1 = vector(1,3);
    epst2 = vector(1,3);
    epst3 = vector(1,3);
    picked_times = vector(1,ntr);
    
    fprintf(FP," ... memory allocation for PE %d was successfull.\n\n", MYID);
    
    
    /* Holberg coefficients for FD operators*/
    hc = holbergcoeff();
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(FORWARD_ONLY==0&&USE_WORKFLOW){
        read_workflow(FILE_WORKFLOW,&workflow, &workflow_lines,workflow_header);
    }
    
    /* create model grids */
    if(L){
        if(!ACOUSTIC){
            if (READMOD){
                readmod(prho,ppi,pu,ptaus,ptaup,peta);
            }else{
                model(prho,ppi,pu,ptaus,ptaup,peta);
            }
        }else{
            if (READMOD){
                readmod_viscac(prho,ppi,ptaup,peta);
            }else{
                model_viscac(prho,ppi,ptaup,peta);
            }
        }
    }else{
        if(!ACOUSTIC){
            if (READMOD){
                readmod_elastic(prho,ppi,pu);
            }else{
                model_elastic(prho,ppi,pu);
            }
        }else{
            if (READMOD){
                readmod_acoustic(prho,ppi);
            }else{
                model_acoustic(prho,ppi);
            }
        }
    }
    
    /* check if the FD run will be stable and free of numerical dispersion */
    checkfd(FP, prho, ppi, pu, ptaus, ptaup, peta, hc, srcpos, nsrc, recpos, ntr_glob);
    
    /* calculate damping coefficients for CPMLs*/
    if(FW>0)
        PML_pro(d_x, K_x, alpha_prime_x, a_x, b_x, d_x_half, K_x_half, alpha_prime_x_half, a_x_half, b_x_half, d_y, K_y, alpha_prime_y, a_y, b_y, d_y_half, K_y_half, alpha_prime_y_half, a_y_half, b_y_half);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    SHOTINC=1;
    RECINC=1;
    
    switch(TIME_FILT){
        case 1: F_LOW_PASS=F_LOW_PASS_START; break;
            /*read frequencies from file*/
        case 2: F_LOW_PASS_EXT=filter_frequencies(&nfrq); F_LOW_PASS=F_LOW_PASS_EXT[FREQ_NR]; break;
    }
    
    /* Save old SOURCE_SHAPE, which is needed for STF */
    SOURCE_SHAPE_OLD = SOURCE_SHAPE;
    if(WAVETYPE==2 || WAVETYPE==3) SOURCE_SHAPE_OLD_SH=SOURCE_SHAPE_SH;
    
    nt_out=10000;
    if(!VERBOSE) nt_out=1e5;
    /*------------------------------------------------------------------------------*/
    /*----------- start fullwaveform iteration loop --------------------------------*/
    /*------------------------------------------------------------------------------*/
    
    for(iter=1;iter<=ITERMAX;iter++){  /* fullwaveform iteration loop */
        
        // At each iteration the workflow is applied
        if(USE_WORKFLOW&&(FORWARD_ONLY==0)){
            
            apply_workflow(workflow,workflow_lines,workflow_header,&iter,&F_LOW_PASS,wavetype_start,&change_wavetype_iter,&LBFGS_iter_start);
            
        }
        
        if(GRAD_METHOD==2&&(FORWARD_ONLY==0)){
            
            /* detect a change in inversion process and restart L-BFGS */
            if(iter==INV_RHO_ITER||iter==INV_VP_ITER||iter==INV_VS_ITER){
                LBFGS_iter_start=iter;
                
                if(WOLFE_CONDITION) {
                    /* Restart Step Length search */
                    alpha_SL_old=1;
                }
                
                /* set values */
                FWI_run=1;
                gradient_optimization=1;
            }
            
            /* restart L-BFGS */
            if(iter==LBFGS_iter_start) {
                lbfgs_reset(iter,N_LBFGS,LBFGS_NPAR,s_LBFGS,y_LBFGS,rho_LBFGS);
                
                /* set values */
                FWI_run=1;
                gradient_optimization=1;
            }
            
            /* Reset fail status of parabolic step length search */
            step3=0;
        }
        
        if (MYID==0){
            time2=MPI_Wtime();
            fprintf(FP,"\n\n\n ------------------------------------------------------------------\n");
            if(FORWARD_ONLY==0) {
                fprintf(FP,"\n\n\n                   TDFWI ITERATION %d \t of %d \n",iter,ITERMAX);
            } else {
                fprintf(FP,"\n\n\n                        FD-SIMULATION \n");
            }
            fprintf(FP,"\n\n\n ------------------------------------------------------------------\n");
        }
        
        countstep=0;
        
        if(GRAD_METHOD==1) {FWI_run=1; steplength_search=0; gradient_optimization=1;}
        
        /*-----------------------------------------------------*/
        /*  While loop for Wolfe step length search            */
        /*-----------------------------------------------------*/
        while(FWI_run || steplength_search || gradient_optimization) {
            
            /*-----------------------------------------------------*/
            /*              Calculate Misfit and gradient          */
            /*-----------------------------------------------------*/
            if(FWI_run){
                /* For the calculation of the material parameters between gridpoints
                 they have to be averaged. For this, values lying at 0 and NX+1,
                 for example, are required on the local grid. These are now copied from the
                 neighbouring grids */
                if (L){
                    if(!ACOUSTIC){
                        matcopy(prho,ppi,pu,ptaus,ptaup);
                    } else {
                        matcopy_viscac(prho,ppi,ptaup);
                    }
                }else{
                    if(!ACOUSTIC){
                        matcopy_elastic(prho, ppi, pu);
                    }else{
                        matcopy_acoustic(prho, ppi);
                    }
                }
                
                MPI_Barrier(MPI_COMM_WORLD);
                
                /* MPI split for processors with ntr>0 */
                int myid_ntr, group_id=0, groupsize;
                MPI_Comm	MPI_COMM_NTR;
                
                if (ntr) group_id = 1;
                else group_id = 0;
                MPI_Comm_split(MPI_COMM_WORLD, group_id, MYID, &MPI_COMM_NTR);
                MPI_Comm_rank(MPI_COMM_NTR, &myid_ntr);
                /* end of MPI split for processors with ntr>0 */
                
                
                if(!ACOUSTIC) av_mue(pu,puipjp,prho);
                av_rho(prho,prip,prjp);
                if (!ACOUSTIC && L) av_tau(ptaus,ptausipjp);
                
                
                /* Preparing memory variables for update_s (viscoelastic) */
                if (L) {
                    if(!ACOUSTIC){
                        prepare_update_s(etajm,etaip,peta,fipjp,pu,puipjp,ppi,prho,ptaus,ptaup,ptausipjp,f,g,bip,bjm,cip,cjm,dip,d,e);
                    } else {
                        prepare_update_p(etajm,peta,ppi,prho,ptaup,g,bjm,cjm,e);
                    }
                }
                
                /* Do some initia calculations */
                if(iter==1){
                    
                    /* Calculationg material parameters according to PARAMETERIZATION */
                    for (j=1;j<=NY;j=j+IDY){
                        for (i=1;i<=NX;i=i+IDX){
                            
                            if(PARAMETERIZATION==1){
                                
                                Vp0[j][i] = ppi[j][i];
                                if(!ACOUSTIC) Vs0[j][i] = pu[j][i];
                                Rho0[j][i] = prho[j][i];}
                            
                            
                            
                            if(PARAMETERIZATION==2){
                                
                                Vp0[j][i] = sqrt((ppi[j][i]+2.0*pu[j][i])/prho[j][i]);
                                Vs0[j][i] = sqrt((pu[j][i])/prho[j][i]);
                                Rho0[j][i] = prho[j][i];
                                
                            }
                            
                            if(PARAMETERIZATION==3){
                                
                                Vp0[j][i] = ppi[j][i];
                                Vs0[j][i] = pu[j][i];
                                Rho0[j][i] = prho[j][i];
                                
                            }
                            
                        }
                    }
                    
                    /* Get average values from material parameters */
                    Vp_avg=average_matrix(ppi);
                    rho_avg=average_matrix(prho);
                    if(!ACOUSTIC) Vs_avg=average_matrix(pu);
        
                    if(!ACOUSTIC) if(VERBOSE) printf("MYID = %d \t Vp_avg = %e \t Vs_avg = %e \t rho_avg = %e \n ",MYID,Vp_avg,Vs_avg,rho_avg);
                    else if(VERBOSE) printf("MYID = %d \t Vp_avg = %e \t rho_avg = %e \n ",MYID,Vp_avg,rho_avg);

                    C_vp = Vp_avg*Vp_avg;
                    if(!ACOUSTIC) C_vs = Vs_avg*Vs_avg;
                    C_rho = rho_avg*rho_avg;
                }
                
                
                /* Seperate PSV and SH logging in case of a joint inversion */
                if(WAVETYPE==3){
                    sprintf(L2_joint_log,"%s_JOINT",MISFIT_LOG_FILE);
                }
                
                /* Open Log File for L2 norm */
                if(!FORWARD_ONLY && MYID==0){
                    
                    if(iter==1){
                        
                        FPL2=fopen(MISFIT_LOG_FILE,"w");
                        
                        /* Write header for misfit log file */
                        if(GRAD_METHOD==1&&VERBOSE) {
                            if (TIME_FILT==0){
                                fprintf(FPL2,"opteps_vp \t epst1[1] \t epst1[2] \t epst1[3] \t L2t[1] \t L2t[2] \t L2t[3] \t L2t[4] \t GAMMA \n");}
                            else{
                                fprintf(FPL2,"opteps_vp \t epst1[1] \t epst1[2] \t epst1[3] \t L2t[1] \t L2t[2] \t L2t[3] \t L2t[4] \t F_LOW_PASS \t GAMMA  \n");
                            }
                        }
                        
                        if(WAVETYPE==3) FPL2_JOINT=fopen(L2_joint_log,"w");
                        
                    } else {
                        
                        FPL2=fopen(MISFIT_LOG_FILE,"a");
                        
                        if(WAVETYPE==3) FPL2_JOINT=fopen(L2_joint_log,"a");
                        
                    }
                }
                
                /* initialization of L2 calculation */
                L2=0.0;
                energy=0.0;
                L2_all_shots=0.0;
                energy_all_shots=0.0;
                killed_traces=0;
                killed_traces_testshots=0;
                
                if(WAVETYPE==2||WAVETYPE==3){
                    L2_SH=0.0;
                    energy_SH=0.0;
                    L2_all_shots_SH=0.0;
                    energy_all_shots_SH=0.0;
                }
                
                EPSILON=0.0;  /* test step length */
                exchange_par();
                
                /* initialize waveconv matrix*/
                if(WAVETYPE==1||WAVETYPE==3){
                    if(FORWARD_ONLY==0){
                        for (j=1;j<=NY;j=j+IDY){
                            for (i=1;i<=NX;i=i+IDX){
                                waveconv[j][i]=0.0;
                                waveconv_rho[j][i]=0.0;
                                if(!ACOUSTIC) waveconv_u[j][i]=0.0;
                            }
                        }
                    }
                }
                /* initialize waveconv matrix*/
                if(WAVETYPE==2||WAVETYPE==3){
                    if(FORWARD_ONLY==0){
                        for (j=1;j<=NY;j=j+IDY){
                            for (i=1;i<=NX;i=i+IDX){
                                waveconv_rho_z[j][i]=0.0;
                                waveconv_u_z[j][i]=0.0;
                                
                            }
                        }
                    }
                }
                
                if((EPRECOND>0)&&(EPRECOND_ITER==iter||(EPRECOND_ITER==0))){
                    for (j=1;j<=NY;j=j+IDY){
                        for (i=1;i<=NX;i=i+IDX){
                            if(WAVETYPE==1||WAVETYPE==3) We_sum[j][i]=0.0;
                            if(WAVETYPE==2||WAVETYPE==3) We_sum_SH[j][i]=0.0;
                        }
                    }
                }
                
                
                
                itestshot=TESTSHOT_START;
                swstestshot=0;
                
                if(INVTYPE==2){
                    if (RUN_MULTIPLE_SHOTS) nshots=nsrc; else nshots=1;
                    
                    /*------------------------------------------------------------------------------*/
                    /*----------- Start of loop over shots -----------------------------------------*/
                    /*------------------------------------------------------------------------------*/
                    
                    for (ishot=1;ishot<=nshots;ishot+=SHOTINC){
     
       			if (SEISMO && READREC==2){
			  if (ntr>0) {
			    dealloc_sections(ntr,ns,recpos_loc,sectionvx,sectionvy,sectionvz,sectionp,sectionpnp1,sectionpn,sectioncurl,sectiondiv,
					     sectionpdata,sectionpdiff,sectionpdiffold,sectionvxdata,sectionvxdiff,sectionvxdiffold,sectionvydata,
					     sectionvydiff,sectionvydiffold,sectionvzdata,sectionvzdiff,sectionvzdiffold);
			  }
			  free_imatrix(recpos,1,3,1,ntr_glob);
			  recpos=receiver(&ntr, srcpos, ishot);
			  recpos_loc = splitrec(recpos,&ntr_loc, ntr, recswitch);
			  ntr_glob=ntr;
			  ntr=ntr_loc;
			  if (ntr>0){
				alloc_sections(ntr,ns,&sectionvx,&sectionvy,&sectionvz,&sectionp,&sectionpnp1,&sectionpn,&sectioncurl,&sectiondiv,
					      &sectionpdata,&sectionpdiff,&sectionpdiffold,&sectionvxdata,&sectionvxdiff,&sectionvxdiffold,&sectionvydata,
					      &sectionvydiff,&sectionvydiffold,&sectionvzdata,&sectionvzdiff,&sectionvzdiffold);
			  }
			  if (ntr) group_id = 1;
			  else group_id = 0;
			  MPI_Comm_split(MPI_COMM_WORLD, group_id, MYID, &MPI_COMM_NTR);
			  MPI_Comm_rank(MPI_COMM_NTR, &myid_ntr);
			}

                        SOURCE_SHAPE = SOURCE_SHAPE_OLD;
                        if(WAVETYPE==2 || WAVETYPE==3) SOURCE_SHAPE_SH=SOURCE_SHAPE_OLD_SH;
                        
                        /*------------------------------------------------------------------------------*/
                        /*----------- Start of inversion of source time function -----------------------*/
                        /*------------------------------------------------------------------------------*/
                        
                        /* Do not Excute STF if this is a step length search run for Wolfe condition
                         * Therefore (gradient_optimization==1) is added.
                         */
                        
                        if(((INV_STF==1)&&( (iter==1) || (do_stf==1) )) && (gradient_optimization==1) ){
                            fprintf(FP,"\n==================================================================================\n");
                            fprintf(FP,"\n MYID=%d *****  Forward simulation for inversion of source time function ******** \n",MYID);
                            fprintf(FP,"\n MYID=%d * Starting simulation (forward model) for shot %d of %d. Iteration %d ** \n",MYID,ishot,nshots,iter);
                            fprintf(FP,"\n==================================================================================\n\n");
                            
                            for (nt=1;nt<=8;nt++) srcpos1[nt][1]=srcpos[nt][ishot];
                            
                            if (RUN_MULTIPLE_SHOTS){
                                /* find this single source positions on subdomains */
                                if (nsrc_loc>0) free_matrix(srcpos_loc,1,8,1,1);
                                srcpos_loc=splitsrc(srcpos1,&nsrc_loc, 1);
                            }else{
                                /* Distribute multiple source positions on subdomains */
                                srcpos_loc = splitsrc(srcpos,&nsrc_loc, nsrc);
                            }
                            
                            if(SOURCE_SHAPE==3) declare_error("SOURCE_SHAPE==3 isn't possible with INV_STF==1");
                            MPI_Barrier(MPI_COMM_WORLD);
                            
                            
                            /*-------------------*/
                            /* calculate wavelet */
                            /*-------------------*/
                            /* calculate wavelet for each source point P SV */
                            if(WAVETYPE==1||WAVETYPE==3){
                                signals=NULL;
                                signals=wavelet(srcpos_loc,nsrc_loc,ishot,0,1);
                            }
                            /* calculate wavelet for each source point SH */
                            if(WAVETYPE==2||WAVETYPE==3){
                                signals_SH=NULL;
                                signals_SH=wavelet(srcpos_loc,nsrc_loc,ishot,1,1);
                            }
                            
                            
                            /* initialize wavefield with zero */
                            if (L){
                                if(!ACOUSTIC) {
                                    zero_fdveps_visc(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,pvz,psxx,psyy,psxy,psxz,psyz,ux,uy,uxy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_sxz_x,psi_vxx,psi_vyx,psi_vzx,psi_syy_y,psi_sxy_y,psi_syz_y,psi_vyy,psi_vxy,psi_vzy,psi_vxxs,pr,pp,pq,pt,po);
                                } else {
                                    zero_fdveps_viscac(-nd+1, NY+nd, -nd+1, NX+nd, pvx, pvy, psp, pvxp1, pvyp1, psi_sxx_x, psi_sxy_x, psi_vxx, psi_vyx, psi_syy_y, psi_sxy_y, psi_vyy, psi_vxy, psi_vxxs, pp);
                                }
                            }else{
                                if(!ACOUSTIC)
                                    zero_fdveps(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,pvz,psxx,psyy,psxy,psxz,psyz,ux,uy,uxy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_sxz_x,psi_vxx,psi_vyx,psi_vzx,psi_syy_y,psi_sxy_y,psi_syz_y,psi_vyy,psi_vxy,psi_vzy,psi_vxxs);
                                else
                                    zero_fdveps_ac(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,psp,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs);
                            }
                            
                            if((!VERBOSE)&&(MYID==0)) fprintf(FP,"\n ****************************************\n ");
                            
                            /*------------------------------------------------------------------------------*/
                            /*----------------------  start loop over timesteps ( STF ) --------------------*/
                            /*------------------------------------------------------------------------------*/
                            
                            lsnap=iround(TSNAP1/DT);
                            lsamp=NDT;
                            nsnap=0;
                            
                            hin=1;
                            hin1=1;
                            
                            imat=1;
                            imat1=1;
                            imat2=1;
                            hi=1;
                            
                            for (nt=1;nt<=NT;nt++){
                                
                                infoout = !(nt%nt_out);
                                if((!VERBOSE)&&(MYID==0)) if(!(nt%(NT/40))) fprintf(FP,"*");
                                
                                /* Check if simulation is still stable P and SV */
                                if (WAVETYPE==1 || WAVETYPE==3) {
                                    if (isnan(pvy[NY/2][NX/2])) {
                                        fprintf(FP,"\n Time step: %d; pvy: %f \n",nt,pvy[NY/2][NX/2]);
                                        declare_error(" Simulation is unstable !");
                                    }
                                }
                                
                                /* Check if simulation is still stable SH */
                                if (WAVETYPE==2 || WAVETYPE==3) {
                                    if (isnan(pvz[NY/2][NX/2])) {
                                        fprintf(FP,"\n Time step: %d; pvy: %f \n",nt,pvy[NY/2][NX/2]);
                                        declare_error(" Simulation is unstable !");
                                    }
                                }
                                
                                if (MYID==0){
                                    if (infoout)  fprintf(FP,"\n Computing timestep %d of %d \n",nt,NT);
                                    time3=MPI_Wtime();
                                }
                                
                                /* update of particle velocities */
                                if(!ACOUSTIC) {
                                    if (WAVETYPE==1 || WAVETYPE==3) {
                                        update_v_PML(1, NX, 1, NY, nt, pvx, pvxp1, pvxm1, pvy, pvyp1, pvym1, uttx, utty, psxx, psyy, psxy, prip, prjp, srcpos_loc,signals,signals,nsrc_loc,absorb_coeff,hc,infoout,0, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_sxx_x, psi_syy_y, psi_sxy_y, psi_sxy_x);
                                    }
                                    
                                    if (WAVETYPE==2 || WAVETYPE==3) {
                                        update_v_PML_SH(1, NX, 1, NY, nt, pvz, pvzp1, pvzm1, psxz, psyz,prjp, srcpos_loc, signals, signals_SH, nsrc_loc, absorb_coeff,hc,infoout,0, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_sxz_x, psi_syz_y);
                                    }
                                } else {
                                    update_v_acoustic_PML(1, NX, 1, NY, nt, pvx, pvxp1, pvxm1, pvy, pvyp1, pvym1, psp, prip, prjp, srcpos_loc,signals,signals,nsrc_loc,absorb_coeff,hc,infoout,0, K_x_half, a_x_half, b_x_half, K_y_half, a_y_half, b_y_half, psi_sxx_x, psi_syy_y);
                                }

                                
                                if (MYID==0){
                                    time4=MPI_Wtime();
                                    time_av_v_update+=(time4-time3);
                                    if (infoout)  fprintf(FP," particle velocity exchange between PEs ...");
                                }
                                
                                /* exchange of particle velocities between PEs */
                                exchange_v(pvx,pvy,pvz, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, req_send, req_rec,wavetype_start);
                                
                                if (MYID==0){
                                    time5=MPI_Wtime();
                                    time_av_v_exchange+=(time5-time4);
                                    if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time5-time4);
                                }
                                
                                if (L) {   /* viscoelastic */
                                    if (WAVETYPE==1 || WAVETYPE==3) {
                                        if(!ACOUSTIC) {
                                            update_s_visc_PML(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppi, pu, puipjp, prho, hc, infoout, pr, pp, pq, fipjp, f, g, bip, bjm, cip, cjm, d, e, dip, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);
                                        }else{
                                            update_p_visc_PML(1, NX, 1, NY, pvx, pvy, psp, ppi, prho, hc, infoout, pp, g, bjm, cjm, e, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);
                                        }
                                    }
                                    if (WAVETYPE==2 || WAVETYPE==3) {
                                        update_s_visc_PML_SH(1, NX, 1, NY, pvz, psxz, psyz, pt, po, bip, bjm, cip, cjm, d, dip,fipjp, f, hc,infoout, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half,psi_vzx, psi_vzy);
                                    }
                                } else {   /* elastic */
                                    if (WAVETYPE==1 || WAVETYPE==3) {
                                        if(!ACOUSTIC) {
                                            update_s_elastic_PML(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppi, pu, puipjp, absorb_coeff, prho, hc, infoout, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);
                                        } else {
                                            update_p_PML(1, NX, 1, NY, pvx, pvy, psp, u, ppi, absorb_coeff, prho, hc, infoout, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);
                                        }
                                    }
                                    if (WAVETYPE==2 || WAVETYPE==3) {
                                        update_s_elastic_PML_SH(1, NX, 1, NY, pvz,psxz,psyz,uxz,uyz,hc,infoout, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half,psi_vzx, psi_vzy,puipjp,pu,prho);
                                    }
                                }
                                
                                /* explosive source */
                                if ((SOURCE_TYPE==1))
                                    psource(nt,psxx,psyy,psp,srcpos_loc,signals,nsrc_loc,0);
                                

                                /* Applying free surface condition */
                                if ((FREE_SURF) && (POS[2]==0)){
                                    if (!ACOUSTIC){
                                        if (L){
                                            /* viscoelastic */
                                            surface_PML(1, pvx, pvy, psxx, psyy, psxy,psyz, pp, pq, ppi, pu, prho, ptaup, ptaus, etajm, peta, hc, K_x, a_x, b_x, psi_vxxs, ux, uy,uxy,uyz,psxz,uxz);
                                        }else{
                                            /* elastic */
                                            surface_elastic_PML(1, pvx, pvy, psxx, psyy, psxy,psyz, ppi, pu, prho, hc, K_x, a_x, b_x, psi_vxxs, ux, uy, uxy,uyz,psxz,uxz);
                                        }
                                    } else {
                                        /* viscoelastic and elastic ACOUSTIC */
                                        surface_acoustic_PML(1, psp);
                                    }
                                }
                                
                                if (MYID==0){
                                    time6=MPI_Wtime();
                                    time_av_s_update+=(time6-time5);
                                    if (infoout)  fprintf(FP," stress exchange between PEs ...");
                                }
                                
                                
                                /* stress exchange between PEs */
                                if(!ACOUSTIC)
                                    exchange_s(psxx,psyy,psxy,psxz,psyz,bufferlef_to_rig, bufferrig_to_lef,buffertop_to_bot, bufferbot_to_top,req_send, req_rec,wavetype_start);
                                else
                                    exchange_p(psp,bufferlef_to_rig, bufferrig_to_lef,buffertop_to_bot, bufferbot_to_top,req_send, req_rec);
                                
                                if (MYID==0){
                                    time7=MPI_Wtime();
                                    time_av_s_exchange+=(time7-time6);
                                    if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time7-time6);
                                }
                                
                                /* store amplitudes at receivers in section-arrays */
                                if (SEISMO){
                                    seismo_ssg(nt, ntr, recpos_loc, sectionvx, sectionvy,sectionvz,sectionp, sectioncurl, sectiondiv,pvx, pvy,pvz, psxx, psyy, psp, ppi, pu, hc);
                                    /*lsamp+=NDT;*/
                                }
                                
                                if (MYID==0){
                                    time8=MPI_Wtime();
                                    time_av_timestep+=(time8-time3);
                                    if (infoout)  fprintf(FP," total real time for timestep %d : %4.2f s.\n",nt,time8-time3);
                                }
                                
                            }
                            
                            /*------------------------------------------------------------------------------*/
                            /*--------------------  End  of loop over timesteps (   STF   ) ----------------*/
                            /*------------------------------------------------------------------------------*/
                            
                            if((!VERBOSE)&&(MYID==0)) fprintf(FP,"\n");
                            
                            // Exchange measured seismogramms and save it to file
                            switch (SEISMO){
                                case 1 : 	/* particle velocities only */
                                    if (WAVETYPE==1 || WAVETYPE==3) {
                                        catseis(sectionvx, fulldata_vx, recswitch, ntr_glob, MPI_COMM_WORLD);
                                        catseis(sectionvy, fulldata_vy, recswitch, ntr_glob, MPI_COMM_WORLD);
                                    }
                                    if (WAVETYPE==2 || WAVETYPE==3) {
                                        catseis(sectionvz, fulldata_vz, recswitch, ntr_glob, MPI_COMM_WORLD);
                                    }
                                    if(LNORM==8){
                                        calc_envelope(fulldata_vy,fulldata_vy,ns,ntr_glob);
                                        calc_envelope(fulldata_vx,fulldata_vx,ns,ntr_glob);}
//                                     if (MYID==0) saveseis_glob(FP,fulldata_vx,fulldata_vy,fulldata_vz,fulldata_p,fulldata_curl,fulldata_div,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter,1);
                                    break;
                                    
                                case 2 :	/* pressure only */
                                    catseis(sectionp, fulldata_p, recswitch, ntr_glob, MPI_COMM_WORLD);
//                                     if (MYID==0) saveseis_glob(FP,fulldata_vx,fulldata_vy,fulldata_vz,fulldata_p,fulldata_curl,fulldata_div,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter,1);
                                    break;
                                    
                                case 3 : 	/* curl and div only */
                                    catseis(sectiondiv, fulldata_div, recswitch, ntr_glob, MPI_COMM_WORLD);
                                    catseis(sectioncurl, fulldata_curl, recswitch, ntr_glob, MPI_COMM_WORLD);
//                                     if (MYID==0) saveseis_glob(FP,fulldata_vx,fulldata_vy,fulldata_vz,fulldata_p,fulldata_curl,fulldata_div,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter,1);
                                    break;
                                    
                                case 4 :	/* everything */
                                    if (WAVETYPE==1 || WAVETYPE==3) {
                                        catseis(sectionvx, fulldata_vx, recswitch, ntr_glob, MPI_COMM_WORLD);
                                        catseis(sectionvy, fulldata_vy, recswitch, ntr_glob, MPI_COMM_WORLD);
                                    }
                                    if (WAVETYPE==2 || WAVETYPE==3) {
                                        catseis(sectionvz, fulldata_vz, recswitch, ntr_glob, MPI_COMM_WORLD);
                                    }
                                    catseis(sectionp, fulldata_p, recswitch, ntr_glob, MPI_COMM_WORLD);
                                    catseis(sectiondiv, fulldata_div, recswitch, ntr_glob, MPI_COMM_WORLD);
                                    catseis(sectioncurl, fulldata_curl, recswitch, ntr_glob, MPI_COMM_WORLD);
//                                     if (MYID==0) saveseis_glob(FP,fulldata_vx,fulldata_vy,fulldata_vz,fulldata_p,fulldata_curl,fulldata_div,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter,1);
                                    break;
                                    
                                case 5 :	/* everything except curl and div*/
                                    if (WAVETYPE==1 || WAVETYPE==3) {
                                        catseis(sectionvx, fulldata_vx, recswitch, ntr_glob, MPI_COMM_WORLD);
                                        catseis(sectionvy, fulldata_vy, recswitch, ntr_glob, MPI_COMM_WORLD);
                                    }
                                    if (WAVETYPE==2 || WAVETYPE==3) {
                                        catseis(sectionvz, fulldata_vz, recswitch, ntr_glob, MPI_COMM_WORLD);
                                    }
                                    catseis(sectionp, fulldata_p, recswitch, ntr_glob, MPI_COMM_WORLD);
//                                     if (MYID==0) saveseis_glob(FP,fulldata_vx,fulldata_vy,fulldata_vz,fulldata_p,fulldata_curl,fulldata_div,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter,1);
                                    break;
                                    
                            } /* end of switch (SEISMO) */
                            
                            /*------------------------------------------------------------------------------*/
                            /*----------- Start of inversion of source time function -----------------------*/
                            /*------------------------------------------------------------------------------*/
                            if((TIME_FILT==1) ||(TIME_FILT==2)){
                                
                                if (!FORWARD_ONLY){
                                    if((INV_STF==1)&&((iter==1)||(do_stf==1))){
                                        
                                        if (nsrc_loc>0){
                                            
                                            /*time domain filtering of the observed data sectionvy_obs */
                                            if(WAVETYPE==1 || WAVETYPE==3){
                                                if ((ADJOINT_TYPE==1)|| (ADJOINT_TYPE==2)){
                                                    inseis(fprec,ishot,sectionvy_obs,ntr_glob,ns,2,iter);
                                                    timedomain_filt(sectionvy_obs,F_LOW_PASS,ORDER,ntr_glob,ns,1);
                                                }
                                                if (ADJOINT_TYPE==4){
                                                    inseis(fprec,ishot,sectionp_obs,ntr_glob,ns,9,iter);
                                                    timedomain_filt(sectionp_obs,F_LOW_PASS,ORDER,ntr_glob,ns,1);
                                                }
                                            }
                                            
                                            if(WAVETYPE==2 || WAVETYPE==3){
                                                inseis(fprec,ishot,sectionvz_obs,ntr_glob,ns,10,iter);
                                                timedomain_filt(sectionvz_obs,F_LOW_PASS,ORDER,ntr_glob,ns,1);
                                            }
                                            
                                            printf("\n ====================================================================================================== \n");
                                            printf("\n Time Domain Filter is used for the inversion: lowpass filter, corner frequency of %.2f Hz, order %d\n",F_LOW_PASS,ORDER);
                                            printf("\n ====================================================================================================== \n");
                                            
                                            if(iter==1){
                                                printf("\n ====================================================== \n");
                                                printf("\n MYID = %d: STF inversion at first iteration \n",MYID);
                                            }
                                            else{
                                                printf("\n ================================================================================================ \n");
                                                printf("\n MYID = %d: STF inversion because of frequency step at the end of the last iteration \n",MYID);
                                            }
                                            
                                            if(WAVETYPE==1 || WAVETYPE==3){
                                                if ((ADJOINT_TYPE==1)|| (ADJOINT_TYPE==2)){
                                                    stf(FP,fulldata_vy,sectionvy_obs,sectionvy_conv,source_time_function,recpos,recpos_loc,ntr_glob,ntr,srcpos,ishot,ns,iter,nsrc_glob,F_LOW_PASS,0,nsrc_glob);
                                                }
                                                if (ADJOINT_TYPE==4){
                                                    stf(FP,fulldata_p,sectionp_obs,sectionp_conv,source_time_function,recpos,recpos_loc,ntr_glob,ntr,srcpos,ishot,ns,iter,nsrc_glob,F_LOW_PASS,0,nsrc_glob);
                                                }
                                            }
                                            
                                            if(WAVETYPE==2 || WAVETYPE==3){
                                                stf(FP,fulldata_vz,sectionvz_obs,sectionvz_conv,source_time_function,recpos,recpos_loc,ntr_glob,ntr,srcpos,ishot,ns,iter,nsrc_glob,F_LOW_PASS,1,nsrc_glob);
                                            }
                                            

                                        }
                                    }
                                }
                            
                            } else {
                                
                                if (FORWARD_ONLY==0){
                                    if((INV_STF==1)&&(iter==N_STF_START)){
                                        
                                        if(ishot==nshots){
                                            N_STF_START=N_STF_START+N_STF;
                                        }
                                        
                                        if (nsrc_loc>0){
                                            printf("\n ====================================================== \n");
                                            printf("\n MYID = %d: STF inversion due to the increment N_STF \n",MYID);
                                            
                                            if(WAVETYPE==1 || WAVETYPE==3){
                                                if ((ADJOINT_TYPE==1)|| (ADJOINT_TYPE==2)){
                                                    inseis(fprec,ishot,sectionvy_obs,ntr_glob,ns,2,iter);
                                                }
                                                if (ADJOINT_TYPE==4){
                                                    inseis(fprec,ishot,sectionp_obs,ntr_glob,ns,9,iter);
                                                }
                                                
                                                if ((ADJOINT_TYPE==1)|| (ADJOINT_TYPE==2)){
                                                    stf(FP,fulldata_vy,sectionvy_obs,sectionvy_conv,source_time_function,recpos,recpos_loc,ntr_glob,ntr,srcpos,ishot,ns,iter,nsrc_glob,F_LOW_PASS,0,nsrc_glob);
                                                }
                                                if (ADJOINT_TYPE==4){
                                                    stf(FP,fulldata_p,sectionp_obs,sectionp_conv,source_time_function,recpos,recpos_loc,ntr_glob,ntr,srcpos,ishot,ns,iter,nsrc_glob,F_LOW_PASS,0,nsrc_glob);
                                                }
                                            }
                                            if(WAVETYPE==2 || WAVETYPE==3){
                                                inseis(fprec,ishot,sectionvz_obs,ntr_glob,ns,10,iter);
                                                stf(FP,fulldata_vz,sectionvz_obs,sectionvz_conv,source_time_function,recpos,recpos_loc,ntr_glob,ntr,srcpos,ishot,ns,iter,nsrc_glob,F_LOW_PASS,1,nsrc_glob);
                                            }
                                        }
                                    }
                                }
                            }
                            
                            
                            MPI_Barrier(MPI_COMM_WORLD);
                            
                            
                        }
                        
                        /*------------------------------------------------------------------------------*/
                        /*----------- End of inversion of source time function -------------------------*/
                        /*------------------------------------------------------------------------------*/
                        
                        fprintf(FP,"\n==================================================================================\n");
                        fprintf(FP,"\n MYID=%d * Starting simulation (forward model) for shot %d of %d. Iteration %d ** \n",MYID,ishot,nshots,iter);
                        fprintf(FP,"\n==================================================================================\n");
                        
                        for (nt=1;nt<=8;nt++) srcpos1[nt][1]=srcpos[nt][ishot];
                        
                        /*-----------------------------------*/
                        /* determine source position on grid */
                        /*-----------------------------------*/
                        if (RUN_MULTIPLE_SHOTS){
                            /* find this single source positions on subdomains */
                            if (nsrc_loc>0) free_matrix(srcpos_loc,1,8,1,1);
                            srcpos_loc=splitsrc(srcpos1,&nsrc_loc, 1);
                        }else{
                            /* Distribute multiple source positions on subdomains */
                            srcpos_loc = splitsrc(srcpos,&nsrc_loc, nsrc);
                        }
                        
                        /*-------------------*/
                        /*  Use STF wavelet  */
                        /*-------------------*/
                        if(INV_STF){
                            
                            SOURCE_SHAPE=7;
                            if( WAVETYPE==1 || WAVETYPE==3 ) fprintf(FP,"\n Using optimized source time function located in %s.shot%d \n",SIGNAL_FILE,ishot);
                            
                            if( WAVETYPE==2 || WAVETYPE==3 ) {
                                SOURCE_SHAPE_SH=7;
                                fprintf(FP,"\n Using optimized source time function located in %s.shot%d  \n",SIGNAL_FILE_SH,ishot);
                            }
                            
                        }
                        
                        MPI_Barrier(MPI_COMM_WORLD);
                        
                        /*-------------------*/
                        /* calculate wavelet */
                        /*-------------------*/
                        /* calculate wavelet for each source point P SV */
                        if(WAVETYPE==1||WAVETYPE==3){
                            signals=NULL;
                            signals=wavelet(srcpos_loc,nsrc_loc,ishot,0,0);
                        }
                        /* calculate wavelet for each source point SH */
                        if(WAVETYPE==2||WAVETYPE==3){
                            signals_SH=NULL;
                            signals_SH=wavelet(srcpos_loc,nsrc_loc,ishot,1,0);
                        }
                        
                        /*------------------------------------------------------------------------------*/
                        /*----------- Start of Time Domain Filtering -----------------------------------*/
                        /*------------------------------------------------------------------------------*/
                        
                        if (((TIME_FILT==1) || (TIME_FILT==2)) && (SOURCE_SHAPE!=6) && (INV_STF==0)){
                            fprintf(FP,"\n Time Domain Filter applied: Lowpass with corner frequency of %.2f Hz, order %d\n",F_LOW_PASS,ORDER);
                            
                            /*time domain filtering of the source signal */
                            if(WAVETYPE==1||WAVETYPE==3) timedomain_filt(signals,F_LOW_PASS,ORDER,nsrc_loc,ns,1);
                            if(WAVETYPE==2||WAVETYPE==3) timedomain_filt(signals_SH,F_LOW_PASS,ORDER,nsrc_loc,ns,1);
                            
                        }
                        /*------------------------------------------------------------------------------*/
                        /*----------- End of Time Domain Filtering -------------------------------------*/
                        /*------------------------------------------------------------------------------*/
                        
                        MPI_Barrier(MPI_COMM_WORLD);
                        
                        /* initialize wavefield with zero */
                        if (L){
                            if(!ACOUSTIC) {
                                zero_fdveps_visc(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,pvz,psxx,psyy,psxy,psxz,psyz,ux,uy,uxy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_sxz_x,psi_vxx,psi_vyx,psi_vzx,psi_syy_y,psi_sxy_y,psi_syz_y,psi_vyy,psi_vxy,psi_vzy,psi_vxxs,pr,pp,pq,pt,po);
                            } else {
                                zero_fdveps_viscac(-nd+1, NY+nd, -nd+1, NX+nd, pvx, pvy, psp, pvxp1, pvyp1, psi_sxx_x, psi_sxy_x, psi_vxx, psi_vyx, psi_syy_y, psi_sxy_y, psi_vyy, psi_vxy, psi_vxxs, pp); }
                        }else{
                            if(!ACOUSTIC)
                                zero_fdveps(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,pvz,psxx,psyy,psxy,psxz,psyz,ux,uy,uxy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_sxz_x,psi_vxx,psi_vyx,psi_vzx,psi_syy_y,psi_sxy_y,psi_syz_y,psi_vyy,psi_vxy,psi_vzy,psi_vxxs);
                            else
                                zero_fdveps_ac(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,psp,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs);
                        }
                        
                        /*initialize gradient matrices for each shot with zeros PSV*/
                        if(WAVETYPE==1 || WAVETYPE==3) {
                            if(FORWARD_ONLY==0){
                                for(j=1;j<=NY;j=j+IDY){
                                    for(i=1;i<=NX;i=i+IDX){
                                        waveconv_shot[j][i]=0.0;
                                        waveconv_rho_shot[j][i]=0.0;
                                    }
                                }
                                if(!ACOUSTIC){
                                    for(j=1;j<=NY;j=j+IDY){
                                        for(i=1;i<=NX;i=i+IDX){
                                            waveconv_u_shot[j][i]=0.0;
                                        }
                                    }
                                }
                            }
                        }
                        /*initialize gradient matrices for each shot with zeros SH*/
                        if(WAVETYPE==2 || WAVETYPE==3){
                            if(FORWARD_ONLY==0){
                                for(j=1;j<=NY;j=j+IDY){
                                    for(i=1;i<=NX;i=i+IDX){
                                        waveconv_rho_shot_z[j][i]=0.0;
                                        waveconv_u_shot_z[j][i]=0.0;
                                    }
                                }
                                
                            }
                        }
                        
                        if((EPRECOND==1)||(EPRECOND==3)&&(EPRECOND_ITER==iter||(EPRECOND_ITER==0))){
                            for(j=1;j<=NY;j=j+IDY){
                                for(i=1;i<=NX;i=i+IDX){
                                    if(WAVETYPE==1 || WAVETYPE==3){
                                        Ws[j][i]=0.0;
                                        Wr[j][i]=0.0;
                                        We[j][i]=0.0;
                                    }
                                    if(WAVETYPE==2 || WAVETYPE==3){
                                        Ws_SH[j][i]=0.0;
                                        Wr_SH[j][i]=0.0;
                                        We_SH[j][i]=0.0;
                                    }
                                }
                            }
                        }
                        
                        
                        lsnap=iround(TSNAP1/DT); lsamp=NDT; nsnap=0;
                        hin=1; hin1=1;
                        imat=1; imat1=1; imat2=1; hi=1;
                        
                        if((!VERBOSE)&&(MYID==0)) fprintf(FP,"\n ****************************************\n ");
                        
                        
                        /*------------------------------------------------------------------------------*/
                        /*----------------------  start loop over timesteps (forward model) ------------*/
                        /*------------------------------------------------------------------------------*/
                        for (nt=1;nt<=NT;nt++){
                            
                            // Ratio to give output to stout
                            infoout = !(nt%nt_out);
                            
                            if((!VERBOSE)&&(MYID==0)) if(!(nt%(NT/40))) fprintf(FP,"*");
                            
                            /* Check if simulation is still stable P and SV */
                            if (WAVETYPE==1 || WAVETYPE==3) {
                                if (isnan(pvy[NY/2][NX/2])) {
                                    fprintf(FP,"\n Time step: %d; pvy: %f \n",nt,pvy[NY/2][NX/2]);
                                    declare_error(" Simulation is unstable !");
                                }
                            }
                            
                            /* Check if simulation is still stable SH */
                            if (WAVETYPE==2 || WAVETYPE==3) {
                                if (isnan(pvz[NY/2][NX/2])) {
                                    fprintf(FP,"\n Time step: %d; pvy: %f \n",nt,pvy[NY/2][NX/2]);
                                    declare_error(" Simulation is unstable !");
                                }
                            }
                            
                            
                            
                            if (MYID==0){
                                if (infoout)  fprintf(FP,"\n Computing timestep %d of %d \n",nt,NT);
                                time3=MPI_Wtime();
                            }
                            
                            
                            /* update of particle velocities */
                            if(!ACOUSTIC) {
                                if (WAVETYPE==1 || WAVETYPE==3) {
                                    update_v_PML(1, NX, 1, NY, nt, pvx, pvxp1, pvxm1, pvy, pvyp1, pvym1, uttx, utty, psxx, psyy, psxy, prip, prjp, srcpos_loc,signals,signals,nsrc_loc,absorb_coeff,hc,infoout,0, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_sxx_x, psi_syy_y, psi_sxy_y, psi_sxy_x);
                                }
                                
                                if (WAVETYPE==2 || WAVETYPE==3) {
                                    update_v_PML_SH(1, NX, 1, NY, nt, pvz, pvzp1, pvzm1, psxz, psyz,prjp, srcpos_loc, signals, signals_SH, nsrc_loc, absorb_coeff,hc,infoout,0, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_sxz_x, psi_syz_y);
                                }
                            } else {
                                update_v_acoustic_PML(1, NX, 1, NY, nt, pvx, pvxp1, pvxm1, pvy, pvyp1, pvym1, psp, prip, prjp, srcpos_loc,signals,signals,nsrc_loc,absorb_coeff,hc,infoout,0, K_x_half, a_x_half, b_x_half, K_y_half, a_y_half, b_y_half, psi_sxx_x, psi_syy_y);
                            }
                            
                            if (MYID==0){
                                time4=MPI_Wtime();
                                time_av_v_update+=(time4-time3);
                                if (infoout)  fprintf(FP," particle velocity exchange between PEs ...");
                            }
                            
                            /* exchange of particle velocities between PEs */
                            exchange_v(pvx,pvy,pvz, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, req_send, req_rec,wavetype_start);
                            
                            if (MYID==0){
                                time5=MPI_Wtime();
                                time_av_v_exchange+=(time5-time4);
                                if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time5-time4);
                            }
                            
                            if (L) {   /* viscoelastic */
                                if (WAVETYPE==1 || WAVETYPE==3) {
                                    if(!ACOUSTIC) {
                                        update_s_visc_PML(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppi, pu, puipjp, prho, hc, infoout, pr, pp, pq, fipjp, f, g, bip, bjm, cip, cjm, d, e, dip, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);
                                    }else{
                                        update_p_visc_PML(1, NX, 1, NY, pvx, pvy, psp, ppi, prho, hc, infoout, pp, g, bjm, cjm, e, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);
                                    }
                                }
                                if (WAVETYPE==2 || WAVETYPE==3) {
                                    update_s_visc_PML_SH(1, NX, 1, NY, pvz, psxz, psyz, pt, po, bip, bjm, cip, cjm, d, dip,fipjp, f, hc,infoout, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half,psi_vzx, psi_vzy);
                                }
                            } else {   /* elastic */
                                if (WAVETYPE==1 || WAVETYPE==3) {
                                    if(!ACOUSTIC) {
                                        update_s_elastic_PML(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppi, pu, puipjp, absorb_coeff, prho, hc, infoout, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);
                                    } else {
                                        update_p_PML(1, NX, 1, NY, pvx, pvy, psp, u, ppi, absorb_coeff, prho, hc, infoout, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);
                                    }
                                }
                                if (WAVETYPE==2 || WAVETYPE==3) {
                                    update_s_elastic_PML_SH(1, NX, 1, NY, pvz,psxz,psyz,uxz,uyz,hc,infoout, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half,psi_vzx, psi_vzy,puipjp,pu,prho);
                                }
                            }
                            
                            
                            
                            /* explosive source */
                            if ((SOURCE_TYPE==1)&&(WAVETYPE==1||WAVETYPE==3))
                                psource(nt,psxx,psyy,psp,srcpos_loc,signals,nsrc_loc,0);
                            
                            /* Applying free surface condition */
                            if ((FREE_SURF) && (POS[2]==0)){
                                if (!ACOUSTIC){
                                    if (L){
                                        /* viscoelastic */
                                        surface_PML(1, pvx, pvy, psxx, psyy, psxy,psyz, pp, pq, ppi, pu, prho, ptaup, ptaus, etajm, peta, hc, K_x, a_x, b_x, psi_vxxs, ux, uy,uxy,uyz,psxz,uxz);
                                    }else{
                                        /* elastic */
                                        surface_elastic_PML(1, pvx, pvy, psxx, psyy, psxy,psyz, ppi, pu, prho, hc, K_x, a_x, b_x, psi_vxxs, ux, uy, uxy,uyz,psxz,uxz);
                                    }
                                } else {
                                    /* viscoelastic and elastic ACOUSTIC */
                                    surface_acoustic_PML(1, psp);
                                }
                            }
                            
                            if (MYID==0){
                                time6=MPI_Wtime();
                                time_av_s_update+=(time6-time5);
                                if (infoout)  fprintf(FP," stress exchange between PEs ...");
                            }
                            
                            /* stress exchange between PEs */
                            if(!ACOUSTIC)
                                exchange_s(psxx,psyy,psxy,psxz,psyz, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, req_send, req_rec,wavetype_start);
                            else
                                exchange_p(psp,bufferlef_to_rig, bufferrig_to_lef,buffertop_to_bot, bufferbot_to_top,req_send, req_rec);
                            
                            
                            if (MYID==0){
                                time7=MPI_Wtime();
                                time_av_s_exchange+=(time7-time6);
                                if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time7-time6);
                            }
                            
                            /* store amplitudes at receivers in section-arrays */
                            if (SEISMO){
                                seismo_ssg(nt, ntr, recpos_loc, sectionvx, sectionvy,sectionvz, sectionp, sectioncurl, sectiondiv, pvx, pvy,pvz, psxx, psyy, psp, ppi, pu, hc);
                            }
                            
                            /* save snapshots from forward model */
                            if(nt==hin1){
                                if(FORWARD_ONLY==0){
                                    
                                    if(WAVETYPE==1||WAVETYPE==3){
                                        for (j=1;j<=NY;j=j+IDYI){
                                            for (i=1;i<=NX;i=i+IDXI){
                                                forward_prop_rho_x[j][i][hin]=pvxp1[j][i];
                                                forward_prop_rho_y[j][i][hin]=pvyp1[j][i];
                                            }
                                        }
                                        
                                        if(!ACOUSTIC){
                                            for (j=1;j<=NY;j=j+IDYI){
                                                for (i=1;i<=NX;i=i+IDXI){
                                                    if(VELOCITY==0){
                                                        forward_prop_x[j][i][hin]=psxx[j][i];
                                                        forward_prop_y[j][i][hin]=psyy[j][i];
                                                    }else{
                                                        forward_prop_x[j][i][hin]=ux[j][i];
                                                        forward_prop_y[j][i][hin]=uy[j][i];
                                                    }
                                                }
                                            }
                                            for (j=1;j<=NY;j=j+IDYI){
                                                for (i=1;i<=NX;i=i+IDXI){
                                                    if(VELOCITY==0){
                                                        forward_prop_u[j][i][hin]=psxy[j][i];
                                                    }else{
                                                        forward_prop_u[j][i][hin]=uxy[j][i];
                                                    }
                                                }
                                            }
                                        }else{
                                            for (j=1;j<=NY;j=j+IDYI){
                                                for (i=1;i<=NX;i=i+IDXI){
                                                    if(VELOCITY==0){
                                                        forward_prop_p[j][i][hin]=psp[j][i];
                                                    }else{
                                                        forward_prop_p[j][i][hin]=u[j][i];
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    if(WAVETYPE==2||WAVETYPE==3){
                                        for (j=1;j<=NY;j=j+IDYI){
                                            for (i=1;i<=NX;i=i+IDXI){
                                                forward_prop_rho_z[j][i][hin]=pvzp1[j][i];
                                            }
                                        }
                                        
                                        if(!ACOUSTIC){
                                            for (j=1;j<=NY;j=j+IDYI){
                                                for (i=1;i<=NX;i=i+IDXI){
                                                    if(VELOCITY==0){
                                                        forward_prop_z_xz[j][i][hin]=psxz[j][i];
                                                        forward_prop_z_yz[j][i][hin]=psyz[j][i];
                                                    }else{
                                                        forward_prop_z_xz[j][i][hin]=uxz[j][i];
                                                        forward_prop_z_yz[j][i][hin]=uyz[j][i];
                                                    }
                                                }
                                            }
                                            
                                        }
                                    }
                                    hin++;
                                    hin1=hin1+DTINV;
                                }
                                DTINV_help[nt]=1;
                            }
                            
                            
                            if((EPRECOND==1)||(EPRECOND==3)&&(EPRECOND_ITER==iter||(EPRECOND_ITER==0))){
                                if (WAVETYPE==1 || WAVETYPE==3) eprecond(Ws,pvx,pvy);
                                if (WAVETYPE==2 || WAVETYPE==3) eprecond_SH(Ws_SH,pvz);
                            }
                            
                            /* WRITE SNAPSHOTS TO DISK */
                            if ((SNAP) && (nt==lsnap) && (nt<=TSNAP2/DT)){
                                ++nsnap;
                                if(WAVETYPE==1||WAVETYPE==3) snap(FP,nt,nsnap,pvx,pvy,psxx,psyy,psp,pu,ppi,hc,ishot);
                                if(WAVETYPE==2||WAVETYPE==3) snap_SH(FP,nt,nsnap,pvz,pu,ppi,hc,ishot);
                                lsnap=lsnap+iround(TSNAPINC/DT);
                            }
                            
                            
                            if (MYID==0){
                                time8=MPI_Wtime();
                                time_av_timestep+=(time8-time3);
                                if (infoout)  fprintf(FP," total real time for timestep %d : %4.2f s.\n",nt,time8-time3);
                            }
                            
                            
                        }

                        /*------------------------------------------------------------------------------*/
                        /*--------------------  End  of loop over timesteps (forward model) ------------*/
                        /*------------------------------------------------------------------------------*/
                        
                        if((!VERBOSE)&&(MYID==0)) fprintf(FP,"\n");
                        
                        // Exchange measured seismogramms and save it to file
                        switch (SEISMO){
                            case 1 : 	/* particle velocities only */
                                if (WAVETYPE==1 || WAVETYPE==3) {
                                    catseis(sectionvx, fulldata_vx, recswitch, ntr_glob, MPI_COMM_WORLD);
                                    catseis(sectionvy, fulldata_vy, recswitch, ntr_glob, MPI_COMM_WORLD);
                                }
                                if (WAVETYPE==2 || WAVETYPE==3) {
                                    catseis(sectionvz, fulldata_vz, recswitch, ntr_glob, MPI_COMM_WORLD);
                                }
                                if(LNORM==8){
                                    calc_envelope(fulldata_vy,fulldata_vy,ns,ntr_glob);
                                    calc_envelope(fulldata_vx,fulldata_vx,ns,ntr_glob);}
                                if (MYID==0) saveseis_glob(FP,fulldata_vx,fulldata_vy,fulldata_vz,fulldata_p,fulldata_curl,fulldata_div,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter,1);
                                break;
                                
                            case 2 :	/* pressure only */
                                catseis(sectionp, fulldata_p, recswitch, ntr_glob, MPI_COMM_WORLD);
                                if (MYID==0) saveseis_glob(FP,fulldata_vx,fulldata_vy,fulldata_vz,fulldata_p,fulldata_curl,fulldata_div,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter,1);
                                break;
                                
                            case 3 : 	/* curl and div only */
                                catseis(sectiondiv, fulldata_div, recswitch, ntr_glob, MPI_COMM_WORLD);
                                catseis(sectioncurl, fulldata_curl, recswitch, ntr_glob, MPI_COMM_WORLD);
                                if (MYID==0) saveseis_glob(FP,fulldata_vx,fulldata_vy,fulldata_vz,fulldata_p,fulldata_curl,fulldata_div,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter,1);
                                break;
                                
                            case 4 :	/* everything */
                                if (WAVETYPE==1 || WAVETYPE==3) {
                                    catseis(sectionvx, fulldata_vx, recswitch, ntr_glob, MPI_COMM_WORLD);
                                    catseis(sectionvy, fulldata_vy, recswitch, ntr_glob, MPI_COMM_WORLD);
                                }
                                if (WAVETYPE==2 || WAVETYPE==3) {
                                    catseis(sectionvz, fulldata_vz, recswitch, ntr_glob, MPI_COMM_WORLD);
                                }
                                catseis(sectionp, fulldata_p, recswitch, ntr_glob, MPI_COMM_WORLD);
                                catseis(sectiondiv, fulldata_div, recswitch, ntr_glob, MPI_COMM_WORLD);
                                catseis(sectioncurl, fulldata_curl, recswitch, ntr_glob, MPI_COMM_WORLD);
                                if (MYID==0) saveseis_glob(FP,fulldata_vx,fulldata_vy,fulldata_vz,fulldata_p,fulldata_curl,fulldata_div,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter,1);
                                break;
                                
                            case 5 :	/* everything except curl and div*/
                                if (WAVETYPE==1 || WAVETYPE==3) {
                                    catseis(sectionvx, fulldata_vx, recswitch, ntr_glob, MPI_COMM_WORLD);
                                    catseis(sectionvy, fulldata_vy, recswitch, ntr_glob, MPI_COMM_WORLD);
                                }
                                if (WAVETYPE==2 || WAVETYPE==3) {
                                    catseis(sectionvz, fulldata_vz, recswitch, ntr_glob, MPI_COMM_WORLD);
                                }
                                catseis(sectionp, fulldata_p, recswitch, ntr_glob, MPI_COMM_WORLD);
                                if (MYID==0) saveseis_glob(FP,fulldata_vx,fulldata_vy,fulldata_vz,fulldata_p,fulldata_curl,fulldata_div,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter,1);
                                break;
                                
                        } /* end of switch (SEISMO) */
                        
                        /*---------------------------------------------------------------*/
                        /*----------------------  start of inversion process ------------*/
                        /*---------------------------------------------------------------*/
                        if(FORWARD_ONLY==0){
                            
                            /*-----------------------------------*/
                            /*------- Calculate residuals -------*/
                            /*-----------------------------------*/
                            
                            if (MYID==0&&VERBOSE){
                                printf("-------------------  \n");
                                printf("Calculate residuals  \n");
                                printf("-------------------  \n");
                            }
                            
                            if ((ntr > 0)){
                                
                                /* calculate L2-Norm and energy ? */
                                if((ishot==itestshot)&&(ishot<=TESTSHOT_END)){swstestshot=1;}
                                
                                if(WAVETYPE==1 || WAVETYPE==3){
                                    /* --------------------------------- */
                                    /* read seismic data from SU file vx */
                                    /* --------------------------------- */
                                    
                                    if((ADJOINT_TYPE==1)||(ADJOINT_TYPE==3)){ /* if ADJOINT_TYPE */
                                        inseis(fprec,ishot,sectionread,ntr_glob,ns,1,iter);
                                        if ((TIME_FILT==1 )|| (TIME_FILT==2)){
                                            timedomain_filt(sectionread,F_LOW_PASS,ORDER,ntr_glob,ns,1);
                                        }
                                        h=1;
                                        for(i=1;i<=ntr;i++){
                                            for(j=1;j<=ns;j++){
                                                sectionvxdata[h][j]=sectionread[recpos_loc[3][i]][j];
                                            }
                                            h++;
                                        }
                                        L2=calc_res(sectionvxdata,sectionvx,sectionvxdiff,sectionvxdiffold,ntr,ns,LNORM,L2,0,1,swstestshot,ntr_glob,recpos_loc,nsrc_glob,ishot,iter,srcpos,recpos);
                                        if(swstestshot==1){energy=calc_energy(sectionvxdata,ntr,ns,energy, ntr_glob, recpos_loc, nsrc_glob, ishot,iter,srcpos,recpos);}
                                        energy_all_shots=calc_energy(sectionvxdata,ntr,ns,energy_all_shots, ntr_glob, recpos_loc, nsrc_glob, ishot,iter,srcpos,recpos);
                                        L2_all_shots=calc_misfit(sectionvxdata,sectionvx,ntr,ns,LNORM,L2_all_shots,0,1,1, ntr_glob, recpos_loc, nsrc_glob, ishot,iter,srcpos,recpos);
                                        /*fprintf(FP,"Energy vxdata for PE %d:   %f\n\n", MYID,energy);*/
                                    } /* end ADJOINT_TYPE */
                                    
                                    /* --------------------------------- */
                                    /* read seismic data from SU file vy */
                                    /* --------------------------------- */
                                    
                                    if((ADJOINT_TYPE==1)||(ADJOINT_TYPE==2)){ /* if ADJOINT_TYPE */
                                        inseis(fprec,ishot,sectionread,ntr_glob,ns,2,iter);
                                        if ((TIME_FILT==1 )|| (TIME_FILT==2)){
                                            timedomain_filt(sectionread,F_LOW_PASS,ORDER,ntr_glob,ns,1);
                                        }
                                        h=1;
                                        for(i=1;i<=ntr;i++){
                                            for(j=1;j<=ns;j++){
                                                sectionvydata[h][j]=sectionread[recpos_loc[3][i]][j];
                                            }
                                            h++;
                                        }
                                        L2=calc_res(sectionvydata,sectionvy,sectionvydiff,sectionvydiffold,ntr,ns,LNORM,L2,0,1,swstestshot,ntr_glob,recpos_loc,nsrc_glob,ishot,iter,srcpos,recpos);
                                        if(swstestshot==1){energy=calc_energy(sectionvydata,ntr,ns,energy, ntr_glob, recpos_loc, nsrc_glob, ishot,iter,srcpos,recpos);}
                                        energy_all_shots=calc_energy(sectionvydata,ntr,ns,energy_all_shots, ntr_glob, recpos_loc, nsrc_glob, ishot,iter,srcpos,recpos);
                                        L2_all_shots=calc_misfit(sectionvydata,sectionvy,ntr,ns,LNORM,L2_all_shots,0,1,1, ntr_glob, recpos_loc, nsrc_glob, ishot,iter,srcpos,recpos);
                                        /*fprintf(FP,"Energy vydata for PE %d:   %f\n\n", MYID,energy);	*/
                                    } /* end ADJOINT_TYPE */
                                    
                                    /* --------------------------------- */
                                    /* read seismic data from SU file p */
                                    /* --------------------------------- */
                                    if(ADJOINT_TYPE==4){ /* if ADJOINT_TYPE */
                                        inseis(fprec,ishot,sectionread,ntr_glob,ns,9,iter);
                                        if ((TIME_FILT==1 )|| (TIME_FILT==2)){
                                            timedomain_filt(sectionread,F_LOW_PASS,ORDER,ntr_glob,ns,1);
                                        }
                                        h=1;
                                        for(i=1;i<=ntr;i++){
                                            for(j=1;j<=ns;j++){
                                                sectionpdata[h][j]=sectionread[recpos_loc[3][i]][j];
                                            }
                                            h++;
                                        }
                                        L2=calc_res(sectionpdata,sectionp,sectionpdiff,sectionpdiffold,ntr,ns,LNORM,L2,0,1,swstestshot,ntr_glob,recpos_loc,nsrc_glob,ishot,iter,srcpos,recpos);
                                        if(swstestshot==1){energy=calc_energy(sectionpdata,ntr,ns,energy, ntr_glob, recpos_loc, nsrc_glob, ishot,iter,srcpos,recpos);}
                                        energy_all_shots=calc_energy(sectionpdata,ntr,ns,energy_all_shots, ntr_glob, recpos_loc, nsrc_glob, ishot,iter,srcpos,recpos);
                                        L2_all_shots=calc_misfit(sectionpdata,sectionp,ntr,ns,LNORM,L2_all_shots,0,1,1, ntr_glob, recpos_loc, nsrc_glob, ishot,iter,srcpos,recpos);
                                    } /* end ADJOINT_TYPE */
                                }
                                
                                /* --------------------------------- */
                                /* read seismic data from SU file vz */
                                /* --------------------------------- */
                                if(WAVETYPE==2 || WAVETYPE==3){
                                    inseis(fprec,ishot,sectionread,ntr_glob,ns,10,iter);
                                    if ((TIME_FILT==1 )|| (TIME_FILT==2)){
                                        timedomain_filt(sectionread,F_LOW_PASS,ORDER,ntr_glob,ns,1);
                                    }
                                    h=1;
                                    for(i=1;i<=ntr;i++){
                                        for(j=1;j<=ns;j++){
                                            sectionvzdata[h][j]=sectionread[recpos_loc[3][i]][j];
                                        }
                                        h++;
                                    }
                                    L2_SH=calc_res(sectionvzdata,sectionvz,sectionvzdiff,sectionvzdiffold,ntr,ns,LNORM,L2_SH,0,1,swstestshot,ntr_glob,recpos_loc,nsrc_glob,ishot,iter,srcpos,recpos);
                                    if(swstestshot==1){energy_SH=calc_energy(sectionvzdata,ntr,ns,energy_SH, ntr_glob, recpos_loc, nsrc_glob, ishot,iter,srcpos,recpos);}
                                    energy_all_shots_SH=calc_energy(sectionvzdata,ntr,ns,energy_all_shots_SH, ntr_glob, recpos_loc, nsrc_glob, ishot,iter,srcpos,recpos);
                                    L2_all_shots_SH=calc_misfit(sectionvzdata,sectionvz,ntr,ns,LNORM,L2_all_shots_SH,0,1,1, ntr_glob, recpos_loc, nsrc_glob, ishot,iter,srcpos,recpos);
                                }
                                
                                // Tracekill
                                if (TRKILL){
                                    count_killed_traces(ntr,swstestshot,ntr_glob,recpos_loc,nsrc_glob,ishot,ptr_killed_traces,ptr_killed_traces_testshots,srcpos,recpos);
                                }
                                                                
                                if((ishot==itestshot)&&(ishot<=TESTSHOT_END)){
                                    swstestshot=0;
                                    itestshot+=TESTSHOT_INCR;
                                }
                                
                                /* Write differences between measured and synthetic seismogramms (adjoint sources) to disk */
                                if (SEISMO&&VERBOSE){
                                    if(WAVETYPE==1 || WAVETYPE==3){
                                        if ((ADJOINT_TYPE==1)|| (ADJOINT_TYPE==3)){
                                        catseis(sectionvxdiff, fulldata_vx, recswitch, ntr_glob, MPI_COMM_NTR);
                                        }
                                        if ((ADJOINT_TYPE==1)|| (ADJOINT_TYPE==2)){
                                        catseis(sectionvydiff, fulldata_vy, recswitch, ntr_glob, MPI_COMM_NTR);
                                        }
                                        if (ADJOINT_TYPE==4){
                                            catseis(sectionpdiff, fulldata_p, recswitch, ntr_glob, MPI_COMM_NTR);
                                    }
                                    }
                                    if(WAVETYPE==2 || WAVETYPE==3){
                                        catseis(sectionvzdiff, fulldata_vz, recswitch, ntr_glob, MPI_COMM_NTR);
                                    }
                                    if(myid_ntr==0){
                                        saveseis_glob(FP,fulldata_vx,fulldata_vy,fulldata_vz,fulldata_p,sectionpdiff,sectionpdiff,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter,2);
                                    }
                                }
                                
                                /* Write measured filtered seismogramms to disk */
                                if (SEISMO && TIME_FILT && WRITE_FILTERED_DATA){
                                    if(WAVETYPE==1 || WAVETYPE==3){
                                        if ((ADJOINT_TYPE==1)|| (ADJOINT_TYPE==3)){
                                        catseis(sectionvxdata, fulldata_vx, recswitch, ntr_glob, MPI_COMM_NTR);
                                        }
                                        if ((ADJOINT_TYPE==1)|| (ADJOINT_TYPE==2)){
                                        catseis(sectionvydata, fulldata_vy, recswitch, ntr_glob, MPI_COMM_NTR);
                                        }
                                        if (ADJOINT_TYPE==4){
                                            catseis(sectionpdata, fulldata_p, recswitch, ntr_glob, MPI_COMM_NTR);
                                    }
                                    }
                                    if(WAVETYPE==2 || WAVETYPE==3){
                                        catseis(sectionvzdata, fulldata_vz, recswitch, ntr_glob, MPI_COMM_NTR);
                                    }
                                    if(myid_ntr==0){
                                        saveseis_glob(FP,fulldata_vx,fulldata_vy,fulldata_vz,fulldata_p,sectionpdiff,sectionpdiff,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter,3);
                                    }
                                }
                                
                                
                                /* Write synthetic filtered seismogramms to disk */
                                if (SEISMO && TIME_FILT && WRITE_FILTERED_DATA==2){
                                    if(WAVETYPE==1 || WAVETYPE==3){
                                        if ((ADJOINT_TYPE==1)|| (ADJOINT_TYPE==3)){
                                            catseis(sectionvx, fulldata_vx, recswitch, ntr_glob, MPI_COMM_NTR);
                                        }
                                        if ((ADJOINT_TYPE==1)|| (ADJOINT_TYPE==2)){
                                            catseis(sectionvy, fulldata_vy, recswitch, ntr_glob, MPI_COMM_NTR);
                                        }
                                        if (ADJOINT_TYPE==4){
                                            catseis(sectionp, fulldata_p, recswitch, ntr_glob, MPI_COMM_NTR);
                                        }
                                    }
                                    if(WAVETYPE==2 || WAVETYPE==3){
                                        catseis(sectionvz, fulldata_vz, recswitch, ntr_glob, MPI_COMM_NTR);
                                    }
                                    if(myid_ntr==0){
                                        saveseis_glob(FP,fulldata_vx,fulldata_vy,fulldata_vz,fulldata_p,sectionpdiff,sectionpdiff,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter,4);
                                    }
                                }
                                
                                
                            }
                            
                            
                            /*----------------------------------------------------------------------------------*/
                            /* ------ start loop over shots at receiver positions (backward model) ------------ */
                            /*----------------------------------------------------------------------------------*/
                            nshots1=1;
                            for (irec=1;irec<=nshots1;irec+=RECINC){ /* loop over shots at receiver positions */
                                
                                hin=1;
                                hin1=1;
                                
                                
                                if(MYID==0){
                                    printf("\n==================================================================================\n");
                                    printf("\n MYID=%d *****  Starting simulation (backward model) for shot %d of %d  ********** \n",MYID,irec,nshots1);
                                    printf("\n==================================================================================\n");
                                }
                                
                                
                                /*------------------------------------------------------------*/
                                /* determine source position out of reciever position on grid */
                                /*------------------------------------------------------------*/
                                
                                /* Distribute multiple source positions on subdomains */
                                /* define source positions at the receivers */
                                srcpos_loc_back = matrix(1,6,1,ntr);
                                for (i=1;i<=ntr;i++){
                                    srcpos_loc_back[1][i] = (recpos_loc[1][i]);
                                    srcpos_loc_back[2][i] = (recpos_loc[2][i]);
                                }
                                ntr1=ntr;
                                
                                /*----------------------------------*/
                                /* initialize wavefield with zero   */
                                /*----------------------------------*/
                                if (L){
                                    if(!ACOUSTIC) {
                                        zero_fdveps_visc(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,pvz,psxx,psyy,psxy,psxz,psyz,ux,uy,uxy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_sxz_x,psi_vxx,psi_vyx,psi_vzx,psi_syy_y,psi_sxy_y,psi_syz_y,psi_vyy,psi_vxy,psi_vzy,psi_vxxs,pr,pp,pq,pt,po);
                                    } else {
                                        zero_fdveps_viscac(-nd+1, NY+nd, -nd+1, NX+nd, pvx, pvy, psp, pvxp1, pvyp1, psi_sxx_x, psi_sxy_x, psi_vxx, psi_vyx, psi_syy_y, psi_sxy_y, psi_vyy, psi_vxy, psi_vxxs, pp);
                                    }
                                }else{
                                    if(!ACOUSTIC) {
                                        zero_fdveps(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,pvz,psxx,psyy,psxy,psxz,psyz,ux,uy,uxy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_sxz_x,psi_vxx,psi_vyx,psi_vzx,psi_syy_y,psi_sxy_y,psi_syz_y,psi_vyy,psi_vxy,psi_vzy,psi_vxxs);
                                    }else{
                                        zero_fdveps_ac(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,psp,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs);
                                    }
                                }
                                
                                lsnap=iround(TSNAP1/DT);
                                lsamp=NDT;
                                nsnap=0;
                                
                                /*--------------------------------------------------------------------------------*/
                                /*---------------------- Start loop over timesteps (backpropagation) -------------*/
                                /*--------------------------------------------------------------------------------*/
                                
                                if((!VERBOSE)&&(MYID==0)) fprintf(FP,"\n ****************************************\n ");
                                
                                for (nt=1;nt<=NT;nt++){
                                    
                                    // Ratio to give output to stout
                                    infoout = !(nt%nt_out);
                                    
                                    if((!VERBOSE)&&(MYID==0)) if(!(nt%(NT/40))) fprintf(FP,"*");
                                    
                                    /* Check if simulation is still stable P and SV */
                                    if (WAVETYPE==1 || WAVETYPE==3) {
                                        if (isnan(pvy[NY/2][NX/2])) {
                                            fprintf(FP,"\n Time step: %d; pvy: %f \n",nt,pvy[NY/2][NX/2]);
                                            declare_error(" Simulation is unstable !");
                                        }
                                    }
                                    /* Check if simulation is still stable SH */
                                    if (WAVETYPE==2 || WAVETYPE==3) {
                                        if (isnan(pvz[NY/2][NX/2])) {
                                            fprintf(FP,"\n Time step: %d; pvy: %f \n",nt,pvy[NY/2][NX/2]);
                                            declare_error(" Simulation is unstable !");
                                        }
                                    }
                                    
                                    if (MYID==0){
                                        if (infoout)  fprintf(FP,"\n Computing timestep %d of %d \n",nt,NT);
                                        time3=MPI_Wtime();
                                    }
                                    
                                    /* update of particle velocities */
                                    if(!ACOUSTIC) {
                                        if (WAVETYPE==1 || WAVETYPE==3) {
                                            update_v_PML(1, NX, 1, NY, nt, pvx, pvxp1, pvxm1, pvy, pvyp1, pvym1, uttx, utty, psxx, psyy, psxy, prip, prjp, srcpos_loc_back,sectionvxdiff,sectionvydiff,ntr1,absorb_coeff,hc,infoout,1, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_sxx_x, psi_syy_y, psi_sxy_y, psi_sxy_x);
                                        }
                                        
                                        if (WAVETYPE==2 || WAVETYPE==3) {
                                            update_v_PML_SH(1, NX, 1, NY, nt, pvz, pvzp1, pvzm1, psxz, psyz,prjp, srcpos_loc_back, sectionvzdiff, sectionvzdiff, ntr1, absorb_coeff,hc,infoout,1, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_sxz_x, psi_syz_y);
                                        }
                                    } else {
                                        update_v_acoustic_PML(1, NX, 1, NY, nt, pvx, pvxp1, pvxm1, pvy, pvyp1, pvym1, psp, prip, prjp, srcpos_loc_back, sectionvxdiff,sectionvydiff,ntr1,absorb_coeff,hc,infoout,1, K_x_half, a_x_half, b_x_half, K_y_half, a_y_half, b_y_half, psi_sxx_x, psi_syy_y);
                                    }
                                    
                                    
                                    if (MYID==0){
                                        time4=MPI_Wtime();
                                        time_av_v_update+=(time4-time3);
                                        if (infoout)  fprintf(FP," particle velocity exchange between PEs ...");
                                    }
                                    
                                    /* exchange of particle velocities between PEs */
                                    exchange_v(pvx,pvy,pvz, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, req_send, req_rec,wavetype_start);
                                    
                                    if (MYID==0){
                                        time5=MPI_Wtime();
                                        time_av_v_exchange+=(time5-time4);
                                        if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time5-time4);
                                    }
                                    
                                    if (L) {
                                        /* viscoelastic */
                                        if (WAVETYPE==1 || WAVETYPE==3) {
                                            if(!ACOUSTIC) {
                                                update_s_visc_PML(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppi, pu, puipjp, prho, hc, infoout,pr, pp, pq, fipjp, f, g, bip, bjm, cip, cjm, d, e, dip,K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);
                                            }else{
                                                update_p_visc_PML(1, NX, 1, NY, pvx, pvy, psp, ppi, prho, hc, infoout, pp, g, bjm, cjm, e, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);
                                            }
                                        }
                                        if (WAVETYPE==2 || WAVETYPE==3) {
                                            update_s_visc_PML_SH(1, NX, 1, NY, pvz, psxz, psyz, pt, po, bip, bjm, cip, cjm, d, dip,fipjp, f, hc,infoout, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half,psi_vzx, psi_vzy);
                                        }
                                    } else{
                                        /* elastic */
                                        if(!ACOUSTIC){
                                            if (WAVETYPE==1 || WAVETYPE==3) {
                                                update_s_elastic_PML(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppi, pu, puipjp, absorb_coeff, prho, hc, infoout, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);
                                            }
                                            if (WAVETYPE==2 || WAVETYPE==3) {
                                                update_s_elastic_PML_SH(1, NX, 1, NY, pvz,psxz,psyz,uxz,uyz,hc,infoout, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half,psi_vzx, psi_vzy,puipjp,pu,prho);
                                            }
                                        } else {
                                            update_p_PML(1, NX, 1, NY, pvx, pvy, psp, u, ppi, absorb_coeff, prho, hc, infoout, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);
                                        }
                                    }
                                    
                                    /* explosive source */
                                    if ((ADJOINT_TYPE==4))
                                        psource(nt,psxx,psyy,psp,srcpos_loc_back,sectionpdiff,ntr1,1);
                                    
                                    /* Applying free surface condition */
                                    if ((FREE_SURF) && (POS[2]==0)){
                                        if (!ACOUSTIC){
                                            if (L){
                                                /* viscoelastic */
                                                surface_PML(1, pvx, pvy, psxx, psyy, psxy,psyz, pp, pq, ppi, pu, prho, ptaup, ptaus, etajm, peta, hc, K_x, a_x, b_x, psi_vxxs, ux, uy,uxy,uyz,psxz,uxz);
                                            }else{
                                                /* elastic */
                                                surface_elastic_PML(1, pvx, pvy, psxx, psyy, psxy,psyz, ppi, pu, prho, hc, K_x, a_x, b_x, psi_vxxs, ux, uy, uxy,uyz,psxz,uxz);
                                            }
                                        } else {
                                            /* viscoelastic and elastic ACOUSTIC */
                                            surface_acoustic_PML(1, psp);
                                        }
                                    }
                                    
                                    if (MYID==0){
                                        time6=MPI_Wtime();
                                        time_av_s_update+=(time6-time5);
                                        if (infoout)  fprintf(FP," stress exchange between PEs ...");
                                    }
                                    
                                    /* stress exchange between PEs */
                                    if(!ACOUSTIC)
                                        exchange_s(psxx,psyy,psxy,psxz,psyz, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, req_send, req_rec,wavetype_start);
                                    else
                                        exchange_p(psp,bufferlef_to_rig, bufferrig_to_lef,buffertop_to_bot, bufferbot_to_top,req_send, req_rec);
                                    
                                    
                                    if (MYID==0){
                                        time7=MPI_Wtime();
                                        time_av_s_exchange+=(time7-time6);
                                        if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time7-time6);
                                    }
                                    
                                    
                                    /*if(nt==hin1){*/
                                    
                                    /*-------------------------------------------------*/
                                    /* Calculate convolution for every DTINV time step */
                                    /*-------------------------------------------------*/
                                    if(DTINV_help[NT-nt+1]==1){
                                        
                                        imat=((nxnyi*(NTDTINV)) - hin*nxnyi)+1;
                                        
                                        if((FORWARD_ONLY==0)){
                                            for (j=1;j<=NY;j=j+IDYI){
                                                for (i=1;i<=NX;i=i+IDXI){
                                                    
                                                    if (WAVETYPE==1 || WAVETYPE==3) {
                                                        waveconv_rho_shot[j][i]+=(pvxp1[j][i]*forward_prop_rho_x[j][i][NTDTINV-hin+1])+(pvyp1[j][i]*forward_prop_rho_y[j][i][NTDTINV-hin+1]);
                                                        
                                                        if(!ACOUSTIC){
                                                            waveconv_shot[j][i]+= (forward_prop_x[j][i][NTDTINV-hin+1]+forward_prop_y[j][i][NTDTINV-hin+1])*(psxx[j][i]+psyy[j][i]);
                                                        }else{
                                                            waveconv_shot[j][i]+= (forward_prop_p[j][i][NTDTINV-hin+1])*(psp[j][i]);
                                                        }
                                                        
                                                        if(!ACOUSTIC){
                                                            muss = prho[j][i] * pu[j][i] * pu[j][i];
                                                            lamss = prho[j][i] * ppi[j][i] * ppi[j][i] - 2.0 * muss;
                                                            
                                                            if(pu[j][i]>0.0){
                                                                waveconv_u_shot[j][i]+= ((1.0/(muss*muss))*(forward_prop_u[j][i][NTDTINV-hin+1] * psxy[j][i]))
                                                                + ((1.0/4.0) * ((forward_prop_x[j][i][NTDTINV-hin+1] + forward_prop_y[j][i][NTDTINV-hin+1]) * (psxx[j][i] + psyy[j][i])) / ((lamss+muss)*(lamss+muss)))
                                                                + ((1.0/4.0) * ((forward_prop_x[j][i][NTDTINV-hin+1] - forward_prop_y[j][i][NTDTINV-hin+1]) * (psxx[j][i] - psyy[j][i])) / (muss*muss));
                                                            }
                                                        }
                                                    }
                                                    
                                                    if (WAVETYPE==2 || WAVETYPE==3) {
                                                        waveconv_rho_shot_z[j][i]+=(pvzp1[j][i]*forward_prop_rho_z[j][i][NTDTINV-hin+1]);
                                                        muss = prho[j][i] * pu[j][i] * pu[j][i];
                                                        waveconv_u_shot_z[j][i]+=(1.0/(muss*muss))*(forward_prop_z_xz[j][i][NTDTINV-hin+1]*psxz[j][i]+forward_prop_z_yz[j][i][NTDTINV-hin+1]*psyz[j][i]);
                                                    }
                                                }
                                            }
                                        }

                                        hin++;
                                    }
                                    
                                    if((EPRECOND==1)&&(EPRECOND_ITER==iter||(EPRECOND_ITER==0))){
                                        if (WAVETYPE==1 || WAVETYPE==3) eprecond(Wr,pvx,pvy);
                                        if (WAVETYPE==2 || WAVETYPE==3) eprecond_SH(Wr_SH,pvz);
                                    }
                                    
                                    /* WRITE SNAPSHOTS TO DISK */
                                    if ((SNAP) && (nt==lsnap) && (nt<=TSNAP2/DT)){
                                        snap(FP,nt,++nsnap,pvx,pvy,psxx,psyy,psp,pu,ppi,hc,ishot);
                                        lsnap=lsnap+iround(TSNAPINC/DT);
                                    }
                                    
                                    if (MYID==0){
                                        time8=MPI_Wtime();
                                        time_av_timestep+=(time8-time3);
                                        if (infoout)  fprintf(FP," total real time for timestep %d : %4.2f s.\n",nt,time8-time3);
                                    }
                                    
                                }
                                
                                if((!VERBOSE)&&(MYID==0)) fprintf(FP,"\n");
                                /*--------------------------------------------------------------------------------*/
                                /*---------------------- End loop over timesteps (backpropagation) -------------*/
                                /*--------------------------------------------------------------------------------*/
                                
                                if(VERBOSE){
                                    /* Output jacobian VS per SHOT  PSV */
                                    if (WAVETYPE==1 || WAVETYPE==3) {
                                        if(!ACOUSTIC){
                                            sprintf(jac,"%s_jacobian_u_shot%i",JACOBIAN,ishot);
                                            write_matrix_disk(waveconv_u_shot, jac);
                                        } else {
                                            sprintf(jac,"%s_jacobian_shot%i",JACOBIAN,ishot);
                                            write_matrix_disk(waveconv_shot, jac);
                                        }
                                    }
                                    
                                    /* Output jacobian VS per SHOT SH */
                                    if (WAVETYPE==2 || WAVETYPE==3) {
                                        sprintf(jac,"%s_jacobian_u_SH_shot%i",JACOBIAN,ishot);
                                        write_matrix_disk(waveconv_u_shot_z, jac);
                                    }
                                }
                                
                            }
                            
                            /*--------------------------------------------------------------------------------*/
                            /* ------ end loop over shots at receiver positions (backward model) ------------ */
                            /*--------------------------------------------------------------------------------*/

                            
                            /* ------------------------------- */
                            /* calculate gradient direction pi */
                            /* ------------------------------- */
                            if((FORWARD_ONLY==0)&&(WAVETYPE==1||WAVETYPE==3)){
                                
                                
                                /* interpolate unknown values */
                                if((IDXI>1)||(IDYI>1)){
                                    interpol(IDXI,IDYI,waveconv_shot,1);
                                }
                                
                                /* calculate complete gradient */
                                for (j=1;j<=NY;j=j+IDY){
                                    for (i=1;i<=NX;i=i+IDX){
                                        
                                        waveconv_lam[j][i] = - DT * waveconv_shot[j][i];
                                        
                                        if(PARAMETERIZATION==1){
                                            if(!ACOUSTIC)
                                                muss = prho[j][i] * pu[j][i] * pu[j][i];
                                            else
                                                muss = 0;
                                            
                                            lamss = prho[j][i] * ppi[j][i] * ppi[j][i] - 2.0 *  muss;
                                            if(!ACOUSTIC)
                                                waveconv_lam[j][i] = (1.0/(4.0 * (lamss+muss) * (lamss+muss))) * waveconv_lam[j][i];
                                            else
                                                waveconv_lam[j][i] = (1.0/((lamss+muss) * (lamss+muss))) * waveconv_lam[j][i];
                                            
                                            /* calculate Vp gradient */
                                            waveconv_shot[j][i] = 2.0 * ppi[j][i] * prho[j][i] * waveconv_lam[j][i];
                                        }
                                        
                                        if(PARAMETERIZATION==2){
                                            /* calculate Zp gradient */
                                            waveconv_shot[j][i] = 2.0 * ppi[j][i] * waveconv_lam[j][i];}
                                        
                                        if(PARAMETERIZATION==3){
                                            waveconv_shot[j][i] = waveconv_lam[j][i];}
                                    }
                                }
                            }
                            
                            /* ---------------------------------- */
                            /* calculate gradient direction u PSV */
                            /* ---------------------------------- */
                            if((WAVETYPE==1 || WAVETYPE==3)&&(FORWARD_ONLY==0)&&(!ACOUSTIC)){
                                
                                
                                /* interpolate unknown values */
                                if((IDXI>1)||(IDYI>1)){
                                    interpol(IDXI,IDYI,waveconv_u_shot,1);
                                }
                                
                                /* calculate complete gradient */
                                for (j=1;j<=NY;j=j+IDY){
                                    for (i=1;i<=NX;i=i+IDX){
                                        
                                        /* calculate mu gradient */
                                        waveconv_mu[j][i] = - DT * waveconv_u_shot[j][i];
                                        
                                        if(PARAMETERIZATION==1){
                                            /* calculate Vs gradient */
                                            waveconv_u_shot[j][i] = (- 4.0 * prho[j][i] * pu[j][i] * waveconv_lam[j][i]) + 2.0 * prho[j][i] * pu[j][i] * waveconv_mu[j][i];
                                        }
                                        
                                        if(PARAMETERIZATION==2){
                                            /* calculate Zs gradient */
                                            waveconv_u_shot[j][i] = (- 4.0 * pu[j][i] * waveconv_lam[j][i]) + (2.0 * pu[j][i] * waveconv_mu[j][i]);}
                                        
                                        if(PARAMETERIZATION==3){
                                            /* calculate u gradient */
                                            waveconv_u_shot[j][i] = waveconv_mu[j][i];
                                        }
                                    }
                                }
                            }
                            
                            /* ---------------------------------- */
                            /* calculate gradient direction u SH */
                            /* ---------------------------------- */
                            if((WAVETYPE==2 || WAVETYPE==3)&&(FORWARD_ONLY==0)&&(!ACOUSTIC)){
                                
                                /* interpolate unknown values */
                                if((IDXI>1)||(IDYI>1)){
                                    interpol(IDXI,IDYI,waveconv_u_shot_z,1);
                                }
                                
                                /* calculate complete gradient */
                                for (j=1;j<=NY;j=j+IDY){
                                    for (i=1;i<=NX;i=i+IDX){
                                        
                                        /* calculate mu gradient */
                                        waveconv_mu_z[j][i] = - DT * waveconv_u_shot_z[j][i];
                                        
                                        if(PARAMETERIZATION==1){
                                            /* calculate Vs gradient */
                                            waveconv_u_shot_z[j][i] = (2.0 * prho[j][i] * pu[j][i] * waveconv_mu_z[j][i]);
                                        }
                                        
                                        if(PARAMETERIZATION==3){
                                            /* calculate u gradient */
                                            waveconv_u_shot_z[j][i] = waveconv_mu_z[j][i];
                                        }
                                    }
                                }
                            }
                            
                            /* ------------------------------------ */
                            /* calculate gradient direction rho PSV */
                            /* ------------------------------------ */
                            if((WAVETYPE==1 || WAVETYPE==3)&&(FORWARD_ONLY==0)){
                                
                                
                                /* interpolate unknown values */
                                if((IDXI>1)||(IDYI>1)){
                                    interpol(IDXI,IDYI,waveconv_rho,1);
                                }
                                
                                /* calculate complete gradient */
                                for (j=1;j<=NY;j=j+IDY){
                                    for (i=1;i<=NX;i=i+IDX){
                                        
                                        /* calculate density gradient rho' */
                                        waveconv_rho_s[j][i]= - DT * waveconv_rho_shot[j][i];
                                        
                                        if(PARAMETERIZATION==1){
                                            /* calculate density gradient */
                                            if(!ACOUSTIC) {
                                                waveconv_rho_shot[j][i] = ((((ppi[j][i] * ppi[j][i])-(2.0 * pu[j][i] * pu[j][i])) * waveconv_lam[j][i]) + (pu[j][i] * pu[j][i] * waveconv_mu[j][i]) + waveconv_rho_s[j][i]);
                                            } else {
                                                waveconv_rho_shot[j][i] = (((ppi[j][i] * ppi[j][i]) * waveconv_lam[j][i]) + waveconv_rho_s[j][i]);
                                            }
                                        }
                                        
                                        if(PARAMETERIZATION==3){
                                            /* calculate density gradient */
                                            waveconv_rho_shot[j][i] = waveconv_rho_s[j][i];
                                        }
                                        
                                    }
                                }
                            }
                            
                            /* ------------------------------------ */
                            /* calculate gradient direction rho SH*/
                            /* ------------------------------------ */
                            if((WAVETYPE==2 || WAVETYPE==3)&&(FORWARD_ONLY==0)){
                                
                                
                                /* interpolate unknown values */
                                if((IDXI>1)||(IDYI>1)){
                                    interpol(IDXI,IDYI,waveconv_rho_shot_z,1);
                                }
                                
                                /* calculate complete gradient */
                                for (j=1;j<=NY;j=j+IDY){
                                    for (i=1;i<=NX;i=i+IDX){
                                        
                                        /* calculate density gradient rho' */
                                        waveconv_rho_s_z[j][i]= - DT * waveconv_rho_shot_z[j][i];
                                        
                                        if(PARAMETERIZATION==1){
                                            /* calculate density gradient */
                                            waveconv_rho_shot_z[j][i] = ( (pu[j][i] * pu[j][i] * waveconv_mu_z[j][i]) + waveconv_rho_s_z[j][i]);
                                            
                                        }
                                        
                                        if(PARAMETERIZATION==3){
                                            /* calculate density gradient */
                                            waveconv_rho_shot_z[j][i] = waveconv_rho_s_z[j][i];
                                        }
                                        
                                    }
                                }
                            }
                            
                            
                            /* -------------------------------------------- */
                            /* calculate and apply energy preconditioning   */
                            /* -------------------------------------------- */
                            if((EPRECOND==1)||(EPRECOND==3)){
                                
                                /* calculate energy weights */
                                if(EPRECOND_ITER==iter||(EPRECOND_ITER==0)) {
                                    fprintf(FP,"\n Calculating approx. Hessian for shot %i. EPRECOND=%i, EPSILON_WE=%f",ishot,EPRECOND,EPSILON_WE);
                                    
                                    if(WAVETYPE==1 || WAVETYPE==3) {
                                        eprecond1(We,Ws,Wr,EPSILON_WE);
                                        if(EPRECOND_PER_SHOT) We_max=global_maximum(We);
                                    }
                                    
                                    if(WAVETYPE==2 || WAVETYPE==3) {
                                        eprecond1(We_SH,Ws_SH,Wr_SH,EPSILON_WE_SH);
                                        if(EPRECOND_PER_SHOT_SH) We_max_SH=global_maximum(We_SH);
                                    }
                                    
                                    if(EPRECOND_PER_SHOT && (WAVETYPE==1 || WAVETYPE==3)){
                                        fprintf(FP,"\n Applying approx. Hessian for shot %i PSV. EPRECOND=%i, EPSILON_WE=%f",ishot,EPRECOND,EPSILON_WE);
                                        for (j=1;j<=NY;j=j+IDY){
                                            for (i=1;i<=NX;i=i+IDX){
                                                We[j][i]=We[j][i]/We_max;
                                                waveconv_shot[j][i] = waveconv_shot[j][i]/(We[j][i]);
                                                if(!ACOUSTIC){
                                                    waveconv_u_shot[j][i] = waveconv_u_shot[j][i]/(We[j][i]);
                                                }
                                                waveconv_rho_shot[j][i] = waveconv_rho_shot[j][i]/(We[j][i]);
                                            }
                                        }
                                    }
                                    
                                    if(EPRECOND_PER_SHOT_SH && (WAVETYPE==2 || WAVETYPE==3)){
                                        fprintf(FP,"\n Applying approx. Hessian for shot %i SH. EPRECOND=%i, EPSILON_WE=%f",ishot,EPRECOND,EPSILON_WE);
                                        for (j=1;j<=NY;j=j+IDY){
                                            for (i=1;i<=NX;i=i+IDX){
                                                We_SH[j][i]=We_SH[j][i]/We_max_SH;
                                                waveconv_u_shot_z[j][i] = waveconv_u_shot_z[j][i]/(We_SH[j][i]);
                                                waveconv_rho_shot_z[j][i] = waveconv_rho_shot_z[j][i]/(We_SH[j][i]);
                                            }
                                        }
                                    }
                                    
                                    for (j=1;j<=NY;j=j+IDY){
                                        for (i=1;i<=NX;i=i+IDX){
                                            if(WAVETYPE==1 || WAVETYPE==3) We_sum[j][i]+=1.0/We[j][i];
                                            if(WAVETYPE==2 || WAVETYPE==3) We_sum_SH[j][i]+=1.0/We_SH[j][i];
                                        }
                                    }
                                }
                            }
                            
                            /* --------------- */
                            /* Apply taper SH  */
                            /* --------------- */
                            if (WAVETYPE==2 || WAVETYPE==3){
                                /* applying a circular taper at the source position to the gradient of each shot */
                                if (SWS_TAPER_CIRCULAR_PER_SHOT){
                                    
                                    /* applying the preconditioning */
                                    taper_grad_shot(waveconv_u_shot_z,taper_coeff,srcpos,nsrc,recpos,ntr_glob,ishot,1);
                                    taper_grad_shot(waveconv_rho_shot_z,taper_coeff,srcpos,nsrc,recpos,ntr_glob,ishot,1);
                                    
                                }
                                /* end of SWS_TAPER_CIRCULAR_PER_SHOT == 1 */
                                
                                /* applying taper file which is read in */
                                if (SWS_TAPER_FILE_PER_SHOT){
                                    
                                    /* applying the preconditioning */
                                    taper_grad_shot(waveconv_u_shot_z,taper_coeff,srcpos,nsrc,recpos,ntr_glob,ishot,3); /*taper vs gradient */
                                    taper_grad_shot(waveconv_rho_shot_z,taper_coeff,srcpos,nsrc,recpos,ntr_glob,ishot,4); /*taper rho gradient */
                                    
                                }
                            }
                            
                            /* --------------- */
                            /* Apply taper PSV */
                            /* --------------- */
                            if (WAVETYPE==1 || WAVETYPE==3) {
                                /* applying a circular taper at the source position to the gradient of each shot */
                                if (SWS_TAPER_CIRCULAR_PER_SHOT){
                                    
                                    /* applying the preconditioning */
                                    taper_grad_shot(waveconv_shot,taper_coeff,srcpos,nsrc,recpos,ntr_glob,ishot,1);
                                    if(!ACOUSTIC){
                                        taper_grad_shot(waveconv_u_shot,taper_coeff,srcpos,nsrc,recpos,ntr_glob,ishot,1);
                                    }
                                    taper_grad_shot(waveconv_rho_shot,taper_coeff,srcpos,nsrc,recpos,ntr_glob,ishot,1);
                                    
                                }
                                /* end of SWS_TAPER_CIRCULAR_PER_SHOT == 1 */
                                
                                /* applying taper file which is read in */
                                if (SWS_TAPER_FILE_PER_SHOT){
                                    
                                    /* applying the preconditioning */
                                    taper_grad_shot(waveconv_shot,taper_coeff,srcpos,nsrc,recpos,ntr_glob,ishot,2);	/* taper vp gradient */
                                    if(!ACOUSTIC){
                                        taper_grad_shot(waveconv_u_shot,taper_coeff,srcpos,nsrc,recpos,ntr_glob,ishot,3); /*taper vs gradient */
                                    }
                                    taper_grad_shot(waveconv_rho_shot,taper_coeff,srcpos,nsrc,recpos,ntr_glob,ishot,4); /*taper rho gradient */
                                    
                                }
                            }
                            
                            /* ----------------------------------------- */
                            /* Summing up the gradient for all shots PSV */
                            /* ----------------------------------------- */
                            if (WAVETYPE==1 || WAVETYPE==3) {
                                for (j=1;j<=NY;j=j+IDY){
                                    for (i=1;i<=NX;i=i+IDX){
                                        waveconv[j][i] += waveconv_shot[j][i];
                                        if(!ACOUSTIC){
                                            waveconv_u[j][i] += waveconv_u_shot[j][i];
                                        }
                                        waveconv_rho[j][i] += waveconv_rho_shot[j][i];
                                    }
                                }
                            }
                            
                            /* ----------------------------------------- */
                            /* Summing up the gradient for all shots SH */
                            /* ----------------------------------------- */
                            if (WAVETYPE==2 || WAVETYPE==3) {
                                for (j=1;j<=NY;j=j+IDY){
                                    for (i=1;i<=NX;i=i+IDX){
                                        waveconv_u_z[j][i] += waveconv_u_shot_z[j][i];
                                        waveconv_rho_z[j][i] += waveconv_rho_shot_z[j][i];
                                    }
                                }
                            }
                            
                        }
                        
                    }
                    
                    /*------------------------------------------------------------------------------*/
                    /*----------- End of loop over shots -----------------------------------------*/
                    /*------------------------------------------------------------------------------*/
                    
                    nsrc_loc=0;
                    
                }
                
                if(FORWARD_ONLY==0){
                    /* ----------------------------------------- */
                    /*     Applying and output approx hessian    */
                    /* ----------------------------------------- */
                    if((EPRECOND==1)||(EPRECOND==3)){
                        
                        if(!EPRECOND_PER_SHOT && (WAVETYPE==1 || WAVETYPE==3)){
                            fprintf(FP,"\n Applying approx. Hessian to summed gradient PSV. EPRECOND=%i, EPSILON_WE=%f",EPRECOND,EPSILON_WE);
                            
                            We_sum_max1=global_maximum(We_sum);
                            for (j=1;j<=NY;j=j+IDY){
                                for (i=1;i<=NX;i=i+IDX){
                                    We_sum[j][i]=We_sum[j][i]/We_sum_max1;
                                    waveconv[j][i] = waveconv[j][i]*We_sum[j][i];
                                    if(!ACOUSTIC){
                                        waveconv_u[j][i] = waveconv_u[j][i]*We_sum[j][i];
                                    }
                                    waveconv_rho[j][i] = waveconv_rho[j][i]*We_sum[j][i];
                                }
                            }
                        }
                        
                        if(!EPRECOND_PER_SHOT_SH && (WAVETYPE==2 || WAVETYPE==3)){
                            fprintf(FP,"\n Applying approx. Hessian to summed gradient SH. EPRECOND=%i, EPSILON_WE=%f",EPRECOND,EPSILON_WE);
                            
                            We_sum_max1=global_maximum(We_sum_SH);
                            for (j=1;j<=NY;j=j+IDY){
                                for (i=1;i<=NX;i=i+IDX){
                                    We_sum_SH[j][i]=We_sum_SH[j][i]/We_sum_max1;
                                    waveconv_u_z[j][i] = waveconv_u_z[j][i]*We_sum_SH[j][i];
                                    waveconv_rho_z[j][i] = waveconv_rho_z[j][i]*We_sum_SH[j][i];
                                }
                            }
                        }
                        
                        if (WAVETYPE==1 || WAVETYPE==3) {
                            sprintf(jac,"%s_approx_hessian_it%i",JACOBIAN,iter);
                            write_matrix_disk(We_sum, jac);
                        }
                        
                        if (WAVETYPE==2 || WAVETYPE==3) {
                            sprintf(jac,"%s_approx_hessian_SH_it%i",JACOBIAN,iter);
                            write_matrix_disk(We_sum_SH, jac);
                        }
                    }
                    
                    /* ----------------------------------------- */
                    /*  Set gradient to zero if no inversion     */
                    /* ----------------------------------------- */
                    if (WAVETYPE==1 || WAVETYPE==3) {
                        for (j=1;j<=NY;j=j+IDY){
                            for (i=1;i<=NX;i=i+IDX){
                                if(iter<INV_VP_ITER) waveconv[j][i] = 0.0;
                                if(iter<INV_VS_ITER && !ACOUSTIC) waveconv_u[j][i] = 0.0;
                                if(iter<INV_RHO_ITER) waveconv_rho[j][i] = 0.0;
                            }
                        }
                    }
                    if (WAVETYPE==2 || WAVETYPE==3) {
                        for (j=1;j<=NY;j=j+IDY){
                            for (i=1;i<=NX;i=i+IDX){
                                if(iter<INV_VS_ITER) waveconv_u_z[j][i] = 0.0;
                                if(iter<INV_RHO_ITER) waveconv_rho_z[j][i] = 0.0;
                            }
                        }
                    }
                    
                    /* ------------------------------*/
                    /* calculate L2 norm of all CPUs */
                    /* ------------------------------*/
                    if(WAVETYPE==1||WAVETYPE==3){
                        L2sum = 0.0;
                        MPI_Allreduce(&L2,&L2sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
                        energy_sum = 0.0;
                        MPI_Allreduce(&energy,&energy_sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
                        L2sum_all_shots = 0.0;
                        MPI_Allreduce(&L2_all_shots,&L2sum_all_shots,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
                        energy_sum_all_shots = 0.0;
                        MPI_Allreduce(&energy_all_shots,&energy_sum_all_shots,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
                        
                        if(MYID==0&&(WAVETYPE==3)) printf("\n\n PSV: L2=%f",L2sum_all_shots/energy_sum_all_shots);
                    }
                    if(WAVETYPE==2||WAVETYPE==3){
                        L2sum_SH = 0.0;
                        MPI_Allreduce(&L2_SH,&L2sum_SH,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
                        energy_sum_SH = 0.0;
                        MPI_Allreduce(&energy_SH,&energy_sum_SH,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
                        L2sum_all_shots_SH = 0.0;
                        MPI_Allreduce(&L2_all_shots_SH,&L2sum_all_shots_SH,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
                        energy_sum_all_shots_SH = 0.0;
                        MPI_Allreduce(&energy_all_shots_SH,&energy_sum_all_shots_SH,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
                        
                        if(MYID==0&&(WAVETYPE==3)) printf("\n  SH: L2=%f",L2sum_all_shots_SH/energy_sum_all_shots_SH);
                    }
                    
                    sum_killed_traces=0;
                    MPI_Allreduce(&killed_traces,&sum_killed_traces,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
                    sum_killed_traces_testshots=0;
                    MPI_Allreduce(&killed_traces_testshots,&sum_killed_traces_testshots,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
                    
                    
                    switch (LNORM){
                            case 2:
                            L2t[1]=0.0; L2t[4]=0.0;
    
                            
                            if(JOINT_EQUAL_WEIGHTING){
                                if(JOINT_EQUAL_new_max){
                                    JOINT_EQUAL_PSV=L2sum/energy_sum;
                                    JOINT_EQUAL_SH=L2sum_SH/energy_sum_SH;
                                    
                                    JOINT_EQUAL_PSV_all=L2sum_all_shots/energy_sum_all_shots;
                                    JOINT_EQUAL_SH_all=L2sum_all_shots_SH/energy_sum_all_shots_SH;
                                    
                                    JOINT_EQUAL_new_max=0;
                                }
                            
                                L2t[1]+=(L2sum/energy_sum)/JOINT_EQUAL_PSV;
                                L2t[4]+=(L2sum_all_shots/energy_sum_all_shots)/JOINT_EQUAL_PSV_all;
                                
                                L2t[1]+=(L2sum_SH/energy_sum_SH)/JOINT_EQUAL_SH;
                                L2t[4]+=(L2sum_all_shots_SH/energy_sum_all_shots_SH)/JOINT_EQUAL_SH_all;
                                
                                break;
                            }
                            
                            if(WAVETYPE==1||WAVETYPE==3){
                                L2t[1]+=L2sum/energy_sum;
                                L2t[4]+=L2sum_all_shots/energy_sum_all_shots;
                            }
                            
                            if(WAVETYPE==2||WAVETYPE==3){
                                L2t[1]+=L2sum_SH/energy_sum_SH;
                                L2t[4]+=L2sum_all_shots_SH/energy_sum_all_shots_SH;
                            }
                            
                        if(MYID==0&&(WAVETYPE==3)) printf("\n Sum: L2=%f",L2t[4]);
                            
                            break;
                            case 7:
                            if (TRKILL){
                                if(ADJOINT_TYPE==1){	/* x and y component are used in the inversion */
                                    L2t[1]=2.0*(1.0+(L2sum/((float)((NO_OF_TESTSHOTS*ntr_glob-sum_killed_traces_testshots)*2.0))));
                                    L2t[4]=2.0*(1.0+(L2sum_all_shots/((float)((nsrc_glob*ntr_glob-sum_killed_traces)*2.0))));
                                    if (MYID==0){
                                        printf("sum_killed_traces_testshots=%d\n",sum_killed_traces_testshots);
                                        printf("sum_killed_traces=%d\n",sum_killed_traces);}}
                                else{
                                    L2t[1]=2.0*(1.0+(L2sum/((float)(NO_OF_TESTSHOTS*ntr_glob-sum_killed_traces_testshots))));
                                    L2t[4]=2.0*(1.0+(L2sum_all_shots/((float)(nsrc_glob*ntr_glob-sum_killed_traces))));
                                    if (MYID==0){
                                        printf("sum_killed_traces_testshots=%d\n",sum_killed_traces_testshots);
                                        printf("sum_killed_traces=%d\n",sum_killed_traces);
                                        printf("ntr_glob=%d\n",ntr_glob);
                                        printf("nsrc_glob=%d\n",nsrc_glob);}}}
                            else{
                                if(ADJOINT_TYPE==1){	/* x and y component are used in the inversion */
                                    L2t[1]=2.0*(1.0+(L2sum/((float)NO_OF_TESTSHOTS*(float)ntr_glob*2.0)));
                                    L2t[4]=2.0*(1.0+(L2sum_all_shots/((float)nsrc_glob*(float)ntr_glob*2.0)));}
                                else{
                                    L2t[1]=2.0*(1.0+(L2sum/((float)NO_OF_TESTSHOTS*(float)ntr_glob)));
                                    L2t[4]=2.0*(1.0+(L2sum_all_shots/((float)nsrc_glob*(float)ntr_glob)));}
                            }
                            break;
                            case 8:
                            L2t[1]=L2sum/energy_sum;
                            L2t[4]=L2sum_all_shots/energy_sum_all_shots;
                            break;
                        default:
                            L2t[1]=L2sum;
                            L2t[4]=L2sum_all_shots;
                            break;
                    }
                    
                    if(!steplength_search) {
                        L2_SL_old=L2t[4];
                    } else {
                        L2_SL_new=L2t[4];
                    }
                    
                    if(MYID==0&&VERBOSE){
                        fprintf(FP,"\n\nL2sum: %f\n", L2sum);
                        fprintf(FP,"energy_sum: %e\n\n", energy_sum);
                        fprintf(FP,"L2sum_all_shots: %f\n", L2sum_all_shots);
                        fprintf(FP,"energy_sum_all_shots: %e\n\n", energy_sum_all_shots);
                    }
                }
                /* Count how often this loop runs */
                if(GRAD_METHOD==2 && iter>LBFGS_iter_start) {
                    wolfe_sum_FWI++;
                }
            }
            
            /*-----------------------------------------------------*/
            /*     Gradient optimization with PCG or L-BFGS        */
            /*-----------------------------------------------------*/
            if(gradient_optimization==1 && FORWARD_ONLY==0) {
                
                /* ----------------------------------------------------------------------*/
                /* ----------- Preconditioned Conjugate Gradient Method (PCG)  ----------*/
                /* ----------------------------------------------------------------------*/
                if(GRAD_METHOD==1){
                    
                    if( (iter-PCG_iter_start) < 2 ) {
                        fprintf(FP,"\n\n ----- Conjugate Gradient Method -----");
                        fprintf(FP,"\n Will not use second last gradient for conjugate gradient");
                        if( (iter-PCG_iter_start) < 1 ) {
                            fprintf(FP,"\n Will not use last gradient for conjugate gradient");
                        }
                    }
                    
                    if (WAVETYPE==1 || WAVETYPE==3) PCG(waveconv, taper_coeff, nsrc, srcpos, recpos, ntr_glob, iter, C_vp, gradp, nfstart_jac, waveconv_u, C_vs, gradp_u, waveconv_rho, C_rho, gradp_rho,Vs_avg,F_LOW_PASS,PCG_iter_start);
                    if (WAVETYPE==2 || WAVETYPE==3) PCG_SH(taper_coeff, nsrc, srcpos, recpos, ntr_glob, iter, nfstart_jac, waveconv_u_z, C_vs, gradp_u_z, waveconv_rho_z, C_rho, gradp_rho_z,Vs_avg,F_LOW_PASS,PCG_iter_start);
                    
                }
                
                /* ---------------------------------------------------------*/
                /* ----------- Beginn Joint Inversion PSV and SH  ----------*/
                /* ---------------------------------------------------------*/
                if(FORWARD_ONLY==0){
                    switch (WAVETYPE) {
                        case 2:
                            /* If only SH is inverted, set these gradients to the actual ones */
                            waveconv_u=waveconv_u_z;
                            waveconv_rho=waveconv_rho_z;
                            //waveconv=NULL; // SH case no VP inversion is possible
                            INV_VP_ITER=iter+10; // Just in case...
                            break;
                            
                        case 3:
                            /* If PSV and SH are simultaneously are inverted, a joint inversion have to be done */
                            fprintf(FP, "\n\n===================================\n");
                            fprintf(FP, " Joint Inversion Process initiated!\n");
                            fprintf(FP, "===================================\n\n");
                            waveconv_u=joint_inversion_grad(waveconv_u,waveconv_u_z,JOINT_INVERSION_PSV_SH_ALPHA_VS,JOINT_INVERSION_PSV_SH_TYPE);
                            waveconv_rho=joint_inversion_grad(waveconv_rho,waveconv_rho_z,JOINT_INVERSION_PSV_SH_ALPHA_RHO,JOINT_INVERSION_PSV_SH_TYPE);
                            
                            /* Output joint gradient to disk */
                            sprintf(jac,"%s_joint_vs_it%i",JACOBIAN,iter);
                            write_matrix_disk(waveconv_u, jac);
                            
                            /* Output joint gradient to disk */
                            sprintf(jac,"%s_joint_rho_it%i",JACOBIAN,iter);
                            write_matrix_disk(waveconv_rho, jac);
                            break;
                            
                        default:
                            /* By default (WAVETYPE=1) no joint inversion */
                            break;
                    }
                }
                /* ---------------------------------------------------------*/
                /* ----------- END Joint Inversion PSV and SH  -------------*/
                /* ---------------------------------------------------------*/
                
                
                /* -------------------------------------------------------------------------*/
                /* ----------- Limited Memory - Broyden-Fletcher-Goldfarb-Shanno  ----------*/
                /* -------------------------------------------------------------------------*/
                if(GRAD_METHOD==2){
                    
                    /*---------------------*/
                    /*         TAPER       */
                    /*---------------------*/
                    
                    if(WAVETYPE==1 || WAVETYPE==3){
                        if (SWS_TAPER_FILE){   /* read taper from BIN-File*/
                            taper_grad(waveconv,taper_coeff,srcpos,nsrc,recpos,ntr_glob,4);
                        }
                        if(GRAD_FILTER==1 && !ACOUSTIC){
                            smooth(waveconv,1,1,Vs_avg,F_LOW_PASS);
                        }else if(GRAD_FILTER==1 && ACOUSTIC){
                            smooth(waveconv,1,1,Vp_avg,F_LOW_PASS);
                        }
                    }
                    
                    if (SWS_TAPER_FILE && !ACOUSTIC){   /* read taper from BIN-File*/
                        taper_grad(waveconv_u,taper_coeff,srcpos,nsrc,recpos,ntr_glob,5);}
                    
                    if (SWS_TAPER_FILE){   /* read taper from BIN-File*/
                        taper_grad(waveconv_rho,taper_coeff,srcpos,nsrc,recpos,ntr_glob,6);}
                    
                    if(GRAD_FILTER==1 && !ACOUSTIC){smooth(waveconv_u,2,1,Vs_avg,F_LOW_PASS);}
                    if(GRAD_FILTER==1 && !ACOUSTIC){
                        smooth(waveconv_rho,3,1,Vs_avg,F_LOW_PASS);
                    }else if(GRAD_FILTER==1 && ACOUSTIC){
                        smooth(waveconv_rho,3,1,Vp_avg,F_LOW_PASS);
                    }
                    
                    if(WOLFE_CONDITION) {
                        for (j=1;j<=NY;j=j+IDY){
                            for (i=1;i<=NX;i=i+IDX){
                                
                                if(WAVETYPE!=2){
                                    waveconv_old[j][i]=waveconv[j][i];
                                    ppinp1[j][i] = ppi[j][i];
                                }
                                
                                if(!ACOUSTIC){
                                    waveconv_u_old[j][i] = waveconv_u[j][i];
                                    punp1[j][i] =pu[j][i];
                                }
                                
                                waveconv_rho_old[j][i] = waveconv_rho[j][i];
                                prhonp1[j][i] = prho[j][i];
                            }
                        }
                    }
                    
                    /*---------------------*/
                    /*       L-BFGS        */
                    /*---------------------*/
                    lbfgs(waveconv_u, waveconv_rho, waveconv,Vs_avg,rho_avg,Vp_avg, rho_LBFGS, s_LBFGS, y_LBFGS, N_LBFGS,LBFGS_NPAR, iter,&LBFGS_iter_start);
                    
                }
                
                if(WOLFE_CONDITION) {
                    for (j=1;j<=NY;j=j+IDY){
                        for (i=1;i<=NX;i=i+IDX){
                            if(WAVETYPE!=2) waveconv_up[j][i]=waveconv[j][i];
                            if(!ACOUSTIC) waveconv_u_up[j][i] = waveconv_u[j][i];
                            waveconv_rho_up[j][i] = waveconv_rho[j][i];
                        }
                    }
                }
                
                /* Set Values */
                FWI_run=0;
                steplength_search=1;
                gradient_optimization=0;
                countstep=0;
            }
            
            /*-----------------------------------------------------*/
            /*   Wolfe condition: Check and step length search     */
            /*-----------------------------------------------------*/
            if(steplength_search==1 && (GRAD_METHOD==2 && iter>LBFGS_iter_start) && FORWARD_ONLY==0 && WOLFE_CONDITION) {
                
                if(countstep>WOLFE_NUM_TEST) {
                    if(wolfe_found_lower_L2) {
                        fprintf(FP,"\n After %d simulations no step length could be found which reduces the misfit and satisfy the wolfe condition.",countstep);
                        fprintf(FP,"\n However, a step length which reduces the mifit could be found and will be used.");
                        fprintf(FP,"\n L2_new=%f and L2_old=%f with step length=%f.",L2_SL_FS,L2_SL_old,alpha_SL_FS);
                        alpha_SL=alpha_SL_FS;
                        L2_SL_new=L2_SL_FS;
                        use_wolfe_failsafe=1;
                        break;
                    } else {
                        fprintf(FP,"\n After %d simulations no step length could be found which reduces the misfit and satisfy the wolfe condition.",countstep);
                        fprintf(FP,"\n Will continue without model update.");
                        wolfe_SLS_failed=1;
                        break;
                    }
                }
                /* Initialisation of calculation of gradient and misfit with first step lenght */
                if(countstep==0) {
                    
                    /* Init test step lengths */
                    alpha_SL=1; alpha_SL_min=0; alpha_SL_max=0;
                    if(iter>2 && WOLFE_TRY_OLD_STEPLENGTH)  alpha_SL=alpha_SL_old;
                    
                    /* Calculate update */
                    calc_mat_change_test(waveconv_up,waveconv_rho_up,waveconv_u_up,prhonp1,prho,ppinp1,ppi,punp1,pu,iter,1,FORWARD_ONLY,alpha_SL,1,nfstart,Vs0,Vp0,Rho0,wavetype_start,s_LBFGS,N_LBFGS,LBFGS_NPAR,Vs_avg,Vp_avg,rho_avg,LBFGS_iter_start);
                    FWI_run=1;
                }
                
                /* Check if current step length satisfy wolfe condition, if not call linesearch for new step length */
                if(countstep>0) {
                    
                    wolfe_status=check_wolfe(alpha_SL, L2_SL_old, L2_SL_new, waveconv_u_old, waveconv_u, waveconv_u_up,waveconv_rho_old, waveconv_rho, waveconv_rho_up,waveconv_old, waveconv, waveconv_up, c1_SL, c2_SL,LBFGS_NPAR);
                    
                    if(wolfe_status==0) {
                        /* Current step length satisfy wolfe condition, abort step length search */
                        steplength_search=0;
                        break;
                    } else {
                        
                        fprintf(FP,"\n FD-Simulation %d out of max. %d",countstep,WOLFE_NUM_TEST);
                        fprintf(FP,"\n Old L2=%f; New L2=%f",L2_SL_old,L2_SL_new);
                        fprintf(FP, "\n Min. steplength=%1.3f; max. steplength=%1.3f;\n Used steplength=%1.3f",alpha_SL_min,alpha_SL_max,alpha_SL);
                        
                        /* Failsafe, if a step length is found which reduce the misfit, save this one */
                        if((countstep==1 && L2_SL_new<L2_SL_old)||(countstep>2 && L2_SL_new<L2_SL_FS)) {
                            wolfe_found_lower_L2=1;
                            alpha_SL_FS=alpha_SL;
                            L2_SL_FS=L2_SL_new;
                        }
                        
                        /* Current step length do not satisfy wolfe condition, try new step length */
                        wolfe_linesearch(wolfe_status, &alpha_SL_min, &alpha_SL_max, &alpha_SL);
                        calc_mat_change_test(waveconv_up,waveconv_rho_up,waveconv_u_up,prhonp1,prho,ppinp1,ppi,punp1,pu,iter,1,FORWARD_ONLY,alpha_SL,1,nfstart,Vs0,Vp0,Rho0,wavetype_start,s_LBFGS,N_LBFGS,LBFGS_NPAR,Vs_avg,Vp_avg,rho_avg,LBFGS_iter_start);
                        fprintf(FP,"; New steplength=%1.3f",alpha_SL);
                        FWI_run=1;
                    }
                    
                }
                
            } else {
                steplength_search=0;
            }
            
            countstep++;
            
            if(FORWARD_ONLY!=0) break;
        }
        
        if(wolfe_SLS_failed) {
            
            if (TIME_FILT==0){
                if(MYID==0) fprintf(FPL2,"%e \t %d \t %d \t %f \t 0 \t %d \t %e \t %e \t  %f\n",0.0,iter,wolfe_sum_FWI,0.0,countstep-1,L2_SL_old,L2_SL_old,GAMMA);}
            else{
                if(MYID==0) fprintf(FPL2,"%e \t %d \t %d \t %f \t 0 \t %d \t %e \t %e \t %f \t %f \n",0.0,iter,wolfe_sum_FWI,0.0,countstep-1,L2_SL_old,L2_SL_old,F_LOW_PASS,GAMMA);
            }
            
            if(WAVETYPE==3 && MYID==0){
                fprintf(FPL2_JOINT,"%d \t %f \t %f\n",iter,L2sum_all_shots/energy_sum_all_shots,L2sum_all_shots_SH/energy_sum_all_shots_SH);
            }
            
            /* No update is done here, however model fils are written to disk for easy post processing */
            alpha_SL=0.0;
            calc_mat_change_test(waveconv_up,waveconv_rho_up,waveconv_u_up,prhonp1,prho,ppinp1,ppi,punp1,pu,iter,1,FORWARD_ONLY,alpha_SL,0,nfstart,Vs0,Vp0,Rho0,wavetype_start,s_LBFGS,N_LBFGS,LBFGS_NPAR,Vs_avg,Vp_avg,rho_avg,LBFGS_iter_start);
            
            alpha_SL_old=1;
            
            /* If minimum number of iterations would be enforced, L-BFGS is likely to crash */
            min_iter_help=0;
        }
        
        /*-----------------------------------------------------*/
        /*       Wolfe condition: Model update                 */
        /*-----------------------------------------------------*/
        if((GRAD_METHOD==2 && iter>LBFGS_iter_start) && FORWARD_ONLY==0 && WOLFE_CONDITION && !wolfe_SLS_failed) {
            
            /* save old step length for next iteration as first guess */
            if(WOLFE_TRY_OLD_STEPLENGTH) alpha_SL_old=alpha_SL;
            
            /* Not sure if this is needed */
            for (j=1;j<=NY;j=j+IDY){
                for (i=1;i<=NX;i=i+IDX){
                    prho[j][i]=prhonp1[j][i];
                    if(WAVETYPE!=2)  ppi[j][i]=ppinp1[j][i];
                    if(!ACOUSTIC) pu[j][i]=punp1[j][i];
                }
            }
            
            /* do the final model update */
            calc_mat_change_test(waveconv_up,waveconv_rho_up,waveconv_u_up,prhonp1,prho,ppinp1,ppi,punp1,pu,iter,1,FORWARD_ONLY,alpha_SL,0,nfstart,Vs0,Vp0,Rho0,wavetype_start,s_LBFGS,N_LBFGS,LBFGS_NPAR,Vs_avg,Vp_avg,rho_avg,LBFGS_iter_start);
            
            L2_hist[iter]=L2t[4];
            
            /* write L2 log file */
            float diff=0.0;
            diff=fabs((L2_hist[iter-2]-L2_hist[iter])/L2_hist[iter-2]);
            if (TIME_FILT==0){
                if(MYID==0) fprintf(FPL2,"%e \t %d \t %d \t %f \t 0 \t %d \t %e \t %e \t %f \n",alpha_SL,iter,wolfe_sum_FWI,diff,countstep-1,L2_SL_old,L2_SL_new,GAMMA);}
            else{
                if(MYID==0) fprintf(FPL2,"%e \t %d \t %d \t %f \t 0 \t %d \t %e \t %e \t %f \t %f\n",alpha_SL,iter,wolfe_sum_FWI,diff,countstep-1,L2_SL_old,L2_SL_new,F_LOW_PASS,GAMMA);
            }
            
            if(WAVETYPE==3 && MYID==0){
                fprintf(FPL2_JOINT,"%d \t %f \t %f\n",iter,L2sum_all_shots/energy_sum_all_shots,L2sum_all_shots_SH/energy_sum_all_shots_SH);
            }
            
            /* initiate variables for next iteration */
            if(use_wolfe_failsafe==1) {
                L2_hist[iter]=L2_SL_new;
                FWI_run=1;
                gradient_optimization=1;
                steplength_search=0;
                
                wolfe_SLS_failed=1;
                alpha_SL_old=1;
            } else {
                FWI_run=0;
                gradient_optimization=1;
                L2_SL_old=L2_SL_new;
            }
            
            /* Reset */
            wolfe_found_lower_L2=0;
            use_wolfe_failsafe=0;
        } else {
            FWI_run=1;
            gradient_optimization=1;
        }
        
        
        opteps_vp=0.0;
        opteps_vs=0.0;
        opteps_rho=0.0;
        
        
        /* ============================================================================================================================*/
        /* =============================================== Step length estimation =====================================================*/
        /* ============================================================================================================================*/
        
        if((FORWARD_ONLY==0)  && (!WOLFE_CONDITION || (WOLFE_CONDITION && GRAD_METHOD==2 && iter==LBFGS_iter_start))){
            
            fprintf(FP,"\n=================================================================================================\n");
            fprintf(FP,"\n *********************** Starting step length estimation at iteration %i ************************\n",iter);
            fprintf(FP,"\n=================================================================================================\n\n");
            
            
            step1=0;
            step2=0;
            
            /* start with first guess for step length alpha */
            eps_scale=EPS_SCALE; /* maximum model change = 1% of the maximum model value */
            
            countstep=0; /* count number of forward calculations */
            
            itests=2;
            iteste=2;
            
            if(GRAD_METHOD==2&&(iter>LBFGS_iter_start)&&LBFGS_STEP_LENGTH) {
                itests=1;
                iteste=1;
                eps_scale=0.1;
            }
            
            /* set min_iter_help to initial global value of MIN_ITER */
            if(iter==1){min_iter_help=MIN_ITER;}
            
            /* ------------------------------------------------------------------------------*/
            /* ----------- Beginn Search three step lengths for parabolic search  -----------*/
            /* ------------------------------------------------------------------------------*/
            while((step2!=1)||(step1!=1)){
                
                /* --------------------- */
                /* calculate 3 L2 values */
                /* --------------------- */
                for (itest=itests;itest<=iteste;itest++){
                    
                    /* calculate change in the material parameters */
                    calc_mat_change_test(waveconv,waveconv_rho,waveconv_u,prho,prhonp1,ppi,ppinp1,pu,punp1,iter,1,FORWARD_ONLY,eps_scale,1,nfstart,Vs0,Vp0,Rho0,wavetype_start,s_LBFGS,N_LBFGS,LBFGS_NPAR,Vs_avg,Vp_avg,rho_avg,LBFGS_iter_start);
                    
                    /* For the calculation of the material parameters beteween gridpoints
                     the have to be averaged. For this, values lying at 0 and NX+1,
                     for example, are required on the local grid. These are now copied from the
                     neighbouring grids */
                    if(!ACOUSTIC) {
                        matcopy_elastic(prhonp1, ppinp1, punp1);	/* no differentiation of elastic and viscoelastic modelling because the viscoelastic parameters did not change during the forward modelling */
                    }else{
                        matcopy_acoustic(prhonp1, ppinp1);
                    }
                    
                    MPI_Barrier(MPI_COMM_WORLD);
                    
                    if(!ACOUSTIC) av_mue(punp1,puipjp,prhonp1);
                    av_rho(prhonp1,prip,prjp);
                    
                    /* Preparing memory variables for update_s (viscoelastic) */
                    if (L) {
                        if(!ACOUSTIC) {
                            prepare_update_s(etajm,etaip,peta,fipjp,punp1,puipjp,ppinp1,prhonp1,ptaus,ptaup,ptausipjp,f,g, bip,bjm,cip,cjm,dip,d,e);
                        }else{
                            prepare_update_p(etajm,peta,ppinp1,prhonp1,ptaup,g,bjm,cjm,e);
                        }
                    }
                    
                    /* initialization of L2 calculation */
                    L2=0.0;
                    L2_SH=0.0;
                    
                    alphanom = 0.0;
                    alphadenom = 0.0;
                    
                    exchange_par();
                    
                    if (RUN_MULTIPLE_SHOTS) nshots=nsrc; else nshots=1;
                    
                    /* ==============================================================================================*/
                    /* ==================================== start of loop over shots (test forward) =================*/
                    /* ==============================================================================================*/
                    for (ishot=TESTSHOT_START;ishot<=TESTSHOT_END;ishot=ishot+TESTSHOT_INCR){
                        /* step length estimation uses not all aviable shots */
                        fprintf(FP,"\n=================================================================================================\n");
                        fprintf(FP,"\n *****  Starting simulation (test-forward model) no. %d for shot %d of %d (rel. step length %.5f) \n",itest,ishot,nshots,eps_scale);
                        fprintf(FP,"\n=================================================================================================\n\n");
                              
			if (SEISMO && READREC==2){
			  if (ntr>0) {
			    dealloc_sections(ntr,ns,recpos_loc,sectionvx,sectionvy,sectionvz,sectionp,sectionpnp1,sectionpn,sectioncurl,sectiondiv,
					     sectionpdata,sectionpdiff,sectionpdiffold,sectionvxdata,sectionvxdiff,sectionvxdiffold,sectionvydata,
					     sectionvydiff,sectionvydiffold,sectionvzdata,sectionvzdiff,sectionvzdiffold);
			  }
			  free_imatrix(recpos,1,3,1,ntr_glob);
			  recpos=receiver(&ntr, srcpos, ishot);
			  recpos_loc = splitrec(recpos,&ntr_loc, ntr, recswitch);
			  ntr_glob=ntr;
			  ntr=ntr_loc;
			  
			  if (ntr>0){
				alloc_sections(ntr,ns,&sectionvx,&sectionvy,&sectionvz,&sectionp,&sectionpnp1,&sectionpn,&sectioncurl,&sectiondiv,
					      &sectionpdata,&sectionpdiff,&sectionpdiffold,&sectionvxdata,&sectionvxdiff,&sectionvxdiffold,&sectionvydata,
					      &sectionvydiff,&sectionvydiffold,&sectionvzdata,&sectionvzdiff,&sectionvzdiffold);
			  }			  
			}
	
                        for (nt=1;nt<=8;nt++) srcpos1[nt][1]=srcpos[nt][ishot];
                        
                        /*-----------------------------------*/
                        /* determine source position on grid */
                        /*-----------------------------------*/
                        if (RUN_MULTIPLE_SHOTS){
                            /* find this single source positions on subdomains */
                            if (nsrc_loc>0) free_matrix(srcpos_loc,1,8,1,1);
                            srcpos_loc=splitsrc(srcpos1,&nsrc_loc, 1);
                        }else{
                            /* Distribute multiple source positions on subdomains */
                            srcpos_loc = splitsrc(srcpos,&nsrc_loc, nsrc);
                        }
                        
                        
                        /*-------------------*/
                        /* calculate wavelet */
                        /*-------------------*/
                        /* calculate wavelet for each source point P SV */
                        if(WAVETYPE==1||WAVETYPE==3){
                            signals=NULL;
                            signals=wavelet(srcpos_loc,nsrc_loc,ishot,0,0);
                        }
                        /* calculate wavelet for each source point SH */
                        if(WAVETYPE==2||WAVETYPE==3){
                            signals_SH=NULL;
                            signals_SH=wavelet(srcpos_loc,nsrc_loc,ishot,1,0);
                        }
                        
                        /*------------------------------------------------------------------------------*/
                        /*----------- Start of Time Domain Filtering -----------------------------------*/
                        /*------------------------------------------------------------------------------*/
                        /*time domain filtering of the source signal */
                        if(WAVETYPE==1||WAVETYPE==3){
                            if(((TIME_FILT==1) || (TIME_FILT==2)) && (INV_STF==0)){
                                timedomain_filt(signals,F_LOW_PASS,ORDER,nsrc_loc,ns,1);
                            }
                            
                        }
                        
                        /*time domain filtering of the source signal */
                        if(WAVETYPE==2||WAVETYPE==3){
                            if(((TIME_FILT==1) || (TIME_FILT==2)) && (INV_STF==0)){
                                timedomain_filt(signals_SH,F_LOW_PASS,ORDER,nsrc_loc,ns,1);
                            }
                            
                        }
                        /*------------------------------------------------------------------------------*/
                        /*----------- End of Time Domain Filtering -------------------------------------*/
                        /*------------------------------------------------------------------------------*/
                        
                        /* initialize wavefield with zero */
                        if (L){
                            if(!ACOUSTIC) {
                                zero_fdveps_visc(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,pvz,psxx,psyy,psxy,psxz,psyz,ux,uy,uxy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_sxz_x,psi_vxx,psi_vyx,psi_vzx,psi_syy_y,psi_sxy_y,psi_syz_y,psi_vyy,psi_vxy,psi_vzy,psi_vxxs,pr,pp,pq,pt,po);
                            }else {
                                zero_fdveps_viscac(-nd+1, NY+nd, -nd+1, NX+nd, pvx, pvy, psp, pvxp1, pvyp1, psi_sxx_x, psi_sxy_x, psi_vxx, psi_vyx, psi_syy_y, psi_sxy_y, psi_vyy, psi_vxy, psi_vxxs, pp);
                            }
                        }else{
                            if(!ACOUSTIC) {
                                zero_fdveps(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,pvz,psxx,psyy,psxy,psxz,psyz,ux,uy,uxy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_sxz_x,psi_vxx,psi_vyx,psi_vzx,psi_syy_y,psi_sxy_y,psi_syz_y,psi_vyy,psi_vxy,psi_vzy,psi_vxxs);
                            }else{
                                zero_fdveps_ac(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,psp,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs);
                            }
                        }
                        
                        
                        lsnap=iround(TSNAP1/DT);
                        lsamp=NDT;
                        nsnap=0;
                        
                        /*-------------------------------------------------------------------------------*/
                        /*----------------------  loop over timesteps (forward model) step length-------*/
                        /*-------------------------------------------------------------------------------*/
                        
                        if((!VERBOSE)&&(MYID==0)) fprintf(FP,"\n ****************************************\n ");
                        
                        for (nt=1;nt<=NT;nt++){
                            
                            // Ratio to give output to stout
                            infoout = !(nt%nt_out);
                            if((!VERBOSE)&&(MYID==0)) if(!(nt%(NT/40))) fprintf(FP,"*");
                            
                            /* Check if simulation is still stable P and SV */
                            if (WAVETYPE==1 || WAVETYPE==3) {
                                if (isnan(pvy[NY/2][NX/2])) {
                                    fprintf(FP,"\n Time step: %d; pvy: %f \n",nt,pvy[NY/2][NX/2]);
                                    declare_error(" Simulation is unstable !");
                                }
                            }
                            /* Check if simulation is still stable SH */
                            if (WAVETYPE==2 || WAVETYPE==3) {
                                if (isnan(pvz[NY/2][NX/2])) {
                                    fprintf(FP,"\n Time step: %d; pvy: %f \n",nt,pvy[NY/2][NX/2]);
                                    declare_error(" Simulation is unstable !");
                                }
                            }
                            
                            if (MYID==0){
                                if (infoout)  fprintf(FP,"\n Computing timestep %d of %d \n",nt,NT);
                                time3=MPI_Wtime();
                            }
                            
                            /* update of particle velocities */
                            if(!ACOUSTIC) {
                                if (WAVETYPE==1 || WAVETYPE==3) {
                                    update_v_PML(1, NX, 1, NY, nt, pvx, pvxp1, pvxm1, pvy, pvyp1, pvym1, uttx, utty, psxx, psyy, psxy, prip, prjp, srcpos_loc,signals,signals,nsrc_loc,absorb_coeff,hc,infoout,0, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_sxx_x, psi_syy_y, psi_sxy_y, psi_sxy_x);
                                }
                                
                                if (WAVETYPE==2 || WAVETYPE==3) {
                                    update_v_PML_SH(1, NX, 1, NY, nt, pvz, pvzp1, pvzm1, psxz, psyz,prjp, srcpos_loc, signals, signals_SH, nsrc_loc, absorb_coeff,hc,infoout,0, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_sxz_x, psi_syz_y);
                                }
                            }
                            else
                                update_v_acoustic_PML(1, NX, 1, NY, nt, pvx, pvxp1, pvxm1, pvy, pvyp1, pvym1, psp, prip, prjp, srcpos_loc,signals,signals,nsrc_loc,absorb_coeff,hc,infoout,0, K_x_half, a_x_half, b_x_half, K_y_half, a_y_half, b_y_half, psi_sxx_x, psi_syy_y);
                            
                            if (MYID==0){
                                time4=MPI_Wtime();
                                time_av_v_update+=(time4-time3);
                                if (infoout)  fprintf(FP," particle velocity exchange between PEs ...");
                            }
                            
                            /* exchange of particle velocities between PEs */
                            exchange_v(pvx, pvy,pvz, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, req_send, req_rec,wavetype_start);
                            
                            if (MYID==0){
                                time5=MPI_Wtime();
                                time_av_v_exchange+=(time5-time4);
                                if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time5-time4);
                            }
                            
                            if (L) {   /* viscoelastic */
                                if (WAVETYPE==1 || WAVETYPE==3) {
                                    if(!ACOUSTIC) {
                                        update_s_visc_PML(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppinp1, punp1, puipjp, prhonp1, hc, infoout, pr, pp, pq, fipjp, f, g, bip, bjm, cip, cjm, d, e, dip, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);
                                    } else {
                                        update_p_visc_PML(1, NX, 1, NY, pvx, pvy, psp, ppinp1, prhonp1, hc, infoout, pp, g, bjm, cjm, e, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);
                                    }
                                }
                                if (WAVETYPE==2 || WAVETYPE==3) {
                                    update_s_visc_PML_SH(1, NX, 1, NY, pvz, psxz, psyz, pt, po, bip, bjm, cip, cjm, d, dip,fipjp, f, hc,infoout, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half,psi_vzx, psi_vzy);
                                }
                            } else {   /* elastic */
                                if(!ACOUSTIC){
                                    if (WAVETYPE==1 || WAVETYPE==3) {
                                        update_s_elastic_PML(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppinp1, punp1, puipjp, absorb_coeff, prhonp1, hc, infoout, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);
                                    }
                                    if (WAVETYPE==2 || WAVETYPE==3) {
                                        update_s_elastic_PML_SH(1, NX, 1, NY, pvz,psxz,psyz,uxz,uyz,hc,infoout, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half,psi_vzx, psi_vzy,puipjp,punp1,prhonp1);
                                    }
                                }
                                else /* acoustic */
                                    update_p_PML(1, NX, 1, NY, pvx, pvy, psp, u, ppinp1, absorb_coeff, prhonp1, hc, infoout, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);
                            }
                            
                            if (MYID==0){
                                time6=MPI_Wtime();
                                time_av_s_update+=(time6-time5);
                                if (infoout)  fprintf(FP," stress exchange between PEs ...");
                            }
                            
                            /* explosive source */
                            if ((SOURCE_TYPE==1))
                                psource(nt,psxx,psyy,psp,srcpos_loc,signals,nsrc_loc,0);
                            
                            /* Applying free surface condition */
                            if ((FREE_SURF) && (POS[2]==0)){
                                if (!ACOUSTIC){
                                    if (L){
                                        /* viscoelastic */
                                        surface_PML(1, pvx, pvy, psxx, psyy, psxy,psyz, pp, pq, ppinp1, punp1, prhonp1, ptaup, ptaus, etajm, peta, hc, K_x, a_x, b_x, psi_vxxs, ux, uy,uxy,uyz,psxz,uxz);
                                    }else{
                                        /* elastic */
                                        surface_elastic_PML(1, pvx, pvy, psxx, psyy, psxy,psyz, ppinp1, punp1, prhonp1, hc, K_x, a_x, b_x, psi_vxxs, ux, uy, uxy,uyz,psxz,uxz);
                                    }
                                } else {
                                    /* viscoelastic and elastic ACOUSTIC */
                                    surface_acoustic_PML(1, psp);
                                }
                            }
                            
                            /* stress exchange between PEs */
                            if(!ACOUSTIC)
                                exchange_s(psxx,psyy,psxy,psxz,psyz,bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, req_send, req_rec,wavetype_start);
                            else
                                exchange_p(psp,bufferlef_to_rig, bufferrig_to_lef,buffertop_to_bot, bufferbot_to_top,req_send, req_rec);
                            
                            if (MYID==0){
                                time7=MPI_Wtime();
                                time_av_s_exchange+=(time7-time6);
                                if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time7-time6);
                            }
                            
                            /* store amplitudes at receivers in section-arrays */
                            if (SEISMO){
                                seismo_ssg(nt, ntr, recpos_loc, sectionvx, sectionvy,sectionvz, sectionp, sectioncurl, sectiondiv, pvx, pvy,pvz, psxx, psyy, psp, ppinp1, punp1, hc);
                                /*lsamp+=NDT;*/
                            }
                            
                            
                            if (MYID==0){
                                time8=MPI_Wtime();
                                time_av_timestep+=(time8-time3);
                                if (infoout)  fprintf(FP," total real time for timestep %d : %4.2f s.\n",nt,time8-time3);
                            }
                            
                        }
                        
                        /*-------------------------------------------------------------------------------*/
                        /*------------------ end loop over timesteps (forward model) step length  -------*/
                        /*-------------------------------------------------------------------------------*/
                        
                        if((!VERBOSE)&&(MYID==0)) fprintf(FP,"\n\n");
                        
                        if (MYID==0&&VERBOSE){
                            printf("Calculate residuals between test forward model m - mu * dm and actual model m \n");
                            printf("----------------------------------------------------------------------------- \n");
                        }
                        
                        /*-----------------------------------*/
                        /*------- Calculate residuals -------*/
                        /*-----------------------------------*/
                        
                        if (ntr > 0){
                            if(WAVETYPE==1 || WAVETYPE==3){
                                /* --------------------------------- */
                                /* read seismic data from SU file vx */
                                /* --------------------------------- */
                                
                                if((ADJOINT_TYPE==1)||(ADJOINT_TYPE==3)){ /* if ADJOINT_TYPE */
                                    inseis(fprec,ishot,sectionread,ntr_glob,ns,1,iter);
                                    if ((TIME_FILT==1 )|| (TIME_FILT==2)){
                                        timedomain_filt(sectionread,F_LOW_PASS,ORDER,ntr_glob,ns,1);
                                    }
                                    h=1;
                                    for(i=1;i<=ntr;i++){
                                        for(j=1;j<=ns;j++){
                                            sectionvxdata[h][j]=sectionread[recpos_loc[3][i]][j];
                                        }
                                        h++;
                                    }
                                    L2=calc_res(sectionvxdata,sectionvx,sectionvxdiff,sectionvxdiffold,ntr,ns,LNORM,L2,0,1,1,ntr_glob,recpos_loc,nsrc_glob,ishot,iter,srcpos,recpos);
                                }
                                
                                /* --------------------------------- */
                                /* read seismic data from SU file vy */
                                /* --------------------------------- */
                                
                                if((ADJOINT_TYPE==1)||(ADJOINT_TYPE==2)){ /* if ADJOINT_TYPE */
                                    inseis(fprec,ishot,sectionread,ntr_glob,ns,2,iter);
                                    if ((TIME_FILT==1 )|| (TIME_FILT==2)){
                                        timedomain_filt(sectionread,F_LOW_PASS,ORDER,ntr_glob,ns,1);
                                    }
                                    h=1;
                                    for(i=1;i<=ntr;i++){
                                        for(j=1;j<=ns;j++){
                                            sectionvydata[h][j]=sectionread[recpos_loc[3][i]][j];
                                        }
                                        h++;
                                    }
                                    L2=calc_res(sectionvydata,sectionvy,sectionvydiff,sectionvydiffold,ntr,ns,LNORM,L2,0,1,1,ntr_glob,recpos_loc,nsrc_glob,ishot,iter,srcpos,recpos);
                                }
                                
                                /* --------------------------------- */
                                /* read seismic data from SU file p */
                                /* --------------------------------- */
                                if(ADJOINT_TYPE==4){ /* if ADJOINT_TYPE */
                                    inseis(fprec,ishot,sectionread,ntr_glob,ns,9,iter);
                                    if ((TIME_FILT==1 )|| (TIME_FILT==2)){
                                        timedomain_filt(sectionread,F_LOW_PASS,ORDER,ntr_glob,ns,1);
                                    }
                                    h=1;
                                    for(i=1;i<=ntr;i++){
                                        for(j=1;j<=ns;j++){
                                            sectionpdata[h][j]=sectionread[recpos_loc[3][i]][j];
                                        }
                                        h++;
                                    }
                                    L2=calc_res(sectionpdata,sectionp,sectionpdiff,sectionpdiffold,ntr,ns,LNORM,L2,0,1,1,ntr_glob,recpos_loc,nsrc_glob,ishot,iter,srcpos,recpos);
                                }
                            }
                            
                            /* --------------------------------- */
                            /* read seismic data from SU file vz */
                            /* --------------------------------- */
                            if(WAVETYPE==2 || WAVETYPE==3){
                                inseis(fprec,ishot,sectionread,ntr_glob,ns,10,iter);
                                if ((TIME_FILT==1 )|| (TIME_FILT==2)){
                                    timedomain_filt(sectionread,F_LOW_PASS,ORDER,ntr_glob,ns,1);
                                }
                                h=1;
                                for(i=1;i<=ntr;i++){
                                    for(j=1;j<=ns;j++){
                                        sectionvzdata[h][j]=sectionread[recpos_loc[3][i]][j];
                                    }
                                    h++;
                                }
                                L2_SH=calc_res(sectionvzdata,sectionvz,sectionvzdiff,sectionvzdiffold,ntr,ns,LNORM,L2_SH,0,1,1,ntr_glob,recpos_loc,nsrc_glob,ishot,iter,srcpos,recpos);
                            }
                        }
                        
                    }
                    /* ==============================================================================================*/
                    /* ==================================== end of loop over shots (test forward) ===================*/
                    /* ==============================================================================================*/
                    epst1[itest]=eps_scale;
                    epst1[1] = 0.0;
                    
                    if(WAVETYPE==1||WAVETYPE==3){
                        L2sum=0.0;
                        MPI_Allreduce(&L2,&L2sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
                    }
                    
                    if(WAVETYPE==2||WAVETYPE==3){
                        L2sum_SH=0.0;
                        MPI_Allreduce(&L2_SH,&L2sum_SH,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
                    }
                    
                    switch (LNORM){
                        case 2:
                            L2t[itest]=0.0;
                            
                            if(JOINT_EQUAL_WEIGHTING){
                            
                                L2t[itest]+=(L2sum/energy_sum)/JOINT_EQUAL_PSV;
                                
                                L2t[itest]+=(L2sum_SH/energy_sum_SH)/JOINT_EQUAL_SH;
                                
                                break;
                            }
                            
                            if(WAVETYPE==1||WAVETYPE==3){
                                L2t[itest]+=L2sum/energy_sum;
                            }
                            
                            if(WAVETYPE==2||WAVETYPE==3){
                                L2t[itest]+=L2sum_SH/energy_sum_SH;
                            }
                            
                            break;
                        case 7:
                            if (TRKILL){
                                if(ADJOINT_TYPE==1){	/* x and y component are used in the inversion */
                                    L2t[itest]=2.0*(1.0+(L2sum/((float)((NO_OF_TESTSHOTS*ntr_glob-sum_killed_traces_testshots)*2.0))));
                                    if (MYID==0){
                                        printf("sum_killed_traces_testshots=%d\n",sum_killed_traces_testshots);}}
                                else{
                                    L2t[itest]=2.0*(1.0+(L2sum/((float)(NO_OF_TESTSHOTS*ntr_glob-sum_killed_traces_testshots))));
                                    if (MYID==0){
                                        printf("sum_killed_traces_testshots=%d\n",sum_killed_traces_testshots);
                                        printf("ntr_glob=%d\n",ntr_glob);
                                        printf("NO_OF_TESTSHOTS=%d\n",NO_OF_TESTSHOTS);}}}
                            else{
                                if(ADJOINT_TYPE==1){	/* x and y component are used in the inversion */
                                    L2t[itest]=2.0*(1.0+(L2sum/((float)NO_OF_TESTSHOTS*(float)ntr_glob*2.0)));}
                                else{
                                    L2t[itest]=2.0*(1.0+(L2sum/((float)NO_OF_TESTSHOTS*(float)ntr_glob)));}
                            }
                            break;
                        case 8:
                            L2t[itest]=L2sum/energy_sum;
                            break;
                        default:
                            L2t[itest]=L2sum;
                            break;
                    }
                    
                }
                /* --------------------- */
                /*    end of L2 test     */
                /* --------------------- */
                
                
                if(GRAD_METHOD==2&&(iter>LBFGS_iter_start)&&LBFGS_STEP_LENGTH) {
                    
                    if(itests==3) {
                        epst1[3]=eps_scale;
                        break;
                    }
                    
                    if(itests==2) {
                        epst1[2]=eps_scale;
                        itests=3;iteste=3;
                        eps_scale=1;
                    }
                    
                    if(itests==1) {
                        epst1[1]=eps_scale;
                        itests=2;iteste=2;
                        eps_scale=0.6;
                    }
                    
                } else {
                    /* step1 search L2t[2] < L2t[1]  */
                    /* step2 search L2t[2] < L2t[3]  */
                    /* -> parabolic fit              */
                    
                    /* Did not find a step size which reduces the misfit function */
                    /* Reduce step length                                         */
                    if((step1==0)&&(L2t[1]<=L2t[2])){
                        eps_scale = eps_scale/scalefac;
                        countstep++;
                    }
                    
                    /* Found a step size with L2t[2] < L2t[3]*/
                    /* -> now parabolic fit                  */
                    if((step1==1)&&(L2t[2]<L2t[3])){
                        epst1[3]=eps_scale;
                        step2=1;
                    }
                    
                    /* Could not find a step size with L2t[2] < L2t[3]*/
                    /* increase step length to find  a larger misfit function than L2t[2]*/
                    if((step1==1)&&(L2t[2]>=L2t[3])){
                        epst1[3]=eps_scale;
                        eps_scale = eps_scale + (eps_scale/scalefac);
                        countstep++;
                    }
                    
                    /* found a step size which reduces the misfit function */
                    /* -> next step2                                       */
                    if((step1==0)&&(L2t[1]>L2t[2])){
                        epst1[2]=eps_scale;
                        step1=1;
                        iteste=3;
                        itests=3;
                        countstep=0;
                        /* find a second step length with a larger misfit function than L2t[2]*/
                        eps_scale = eps_scale + (eps_scale/scalefac);
                    }
                    
                    step3=0;
                    
                    if((step1==0)&&(countstep>stepmax)){
                        if(MYID==0){
                            printf(" Steplength estimation failed!\n\n");}
                        step3=1;
                        break;
                    }
                    
                    if((step1==1)&&(countstep>stepmax)){
                        if(MYID==0){
                            printf("Could not find a proper 3rd step length which brackets the minimum\n");}
                        step1=1;
                        step2=1;
                    }
                }
                if((MYID==0)){
                    fprintf(FP,"\n=================================================================================================\n");
                    printf("Step length estimation status at iteration %i\n",iter);
                    if(VERBOSE){printf("iteste = %d \t itests = %d \t step1 = %d \t step2 = %d \t eps_scale = %e \t countstep = %d \t stepmax= %d \t scalefac = %e \t MYID = %d \t L2t[1] = %e \t L2t[2] = %e \t L2t[3] = %e \n",iteste,itests,step1,step2,eps_scale,countstep,stepmax,scalefac,MYID,L2t[1],L2t[2],L2t[3]);}
                    if((!VERBOSE)){printf("FD-Simulation %d of %d\n",countstep+1,stepmax+2);}
                    if((!VERBOSE)){printf("L2-Norm[1] = %e\n",L2t[1]);}
                    if((!VERBOSE)){printf("L2-Norm[2] = %e\n",L2t[2]);}
                    if((!VERBOSE)){printf("L2-Norm[3] = %e\n",L2t[3]);}
                    if((VERBOSE)){printf("iteste = %d \t itests = %d \t step1 = %d \t step2 = %d \t eps_scale = %e \n",iteste,itests,step1,step2,eps_scale);}
                }
                
            }
            /* ------------------------------------------------------------------------------*/
            /* ----------- End Search three step lengths for parabolic search  -----------*/
            /* ------------------------------------------------------------------------------*/
            
            if(GRAD_METHOD==2&&(iter>LBFGS_iter_start)&&LBFGS_STEP_LENGTH) {
                if(L2t[1]>L2t[2] && L2t[3]>L2t[2]) {
                    step1=1;
                }
                
                if(L2t[1]>L2t[2] && L2t[2]>L2t[3]) {
                    step1=1;
                }
                
                if(L2t[1]<L2t[2] && L2t[2]<L2t[3]) {
                    step1=1;
                }
                
                if(L2t[1]>L2t[2] && L2t[2]>L2t[3]) {
                    step1=1;
                }
            }
            
            if(step1==1){ /* only find an optimal step length if step1==1 */
                /* calculate optimal step length epsilon for Vp and Vs*/
                if(MYID==0){
                    printf("======================================================================================= \n");
                    printf("calculate optimal step length epsilon for Vp, Vs and density (step1 == %d, step3 == %d) \n",step1,step3);
                    printf("======================================================================================= \n");
                }
                opteps_vp=calc_opt_step(L2t,waveconv,gradg,epst1,1,C_vp);
                eps_scale = opteps_vp;
            }else{
                if (iter < min_iter_help){
                    eps_scale=EPS_SCALE;
                    opteps_vp=EPS_SCALE;
                }
                else{
                    eps_scale=0.0;
                    opteps_vp=0.0;}
                
                if(GRAD_METHOD==2&&(iter>LBFGS_iter_start)&&LBFGS_STEP_LENGTH) {
                    eps_scale=0.05;
                    opteps_vp=eps_scale;
                    
                    if(MYID==0) printf("\n\n=========================================");
                    if(MYID==0) printf("\n L-BFGS step length search failed");
                    if(MYID==0) printf("\n eps_scale is set to %f",eps_scale);
                    if(MYID==0) printf("\n=========================================\n\n");
                }
                
                if(MYID==0){
                    printf("================================================================================== \n");
                    printf("using EPS_SCALE for the step length (step1 == %d, step3 == %d, EPS_SCALE == %4.2f) \n",step1,step3,EPS_SCALE);
                    printf("================================================================================== \n");
                }
            }
            if(MYID==0) fprintf(FP,"\n=================================================================================================\n");
            if(MYID==0){printf("Step length search final status:\n");}
            if(MYID==0){
                printf("MYID = %d \t opteps_vp = %e \t opteps_vs = %e \t opteps_rho = %e \n",MYID,opteps_vp,opteps_vs,opteps_rho);
                if(VERBOSE) printf("MYID = %d \t L2t[1] = %e \t L2t[2] = %e \t L2t[3] = %e \t Final-L2= %e \n",MYID,L2t[1],L2t[2],L2t[3],L2t[4]);
                
                if(MYID==0&&(!VERBOSE)){printf("L2-Norm[1] = %e\n",L2t[1]);}
                if(MYID==0&&(!VERBOSE)){printf("L2-Norm[2] = %e\n",L2t[2]);}
                if(MYID==0&&(!VERBOSE)){printf("L2-Norm[3] = %e\n",L2t[3]);}
                if(MYID==0&&(!VERBOSE)){printf("L2-Norm[4] = %e (final L2 of all shots)\n",L2t[4]);}
                printf("MYID = %d \t epst1[1] = %e \t epst1[2] = %e \t epst1[3] = %e \n",MYID,epst1[1],epst1[2],epst1[3]);
                if (TIME_FILT==0){
                    fprintf(FPL2,"%e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %f \n",opteps_vp,epst1[1],epst1[2],epst1[3],L2t[1],L2t[2],L2t[3],L2t[4],GAMMA);}
                else{
                    fprintf(FPL2,"%e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %f \t %f \n",opteps_vp,epst1[1],epst1[2],epst1[3],L2t[1],L2t[2],L2t[3],L2t[4],F_LOW_PASS,GAMMA);}
                if(WAVETYPE==3 && MYID==0){
                    fprintf(FPL2_JOINT,"%d \t %f \t %f\n",iter,L2sum_all_shots/energy_sum_all_shots,L2sum_all_shots_SH/energy_sum_all_shots_SH);
                }
            }
            
            /* saving history of final L2*/
            L2_hist[iter]=L2t[4];
            do_stf=0;
            
            
            /* -----------------------------------------------------------------------*/
            /* ----------- Do the actual update to the material parameters -----------*/
            /* -----------------------------------------------------------------------*/
            calc_mat_change_test(waveconv,waveconv_rho,waveconv_u,prho,prhonp1,ppi,ppinp1,pu,punp1,iter,1,FORWARD_ONLY,eps_scale,0,nfstart,Vs0,Vp0,Rho0,wavetype_start,s_LBFGS,N_LBFGS,LBFGS_NPAR,Vs_avg,Vp_avg,rho_avg,LBFGS_iter_start);
            fprintf(FP,"=================================================================================================\n");
            
        } /* end of if(FORWARD_ONLY!=4) */
        
        /* ------------------------------------*/
        /* smoothing the models vp, vs and rho */
        /* ------------------------------------*/
        if (FORWARD_ONLY==0 && (opteps_vp>0.0 || WOLFE_CONDITION)){
            if(!ACOUSTIC){
                if(WAVETYPE==1||WAVETYPE==3) if(MODEL_FILTER)smooth(ppi,4,2,Vs_avg,F_LOW_PASS);
                if(MODEL_FILTER)smooth(pu,5,2,Vs_avg,F_LOW_PASS);
                if(MODEL_FILTER)smooth(prho,6,2,Vs_avg,F_LOW_PASS);
            }else{
                if(WAVETYPE==1||WAVETYPE==3) if(MODEL_FILTER)smooth(ppi,4,2,Vp_avg,F_LOW_PASS);
                if(MODEL_FILTER)smooth(prho,6,2,Vp_avg,F_LOW_PASS);
            }
        }
        
        
        if(FORWARD_ONLY!=1){
            if(MYID==0){
                fclose(FPL2);
            }
            if(WAVETYPE==3 && MYID==0) {
                fclose(FPL2_JOINT);
            }
        }
        
        if(iter==nfstart){
            nfstart = nfstart + nf;
        }
        
        if(iter==nfstart_jac){
            nfstart_jac = nfstart_jac + nf_jac;
        }
        
        if (MYID==0&&VERBOSE){
            fprintf(FP,"\n **Info from main (written by PE %d): \n",MYID);
            fprintf(FP," CPU time of program per PE: %li seconds.\n",clock()/CLOCKS_PER_SEC);
            time8=MPI_Wtime();
            fprintf(FP," Total real time of program: %4.2f seconds.\n",time8-time1);
            time_av_v_update=time_av_v_update/(double)NT;
            time_av_s_update=time_av_s_update/(double)NT;
            time_av_v_exchange=time_av_v_exchange/(double)NT;
            time_av_s_exchange=time_av_s_exchange/(double)NT;
            time_av_timestep=time_av_timestep/(double)NT;
            fprintf(FP," Average times for \n");
            fprintf(FP," velocity update:  \t %5.3f seconds  \n",time_av_v_update);
            fprintf(FP," stress update:  \t %5.3f seconds  \n",time_av_s_update);
            fprintf(FP," velocity exchange:  \t %5.3f seconds  \n",time_av_v_exchange);
            fprintf(FP," stress exchange:  \t %5.3f seconds  \n",time_av_s_exchange);
            fprintf(FP," timestep:  \t %5.3f seconds  \n",time_av_timestep);
            
        }
        
        
        /* ----------------------------------------------*/
        /* ----------- Check abort criteriums -----------*/
        /* ----------------------------------------------*/
        if (iter>min_iter_help){
            
            float diff=0.0, pro=PRO;
            
            /* calculating differnce of the actual L2 and before two iterations, dividing with L2_hist[iter-2] provide changing in procent*/
            diff=fabs((L2_hist[iter-2]-L2_hist[iter])/L2_hist[iter-2]);
            
            /* abort criterion: if diff is smaller than pro (1% ?? is this reasonable?) than the inversion abort or switch to another frequency range*/
            
            
            /* ------------------------------------------- */
            /*  Check when NO workflow and TIME_FILT==0    */
            /* ------------------------------------------- */
            if((diff<=pro)&&(TIME_FILT==0)&&USE_WORKFLOW==0){
                if(MYID==0){
                    printf("\n Reached the abort criterion of pro = %4.2f: diff = %4.2f \n",pro,diff);
                }
                
                break;
            }
            
            /* abort criterion: did not found a step length which decreases the misfit*/
            if((step3==1||wolfe_SLS_failed)&&(TIME_FILT==0&&USE_WORKFLOW==0)){
                if(MYID==0){
                    printf("\n Did not find a step length which decreases the misfit.\n");
                }
                step3=0;
                break;
            }
            
            /* ------------------------------------------- */
            /*       Check when Workflow is used           */
            /* ------------------------------------------- */
            if(USE_WORKFLOW && ( diff<=pro || wolfe_SLS_failed || step3==1 )){
                
                if(workflow_lines==WORKFLOW_STAGE){
                    fprintf(FP,"\n Reached the abort criterion of pro = %4.2f: diff = %4.2f \n",pro,diff);
                    fprintf(FP,"\n No new line in workflow file");
                    break;
                }
                
                if(diff<=pro) fprintf(FP,"\n Reached the abort criterion of pro = %4.2f: diff = %4.2f \n",pro,diff);
                if(wolfe_SLS_failed) fprintf(FP,"\n Wolfe step length search failed \n");
                fprintf(FP,"\n Switching to next line in workflow");
                
                WORKFLOW_STAGE++;
                do_stf=1;
                min_iter_help=0;
                min_iter_help=iter+MIN_ITER;
                
                /* Sync WORKFLOW_STAGE on all PEs */
                buf1=WORKFLOW_STAGE;
                buf2=0;
                MPI_Allreduce(&buf1,&buf2, 1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
                WORKFLOW_STAGE=buf2;
                
                alpha_SL_old=1;
                
                /* Restart L-BFGS at next iteration */
                LBFGS_iter_start=iter+1;
                
                /* Restart conjugate gradient at next iteration */
                PCG_iter_start=iter+1;
                
                wolfe_SLS_failed=0;
                
                step3=0;
                JOINT_EQUAL_new_max=1;
            }
            
            /* ------------------------------------------------- */
            /* Check when Workflow is NOT used and TIME_FILT==1  */
            /* ------------------------------------------------- */
            if( (TIME_FILT==1) && (!USE_WORKFLOW) && ( (diff<=pro) || (step3==1) || wolfe_SLS_failed ) ){
                if(MYID==0){
                    if (diff<=pro){
                        printf("\n Reached the abort criterion of pro = %4.2f: diff = %4.2f \n",pro,diff);}
                    if(step3==1){
                        printf("\n Did not find a step length which decreases the misfit.\n");}
                }
                
                F_LOW_PASS=F_LOW_PASS+F_LOW_PASS_INCR;
                do_stf=1;
                min_iter_help=0;
                min_iter_help=iter+MIN_ITER;
                
                if(F_LOW_PASS>F_LOW_PASS_END){
                    if(MYID==0){
                        printf("\n Reached the maximum frequency of %4.2f Hz \n",F_LOW_PASS);
                    }
                    break;
                }
                
                if(MYID==0) printf("\n Changing to corner frequency of %4.2f Hz \n",F_LOW_PASS);
                
                /* Restart L-BFGS at next iteration */
                LBFGS_iter_start=iter+1;
                
                /* Restart conjugate gradient at next iteration */
                PCG_iter_start=iter+1;
                
                wolfe_SLS_failed=0;
                alpha_SL_old=1;
                
                step3=0;
                JOINT_EQUAL_new_max=1;
            }
            
            /* ------------------------------------------------- */
            /* Check when Workflow is NOT used and TIME_FILT==2  */
            /* ------------------------------------------------- */
            if( (TIME_FILT==2) && (!USE_WORKFLOW) && ( (diff<=pro) || (step3==1) || wolfe_SLS_failed ) ){
                if(MYID==0){
                    if (diff<=pro){
                        printf("\n Reached the abort criterion of pro = %4.2f: diff = %4.2f \n",pro,diff);}
                    if(step3==1){
                        printf("\n Did not find a step length which decreases the misfit.\n");}
                }
                
                if(FREQ_NR==nfrq){
                    if(MYID==0){
                        printf("\n Finished at the maximum frequency of %4.2f Hz \n",F_LOW_PASS);
                    }
                    break;
                }
                
                FREQ_NR=FREQ_NR+1;
                F_LOW_PASS=F_LOW_PASS_EXT[FREQ_NR];
                do_stf=1;
                min_iter_help=0;
                min_iter_help=iter+MIN_ITER;
                if(MYID==0) printf("\n Changing to corner frequency of %4.2f Hz \n",F_LOW_PASS);
                
                /* Restart L-BFGS at next iteration */
                LBFGS_iter_start=iter+1;
                
                /* Restart conjungate gradient at next iteration */
                PCG_iter_start=iter+1;
                
                wolfe_SLS_failed=0;
                alpha_SL_old=1;
                
                step3=0;
                JOINT_EQUAL_new_max=1;
            }
            
        }
        

        
    }
    /*------------------------------------------------------------------------------*/
    /*----------- End fullwaveform iteration loop ----------------------------------*/
    /*------------------------------------------------------------------------------*/
    

    /* ====================================== */
    /* ====== deallocation of memory =========*/
    /* ====================================== */
    if (SEISMO){
        free_imatrix(recpos,1,3,1,ntr_glob);
        free_ivector(recswitch,1,ntr_glob);
    }
    
    /* free memory for abort criterion */
    free_vector(L2_hist,1,1000);
    
    if(INV_STF) free_matrix(fulldata,1,ntr_glob,1,NT);
    
    free_ivector(DTINV_help,1,NT);
    
    if(!ACOUSTIC){
        if (WAVETYPE==1 || WAVETYPE==3) {
            free_matrix(psxx,-nd+1,NY+nd,-nd+1,NX+nd);
            free_matrix(psxy,-nd+1,NY+nd,-nd+1,NX+nd);
            free_matrix(psyy,-nd+1,NY+nd,-nd+1,NX+nd);
        }
        if (WAVETYPE==2 || WAVETYPE==3) {
            free_matrix(psyz,-nd+1,NY+nd,-nd+1,NX+nd);
            free_matrix(psxz,-nd+1,NY+nd,-nd+1,NX+nd);
        }
    }else{
        free_matrix(psp,-nd+1,NY+nd,-nd+1,NX+nd);
    }
    if (WAVETYPE==1 || WAVETYPE==3) {
        free_matrix(pvx,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(pvy,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(pvxp1,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(pvyp1,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(pvxm1,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(pvym1,-nd+1,NY+nd,-nd+1,NX+nd);
    }
    if (WAVETYPE==2 || WAVETYPE==3) {
        free_matrix(pvz,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(pvzp1,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(pvzm1,-nd+1,NY+nd,-nd+1,NX+nd);
    }
    if(!ACOUSTIC){
        if (WAVETYPE==1 || WAVETYPE==3) {
            free_matrix(ux,-nd+1,NY+nd,-nd+1,NX+nd);
            free_matrix(uy,-nd+1,NY+nd,-nd+1,NX+nd);
            free_matrix(uxy,-nd+1,NY+nd,-nd+1,NX+nd);
            free_matrix(uyx,-nd+1,NY+nd,-nd+1,NX+nd);
            free_matrix(uttx,-nd+1,NY+nd,-nd+1,NX+nd);
            free_matrix(utty,-nd+1,NY+nd,-nd+1,NX+nd);
        }
    }
    free_matrix(Vp0,-nd+1,NY+nd,-nd+1,NX+nd);
    if(!ACOUSTIC)
        free_matrix(Vs0,-nd+1,NY+nd,-nd+1,NX+nd);
    free_matrix(Rho0,-nd+1,NY+nd,-nd+1,NX+nd);
    
    free_matrix(prho,-nd+1,NY+nd,-nd+1,NX+nd);
    free_matrix(prhonp1,-nd+1,NY+nd,-nd+1,NX+nd);
    free_matrix(prip,-nd+1,NY+nd,-nd+1,NX+nd);
    free_matrix(prjp,-nd+1,NY+nd,-nd+1,NX+nd);
    free_matrix(pripnp1,-nd+1,NY+nd,-nd+1,NX+nd);
    free_matrix(prjpnp1,-nd+1,NY+nd,-nd+1,NX+nd);
    free_matrix(ppi,-nd+1,NY+nd,-nd+1,NX+nd);
    free_matrix(ppinp1,-nd+1,NY+nd,-nd+1,NX+nd);
    if(!ACOUSTIC){
        free_matrix(pu,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(punp1,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(puipjp,-nd+1,NY+nd,-nd+1,NX+nd);
    }
    free_matrix(vpmat,-nd+1,NY+nd,-nd+1,NX+nd);
    
    /* free memory for viscoelastic modeling variables */
    if (L) {
        free_f3tensor(pr,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        free_f3tensor(pp,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        free_f3tensor(pq,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        free_f3tensor(dip,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        free_f3tensor(d,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        free_f3tensor(e,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        free_matrix(ptaus,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(ptausipjp,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(ptaup,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(fipjp,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(f,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(g,-nd+1,NY+nd,-nd+1,NX+nd);
        free_vector(peta,1,L);
        free_vector(etaip,1,L);
        free_vector(etajm,1,L);
        free_vector(bip,1,L);
        free_vector(bjm,1,L);
        free_vector(cip,1,L);
        free_vector(cjm,1,L);
    }
    if(FORWARD_ONLY==0){
        free_matrix(waveconv,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(waveconv_lam,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(waveconv_shot,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(waveconvtmp,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(wcpart,1,3,1,3);
        free_matrix(wavejac,-nd+1,NY+nd,-nd+1,NX+nd);
        if(!ACOUSTIC){
            free_f3tensor(forward_prop_x,-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
        }else{
            free_f3tensor(forward_prop_p,-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
        }
        free_matrix(gradg,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(gradp,-nd+1,NY+nd,-nd+1,NX+nd);
        
        if(WAVETYPE==1 || WAVETYPE==3){
            free_f3tensor(forward_prop_rho_x,-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
            free_f3tensor(forward_prop_rho_y,-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
        }
        
        free_matrix(gradg_rho,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(gradp_rho,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(waveconv_rho,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(waveconv_rho_s,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(waveconv_rho_shot,-nd+1,NY+nd,-nd+1,NX+nd);
        
        if(!ACOUSTIC){
            free_f3tensor(forward_prop_u,-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
            free_matrix(gradg_u,-nd+1,NY+nd,-nd+1,NX+nd);
            free_matrix(gradp_u,-nd+1,NY+nd,-nd+1,NX+nd);
            free_matrix(waveconv_u,-nd+1,NY+nd,-nd+1,NX+nd);
            free_matrix(waveconv_mu,-nd+1,NY+nd,-nd+1,NX+nd);
            free_matrix(waveconv_u_shot,-nd+1,NY+nd,-nd+1,NX+nd);
        }
        
    }
    
    if(FW>0){
        free_vector(d_x,1,2*FW);
        free_vector(K_x,1,2*FW);
        free_vector(alpha_prime_x,1,2*FW);
        free_vector(a_x,1,2*FW);
        free_vector(b_x,1,2*FW);
        
        free_vector(d_x_half,1,2*FW);
        free_vector(K_x_half,1,2*FW);
        free_vector(alpha_prime_x_half,1,2*FW);
        free_vector(a_x_half,1,2*FW);
        free_vector(b_x_half,1,2*FW);
        
        free_vector(d_y,1,2*FW);
        free_vector(K_y,1,2*FW);
        free_vector(alpha_prime_y,1,2*FW);
        free_vector(a_y,1,2*FW);
        free_vector(b_y,1,2*FW);
        
        free_vector(d_y_half,1,2*FW);
        free_vector(K_y_half,1,2*FW);
        free_vector(alpha_prime_y_half,1,2*FW);
        free_vector(a_y_half,1,2*FW);
        free_vector(b_y_half,1,2*FW);
        if (WAVETYPE==1||WAVETYPE==3){
            free_matrix(psi_sxx_x,1,NY,1,2*FW);
            free_matrix(psi_syy_y,1,2*FW,1,NX);
            free_matrix(psi_sxy_x,1,NY,1,2*FW);
            free_matrix(psi_sxy_y,1,2*FW,1,NX);
            free_matrix(psi_vxx,1,NY,1,2*FW);
            free_matrix(psi_vxxs,1,NY,1,2*FW);
            free_matrix(psi_vyy,1,2*FW,1,NX);
            free_matrix(psi_vxy,1,2*FW,1,NX);
            free_matrix(psi_vyx,1,NY,1,2*FW);
        }
        if (WAVETYPE==2||WAVETYPE==3){
            free_matrix(psi_sxz_x,1,NY,1,2*FW);
            free_matrix(psi_syz_y,1,2*FW,1,NX);
            free_matrix(psi_vzx,1,NY,1,2*FW);
            free_matrix(psi_vzy,1,2*FW,1,NX);
        }
    }
    
    free_matrix(taper_coeff,1,NY,1,NX);
    free_matrix(bufferlef_to_rig,1,NY,1,fdo3);
    free_matrix(bufferrig_to_lef,1,NY,1,fdo3);
    free_matrix(buffertop_to_bot,1,NX,1,fdo3);
    free_matrix(bufferbot_to_top,1,NX,1,fdo3);
    switch (SEISMO){
        case 1 : /* particle velocities only */
            if (WAVETYPE==1 || WAVETYPE==3) {
                free_matrix(fulldata_vx,1,ntr_glob,1,NT);
                free_matrix(fulldata_vy,1,ntr_glob,1,NT);
            }
            if (WAVETYPE==2 || WAVETYPE==3) {
                free_matrix(fulldata_vz,1,ntr_glob,1,NT);
            }
            break;
        case 2 : /* pressure only */
            if (WAVETYPE==1 || WAVETYPE==3) {
                free_matrix(fulldata_p,1,ntr_glob,1,NT);
            }
            break;
        case 3 : /* curl and div only */
            if (WAVETYPE==1 || WAVETYPE==3) {
                free_matrix(fulldata_div,1,ntr_glob,1,NT);
                free_matrix(fulldata_curl,1,ntr_glob,1,NT);
            }
            break;
        case 4 : /* everything */
            if (WAVETYPE==1 || WAVETYPE==3) {
                free_matrix(fulldata_vx,1,ntr_glob,1,NT);
                free_matrix(fulldata_vy,1,ntr_glob,1,NT);
                free_matrix(fulldata_p,1,ntr_glob,1,NT);
                free_matrix(fulldata_div,1,ntr_glob,1,NT);
                free_matrix(fulldata_curl,1,ntr_glob,1,NT);
            }
            if (WAVETYPE==2 || WAVETYPE==3) {
                free_matrix(fulldata_vz,1,ntr_glob,1,NT);
            }
            break;
        case 5 : /* everything except curl and div */
            if (WAVETYPE==1 || WAVETYPE==3) {
                free_matrix(fulldata_vx,1,ntr_glob,1,NT);
                free_matrix(fulldata_vy,1,ntr_glob,1,NT);
                free_matrix(fulldata_p,1,ntr_glob,1,NT);
            }
            if (WAVETYPE==2 || WAVETYPE==3) {
                free_matrix(fulldata_vz,1,ntr_glob,1,NT);
            }
            break;
    }
    if ((ntr>0) && (SEISMO)){
      dealloc_sections(ntr,ns,recpos_loc,sectionvx,sectionvy,sectionvz,sectionp,sectionpnp1,sectionpn,sectioncurl,sectiondiv,
	sectionpdata,sectionpdiff,sectionpdiffold,sectionvxdata,sectionvxdiff,sectionvxdiffold,sectionvydata,
	sectionvydiff,sectionvydiffold,sectionvzdata,sectionvzdiff,sectionvzdiffold);
    }
    
    
    if((EPRECOND==1)||(EPRECOND==3)){
        if(WAVETYPE==1 || WAVETYPE==3) {
            free_matrix(We_sum,-nd+1,NY+nd,-nd+1,NX+nd);
            free_matrix(Ws,-nd+1,NY+nd,-nd+1,NX+nd);
            free_matrix(Wr,-nd+1,NY+nd,-nd+1,NX+nd);
            free_matrix(We,-nd+1,NY+nd,-nd+1,NX+nd);
        }
        if(WAVETYPE==2 || WAVETYPE==3) {
            free_matrix(We_sum_SH,-nd+1,NY+nd,-nd+1,NX+nd);
            free_matrix(Ws_SH,-nd+1,NY+nd,-nd+1,NX+nd);
            free_matrix(Wr_SH,-nd+1,NY+nd,-nd+1,NX+nd);
            free_matrix(We_SH,-nd+1,NY+nd,-nd+1,NX+nd);
        }
    }

    free_matrix(sectionread,1,ntr_glob,1,ns);
    
    
    if((INV_STF==1)||(TIME_FILT==1) || (TIME_FILT==2)) {
        /* free memory for inversion of source time function */
        if (WAVETYPE==1 || WAVETYPE==3) {
            free_matrix(sectionvy_conv,1,ntr_glob,1,NT);
            free_matrix(sectionvx_conv,1,ntr_glob,1,NT);
            free_matrix(sectionp_conv,1,ntr_glob,1,NT);
            free_matrix(sectionvy_obs,1,ntr_glob,1,NT);
            free_matrix(sectionvx_obs,1,ntr_glob,1,NT);
            free_matrix(sectionp_obs,1,ntr_glob,1,NT);
        }
        if (WAVETYPE==2 || WAVETYPE==3) {
            free_matrix(sectionvz_conv,1,ntr_glob,1,NT);
            free_matrix(sectionvz_obs,1,ntr_glob,1,NT);
        }
        free_vector(source_time_function,1,NT);
    }
    
    if(WOLFE_CONDITION){
        free_matrix(waveconv_old,-nd+1,NY+nd,-nd+1,NX+nd);
        if(!ACOUSTIC) free_matrix(waveconv_u_old,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(waveconv_rho_old,-nd+1,NY+nd,-nd+1,NX+nd);
        
        free_matrix(waveconv_up,-nd+1,NY+nd,-nd+1,NX+nd);
        if(!ACOUSTIC) free_matrix(waveconv_u_up,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(waveconv_rho_up,-nd+1,NY+nd,-nd+1,NX+nd);
    }
    
    if((WAVETYPE==2 || WAVETYPE==3) && (FORWARD_ONLY==0)){
        free_f3tensor(forward_prop_rho_z,-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
        free_f3tensor(forward_prop_z_xz,-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
        free_f3tensor(forward_prop_z_yz,-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
        free_matrix(waveconv_rho_shot_z,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(waveconv_u_shot_z,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(waveconv_mu_z,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(waveconv_rho_s_z,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(waveconv_u_z,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(waveconv_rho_z,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(gradp_u_z,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(gradp_rho_z,-nd+1,NY+nd,-nd+1,NX+nd);
        if(L){
            free_f3tensor(pt,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
            free_f3tensor(po,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        }
    }
    
    
    /* free memory for source position definition */
    free_matrix(srcpos1,1,8,1,1);
    
    free_vector(L2t,1,4);
    free_vector(epst1,1,3);
    free_vector(epst2,1,3);
    free_vector(epst3,1,3);
    free_vector(picked_times,1,ntr);
    
    free_vector(hc,0,6);
    
    /* free memory for global source positions */
    free_matrix(srcpos,1,8,1,nsrc);
    

    
    if(TIME_FILT==2){
        free_vector(F_LOW_PASS_EXT,1,nfrq);
    }
    
    if (nsrc_loc>0){
        free_matrix(signals,1,nsrc_loc,1,NT);
        free_matrix(srcpos_loc,1,8,1,nsrc_loc);
        free_matrix(srcpos_loc_back,1,6,1,nsrc_loc);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (MYID==0){
        if(VERBOSE){
            fprintf(FP,"\n **Info from main (written by PE %d): \n",MYID);
            fprintf(FP," CPU time of program per PE: %li seconds.\n",clock()/CLOCKS_PER_SEC);
            time8=MPI_Wtime();
            fprintf(FP," Total real time of program: %4.2f seconds.\n",time8-time1);
            time_av_v_update=time_av_v_update/(double)NT;
            time_av_s_update=time_av_s_update/(double)NT;
            time_av_v_exchange=time_av_v_exchange/(double)NT;
            time_av_s_exchange=time_av_s_exchange/(double)NT;
            time_av_timestep=time_av_timestep/(double)NT;
            fprintf(FP," Average times for \n");
            fprintf(FP," velocity update:  \t %5.3f seconds  \n",time_av_v_update);
            fprintf(FP," stress update:  \t %5.3f seconds  \n",time_av_s_update);
            fprintf(FP," velocity exchange:  \t %5.3f seconds  \n",time_av_v_exchange);
            fprintf(FP," stress exchange:  \t %5.3f seconds  \n",time_av_s_exchange);
            fprintf(FP," timestep:  \t %5.3f seconds  \n",time_av_timestep);
        }
        time8=MPI_Wtime();
        fprintf(FP," Total real time of program: %4.2f seconds.\n",time8-time1);
        if(FORWARD_ONLY==0) {
            printf("\n Inversion finished after %d iterations. \n\n",iter);
        } else {
            printf("\n Forward calculation finished. \n\n");
        }
    }
    
    /* de-allocate buffer for messages */
    MPI_Buffer_detach(buff_addr,&buffsize);
    
    fclose(FP);
    
    MPI_Finalize();
    return 0;
    
}/*main*/
