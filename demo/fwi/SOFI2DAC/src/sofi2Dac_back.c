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

int main(int argc, char **argv){
    /* variables in main */
    int ns, nseismograms=0, nt, nd, fdo3, j, i, infoout;
    int lsnap, incsnap, jsnap, nsnap=0, lsamp=0, buffsize;
    int ntr=0, ntr_loc=0, ntr_glob=0, nsrc=0, nsrc_glob=0, ishot, nshots=0;
    /*Limits for local grids defined in subgrid_bounds.c */
    int * gx=NULL, * gy=NULL;
    
    /* General parameters */
    int nt_out;
    double time1, time2, time3, time4, time5, time6, time7, time8,
    time_av_v_update=0.0, time_av_s_update=0.0, time_av_v_exchange=0.0,
    time_av_s_exchange=0.0, time_av_timestep=0.0;
    float memdyn, memmodel, memseismograms, membuffer, memtotal;
    float fac1, fac2;
    char *buff_addr, ext[10], *fileinp, fsrc[STRING_SIZE];
    
    // Pointer for dynamic wavefields:
    float  **  psxx, **psyy, **psp, **u;
    float  **  pvx, **  pvy, **  uvx, **  uvy;
    float  **  pvxp1, **  pvyp1, **  pvxm1, **  pvym1;
    float  **  prho, **prip=NULL, **prjp=NULL, **ppi, **pu;
    float ** bufferlef_to_rig,  ** bufferrig_to_lef, ** buffertop_to_bot, ** bufferbot_to_top;
    
    /* PML variables */
    float * d_x, * K_x, * alpha_prime_x, * a_x, * b_x, * d_x_half, * K_x_half, * alpha_prime_x_half, * a_x_half, * b_x_half, * d_y, * K_y, * alpha_prime_y, * a_y, * b_y, * d_y_half, * K_y_half, * alpha_prime_y_half, * a_y_half, * b_y_half;
    float ** psi_sxx_x, ** psi_syy_y, ** psi_sxy_y, ** psi_sxy_x, ** psi_vxx, ** psi_vyy, ** psi_vxy, ** psi_vyx, ** psi_vxxs;
    
    /* Variables for viscoelastic modeling */
    float **ptaus=NULL, **ptaup=NULL, *etajm=NULL, *peta=NULL, *bjm=NULL;
    float *cjm=NULL, ***e=NULL, ***pp=NULL, **g=NULL;
    
    float  ** sectionvx=NULL, ** sectionvy=NULL, ** sectionp=NULL, ** sectioncurl=NULL, ** sectiondiv=NULL;
    float ** seismo_fulldata=NULL;
    float  **srcpos = NULL, *hc=NULL;
    int    ** recpos=NULL, ** recpos_loc=NULL;
    int    *  recswitch=NULL;
    
    FILE *fpinp;
    
    MPI_Request *req_send, *req_rec;
    MPI_Status  *send_statuses, *rec_statuses;
    
    /*******************************************************/
    /* modified by Tao, in order to insert different source 
    * for FWI, RTM or some other cases */
    /* reads source information, including positon, signals and so on */
    /* source and boundary structure */
    int ntr_bdr=0, ntr_glob_bdr=0;
    float ** sectionvx_bdr=NULL, ** sectionvy_bdr=NULL, ** sectionp_bdr=NULL;
    float ** seismo_fulldata_bdr=NULL;
    int ** recpos_bdr=NULL, ** recpos_loc_bdr=NULL;
    int *  recswitch_bdr=NULL;
    st_acquisition  *acq =NULL, *acqx =NULL, *acqy =NULL;
    st_wavefildfreq   *wfd_f=NULL;
    /* flag = 1.0 : forward-propagation in clockwise direction;
       flag =-1.0 : backward-propagation in anticlockwise direction. */
    float flag=-1.0;
    float **hc_all = NULL, **eng = NULL;
    /*******************************************************/
    
    /* Initialize MPI environment */
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&NP);
    MPI_Comm_rank(MPI_COMM_WORLD,&MYID);
    
    acq = (st_acquisition *) malloc(sizeof(st_acquisition));
    acqx = (st_acquisition *) malloc(sizeof(st_acquisition));
    acqy = (st_acquisition *) malloc(sizeof(st_acquisition));
    wfd_f = (st_wavefildfreq *) malloc(sizeof(st_wavefildfreq));
    
    setvbuf(stdout, NULL, _IONBF, 0);
    if (MYID == 0){
        time1=MPI_Wtime();
        clock();
    }
    
    /* print program name, version etc to stdout*/
    if (MYID == 0) info(stdout);
    
    /* read parameters from parameter-file (stdin) */
    fileinp=argv[1];
    fpinp=fopen(fileinp,"r");
    if(fpinp==NULL) {
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
    
    /* open log-file (each PE is using different file) */
    /*	fp=stdout; */
    sprintf(ext,".%i",MYID);
    strcat(LOG_FILE,ext);
    
    /* If Verbose==0, no PE will write a log file */
//     if(!VERBOSE) sprintf(LOG_FILE,"/dev/null");
    
    if ((MYID==0)) FP=stdout;
    else FP=fopen(LOG_FILE,"w");
    fprintf(FP," This is the log-file generated by PE %d \n\n",MYID);
    
    if (MYID == 0) note(stdout);
    
    /* domain decomposition */
    initproc();
    
    NT=iround(TIME/DT);  	  /* number of timesteps */
    ns=iround(NT/NDT);           /* number of samples per trace */
    lsnap = iround ( TSNAP1/DT ); /* first snapshot at this timestep */
    incsnap = iround ( TSNAPINC/DT );
    jsnap  = iround ( TSNAP2/DT );
    lsamp = NDT;
    
    /* output of parameters to log-file or stdout */
    if (MYID==0) write_par(FP);
    
    /* NXG, NYG denote size of the entire (global) grid */
    NXG=NX;
    NYG=NY;
    
    /* In the following, NX and NY denote size of the local grid ! */
    NX = IENDX;
    NY = IENDY;
    calc_grid_for_inv();
    
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
    
    nd = FDORDER/2;
    fdo3 = 2 * nd;
    
    /*allocate memory for dynamic, static and buffer arrays */
    fac1=(NX+FDORDER)*(NY+FDORDER);
    fac2=sizeof(float)*pow(2.0,-20.0);
    
    if (L){
        memdyn=(5.0+3.0*(float)L)*fac1*fac2;
        memmodel=(12.0+3.0*(float)L)*fac1*fac2;
        
    } else {
        memdyn=5.0*fac1*fac2;
        memmodel=6.0*fac1*fac2;
    }
    memseismograms=nseismograms*ntr*ns*fac2;
    
    membuffer=2.0*fdo3*(NY+NX)*fac2;
    buffsize=2.0*2.0*fdo3*(NX+NY)*sizeof(MPI_FLOAT);
    memtotal=memdyn+memmodel+memseismograms+membuffer+(buffsize*pow(2.0,-20.0));
    
    if (MYID==0){
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
//     send_statuses=(MPI_Status *)malloc(REQUEST_COUNT*sizeof(MPI_Status));
//     rec_statuses=(MPI_Status *)malloc(REQUEST_COUNT*sizeof(MPI_Status));
    
    /* ------------ memory allocation for arrays ------------- */
    /* subgrid arrays*/
    gy = ivector ( 1,4 );
    gx = ivector ( 1,4 );
    
    /* memory allocation for dynamic (wavefield) arrays */
    psp  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    
    pvx  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    pvy  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    uvx  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    uvy  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    eng  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    pvxp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    pvyp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    pvxm1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    pvym1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    
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
    
    /* static (model) arrays (viscoelastic) */
    if (L) {
        /* dynamic (wavefield) arrays for viscoelastic modeling */
        pp = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        /* memory allocation for static arrays for viscoelastic modeling */
        e =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        ptaus =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        ptaup =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        g =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        peta =  vector(1,L);
        etajm =  vector(1,L);
        bjm =  vector(1,L);
        cjm =  vector(1,L);
    }
    
    /* memory allocation for static (model) arrays */
    prho =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    prip =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    prjp =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    ppi  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    
    /* memory allocation for buffer arrays in which the wavefield
     information which is exchanged between neighbouring PEs is stored */
    bufferlef_to_rig = matrix(1,NY,1,fdo3);
    bufferrig_to_lef = matrix(1,NY,1,fdo3);
    buffertop_to_bot = matrix(1,NX,1,fdo3);
    bufferbot_to_top = matrix(1,NX,1,fdo3);
    
    fprintf(FP," ... memory allocation for PE %d was successfull.\n\n", MYID);
    
    /* Holberg coefficients for FD operators*/
    hc = holbergcoeff();
    hc_all = holbergcoeff_all();
    
    MPI_Barrier(MPI_COMM_WORLD);
    /* create model grids */
    if(L){
        readmod_viscac(prho,ppi,ptaup,peta);
    }else{
        readmod_acoustic(prho,ppi);
    }
    /* For the calculation of the material parameters between gridpoints
        they have to be averaged. For this, values lying at 0 and NX+1,
        for example, are required on the local grid. These are now copied from the
        neighbouring grids */
    if (L){
        matcopy_viscac(prho,ppi,ptaup);
    }else{
        matcopy_acoustic(prho, ppi);
    }
    av_rho(prho,prip,prjp);
    
    /* check if the FD run will be stable and free of numerical dispersion */
    checkfd(FP, prho, ppi, ptaup, peta, hc);
    
    /* calculate damping coefficients for CPMLs*/
    if(FW>0)
        PML_pro(d_x, K_x, alpha_prime_x, a_x, b_x, d_x_half, K_x_half, alpha_prime_x_half, a_x_half, b_x_half, d_y, K_y, alpha_prime_y, a_y, b_y, d_y_half, K_y_half, alpha_prime_y_half, a_y_half, b_y_half);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    /* Preparing memory variables for update_s (viscoelastic) */
    if (L) prepare_update_p(etajm,peta,ppi,prho,ptaup,g,bjm,cjm,e);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    time2 = MPI_Wtime();
    fprintf ( FP, "\n\n\n **************************************************\n" );
    fprintf ( FP, " *********** STARTING TIME STEPPING ***************\n" );
    fprintf ( FP, " **************************************************\n\n" );
    if(MYID==0) fprintf(FP, " real time before starting time loop: %4.2f s.\n", time2-time1);
    
    /*******************************************************/
    /*******************************************************/
    /*wavefield in freq. domain*/
    if(SNAP==6) allocate_mem_wfdf(1, wfd_f);
    ishot=1;
    NSHOT1=1;
    fprintf (FP, "\n==================================================================================\n" );
    fprintf (FP, "   MYID=%d *****  Starting simulation for shot %d of 1  ********** \n", MYID, ishot);
    fprintf (FP,   "==================================================================================\n\n" );
    MPI_Barrier ( MPI_COMM_WORLD );
    
    /*******************************************************/
    /*******************************************************/
    /* reads source information, including positon, signals and so on */
    /* Source part */
    sprintf(fsrc,"bdr_p_shot%i.su",ishot);
    extract_current_sources(fsrc, acq);
    sprintf(fsrc,"bdr_vx_shot%i.su",ishot);
    extract_current_sources(fsrc, acqx);
    sprintf(fsrc,"bdr_vy_shot%i.su",ishot);
    extract_current_sources(fsrc, acqy);
    /* Receiver part */
    if(SEISMO){
        recpos = receiver ( FP, &ntr_glob );
        recswitch = ivector ( 1, ntr_glob );
        /* allocate buffer for seismogram output, merged seismogram section of all PEs */
        seismo_fulldata=matrix ( 1,ntr_glob,1,ns );
        recpos_loc = splitrec ( recpos, &ntr, ntr_glob, recswitch );
        if(ntr>0) allocate_mem_seis(1, &sectionvx, &sectionvy, &sectionp, &sectiondiv, ntr, ns);
    }
    /* boundary saving strategy part */
    SEISMO2=0;
    if(SEISMO2){
        recpos_bdr = receiver_boundary ( FP, &ntr_glob_bdr );
        recswitch_bdr = ivector ( 1, ntr_glob_bdr );
        /* allocate buffer for seismogram output, merged seismogram section of all PEs */
        seismo_fulldata_bdr=matrix ( 1,ntr_glob_bdr,1,ns );
        recpos_loc_bdr = splitrec ( recpos_bdr, &ntr_bdr, ntr_glob_bdr, recswitch_bdr );
        if(ntr_bdr>0) allocate_mem_seis2(1, &sectionvx_bdr, &sectionvy_bdr, &sectionp_bdr, ntr_bdr, ns);
    }
    
    if(EPRECOND) M_init(-nd+1, NX+nd, -nd+1, NY+nd, eng);
    /*******************************************************/
    /*******************************************************/
    /* initialize wavefield with zero */
    if (L){
        zero_fdveps_viscac(-nd+1, NY+nd, -nd+1, NX+nd, pvx, pvy, psp, pvxp1, pvyp1, psi_sxx_x, psi_sxy_x, psi_vxx, psi_vyx, psi_syy_y, psi_sxy_y, psi_vyy, psi_vxy, psi_vxxs, pp); 
    }else{
        zero_fdveps_ac(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,psp,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x, psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y, psi_vyy,psi_vxy,psi_vxxs);
    }
    
    /* Reseting lsmap to NDT for saving seismograms  */
    lsamp = NDT;
    
    subgrid_bounds ( 1, NX, 1, NY, gx, gy );
    /*---------------------------------------------------------------*/
    /*----------------------  loop over timesteps  ------------------*/
    /*---------------------------------------------------------------*/
    for (nt=NT;nt>=1;nt--){
        
        /* before update particle velocity */
        if(MYID == 0&&nt%10==0) time3 = MPI_Wtime();
        /* Check if simulation is still stable P and SV */
        if (isnan(pvy[NY/2][NX/2])) {
            fprintf(FP,"\n Time step: %d; pvy: %f \n",nt,pvy[NY/2][NX/2]);
            declare_error(" Simulation is unstable !");
        }
        
        /* update stress */
        if (L) {   /* viscoelastic */
                update_p_visc_PML(1, NX, 1, NY, pvx, pvy, psp, ppi, prho, hc, infoout, pp, g, bjm, cjm, e, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);
        } else {   /* elastic */
                update_p_interior(gx, gy, prho, ppi, hc, pvx, pvy, psp, u, flag);
                //update_p_PML(gx, gy, prho, ppi, hc, pvx, pvy, psp, u, psi_vxx, psi_vyy, K_x, a_x, b_x, K_y, a_y, b_y, NX, NY, flag);
        }
        
        /* explosive source */
        source_boundary(psp, acq->srcpos_loc, acq->signals_loc, nt, acq->nsrc_loc);
        if(nt==NT) source_last_snap(ishot, psp, 3);
        
        /* Applying free surface condition */
        if ((FREE_SURF) && (POS[2]==0)) surface_acoustic_PML(1, psp);
        if(MYID == 0&&nt%10==0) time_av_s_update+=time_count(time5,&time6);
        
        /* stress exchange between PEs */
        exchange_p (nd, psp, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, req_send, req_rec );
        if(MYID == 0&&nt%10==0) time_av_s_exchange+=time_count(time6,&time7);
        
        /* update of particle velocities */
        update_v_interior(gx, gy, prip, prjp, hc, pvx, pvy, uvx, uvy, psp, flag);
        //update_v_PML(gx, gy, prip, prjp, hc, pvx, pvy, uvx, uvy, psp, psi_sxx_x, psi_syy_y, K_x_half, a_x_half, b_x_half, 
                    //K_y_half, a_y_half, b_y_half, NX, NY, flag);
        if(MYID == 0&&nt%10==0) time_av_v_update+=time_count(time3,&time4);
        
        source_boundary(pvx, acqx->srcpos_loc, acqx->signals_loc, nt, acqx->nsrc_loc);
        source_boundary(pvy, acqy->srcpos_loc, acqy->signals_loc, nt, acqy->nsrc_loc);
        if(nt==NT){
            source_last_snap(ishot, pvx, 1);
            source_last_snap(ishot, pvy, 2);
        }
        
        /* exchange of particle velocities between PEs */
        exchange_v ( nd, pvx, pvy, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, req_send, req_rec );        
        if(MYID == 0&&nt%10==0) time_av_v_exchange+=time_count(time4,&time5);
        
        /* store amplitudes at receivers in section-arrays */
        if ( ( SEISMO ) && ( nt == lsamp ) && ( nt < NT ) ) {
            seismo_ssg ( lsamp, ntr, recpos_loc, sectionvx, sectionvy, sectionp, sectiondiv, pvx, pvy, psp, ppi, hc );
            lsamp += NDT;//!!!no output here
        }
        //seismo_ssg ( nt, ntr, recpos_loc, sectionvx, sectionvy, sectionp, sectiondiv, pvx, pvy, psp, ppi, hc );
        if ( SEISMO2){
            seismo_ssg2(nt, ntr_bdr, recpos_loc_bdr, sectionvx_bdr,sectionvy_bdr, sectionp_bdr,pvx, pvy, psp);
            if (nt==NT) snap_last(FP, pvx, pvy,psp, ishot);
        }
        
        /*******************************************************/
        /* WRITE SNAPSHOTS TO DISK */
        if (SNAP){
            if ((nt-lsnap)%incsnap==0 &&(nt>=lsnap&&nt<=jsnap)) {
                //nsnap=(nt-lsnap)/incsnap+1;
                snap(FP, nt, ++nsnap, pvx, pvy, psxx,psyy, psp, pu, ppi, hc, ishot, wfd_f);
            }
        }
        if(EPRECOND) eprecond(eng, pvx, pvy);
        /*******************************************************/
        /* update in one step  */
        if(MYID == 0&&nt%10==0) time_av_timestep+=time_count(time3,&time8);
        
    }
    /*------------------------------------------------------------------------------*/
    /*--------------------  End  of loop over timesteps (forward model) ------------*/
    /*------------------------------------------------------------------------------*/
    
    fprintf ( FP, "\n\n *********** Finish TIME STEPPING ****************\n" );
    fprintf ( FP, " **************************************************\n\n" );
    
    /* write seismograms to file(s) */
    if ( SEISMO ) {
        save_seis(FP, sectionvx, sectionvy, sectionp, sectiondiv, seismo_fulldata, 
                    recswitch, recpos, recpos_loc, ntr_glob, ns, ishot, acq->Sx, acq->Sy);
        fprintf ( FP, "\n\n" );
    }
    if ( SEISMO2 ) {
        save_seis2(FP, sectionvx_bdr, sectionvy_bdr, sectionp_bdr, seismo_fulldata_bdr, 
                    recswitch_bdr, recpos_bdr, recpos_loc_bdr, ntr_glob_bdr, ns, ishot, acq->Sx, acq->Sy);
        fprintf ( FP, "\n\n" );
    }

    /*******************************************************/
    if(SNAP==6) save_wavefields_freq(FP, wfd_f, ishot);
    if(EPRECOND){
        sprintf(fsrc,"eng_shot%i.bin",ishot);
        savematIDXY(eng, fsrc);
    }
    /*******************************************************/
                
    /* ====================================== */
    /* ====== deallocation of memory =========*/
    /* ====================================== */
    /* Source part*/
    if(acq->nsrc_loc>0){
        free_imatrix(acq->srcpos_loc, 1,3,1,acq->nsrc_loc);
        free_matrix(acq->signals_loc, 1,acq->nsrc_loc,1,NT);
        free_imatrix(acqx->srcpos_loc, 1,3,1,acqx->nsrc_loc);
        free_matrix(acqx->signals_loc, 1,acqx->nsrc_loc,1,NT);
        free_imatrix(acqy->srcpos_loc, 1,3,1,acqy->nsrc_loc);
        free_matrix(acqy->signals_loc, 1,acqy->nsrc_loc,1,NT);
    }
    /* Receiver part */
    if(SEISMO){
        /* Receiver part */
        free_imatrix ( recpos, 1, 3, 1, ntr_glob );
        free_ivector(recswitch, 1, ntr_glob );
        free_matrix ( seismo_fulldata,1,ntr_glob,1,ns );
        if(ntr>0){
            free_imatrix ( recpos_loc, 1, 3, 1, ntr );
            allocate_mem_seis(0, &sectionvx, &sectionvy, &sectionp, &sectiondiv,ntr, ns);
        }
    }
    if(SEISMO2){
        /* Receiver part */
        free_imatrix ( recpos_bdr, 1, 3, 1, ntr_glob_bdr );
        free_ivector(recswitch_bdr, 1, ntr_glob_bdr );
        free_matrix ( seismo_fulldata_bdr,1,ntr_glob_bdr,1,ns );
        if(ntr_bdr>0){
            free_imatrix ( recpos_loc_bdr, 1, 3, 1, ntr_bdr );
            allocate_mem_seis2(0, &sectionvx_bdr, &sectionvy_bdr, &sectionp_bdr, ntr_bdr, ns);
        }
    }
    /*wavefield in freq. domain*/
    if(SNAP==6) allocate_mem_wfdf(0, wfd_f);
    /*******************************************************/
    
    free_matrix(psp,-nd+1,NY+nd,-nd+1,NX+nd);
    free_matrix(pvx,-nd+1,NY+nd,-nd+1,NX+nd);
    free_matrix(pvy,-nd+1,NY+nd,-nd+1,NX+nd);
    free_matrix(u,-nd+1,NY+nd,-nd+1,NX+nd);
    free_matrix(uvx,-nd+1,NY+nd,-nd+1,NX+nd);
    free_matrix(uvy,-nd+1,NY+nd,-nd+1,NX+nd);
    free_matrix(eng,-nd+1,NY+nd,-nd+1,NX+nd);
    free_matrix(pvxp1,-nd+1,NY+nd,-nd+1,NX+nd);
    free_matrix(pvyp1,-nd+1,NY+nd,-nd+1,NX+nd);
    free_matrix(pvxm1,-nd+1,NY+nd,-nd+1,NX+nd);
    free_matrix(pvym1,-nd+1,NY+nd,-nd+1,NX+nd);
    
    free_matrix(prho,-nd+1,NY+nd,-nd+1,NX+nd);
    free_matrix(prip,-nd+1,NY+nd,-nd+1,NX+nd);
    free_matrix(prjp,-nd+1,NY+nd,-nd+1,NX+nd);
    free_matrix(ppi,-nd+1,NY+nd,-nd+1,NX+nd);
    
    free_ivector ( gx,1,4 );
    free_ivector ( gy,1,4 );
    
    /* free memory for viscoelastic modeling variables */
    if (L) {
        free_f3tensor(pp,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        free_f3tensor(e,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        free_matrix(ptaus,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(ptaup,-nd+1,NY+nd,-nd+1,NX+nd);
        free_matrix(g,-nd+1,NY+nd,-nd+1,NX+nd);
        free_vector(peta,1,L);
        free_vector(etajm,1,L);
        free_vector(bjm,1,L);
        free_vector(cjm,1,L);
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
        free_matrix(psi_sxx_x,1,NY,1,2*FW);
        free_matrix(psi_syy_y,1,2*FW,1,NX);
        free_matrix(psi_vxx,1,NY,1,2*FW);
        free_matrix(psi_vyy,1,2*FW,1,NX);
        free_matrix(psi_vxy,1,2*FW,1,NX);
        free_matrix(psi_vyx,1,NY,1,2*FW);
    }
    
    free_matrix(bufferlef_to_rig,1,NY,1,fdo3);
    free_matrix(bufferrig_to_lef,1,NY,1,fdo3);
    free_matrix(buffertop_to_bot,1,NX,1,fdo3);
    free_matrix(bufferbot_to_top,1,NX,1,fdo3);
    
    free_vector(hc,0,6);
    
    /* de-allocate buffer for messages */
    MPI_Buffer_detach(buff_addr,&buffsize);
    MPI_Barrier(MPI_COMM_WORLD);
    
    if ( MYID == 0 ) {
        fprintf ( FP, "\n **Info from main (written by PE %d): \n", MYID );
        
        time_av_v_update = time_av_v_update*10.0 / ( double ) NT;
        time_av_s_update = time_av_s_update*10.0 / ( double ) NT;
        time_av_v_exchange = time_av_v_exchange*10.0 / ( double ) NT;
        time_av_s_exchange = time_av_s_exchange*10.0 / ( double ) NT;
        time_av_timestep = time_av_timestep*10.0 / ( double ) NT;
        fprintf ( FP, " Average times for \n" );
        fprintf ( FP, "   velocity update:  \t %5.6f seconds  \n",
                 time_av_v_update );
        fprintf ( FP, "   stress update:  \t %5.6f seconds  \n",
                 time_av_s_update );
        fprintf ( FP, "   velocity exchange:  \t %5.6f seconds  \n",
                 time_av_v_exchange );
        fprintf ( FP, "   stress exchange:  \t %5.6f seconds  \n",
                 time_av_s_exchange );
        fprintf ( FP, "   timestep:  \t\t %5.6f seconds  \n\n", time_av_timestep );
        
        fprintf ( FP, "   time consuming percent:  \n\n");
        fprintf ( FP, "   velocity update:  \t %5.2f %%\n",
                 time_av_v_update/time_av_timestep*100.0 );
        fprintf ( FP, "   stress update:  \t %5.2f %%\n",
                 time_av_s_update/time_av_timestep*100.0  );
        fprintf ( FP, "   velocity exchange:  \t %5.2f %%\n",
                 time_av_v_exchange/time_av_timestep*100.0  );
        fprintf ( FP, "   stress exchange:  \t %5.2f %%\n",
                 time_av_s_exchange/time_av_timestep*100.0  );
        
        fprintf ( FP, " CPU time of program per PE: %li seconds.\n",
                 clock() / CLOCKS_PER_SEC );
        time8 = MPI_Wtime();
        fprintf ( FP, " Total real time of program: %4.3f seconds.\n\n",
                 time8 - time1 );
        fprintf ( FP," ******************************************************\n" );
        fprintf ( FP," **************** SOFI2D_AC has finished *****************\n" );
        fprintf ( FP," ******************************************************\n\n" );
        
    }
    
    fclose(FP);
    MPI_Finalize();
    return 0;
    
}/*main*/
