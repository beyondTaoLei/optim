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
 *   Write 2D snapshot for current timestep  to file                                   
 *
 *  See COPYING file for copying and redistribution conditions.
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void snap(FILE *fp,int nt, int nsnap, float **vx, float **vy, float **sxx,
          float **syy, float **sp, float **u, float **pi, float *hc, int ishot, st_wavefildfreq *wfd_f){

    /* 
    different data formats of output:
    SNAP_FORMAT=1  :  SU (IEEE)
    SNAP_FORMAT=2  :  ASCII
    SNAP_FORMAT=3  :  BINARY (IEEE)

    different types:
    SNAP=1 : values in vx and vy
    SNAP=2 : -(vx+vy) (pressure field)
    SNAP=3 : divergence of vx and vy (energy of compressional waves)
    and curl of vx and vy (energy of shear waves)
    SNAP=4 : both particle velocities (type=1) and energy (type=3)
    */


    int i,j, m, fdoh, nd;
    float amp, vyx, vxy, vxx, vyy, dhi;
    float **divfield, **curlfield;
    char snapfile_x[STRING_SIZE], snapfile_y[STRING_SIZE], snapfile_div[STRING_SIZE];
    char snapfile_rot[STRING_SIZE], snapfile_p[STRING_SIZE], ext[8], wm[2];
    FILE *fpx1, *fpy1, *fpx2, *fpy2;

    extern float DH, DT;
    extern char SNAP_FILE[STRING_SIZE];
    extern int NX, NY,  SNAP_FORMAT, SNAP, FDORDER;
    extern int MYID, POS[3], IDX, IDY,VERBOSE;
    extern int INV_RHO;

    dhi = 1.0/DH;
    fdoh = FDORDER/2;

    switch(SNAP_FORMAT){
    case 1:
        sprintf(ext,".su");
        break;
    case 2:
        sprintf(ext,".asc");
        break;
    case 3:
        sprintf(ext,".bin");
        break;
    }

    sprintf(snapfile_x,"%s%s.x.shot%i.%i.%i",SNAP_FILE,ext,ishot,POS[1],POS[2]);
    sprintf(snapfile_y,"%s%s.y.shot%i.%i.%i",SNAP_FILE,ext,ishot,POS[1],POS[2]);
    sprintf(snapfile_div,"%s%s.div.shot%i.%i.%i",SNAP_FILE,ext,ishot,POS[1],POS[2]);
    sprintf(snapfile_rot,"%s%s.rot.shot%i.%i.%i",SNAP_FILE,ext,ishot,POS[1],POS[2]);
    sprintf(snapfile_p,"%s%s.p.shot%i.%i.%i",SNAP_FILE,ext,ishot,POS[1],POS[2]);

    fprintf(fp,"\n\n PE %d is writing snapshot-data at T=%fs to \n",MYID,nt*DT);
    //printf("nsnap=%d\n",nsnap);
    if (nsnap==1)
        sprintf(wm,"w");
    else 
        sprintf(wm,"a");

    switch(SNAP){
    case 1 :
        fprintf(fp,"%s\n", snapfile_x);
        fprintf(fp,"%s\n\n", snapfile_y);
        fpx1=fopen(snapfile_x,wm);
        fpy1=fopen(snapfile_y,wm);
        for (i=1;i<=NX;i+=IDX)
        for (j=1;j<=NY;j+=IDY){
            writedsk(fpx1,vx[j][i],SNAP_FORMAT);
            writedsk(fpy1,vy[j][i],SNAP_FORMAT);
        }
        fclose(fpx1);
        fclose(fpy1);
        break;


    case 2 :
        fprintf(fp,"%s\n\n",snapfile_p);
        fpx1=fopen(snapfile_p,wm);

        for (i=1;i<=NX;i+=IDX)
        for (j=1;j<=NY;j+=IDY){
//             if(!ACOUSTIC)
//                 amp=-sxx[j][i]-syy[j][i];
//             else
            amp=-sp[j][i];
            writedsk(fpx1,amp,SNAP_FORMAT);
        }
        fclose(fpx1);
        break;

    case 4 :
        fprintf(fp,"%s\n", snapfile_x);
        fprintf(fp,"%s\n", snapfile_y);
        fpx1=fopen(snapfile_x,wm);
        fpy1=fopen(snapfile_y,wm);
        for (i=1;i<=NX;i+=IDX)
        for (j=1;j<=NY;j+=IDY){
            writedsk(fpx1,vx[j][i],SNAP_FORMAT);
            writedsk(fpy1,vy[j][i],SNAP_FORMAT);
        }
        fclose(fpx1);
        fclose(fpy1);

        fprintf(fp,"%s\n\n",snapfile_p);
        fpx1=fopen(snapfile_p,wm);

        for (i=1;i<=NX;i+=IDX)
        for (j=1;j<=NY;j+=IDY){
//             if(!ACOUSTIC)
//                 amp=-sxx[j][i]-syy[j][i];
//             else
            amp=-sp[j][i];
            writedsk(fpx1,amp,SNAP_FORMAT);
        }
        fclose(fpx1);
    case 3 :
        /* output of the curl of the velocity field according to Dougherty and
        Stephen (PAGEOPH, 1988) */
        /*if (NY1<=2) error("NY1 must be greater than 2.");*/
        fprintf(fp,"%s\n", snapfile_div);
        fprintf(fp,"%s\n\n", snapfile_rot);
        fpx2=fopen(snapfile_div,wm);
        fpy2=fopen(snapfile_rot,wm);

        nd = FDORDER/2+1;
        curlfield  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);


        for (i=1;i<=NX;i+=IDX)
        for (j=1;j<=NY;j+=IDY){

            /* spatial derivatives using Holberg coefficients */
            vyx = 0;
            vxy = 0;
            for (m=1; m<=fdoh; m++) {
                vyx += hc[m]*(vy[j][i+m] - vy[j][i-m+1]);
                vxy += hc[m]*(vx[j+m][i] - vx[j-m+1][i]);
            }
            vyx *= dhi;
            vxy *= dhi;

//             if(!ACOUSTIC)
//                 curlfield[j][i]=(vxy-vyx)*sqrt(u[j][i]);
//             else
            curlfield[j][i]=0;
        }


        for (i=1;i<=NX;i+=IDX)
        for (j=1;j<=NY;j+=IDY){
            writedsk(fpy2,curlfield[j][i],SNAP_FORMAT);
        }
        free_matrix(curlfield,-nd+1,NY+nd,-nd+1,NX+nd);

        /* output of the divergence of the velocity field according to Dougherty and
        Stephen (PAGEOPH, 1988) */
        divfield  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        for (i=1;i<=NX;i+=IDX)
        for (j=1;j<=NY;j+=IDY){

            /* spatial derivatives using Holberg coefficients */
            vxx = 0;
            vyy = 0;
            for (m=1; m<=fdoh; m++) {
                vxx += hc[m]*(vx[j][i+m-1] - vx[j][i-m]  );
                vyy += hc[m]*(vy[j+m-1][i] - vy[j-m][i]  );
            }
            vxx *= dhi;
            vyy *= dhi;

            divfield[j][i]=(vxx+vyy)*sqrt(pi[j][i]);
        }


        for (i=1;i<=NX;i+=IDX)
        for (j=1;j<=NY;j+=IDY){
            writedsk(fpx2,divfield[j][i],SNAP_FORMAT);
        }

        free_matrix(divfield,-nd+1,NY+nd,-nd+1,NX+nd);
        fclose(fpx2);
        fclose(fpy2);
        break;
    case 5://for fwi
        fprintf(fp,"%s\n\n",snapfile_p);
        fpx1=fopen(snapfile_p,wm);

        for (i=1;i<=NX;i+=IDX){
        for (j=1;j<=NY;j+=IDY){
            amp=sp[j][i];
            writedsk(fpx1,amp,SNAP_FORMAT);
        }}
        fclose(fpx1);
        break;
    case 6://for gradient calculation in frequency domain
        /* transform wavefields from time domian into frequency domain */
        if(INV_RHO==1){
            dft_forward2(nt, wfd_f->num_freqs, wfd_f->freqs, vx, wfd_f->pvx_f);
            dft_forward2(nt, wfd_f->num_freqs, wfd_f->freqs, vy, wfd_f->pvy_f);
        }
        dft_forward2(nt, wfd_f->num_freqs, wfd_f->freqs, sp, wfd_f->psp_f);
        break;
        
    }

}

void snap_last(FILE *fp, float **vx, float **vy,float **sp, int ishot){

    /* 
    different data formats of output:
    SNAP_FORMAT=1  :  SU (IEEE)
    SNAP_FORMAT=2  :  ASCII
    SNAP_FORMAT=3  :  BINARY (IEEE)

    different types:
    SNAP=1 : values in vx and vy
    SNAP=2 : -(vx+vy) (pressure field)
    SNAP=3 : divergence of vx and vy (energy of compressional waves)
    and curl of vx and vy (energy of shear waves)
    SNAP=4 : both particle velocities (type=1) and energy (type=3)
    */


    int i,j;
    char snapfile_x[STRING_SIZE], snapfile_y[STRING_SIZE];
    char snapfile_p[STRING_SIZE], ext[8], wm[2];
    FILE *fpx1, *fpy1, *fpx2, *fpy2;

    extern float DH, DT;
    extern char SNAP_FILE[STRING_SIZE];
    extern int NX, NY,  SNAP_FORMAT, SNAP, FDORDER;
    extern int MYID, POS[3];


    sprintf(ext,".bin");
    sprintf(snapfile_x,"snapl%s.x.shot%i.%i.%i",ext,ishot,POS[1],POS[2]);
    sprintf(snapfile_y,"snapl%s.y.shot%i.%i.%i",ext,ishot,POS[1],POS[2]);
    sprintf(snapfile_p,"snapl%s.p.shot%i.%i.%i",ext,ishot,POS[1],POS[2]);

    fprintf(fp,"\n\n PE %d is writing snapshot-data at T=nt*DT s to \n",MYID);


    //save the last snapshot for boundary saving strategy
    fprintf(fp,"%s\n", snapfile_x);
    fprintf(fp,"%s\n\n", snapfile_y);
    fprintf(fp,"%s\n\n", snapfile_p);
    fpx1=fopen(snapfile_x,"w");
    fpy1=fopen(snapfile_y,"w");
    fpx2=fopen(snapfile_p,"w");
    for (i=1;i<=NX;i++)
    for (j=1;j<=NY;j++){
        writedsk(fpx1,vx[j][i],3);
        writedsk(fpy1,vy[j][i],3);
        writedsk(fpx2,sp[j][i],3);
    }
    fclose(fpx1);
    fclose(fpy1);
    fclose(fpx2);
    
}