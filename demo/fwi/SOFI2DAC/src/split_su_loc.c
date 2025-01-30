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

#include "fd.h"

void split_su_pos_loc(char srcname[256], int myid, int myid1, 
                      MPI_Comm mpi_comm_no, float dh, float *refrec, 
                      int iendx, int iendy, int *POS1, int nt_in, 
                      float dt_in, int ***srcpos_loc, int *nsrc_loc){        
    /*so far this software just supports the SU-format file.*/
    /* INPUT */
    /* 
     * char srcname[256]
     * int myid, myid1, 
     * MPI_Comm mpi_comm_no
     * float dh 
     * float *refrec
     * int iendx, iendy  
     * int *POS1 
     * int nt_in
     * float dt_in
     */
    
    /* OUTPUT */
    /* 
     * int ***srcpos_loc
     * int *nsrc_loc
     */
    
    /* local variables */
    int nsrc_glob, nt, nsrc_tmp2;
    float dt;
    float **srcpos_tmp3d=NULL;
    int *repeated_tmp=NULL;
    int **srcpos;
    
    if(myid==myid1){
        read_su_pos(srcname, dh, refrec, nt_in, dt_in, &srcpos, &nsrc_glob);
    }
    MPI_Barrier(mpi_comm_no);
    MPI_Bcast(&nsrc_glob,1,MPI_INT,0,mpi_comm_no);
    if(myid!=myid1) srcpos = imatrix(1,3,1,nsrc_glob);
    MPI_Bcast(&srcpos[1][1],nsrc_glob*3,MPI_INT,0,mpi_comm_no);
    
    /***************************************************************/
    /*From the above, we just get srcpos                           */
    /***************************************************************/
    
    /* find this single source positions on subdomains */ 
    *srcpos_loc = splitpos(iendx, iendy, POS1, srcpos, nsrc_loc, nsrc_glob);
    free_imatrix(srcpos,1,3,1,nsrc_glob);
}

void split_su_sig_loc(char srcname[256], int myid, int myid1, 
                      MPI_Comm mpi_comm_no, int **srcpos_loc, 
                      int ntr_loc, float ***signals){        
/* reads the SU-format file to local distributed memory */
    /* INPUT */
    /* 
     * char srcname[256]
     * int myid, myid1, 
     * MPI_Comm mpi_comm_no
     * int ***srcpos_loc
     * int ntr_loc
     */
    
    /* OUTPUT */
    /* float ***signals
     * int *nsrc_loc
     */

    int nsrc_glob,nt;
    float dt;
    float **signals_all;
    int itr,it;
    
    /*line sources*/
    /*so far this software just supports the SU-format file.*/
    if(myid==myid1){
        signals_all = read_su(srcname, &nsrc_glob, &nt);
    }
    MPI_Barrier(mpi_comm_no);
    MPI_Bcast(&nsrc_glob,1,MPI_INT,0,mpi_comm_no);
    MPI_Bcast(&nt,1,MPI_INT,0,mpi_comm_no);
    if(myid!=myid1){
        signals_all = matrix(1,nsrc_glob,1,nt);
    }
    MPI_Bcast(&signals_all[1][1],nsrc_glob*nt,MPI_FLOAT,0,mpi_comm_no);

    /***************************************************************/
    /*From the above, we just get signals_all                      */
    /***************************************************************/
    
    /* find this single source positions on subdomains */ 
    *signals = matrix(1,ntr_loc,1,nt);
    seis_glob2loc(ntr_loc, nt, signals_all, *signals, srcpos_loc);
    free_matrix(signals_all,1,nsrc_glob,1,nt);
}