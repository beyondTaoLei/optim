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

void split_wfd_pos_loc(int myid, int myid1, MPI_Comm mpi_comm_no, 
                        int iendx, int iendy, int *POS1, float dh, 
                        float *refrec, int fdorder, int free_surf, 
                        int fw, int nxg, int nyg, float memory_depth, 
                        int source_snap_type, char snap_coords_file[256],
                        int *nsrc_glob, int ***srcpos_loc, int *nsrc_loc){
    /* If the sources are added as snapshots, which are useful for some 
     * algorithms, such as truncated Newton, boundary saving strategy. 
     * One should read the signals from BIN-format file.
     */
    
    /* local variables */
    int num;
    int **grids=NULL;
    
    if(myid==myid1){
        if(source_snap_type==3){
            grids = read_coordinates(&num, dh, refrec, snap_coords_file);
        }else{
            grids = evaluate_wfd_coords(dh, refrec, fdorder, free_surf,
                    fw, nxg, nyg, memory_depth, source_snap_type, &num);
        }
        //printf("num=%d %d",source_snap_type,num);
    }
    MPI_Barrier(mpi_comm_no);
    MPI_Bcast(&num,1,MPI_INT,0,mpi_comm_no);
    if (myid!=myid1) grids=imatrix(1,3,1,num);
    MPI_Bcast(&grids[1][1],num*3,MPI_INT,0,mpi_comm_no);
    
    *srcpos_loc = splitpos(iendx, iendy, POS1, grids, nsrc_loc, num);
    free_imatrix(grids,1,3,1,num);
    
    *nsrc_glob = num;

}

void split_wfd_sig_loc(FILE *fp_sp, int myid, int myid1, 
                       MPI_Comm mpi_comm_no, float *snap_all, 
                       float *snap_loc, int ** recpos_loc, 
                       int nsrc_glob, int nsrc_loc, int it){
    /* extracts current signals(one snapshot) from local source 
     * position(a time series) in current node */
    
    /*if(myid==myid1){
        fseek(fp_sp, it*nsrc_glob*sizeof(float),SEEK_SET);
        fread(&snap_all[1],sizeof(float),nsrc_glob,fp_sp);
    }*/
    if(myid==myid1){
        read_wfd_file(fp_sp, &snap_all[1], nsrc_glob, it);
    }
    MPI_Barrier(mpi_comm_no);
    MPI_Bcast(&snap_all[1],nsrc_glob,MPI_FLOAT,0,mpi_comm_no);
    /* splits snap_all into snap_sp in the local node */
    seis_glob2loc_vec(nsrc_loc, snap_all, snap_loc, recpos_loc);

}



                      
                         