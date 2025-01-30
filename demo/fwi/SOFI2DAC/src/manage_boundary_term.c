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
void prepare_boundary_saving(int ishot, int nshot1, st_bdr *bdr){
    
    extern int MYID, IENDX, IENDY, NXG, NYG;
    extern int FREE_SURF, FW, FDORDER, FDORDER_TIME, POS[3];
    extern int FLAG_SU_BIN, SRC_SNAP_TYPE;
    
    extern float DH, MEMORY_DEPTH, REFREC[4];
    extern char SNAP_COORDS_FILE[STRING_SIZE];
    extern char SNAP_FILE[STRING_SIZE];
    
    /* local variables */
    char srcname[256];
    
    /* The following is for saving the snapshot. 
     * in this case, the source is excited from SU-format file,
     * and, one will save the wavefield on the boundaries for 
     * boundary saving strategy, RWI, truncated Newton method and 
     * so on. and for which target, it depends on SRC_SNAP_TYPE
     */
    //if(FLAG_SU_BIN!=1) declare_error ( "Please be sure the source type and snapshot type!!" );
    if(ishot == nshot1){
        split_wfd_pos_loc(MYID, 0, MPI_COMM_WORLD, IENDX, IENDY, POS, DH, REFREC, FDORDER, 
                            FREE_SURF,FW, NXG, NYG, MEMORY_DEPTH, SRC_SNAP_TYPE, SNAP_COORDS_FILE, 
                            &(bdr->nsnap_glob), &(bdr->pos_loc), &(bdr->nsnap_loc));
    }
    if(bdr->nsnap_loc>0){
        /* Sp source */
        sprintf(srcname,"%s.%s.shot%i.%i.%i.bin",SNAP_FILE,"sp",ishot,POS[1],POS[2]);
        bdr->fp_sp = fopen(srcname, "w");
        //bdr->fp_sp = open_wfd_file_cmp(SNAP_FILE, ishot, POS, "sp", "w");
        /* Vx source */
        sprintf(srcname,"%s.%s.shot%i.%i.%i.bin",SNAP_FILE,"vx",ishot,POS[1],POS[2]);
        bdr->fp_vx = fopen(srcname, "w");
        //bdr->fp_vx  = open_wfd_file_cmp(SNAP_FILE, ishot, POS, "vx", "w");
        /* Vy source */
        sprintf(srcname,"%s.%s.shot%i.%i.%i.bin",SNAP_FILE,"vy",ishot,POS[1],POS[2]);
        bdr->fp_vy = fopen(srcname, "w");
        //bdr->fp_vy  = open_wfd_file_cmp(SNAP_FILE, ishot, POS, "vy", "w");
    }
    
}

void write_snap_boundary_file(int ishot, int nt, float **psp, float **pvx, float **pvy, st_bdr *bdr){
    
    extern int NT, NX, NY, SAVE_LAST_SNAP, POS[3];
    extern char SNAP_FILE[STRING_SIZE];
    
    /* local variables */
    int **coord_loc = NULL;
    char srcname[256], str_s[256];
    FILE *fp_tmp=NULL;
    int i,j,l;
    
    if(bdr->nsnap_loc>0){
        //printf("%d %d %d ==\n",MYID, bdr->nsnap_loc, nt);
        write_snap_boundary(bdr->fp_sp,  psp, bdr->pos_loc, bdr->nsnap_loc, nt-1);
        write_snap_boundary(bdr->fp_vx,  pvx,  bdr->pos_loc, bdr->nsnap_loc, nt-1);
        write_snap_boundary(bdr->fp_vy,  pvy,  bdr->pos_loc, bdr->nsnap_loc, nt-1);
    }
    //printf("nt,NT,SAVE_LAST_SNAP=%d %d %d\n",nt,NT,SAVE_LAST_SNAP);
    if(nt==NT&&SAVE_LAST_SNAP==1){
        coord_loc=imatrix(1,3,1,NX*NY);
        l=1;
        for(j=1;j<=NY;j++){
        for(i=1;i<=NX;i++){ 
            coord_loc[1][l]=i;
            coord_loc[2][l]=j;
            l++;
        }}
        sprintf(str_s,"%s_lst",SNAP_FILE);
        sprintf(srcname,"%s.%s.shot%i.%i.%i.bin",str_s,"sp",ishot,POS[1],POS[2]);
        fp_tmp = fopen(srcname, "w");
        //fp_tmp = open_wfd_file_cmp(str_s, ishot, POS, "sp", "w");
        write_snap_boundary(fp_tmp, psp, coord_loc, NX*NY, 0);
        fclose(fp_tmp);
        sprintf(srcname,"%s.%s.shot%i.%i.%i.bin",str_s,"vx",ishot,POS[1],POS[2]);
        fp_tmp = fopen(srcname, "w");
        //fp_tmp = open_wfd_file_cmp(str_s, ishot, POS, "vx", "w");
        write_snap_boundary(fp_tmp, pvx, coord_loc, NX*NY, 0);
        fclose(fp_tmp);
        sprintf(srcname,"%s.%s.shot%i.%i.%i.bin",str_s,"vy",ishot,POS[1],POS[2]);
        fp_tmp = fopen(srcname, "w");
        //fp_tmp = open_wfd_file_cmp(str_s, ishot, POS, "vy", "w");
        write_snap_boundary(fp_tmp, pvy, coord_loc, NX*NY, 0);
        fclose(fp_tmp);
        free_imatrix(coord_loc,1,3,1,NX*NY);
    }
}

void finish_boundary_file(st_bdr *bdr){
    
    if(bdr->nsnap_loc>0){
        /* Sp source */
        fclose(bdr->fp_sp);
        /* Vx source */
        fclose(bdr->fp_vx);
        /* Vy source */
        fclose(bdr->fp_vy);
    }
}

void free_boundary_struct_out_shot(st_bdr *bdr){
    
    if(bdr->nsnap_loc>0) free_imatrix(bdr->pos_loc,1,3,1,bdr->nsnap_loc);
}
