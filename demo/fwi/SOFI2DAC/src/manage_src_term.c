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
void init_struct_src_elmts(int flag, st_source *geo_src){
    /* call this function when FLAG_SU_BIN = 1 or 2 */
    extern int SRC_SNAP_TYPE, READ_LAST_SNAP;
    extern int ADD_REPLACE_VX, ADD_REPLACE_VY, ADD_REPLACE_SP, ADD_REPLACE_SNAP;
    extern float MEMORY_DEPTH;
    extern char SNAP_COORDS_FILE[STRING_SIZE];
    extern char SRC_F_VX[STRING_SIZE], SRC_F_VY[STRING_SIZE];
    extern char SRC_F_SP[STRING_SIZE], SRC_F_SNAP[STRING_SIZE];
        
    /* copies the values in json file into the elements in structure *geo_src */
    /* incident wavefields */
//     geo_srci->flag_su_bin = 2;
    /* adjoint wavefields */
//     geo_srca->flag_su_bin = 1;
    
    geo_src->flag_su_bin = flag;
    
    geo_src->src_snap_type = SRC_SNAP_TYPE;
    geo_src->read_last_snap = READ_LAST_SNAP;
    geo_src->add_replace_vx = ADD_REPLACE_VX;
    geo_src->add_replace_vy = ADD_REPLACE_VY;
    geo_src->add_replace_sp = ADD_REPLACE_SP;
    geo_src->add_replace_snap = ADD_REPLACE_SNAP;
    geo_src->memory_depth = MEMORY_DEPTH;
    sprintf(geo_src->src_f_vx, "%s", SRC_F_VX);
    sprintf(geo_src->src_f_vy, "%s", SRC_F_VY);
    sprintf(geo_src->src_f_sp, "%s", SRC_F_SP);
    sprintf(geo_src->src_f_snap, "%s", SRC_F_SNAP);
    sprintf(geo_src->snap_coords_file, "%s", SNAP_COORDS_FILE);
      
    /* source without physical meaning, just for boundary saving strategy */
    if(geo_src->flag_su_bin==2){
        geo_src->add_replace_vx  = geo_src->add_replace_snap;
        geo_src->add_replace_vy  = geo_src->add_replace_snap;
        geo_src->add_replace_sp  = geo_src->add_replace_snap;
    }
    
}

void prepare_src_position(int ishot, int nshot1, st_source *geo_src){
    
    extern int MYID, IENDX, IENDY,NT, NXG, NYG;
    extern int FREE_SURF, FW, FDORDER, FDORDER_TIME, POS[3];
    extern float DT, DH, REFREC[4];
    extern char SRC_POS_FILE[STRING_SIZE];
    
    /* local variables */
    char srcname[256];
    /* Source part */
    if(geo_src->flag_su_bin==1){/* SU-format file */
        /* source positions */
        if(geo_src->add_replace_sp!=0){
            update_shot_num(geo_src->src_f_sp, ishot, srcname);
        }else if(geo_src->add_replace_vx!=0){
            update_shot_num(geo_src->src_f_vx, ishot, srcname);
        }else if(geo_src->add_replace_vy!=0){
            update_shot_num(geo_src->src_f_vy, ishot, srcname);
        }
        if(ishot == nshot1) update_shot_num(srcname, ishot, SRC_POS_FILE);
        split_su_pos_loc(srcname, MYID, 0, MPI_COMM_WORLD, DH, REFREC, 
                         IENDX, IENDY, POS, NT, DT, &(geo_src->srcpos_loc), 
                         &(geo_src->nsrc_loc));
    }else if(geo_src->flag_su_bin==2){/* BIN-format file */ /* (FDORDER_TIME??????????)*/
        /* source positions */
        if(ishot == nshot1){
            split_wfd_pos_loc(MYID, 0, MPI_COMM_WORLD, IENDX, IENDY, 
                              POS, DH, REFREC, FDORDER, FREE_SURF,FW, 
                              NXG, NYG, geo_src->memory_depth, 
                              geo_src->src_snap_type, geo_src->snap_coords_file, 
                              &(geo_src->nsrc_glob), &(geo_src->srcpos_loc), 
                              &(geo_src->nsrc_loc));
            
        }
    }
}

void prepare_sig_file(int ishot, st_source *geo_src){
    
    extern int MYID, NX, NY;
    extern int POS[3];
        
    /* local variables */
    char srcname[256];
    /* Source part */
    if(geo_src->flag_su_bin==1){/* SU-format file */
        
        /* reads seismograms file */
        /* P source */
        if(geo_src->add_replace_sp!=0){
            update_shot_num(geo_src->src_f_sp, ishot, srcname);
            split_su_sig_loc(srcname, MYID, 0, MPI_COMM_WORLD, 
                             geo_src->srcpos_loc, geo_src->nsrc_loc, 
                             &(geo_src->signal_sp));
            geo_src->snap_sp = vector(1, geo_src->nsrc_loc);
        }
        
        /* Vx source */
        if(geo_src->add_replace_vx!=0){
            update_shot_num(geo_src->src_f_vx, ishot, srcname);
            split_su_sig_loc(srcname, MYID, 0, MPI_COMM_WORLD, 
                             geo_src->srcpos_loc, geo_src->nsrc_loc, 
                             &(geo_src->signal_vx));
            geo_src->snap_vx = vector(1, geo_src->nsrc_loc);
        }
        
        /* Vy source */
        if(geo_src->add_replace_vy!=0){
            update_shot_num(geo_src->src_f_vy, ishot, srcname);
            split_su_sig_loc(srcname, MYID, 0, MPI_COMM_WORLD, 
                             geo_src->srcpos_loc, geo_src->nsrc_loc, 
                             &(geo_src->signal_vy));
            geo_src->snap_vy = vector(1, geo_src->nsrc_loc);
        }
    }else if(geo_src->flag_su_bin==2){/* BIN-format file */ /* (FDORDER_TIME??????????)*/
        
        /* opens wavefield file */
        if(geo_src->nsrc_loc>0){
            /* Sp source */
            sprintf(srcname,"%s.%s.shot%i.%i.%i.bin",geo_src->src_f_snap,"sp",ishot,POS[1],POS[2]);
            geo_src->fp_sp = fopen(srcname, "r");
            /* Vx source */
            sprintf(srcname,"%s.%s.shot%i.%i.%i.bin",geo_src->src_f_snap,"vx",ishot,POS[1],POS[2]);
            geo_src->fp_vx = fopen(srcname, "r");
            /* Vy source */
            sprintf(srcname,"%s.%s.shot%i.%i.%i.bin",geo_src->src_f_snap,"vy",ishot,POS[1],POS[2]);
            geo_src->fp_vy = fopen(srcname, "r");
        }
        geo_src->snap_sp = vector(1, geo_src->nsrc_loc);
        geo_src->snap_vx  = vector(1, geo_src->nsrc_loc);
        geo_src->snap_vy  = vector(1, geo_src->nsrc_loc);
        geo_src->snap_all = vector(1, NX*NY);
    }
}

void extract_snap_time(int nt, int nt_snap, st_source *geo_src){
    
    
    /*******************************************************/
    /* modified by Tao, in order to insert different source 
    * for FWI, RTM or some other cases */
    if(geo_src->flag_su_bin==1){/* SU-format file */
        if(geo_src->add_replace_sp!=0) 
            extract_snap_from_su(geo_src->signal_sp, geo_src->snap_sp, geo_src->nsrc_loc, nt);
        if(geo_src->add_replace_vx!=0) 
            extract_snap_from_su(geo_src->signal_vx, geo_src->snap_vx, geo_src->nsrc_loc, nt);
        if(geo_src->add_replace_vy!=0) 
            extract_snap_from_su(geo_src->signal_vy, geo_src->snap_vy, geo_src->nsrc_loc, nt);
    }else if(geo_src->flag_su_bin==2){/* BIN-format file */ /* (FDORDER_TIME??????????)*/
        if(geo_src->nsrc_loc>0){
            read_wfd_file(geo_src->fp_sp, &geo_src->snap_sp[1], geo_src->nsrc_loc, nt_snap);
            read_wfd_file(geo_src->fp_vx,  &geo_src->snap_vx[1],  geo_src->nsrc_loc, nt_snap);
            read_wfd_file(geo_src->fp_vy,  &geo_src->snap_vy[1],  geo_src->nsrc_loc, nt_snap);
        }
    }
    /*******************************************************/
}

void insert_particle_velocity(float **pvx, float **pvy, float  **prip, float **prjp, st_source *geo_src){
    
    extern float DT, DH;
    
    /* local variables */
    int **srcpos_loc = NULL;
    int i,j,l;
    
    /*******************************************************/
    /* modified by Tao, in order to insert different source 
    * for FWI, RTM or some other cases */
    /* reads source information, including positon, signals and so on */
    srcpos_loc = geo_src->srcpos_loc;
    /* vx */
    if(geo_src->add_replace_vx){
        for (l=1;l<=geo_src->nsrc_loc;l++){
            i=(int)srcpos_loc[1][l];
            j=(int)srcpos_loc[2][l];
            pvx[j][i] = insert_src(pvx[j][i], geo_src->snap_vx[l], prip[j][i]*DT/( DH*DH ), geo_src->add_replace_vx);
        }
    }
    /* vy */
    if(geo_src->add_replace_vy){
        for (l=1;l<=geo_src->nsrc_loc;l++){
            i=(int)srcpos_loc[1][l];
            j=(int)srcpos_loc[2][l];
            pvy[j][i] = insert_src(pvy[j][i], geo_src->snap_vy[l], prjp[j][i]*DT/( DH*DH ), geo_src->add_replace_vy);
        }
    }
    /*******************************************************/
}

void insert_velocity_last_snap(int ishot, float **pvx, float **pvy, st_source *geo_src){
    
    extern int NX, NY, POS[3];
    
    /* local variables */
    int i,j,l;
    char srcname[256], str_s[256];
    FILE *fp_tmp=NULL;
    
    sprintf(str_s,"%s_lst",geo_src->src_f_snap);
    sprintf(srcname,"%s.%s.shot%i.%i.%i.bin",str_s,"vx",ishot,POS[1],POS[2]);
    fp_tmp = fopen(srcname, "r");
    read_wfd_file(fp_tmp,  &geo_src->snap_all[1],  NX*NY, 0);
    fclose(fp_tmp);
    l=1;
    for(j=1;j<=NY;j++){
    for(i=1;i<=NX;i++){                        
        pvx[j][i]=geo_src->snap_all[l];
        l++;
    }}
    
    sprintf(srcname,"%s.%s.shot%i.%i.%i.bin",str_s,"vy",ishot,POS[1],POS[2]);
    fp_tmp = fopen(srcname, "r");
    read_wfd_file(fp_tmp,  &geo_src->snap_all[1],  NX*NY, 0);
    fclose(fp_tmp);
    l=1;
    for(j=1;j<=NY;j++){
    for(i=1;i<=NX;i++){                        
        pvy[j][i]=geo_src->snap_all[l];
        l++;
    }}
    
}

void insert_stress(float **psp, st_source *geo_src){
    
    extern float DH;
    
    /* local variables */
    int **srcpos_loc = NULL;
    int i,j,l;
    
    /*******************************************************/
    /* modified by Tao, in order to insert different source 
    * for FWI, RTM or some other cases */
    /* reads source information, including positon, signals and so on */
    /* explosive source */
    //if ( SOURCE_TYPE == 1 )
    //    psource ( nt, psxx, psyy, srcpos_loc, signals, nsrc_loc );
    srcpos_loc = geo_src->srcpos_loc;
    /* sp */
    if(geo_src->add_replace_sp){
        for (l=1;l<=geo_src->nsrc_loc;l++){
            i=(int)srcpos_loc[1][l];
            j=(int)srcpos_loc[2][l];
            psp[j][i] = insert_src(psp[j][i],geo_src->snap_sp[l],1.0/(DH*DH),geo_src->add_replace_sp);
        }
    }
    /*******************************************************/
}

void insert_stress_last_snap(int ishot, float **psp, st_source *geo_src){
    
    extern int NX, NY, POS[3];
    
    /* local variables */
    int i,j,l;
    char srcname[256], str_s[256];
    FILE *fp_tmp=NULL;
    
    sprintf(str_s,"%s_lst",geo_src->src_f_snap);
    sprintf(srcname,"%s.%s.shot%i.%i.%i.bin",str_s,"sp",ishot,POS[1],POS[2]);
    fp_tmp = fopen(srcname, "r");
    read_wfd_file(fp_tmp,  &geo_src->snap_all[1],  NX*NY, 0);
    fclose(fp_tmp);
    l=1;
    for(j=1;j<=NY;j++){
    for(i=1;i<=NX;i++){                        
        psp[j][i]=geo_src->snap_all[l];
        l++;
    }}
    
}
    
void finish_src_file(int ns, st_source *geo_src){
    
    /* Source part */
    if(geo_src->flag_su_bin==1){/* SU-format file */
        if(geo_src->nsrc_loc>0){
            free_imatrix(geo_src->srcpos_loc,1,3,1,geo_src->nsrc_loc);
            /* reads seismograms file */
            /* P source */
            if(geo_src->add_replace_sp!=0){
                free_matrix(geo_src->signal_sp, 1, geo_src->nsrc_loc, 1, ns);
                free_vector(geo_src->snap_sp, 1, geo_src->nsrc_loc);
            }
            
            /* Vx source */
            if(geo_src->add_replace_vx!=0){
                free_matrix(geo_src->signal_vx, 1, geo_src->nsrc_loc, 1, ns);
                free_vector(geo_src->snap_vx, 1, geo_src->nsrc_loc);
            }
            
            /* Vy source */
            if(geo_src->add_replace_vy!=0){
                free_matrix(geo_src->signal_vy, 1, geo_src->nsrc_loc, 1, ns);
                free_vector(geo_src->snap_vy, 1, geo_src->nsrc_loc);
            }
        }
    }else if(geo_src->flag_su_bin==2){/* BIN-format file */
        /* opens wavefield file */
        if(geo_src->nsrc_loc>0){
            /* Sp source */
            fclose(geo_src->fp_sp);
            /* Vx source */
            fclose(geo_src->fp_vx);
            /* Vy source */
            fclose(geo_src->fp_vy);
        }
    }
}

void free_src_struct_out_shot(int nx, int ny, st_source *geo_src){
    
    /* Source part from BIN-format file flag_su_bin==2 */
    free_imatrix(geo_src->srcpos_loc,1,3,1,geo_src->nsrc_loc);
    free_vector(geo_src->snap_sp,1, geo_src->nsrc_loc);
    free_vector(geo_src->snap_vx ,1, geo_src->nsrc_loc);
    free_vector(geo_src->snap_vy ,1, geo_src->nsrc_loc);
    free_vector(geo_src->snap_all,1, nx*ny);
}