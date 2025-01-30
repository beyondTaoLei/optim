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

void allocate_mem_seis(int flag, float ***sectionvx, float ***sectionvy, float ***sectionp, 
                       float ***sectiondiv, int ntr, int ns){
    
    extern int SEISMO;
    
    if(flag==1){/* allocate memory */
        switch ( SEISMO ) {
            case 1: /* particle velocities only */
                *sectionvx = matrix ( 1, ntr, 1, ns );
                *sectionvy = matrix ( 1, ntr, 1, ns );
                break;
            case 2: /* pressure only */
                *sectionp = matrix ( 1, ntr, 1, ns );
                break;
            case 3: /* curl and div only */
//                 *sectioncurl = matrix ( 1, ntr, 1, ns );
                *sectiondiv = matrix ( 1, ntr, 1, ns );
                break;
            case 4: /* everything */
                *sectionvx = matrix ( 1, ntr, 1, ns );
                *sectionvy = matrix ( 1, ntr, 1, ns );
//                 *sectioncurl = matrix ( 1, ntr, 1, ns );
                *sectiondiv = matrix ( 1, ntr, 1, ns );
                *sectionp = matrix ( 1, ntr, 1, ns );
                break;
        }
    }else{ /* deallocate memory */
        switch ( SEISMO ) {
            case 1: /* particle velocities only */
                free_matrix(*sectionvx, 1, ntr, 1, ns );
                free_matrix(*sectionvy, 1, ntr, 1, ns );
                break;
            case 2: /* pressure only */
                free_matrix(*sectionp, 1, ntr, 1, ns );
                break;
            case 3: /* curl and div only */
//                 free_matrix(*sectioncurl, 1, ntr, 1, ns );
                free_matrix(*sectiondiv, 1, ntr, 1, ns );
                break;
            case 4: /* everything */
                free_matrix(*sectionvx, 1, ntr, 1, ns );
                free_matrix(*sectionvy, 1, ntr, 1, ns );
//                 free_matrix(*sectioncurl, 1, ntr, 1, ns );
                free_matrix(*sectiondiv, 1, ntr, 1, ns );
                free_matrix(*sectionp, 1, ntr, 1, ns );
                break;
        }
    }
    
}

void allocate_mem_seis2(int flag, float ***sectionvx, float ***sectionvy, float ***sectionp, 
                       int ntr, int ns){
    
    if(flag==1){/* allocate memory */
        *sectionvx = matrix ( 1, ntr, 1, ns );
        *sectionvy = matrix ( 1, ntr, 1, ns );
        *sectionp  = matrix ( 1, ntr, 1, ns );
    }else{ /* deallocate memory */
        free_matrix(*sectionvx, 1, ntr, 1, ns );
        free_matrix(*sectionvy, 1, ntr, 1, ns );
        free_matrix(*sectionp, 1, ntr, 1, ns );
    }     
}

void save_seis(FILE *fp, float **sectionvx, float **sectionvy, 
               float **sectionp, float **sectiondiv,
               float ** seismo_fulldata, int * recswitch, int **recpos, 
               int **recpos_loc, int ntr_glob, int ns, int ishot,
               float Sx, float Sy){
    
    extern int SEISMO, MYID;
    
    /* merge of seismogram data from all PE and output data collectively */
//     printf("leit=%d %d %d %d %d %f %f\n",MYID, SEISMO,ntr_glob,ns,ishot,Sx, Sy);
    switch ( SEISMO ) {
        case 1 : /* particle velocities only */
            catseis ( sectionvx, seismo_fulldata, recswitch, ntr_glob,ns );
            if ( MYID==0 ) saveseis_glob ( fp,seismo_fulldata,recpos,recpos_loc,ntr_glob,Sx, Sy,ishot,ns,1 );
            catseis ( sectionvy, seismo_fulldata, recswitch, ntr_glob,ns );
            if ( MYID==0 ) saveseis_glob ( fp,seismo_fulldata,recpos,recpos_loc,ntr_glob,Sx, Sy,ishot,ns,2 );
            
            break;
        case 2 : /* pressure only */
            catseis ( sectionp, seismo_fulldata, recswitch, ntr_glob,ns );
            if ( MYID==0 ) saveseis_glob ( fp,seismo_fulldata,recpos,recpos_loc,ntr_glob,Sx, Sy,ishot,ns,4 );
            
            break;
        case 3 : /* curl and div only */
            catseis ( sectiondiv, seismo_fulldata, recswitch, ntr_glob,ns );
            if ( MYID==0 ) saveseis_glob ( fp,seismo_fulldata,recpos,recpos_loc,ntr_glob,Sx, Sy,ishot,ns,5 );
//             catseis ( sectioncurl, seismo_fulldata, recswitch, ntr_glob,ns );
//             if ( MYID==0 ) saveseis_glob ( fp,seismo_fulldata,recpos,recpos_loc,ntr_glob,Sx, Sy,ishot,ns,6 );
            
            break;
        case 4 : /* everything */
            /*fprintf(fp," start merging, ntr= %d : \n",ntr_glob);
                fprintf(stdout,"Message from PE %d\n",MYID);*/
            catseis ( sectionvx, seismo_fulldata, recswitch, ntr_glob,ns );
            if ( MYID==0 ) saveseis_glob ( fp,seismo_fulldata,recpos,recpos_loc,ntr_glob,Sx, Sy,ishot,ns,1 );
            catseis ( sectionvy, seismo_fulldata, recswitch, ntr_glob,ns );
            if ( MYID==0 ) saveseis_glob ( fp,seismo_fulldata,recpos,recpos_loc,ntr_glob,Sx, Sy,ishot,ns,2 );
            catseis ( sectionp, seismo_fulldata, recswitch, ntr_glob,ns );
            if ( MYID==0 ) saveseis_glob ( fp,seismo_fulldata,recpos,recpos_loc,ntr_glob,Sx, Sy,ishot,ns,4 );
            catseis ( sectiondiv, seismo_fulldata, recswitch, ntr_glob,ns );
            if ( MYID==0 ) saveseis_glob ( fp,seismo_fulldata,recpos,recpos_loc,ntr_glob,Sx, Sy,ishot,ns,5 );
//             catseis ( sectioncurl, seismo_fulldata, recswitch, ntr_glob,ns );
//             if ( MYID==0 ) saveseis_glob ( fp,seismo_fulldata,recpos,recpos_loc,ntr_glob,Sx, Sy,ishot,ns,6 );
            
            break;
        default :
            break;
            
    }
    
}

void save_seis2(FILE *fp, float **sectionvx, float **sectionvy, float **sectionp, 
                float ** seismo_fulldata, int * recswitch, int **recpos, 
                int **recpos_loc, int ntr_glob, int ns, int ishot, float Sx, float Sy){
    
    extern int MYID;
    
    /* merge of seismogram data from all PE and output data collectively */
    catseis ( sectionvx, seismo_fulldata, recswitch, ntr_glob,ns );
    if ( MYID==0 ) saveseis_glob2 ( fp,seismo_fulldata,recpos,recpos_loc,ntr_glob,Sx, Sy,ishot,ns,1 );
    catseis ( sectionvy, seismo_fulldata, recswitch, ntr_glob,ns );
    if ( MYID==0 ) saveseis_glob2 ( fp,seismo_fulldata,recpos,recpos_loc,ntr_glob,Sx, Sy,ishot,ns,2 );
    catseis ( sectionp, seismo_fulldata, recswitch, ntr_glob,ns );
    if ( MYID==0 ) saveseis_glob2 ( fp,seismo_fulldata,recpos,recpos_loc,ntr_glob,Sx, Sy,ishot,ns,4 );
                   
}
            