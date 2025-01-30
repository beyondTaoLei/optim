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
#include "globvar.h"      /* definition of global variables  */

char * trim(char * src){
    int i = 0;
    char *begin = src;
    while(src[i] != '\0'){
        if(src[i] != ' '){
                break;
        }else{
                begin++;
        }
        i++;
    }
    for(i = strlen(src)-1; i >= 0;  i--){
        if(src[i] != ' '){
                break;
        }else{
                src[i] = '\0';
        }
    }
    return begin;
}

int main(int argc, char **argv){
    
    
    char *fileinp="";
    char *cmpnt, *str_shot;
    
    int num, nsrc_loc,myid,it,l,hin,i;
    int NXG,NYG;
    int **grids=NULL;
    int   POS1[3];
    int **srcpos_loc;
    float wfd;
    char snap_in[256], snap_out[256];
    char str_grid_file[256];
    FILE *fp2, *fp1;
    
    

    /* ============================================== */
    /* Open parameter-file to check if auto mode or not*/
    fileinp = argv[1];
    cmpnt = argv[2];
    str_shot = argv[3];
    
    printf(" ***********************************************************\n");
    printf(" This is program SNAPMERGE. \n");
    printf(" Merge of snapshot files from the parallel  \n 2-D Viscoelastic Finite Difference Modeling      \n");
    printf("                                                            \n");
    printf(" written by  T. Bohlen                          \n");
    printf(" Geophysical Institute, Department of Physics,         \n");
    printf(" Institute of Technology, Karlsruhe, Germany         \n");
    printf(" http://www.gpi.kit.edu \n");
    printf(" ***********************************************************\n");
    printf("\n");
    printf(" Syntax example if excecuted from ./par directory: ../bin/snapmerge in_and_out/sofi2D.json \n");
    printf(" Input file for the snapmerge process from command line : %s \n",fileinp);

    //FP = fopen(fileinp,"r");
    if ((FP=fopen(fileinp,"r"))==NULL){
        printf(" Opening input file failed.");
        exit(0);
    }
    else printf(" Opening input file was successful.\n\n");
    fclose(FP);
    //fscanf(FP, "%s %s = %i", infostr, modestr, &RUNMODE);

    /* =================================================== */
    
    /* read standard input file */
    if (strstr(fileinp,".json"))
            //read json formated input file
            read_par_json(stdout, fileinp);
    /* =================================================== */
    printf(" \n\n\n");
    printf(" Please make sure SNAP =5 and SNAP_FILE is correct,\n");
    printf(" then type command as < ../bin/snapmerge2  in_and_out/sofi2D_new.json vx 1 >");
    printf(" \n\n\n");

    NXG=NX;
    NYG=NY;	
    NX = NXG/NPROCX;
    NY = NYG/NPROCY;
    NT = iround ( TIME/DT );
    
    if(SRC_SNAP_TYPE==3){
        grids = read_coordinates(&num, DH, REFREC, SNAP_COORDS_FILE);
    }else{
        grids = evaluate_wfd_coords(DH, REFREC, FDORDER, FREE_SURF,
                FW, NXG, NYG, MEMORY_DEPTH, SRC_SNAP_TYPE, &num);
    }
    
    /* write the grid points into the disk */
    sprintf(str_grid_file,"%s.grids.dat",SNAP_FILE);
    fp2=fopen(str_grid_file,"w");
    for(l=1;l<=num;l++){
        fprintf(fp2, "%d %d\n",grids[1][l],grids[2][l]);}
    fclose(fp2);
    printf("You will save the wavefields from %d points.\n",num);
    
    sprintf(snap_out,"%s.%s.shot%s.bin",SNAP_FILE,trim(cmpnt),str_shot);
    fp2 = fopen(snap_out, "w");
    printf("NT=%d\n",NT);
    for(myid=0;myid<NPROCX*NPROCY;myid++){
        POS1[1] = myid % NPROCX;
        POS1[2] = (myid/NPROCX);
        
        srcpos_loc = splitpos(NX, NY, POS1, grids, &nsrc_loc, num);
        //printf("myid=%d,it=%d,POS1=%d %d,nsrc_loc=%d\n",myid,it,POS1[1],POS1[2],nsrc_loc);
        if(nsrc_loc>0){
            sprintf(snap_in,"%s.%s.shot%s.%i.%i.bin",SNAP_FILE,
                trim(cmpnt),str_shot,POS1[1],POS1[2]);
            printf("%s\n",snap_in);
            fp1 = fopen(snap_in, "r");
            for(it=1;it<=NT;it++){
                if(it%100==0) printf("it=%d of %d\n",it,NT);
                hin=it-1;
                fseek(fp1, hin*nsrc_loc*sizeof(float),SEEK_SET);
                for(l=1;l<=nsrc_loc;l++){
                    fread(&wfd,sizeof(float),1,fp1);
                    fseek(fp2, (hin*num+srcpos_loc[3][l]-1)*sizeof(float),SEEK_SET);
                    fwrite(&wfd,sizeof(float),1,fp2);
                }
            }
            fclose(fp1);
            printf("myid=%d\n\n",myid);
            for(l=1;l<=nsrc_loc;l++){
                i=srcpos_loc[3][l];
                printf("%d %d %d %d %d\n",i,srcpos_loc[1][l],srcpos_loc[2][l],grids[1][i],grids[2][i]);
            }
            free_imatrix(srcpos_loc,1,3,1,nsrc_loc);
        }
        
        
    }
    
    fclose(fp2);
    
    free_imatrix(grids,1,3,1,num);
    printf("You will save the wavefields from %d points.\n",num);
    
}