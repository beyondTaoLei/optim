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
@brief find the intersection of trace indexes between two measurements.

@details For two shot gathers, suppose one is observed data, another 
is synthetic data, one plans to measure the traveltime misfit, and the 
first-breaks picking time maybe in different receivers' position, which 
means the number and the index of receivers can be different. In this 
case, we should find the common receivers.

@file generate_par_lsqr.c
@author Tao Lei (leit@sustech.edu.cn), Wei Zhang
@date 20 Jul 2019
@param[in]  verbose \b 1: print the details; \b 0: run in quiet pattern, default value;
@param[in]  fobsf    the folder of the first-break picking time of observed data;
@param[in]  fsynf    the folder of the first-break picking time of synthetic data;
@param[in]  flist    the list of events which are your targeted data;
@param[in]  fpos_s   the relation between local index and global index for sources;
@param[in]  fpos_r   the relation between local index and global index for receivers;
@param[in,out]  fout_k   the list of kernels which will be calculated in the 
            next step (<em>output for PLSQR3 software</em>);
@param[in,out]  fout_t   the first-break picking time difference 
            (<em>diff=Tobs-Tsyn, output for PLSQR3 software</em>);

@code  
generate_par_lsqr fobsf=picked_times/obs fsynf=picked_times/syn flist=list/event.list fpos_s=list/S_POS.dat fpos_r=list/R_POS.dat fout_k=list/ker.list fout_t=list/measure.list
@endcode
*/
#include "comm.h"
#include "par.h"

    
/* Global variables */
int xargc; 
char **xargv;

void find_nums(char *kername, int *n1, int *n2);
int main(int argc, char *argv[]) {
    
    /* IN and OUT */
    int verbose;
    char *fobsf, *fsynf, *flist, *fout_k, *fout_t, *fpos_s, *fpos_r;

    /* Local variables */
    int ns, ntr, ishot, itr, i, j, nrow, nrow2, icount, num;
    float **pos_s, **pos_r, **tobs, **tsyn, **tdiff;
    char filename[256], fobs[256], fsyn[256];
    FILE *fp= NULL, *fp_ker= NULL, *fp_time= NULL;
    
    /* initialize getpar */
    initargs(argc,argv);
    
    /* 1. read formatted input list */
    verbose=0; getparint("verbose",&verbose);
    if(!getparstring("fobsf", &fobsf)) throw_error("fobsf must be specified!\n");
    if(!getparstring("fsynf", &fsynf)) throw_error("fsynf must be specified!\n");
    if(!getparstring("flist", &flist)) throw_error("flist must be specified!\n");
    if(!getparstring("fpos_s", &fpos_s)) throw_error("fpos_s must be specified!\n");
    if(!getparstring("fpos_r", &fpos_r)) throw_error("fpos_r must be specified!\n");
    if(!getparstring("fout_k", &fout_k)) throw_error("fout_k must be specified!\n");
    if(!getparstring("fout_t", &fout_t)) throw_error("fout_t must be specified!\n");
    
    if(verbose){
        printf("\n******\n");
        printf("%-18s| %s\n","fobsf", fobsf);
        printf("%-18s| %s\n","fsynf",fsynf);
        printf("%-18s| %s\n","flist",flist);
        printf("%-18s| %s\n","fpos_s",fpos_s);
        printf("%-18s| %s\n","fpos_r",fpos_r);
        printf("%-18s| %s\n","fout_k",fout_k);
        printf("%-18s| %s\n","fout_t",fout_t);
        printf("******\n");
    }
    
    /* 2. read ASCII-format file */ 
    pos_s=read_text(2, 2, &ns,fpos_s);
    pos_r=read_text(2, 2, &ntr,fpos_r);
    tdiff=matrix(1, ns*ntr, 1, 3);
    
    fp = fopen(flist, "r");
    if (fp == NULL) {
        printf("Can't open input list file!\n");
        exit(1);
    }
    icount=0;
    num=0;
    while (!feof(fp)) {
        if((fscanf(fp,"%s",filename)!=1)||(strlen(filename)<3)){
            break;
        }
        icount++;
        if(verbose) printf("%s\n", filename);
        sprintf(fobs,"%s/%s", fobsf, filename);
        sprintf(fsyn,"%s/%s", fsynf, filename);
        ishot = StringToInteger(filename);
        
        tobs=read_text(2, 2, &nrow, fobs);
        tsyn=read_text(2, 2, &nrow2,fsyn);
        if(nrow!=nrow2||ntr!=nrow){
            printf("fobs= %s\n",fobs);
            printf("fsyn= %s\n",fsyn);
            printf("fpos_r= %s\n",fpos_r);
            throw_error("Please check the above 3 files!\n");
        }
        for(itr=1;itr<=ntr;itr++){
            if(tobs[itr][2]==-1.0||tsyn[itr][2]==-1.0){
                continue;
            }
            num++;
            tdiff[num][1]=pos_s[ishot][2];
            tdiff[num][2]=pos_r[itr][2];
            tdiff[num][3]=tobs[itr][2]-tsyn[itr][2];       
        }
    }
    
    /*remove the duplication, such as ker46_150.bin and ker150_46.bin*/
    for(i=1;i<num;i++){
        for(j=i+1;j<=num;j++){
            if((int)tdiff[i][1]==(int)tdiff[j][2]&&
                (int)tdiff[i][2]==(int)tdiff[j][1]){
                /*printf("%d %d  %d %d\n",(int)tdiff[i][1],(int)tdiff[j][2],\
                (int)tdiff[i][2],(int)tdiff[j][1]);*/
                if(tdiff[i][3]*tdiff[j][3]>0){
                    tdiff[i][3]=(tdiff[i][3]+tdiff[j][3])/2.0;
                    tdiff[j][3]=-100.;
                }else{
                    tdiff[i][3]=-100.;
                    tdiff[j][3]=-100.;
                }
                break;
            }
        }
    }
    /* write in the disk */
    fp_time = fopen(fout_t, "w");
    fp_ker = fopen(fout_k, "w");
    for(i=1;i<num;i++){
        if(tdiff[i][3]!=-100.){
            fprintf(fp_time,"%8.3e\n",tdiff[i][3]);
            fprintf(fp_ker,"ker%d_%d.bin\n",(int)tdiff[i][1],(int)tdiff[i][2]);
        }
    }
    fclose(fp_ker);
    fclose(fp);
    
    if(verbose) printf("There are %d shots\n",icount);
    free_matrix(tdiff, 1, ns*ntr, 1, 3);
    free_matrix(pos_s, 1, ns, 1, 2);
    free_matrix(pos_r, 1, ntr, 1, 2);
    free_matrix(tobs, 1, ntr, 1, 2);
    free_matrix(tsyn, 1, ntr, 1, 2);
    
    if(verbose) printf("finished...\n");
    return 0;
}

void find_nums(char *kername, int *n1, int *n2){
    /*kername is as ker46_150.bin, here n1=46, n2=150 */
    char str_s[40];
    char *pline,*pdot;
    int num_l;
    
    memset(str_s, '\0', sizeof(str_s)); 
    pline = strstr(kername, "_");
    pdot = strstr(kername, ".");
    if(pline!=NULL){
        num_l=pline-kername-3;
        strncpy(str_s,kername+3, num_l);
        *n1 = atoi(str_s);
        num_l=pdot-pline+1;
        strncpy(str_s,pline+1, num_l);
        *n2 = atoi(str_s);
        //printf("is=%d, ir=%d\n",*n1, *n2);
    }else{
        printf("Can't extract the number from string %s",kername);
        exit(1);
    }
    return;
    
}
