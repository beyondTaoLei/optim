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
/*
 * This function provides several matrix operations
 */

#include "fd.h"

void write_matrix_disk(float ** local_matrix,char path_name[STRING_SIZE]){
    char joint[225];
    FILE *FPjoint;
    extern int POS[3],MYID;
    extern int IDX, IDY,NX,NY;
    int i,j;
    sprintf(joint,"%s.bin.%i.%i",path_name,POS[1],POS[2]);
    FPjoint=fopen(joint,"wb");
    
    for (i=1;i<=NX;i=i+IDX){
        for (j=1;j<=NY;j=j+IDY){
            fwrite(&local_matrix[j][i],sizeof(float),1,FPjoint);
        }
    }
    
    fclose(FPjoint);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    /* merge gradient file */
    sprintf(joint,"%s.bin",path_name);
    if (MYID==0) mergemod(joint,3);
    
    MPI_Barrier(MPI_COMM_WORLD);
    sprintf(joint,"%s.bin.%i.%i",path_name,POS[1],POS[2]);
    remove(joint);
};

/* Return global maximum of Matrix */
float global_maximum(float ** gradiant_1) {
    int i, j;
    float max=0.0,max1=0.0;
    extern int NX, NY;
    for (j=1;j<=NY;j++){
        for (i=1;i<=NX;i++){
            if((i*j == 1) || (max == 0.0)) {
                max = fabs(gradiant_1[j][i]);
            } else {
                if(fabs(gradiant_1[j][i]) > max) {
                    max = fabs(gradiant_1[j][i]);
                }
            }
        }
    }
    MPI_Allreduce(&max,  &max1,  1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
    return max1;
};

/* Return average  of Matrix */
float average_matrix(float ** matrix){
    
    extern int POS[3],MYID;
    extern int IDX, IDY,NX,NY;
    extern int NXG, NYG;
    
    int i,j;
    float local_sum=0;
    float global_sum;
    float buf1=0, buf2=0;
    float average;
    
    for (j=1;j<=NY;j=j+IDY){
        for (i=1;i<=NX;i=i+IDX){
            local_sum+=matrix[j][i];
        }
    }
    
    buf1=local_sum;
    buf2=0;
    MPI_Allreduce(&buf1,&buf2, 1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    global_sum=buf2;
    
    average=global_sum/(NXG*NYG);
    
    MPI_Bcast(&average,1,MPI_FLOAT,0,MPI_COMM_WORLD);
    
    return average;
};

float matrix_product(float ** matrix1, float **matrix2) {
    
    extern int IDX, IDY,NX,NY;
    
    float local_sum=0.0;
    float global_sum=0.0;
    float buf1=0, buf2=0;
    int i,j;
    
    for (j=1;j<=NY;j=j+IDY){
        for (i=1;i<=NX;i=i+IDX){
            local_sum+=(matrix1[j][i]*matrix2[j][i]);
        }
    }
    
    buf1=local_sum; buf2=0;
    MPI_Allreduce(&buf1,&buf2, 1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    global_sum=buf2;
    
    return global_sum;
}

float ** get_global_from_local_matrix(float ** local_matrix) {
    
    extern int NXG, NYG;
    extern int NX,NY;
    extern int POS[3];
    
    
    float ** global_matrix=NULL,** global_matrix_temp=NULL;
    int i=0,j=0;
    int ii=0, jj=0;
    
    /* Allocate global matrix temp */
    global_matrix_temp=matrix(1,NYG,1,NXG);
    if(global_matrix_temp==NULL) {
        declare_error("Allocation of global_matrix_temp in get_global_from_local_matrix failed!");
    }
    
    /* Allocate global matrix */
    /* You have to deallocate this matrix on our own */
    global_matrix=matrix(1,NYG,1,NXG);
    if(global_matrix==NULL) {
        declare_error("Allocation of global_matrix in get_global_from_local_matrix failed!");
    }
    
    /* Store local matrix in global matrix */
    for (i=1;i<=NXG;i++){
        for (j=1;j<=NYG;j++){
            
            if ( (POS[1]==((i-1)/NX)) && (POS[2]==((j-1)/NY)) ) {
                ii=i-POS[1]*NX;
                jj=j-POS[2]*NY;
                
                global_matrix_temp[j][i]=local_matrix[jj][ii];
            }
        }
    }
    
    MPI_Allreduce(&global_matrix_temp[1][1],&global_matrix[1][1],NXG*NYG,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    
    free_matrix(global_matrix_temp,1,NYG,1,NXG);
    
    return global_matrix;
}

void get_local_from_global_matrix(float ** global_matrix,float ** local_matrix) {
    
    extern int NXG, NYG;
    extern int NX,NY;
    extern int POS[3];
    
    int i=0,j=0;
    int ii=0, jj=0;
    /* Store local matrix in global matrix */
    for (i=1;i<=NXG;i++){
        for (j=1;j<=NYG;j++){
            
            if ( (POS[1]==((i-1)/NX)) && (POS[2]==((j-1)/NY)) ) {
                ii=i-POS[1]*NX;
                jj=j-POS[2]*NY;
                
                local_matrix[jj][ii]=global_matrix[j][i];
            }
        }
    }
    
}

