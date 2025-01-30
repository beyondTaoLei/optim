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
@brief Data format conversion in order to use PLSQR3 software
@details This program is meant to be run to convert the raw binary 
file to non-zero tuple-struct binary file which is defined by PLSQR3 software 
[v1, v2, ...., vn] ---> [(ix1, idx1, v1), (ix2, idx2, v2), ..., 
(ixn, idxn, vn)], which is modified based on the file 
PLSQR3/PLSQR3_tools/kernel_format/ker_ascii2bin.c. We usually calculate 
the kernel from other software and just save the value in raw BIN-format, 
such as SOFI2D and so on. we can convert the BIN-format file to the 
tuple binary file(defined by PLSQR3).

<em>
struct tuple{

    int iX; //not useful here

    int col_idx; //col_idx=iX+(iY-1)*nx+(iZ-1)*nx*ny;

    double val; //val is the kernel value in coordinate of (iX, iY, iZ)
}; 
</em>
@file mpi_ker_bin2nzbin.c
@author Tao Lei (leit@sustech.edu.cn), Wei Zhang
@date 16 Sep 2019
@param[in]  argv[1] the list of kernel in the raw binary file;
@param[in]  argv[2] the directory of kernel in the raw binary file
@param[in]  argv[3] the directory of kernel in the tuple-struct binary file
@param[in]  argv[4] NX: the global number of model grids in X direction (2D);
@param[in]  argv[5] NY: the global number of model grids in Y direction (2D);
@param[in]  argv[6] P1: the percentile(0~1.0) of the absolute of kernel values, 
            except for the zero value. 
@code
mpicc -o mpi_ker_bin2nzbin mpi_ker_bin2nzbin.c -I ../include/ -L ../lib/ -lIN_OUT -lm

mpirun -np 1 ./mpi_ker_bin2nzbin ker.list . aa 200 200 0.1
@endcode
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "comm.h"

/* Feb 9, 2013
 * a.out input_file input_dir output_dir
 *
 * input_file(list: first line: number of rows; second to last line: row file name)
 *
 *
 * mpiexec -np 8 ./ker_ascii_bin_intintdouble_MPI row_list_1.txt /seismic10/hhuang1/elee8_20130204/test_seismicmst /seismic10/hhuang1/elee8_20130204/test_seismicmst
 * */

#define NLEN 1024

int main(int argc, char *argv[]) {

	int rank, mpisize;
	FILE *ifp, *ifpb, *fp_tmp;
	FILE *infile;
	FILE *outfile;
	char kerout[NLEN];
	char kerfnam[NLEN];
        char kerout_bin[NLEN];
        char kerfnam_fullname[NLEN];
        char kerout_fullname[NLEN];
        char filein[NLEN];
	int ker_num, mpi_count;
	int* rank_id;
	int ik, ix, iy, iz, ic, inv, col, col_prev, inv_arr[3], jj;
	double ker;
        
        int size, size2;
        int NX, NY, NZ, nxy,nxyz;
        int INV_VP, INV_VS, INV_RHO, NO_CLASS;
        float btmp, P1=0.1;
        float *vec=NULL,*vec_abs=NULL, *vec_all=NULL;
        

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

        if(argc !=11){
            printf("The comamnd should be as: mpirun -np 1 ./mpi_ker_bin2nzbin ker.list . aa/ 200 200 0.5\n");
            exit(1);
        }else{
            P1 = atof(argv[4]);
            NX = atoi(argv[5]);
            NY = atoi(argv[6]);
            NZ = atoi(argv[7]);
            INV_VP = atoi(argv[8]);
            INV_VS = atoi(argv[9]);
            INV_RHO= atoi(argv[10]);
            
            inv_arr[0]=INV_VP;
            inv_arr[1]=INV_VS;
            inv_arr[2]=INV_RHO;
            NO_CLASS =INV_VP +INV_VS +INV_RHO;
            nxy=NX*NY;
            nxyz=NX*NY*NZ;
            if(rank==0){
                printf("P1: %.3f\n", P1);
                printf("NX, NY, NZ: %d %d %d\n", NX, NY, NZ);
                printf("INV_VP, INV_VS, INV_RHO: %d %d %d\n", INV_VP, INV_VS, INV_RHO);
            }
            
        }
        
        ifp = fopen(argv[1], "r");
        if (ifp == NULL ) {
            fprintf(stderr, "Can't open input file!\n");
            exit(1);
        }
        fscanf(ifp, "%d", &ker_num);
        if(rank==0) printf("%d cores and %d files\n", mpisize, ker_num);

        //***** assign rank ID
        rank_id = (int*) malloc(sizeof(int) * ker_num);
        mpi_count = 0;
        for (ik = 0; ik < ker_num; ik++) {
            rank_id[ik] = mpi_count;
            mpi_count++;
            if (mpi_count == mpisize)
                    mpi_count = 0;
        }

	//  ****** convert ascii to bin ******
	vec=vector(1, nxyz);
        vec_abs=vector(1, nxyz);
        vec_all=vector(1, nxyz*NO_CLASS);
	for (ik = 0; ik < ker_num; ik++) {
            fscanf(ifp, "%s", kerfnam);
            if (rank == rank_id[ik]) {
                //sprintf(kerout, "%s.bin", kerfnam);
                sprintf(kerout, "%s", kerfnam);
                //printf("rank %d : convert %s to %s\n", rank, kerfnam, kerout);

                //full name
                sprintf(kerfnam_fullname, "%s/%s", argv[2], kerfnam);
                sprintf(kerout_fullname, "%s/%s", argv[3], kerout);
                //printf("input: %s\n", kerfnam_fullname);
                printf("output: %s\n", kerout_fullname);
                outfile = fopen(kerout_fullname, "wb");
                if(!outfile){
                    fprintf(stderr, "cannot open %s \n", kerout_fullname);
                }
                ////////////////
                
                infile = fopen(kerfnam_fullname, "r");
                if(!infile){
                    fprintf(stderr, "cannot open %s \n", kerfnam_fullname);
                }else
                    fclose(infile);
                
                /*1. firstly find the appropriate meanful value
                * (the mininum positive value in kernel matrix) */
                V_read(kerfnam_fullname, &vec_all[1], nxyz*NO_CLASS);
                ////////////////
                ic=0;
                for(inv=0;inv<3;inv++){
                    if(inv_arr[inv]==1){/* vp,vs,rho in sequence*/
                        ic=ic+1;
//                         if(inv==0){
//                             strrpc2(kerfnam_fullname,".bin",".vp.bin",filein);
//                         }else if(inv==1){
//                             strrpc2(kerfnam_fullname,".bin",".vs.bin",filein);
//                         }else if(inv==2){
//                             strrpc2(kerfnam_fullname,".bin",".rho.bin",filein);
//                         }
//                         infile = fopen(filein, "r");
//                         if(!infile){
//                             fprintf(stderr, "cannot open %s \n", filein);
//                         }else
//                             fclose(infile);
//                         
//                         /*1. firstly find the appropriate meanful value
//                         * (the mininum positive value in kernel matrix) */
//                         V_read(filein, &vec[1], nxyz);
                        for(jj=1;jj<=nxyz;jj++) vec[jj]=vec_all[jj+(ic-1)*nxyz];
                        V_abs(1, nxyz, vec, vec_abs);
                        btmp=global_percentile(vec_abs, 1, nxyz, P1);
                        //printf("MYID, btmp= %d %e\n", rank, btmp);
                        
                        //=== convert to bin
                        /* For SOFI2D software, Firstly we save the 2D array 
                        * in Y direction, then in X, Z direction, which is 
                        * different from the requirements by PLSQR3
                        */
                        for (iz=1;iz<=NZ;iz++){ 
                        for (ix=1;ix<=NX;ix++){
                        for (iy=1;iy<=NY;iy++){
                            /*1. calculate the index in SOFI2D */
                            /* ii=iy+(ix-1)*NY, which is 1-based index */
                            col_prev=iy+(ix-1)*NY+(iz-1)*nxy;
                            /*2. here col is defined according to the PLSQR3 */
                            col=ix+(iy-1)*NX+(iz-1)*nxy+(ic-1)*nxyz;
                            /*3. write the value */
                            if(fabs(vec[col_prev])>btmp){
                                ker=(double)vec[col_prev];
                                fwrite(&ix, sizeof(int), 1, outfile);/* actually you can set 1 to ix */
                                fwrite(&col, sizeof(int), 1, outfile);
                                fwrite(&ker, sizeof(double), 1, outfile);
                                //printf("%7d %7d %5d %5d %e\n",ii, col, iy, ix, tmp);
                                //printf("%7d %5d %5d %5d %e\n",col, ix, iy, 1, tmp);
                            }
                        }}}
                    }
                }
                fclose(outfile);
                    
            }	  //if
	}	  //for
        free_vector(vec, 1, nxyz);
        free_vector(vec_abs, 1, nxyz);
        free_vector(vec_all, 1, nxyz*NO_CLASS);
	printf("rank %d : Finish!!\n", rank);
	free(rank_id);
        fclose(ifp);
        
        /* generate the ker_bin.list file */
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank ==0){
            sprintf(kerout_bin, "%s/ker_bin.list", argv[3]);
            ifpb= fopen(kerout_bin, "w");
            ifp = fopen(argv[1], "r");
            fscanf(ifp, "%d", &ker_num);
            for (ik = 0; ik < ker_num; ik++) {
                fscanf(ifp, "%s", kerout);
                //full name
                
                sprintf(kerout_fullname, "%s/%s", argv[3], kerout);
                
                fp_tmp = fopen(kerout_fullname,"r");
                if (fp_tmp==NULL) {
                    printf(" Can't read file: %s\n",kerout_fullname);
                    exit(1);
                }
                fseek(fp_tmp,0L,SEEK_END);
                size=ftell(fp_tmp);
                fclose(fp_tmp);
                if(size%16==0 &&size!=0){
                    size2=size/16;
                    fprintf(ifpb,"%s %d\n",kerout,size2);
                    printf("%s %d\n",kerout,size2);
                }else{
                    printf("Non-zero file: %s isn't generated, check para P1\n",kerout_fullname);
                    exit(1);
                }
            } 
            fclose(ifp);
            fclose(ifpb);
        }
        
	MPI_Finalize();
        return 0;
}	  // main

//function: replaces oldstr with newstr in [str] string
void strrpc2(char str[256],char *oldstr,char *newstr, char str_out[256]){
    int i;
    char bstr[strlen(str)];//buffer
    memset(bstr,0,sizeof(bstr));
    memset(str_out,0,strlen(str_out));
 
    for(i = 0; i < strlen(str); i++){
        if(!strncmp(str+i,oldstr,strlen(oldstr))){//search oldstr
            strcat(bstr,newstr);
            i += strlen(oldstr) - 1;
        }else{//save a char into buffer
            strncat(bstr,str + i,1);
        }
    }
    strcpy(str_out, bstr);
}
