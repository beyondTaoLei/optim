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
@brief Reverse the snapshot binary file
@details This program is meant to be run to reverse the snapshot 
        (1:1:NT <--->NT:-1:1) when one calculates the adjoint wavefields 
        with clockwise TDFD operator.
@file mpi_rtime_snap.c
@author Tao Lei (leit@sustech.edu.cn), Wei Zhang
@date 18 Sep 2019
@param[in]  verbose \b 1: print the details; \b 0(\e default): run in quiet pattern;
@param[in]  n1 the size of one frame;
@param[in]  nt the number of frames (\em optional);
@param[in]  fin the snapshot file;
@param[out] fout the reversed snapshot file (\em optional).
@code 
mpirun -np 8 mpi_rtime_snap verbose=1 n1=40000 fin=./snap/inc/test.sp.shot1
@endcode
*/

/*
for example: 
the file list is as,
    test.sp.shot1.0.0.bin (Nodel 0, N0)
    test.sp.shot1.1.0.bin (Nodel 1, N1)
    test.sp.shot1.0.1.bin (Nodel 2, N2)
    test.sp.shot1.1.1.bin (Nodel 3, N3)
then, merge the above files into:
    test.sp.shot1.bin
    
 * In the visualisation below, the rightmost dim (dim[1]) has been chosen as
 * being the one having consecutive elements in memory.
 *
 *                                                 The subarray we
 *               The full array                      want to read
 *
 *       +---------- dim[0] X------>         +--------- dim[0] X------->
 *       | +-----+-----+-----+-----+         | +-----+-----+-----+-----+ 
 *       | |  0  |  4  |   8 |  12 |         | |  0     4 ||  8     12 | 
 *       | +-----+-----+-----+-----+         | |     N0   ||     N1    | 
 *       | |  1  |  5  |   9 |  13 |         | |  1     5 ||  9     13 | 
 * dim[1]| +-----+-----+-----+-----+   dim[1]| +++++++++++++++++++++++++
 *   Y   | |  2  |  6  |  10 |  14 |     Y   | |  2     6 ||  10    14 | 
 *       | +-----+-----+-----+-----+         | |     N2   ||     N3    | 
 *       | |  3  |  7  |  11 |  15 |         | |  3     7 ||  11    15 | 
 *       V +-----+-----+-----+-----+         V +-----------------------+
 *                                            
 *
*/

#include "par.h"
#include "comm.h"
#include"mpi.h"

/* Global variables */
int xargc; 
char **xargv;

int main(int argc, char* argv[]){
/*

*/
    
    /* IN and OUT */
    int verbose, n1, nt;
    char *fin, *fout;
    
    /* Local variables */
    int NP, MYID;
    int size1, it, view_offset, view_offset2;
 
    float *data=NULL;
    
    MPI_File fh_in, fh_out;
    MPI_Status status;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&MYID);
    MPI_Comm_size(MPI_COMM_WORLD,&NP);
    
    /* 1. read formatted input list */
    /* initialize getpar */
    initargs(argc,argv);
    
    verbose=0; getparint("verbose",&verbose);
    if(!getparstring("fin", &fin)) throw_error("fin must be specified!\n");
    if(!getparstring("fout", &fout))  fout= "out.bin";
    if(!getparint("n1",&n1)) throw_error("n1 must be specified!\n");
    if(!getparint("nt",&nt)){
        size1=get_file_size(fin);
        if(size1%(4*n1)==0){
            nt=size1/(4*n1);
        }else{
            throw_error("Make sure the size of snapshot file is correct!\n");
        }
    }
    
    if(MYID==1||NP==1){
        if(verbose){
            printf("\n******\n");
            printf("n1      \t| %d\n", n1);
            printf("nt      \t| %d\n", nt);
            printf("fin     \t| %s\n", fin);
            printf("fout    \t| %s\n", fout);
            printf("******\n");
        }
    }
    
    /* allocate the memory */
    data = vector(1, n1);
    
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_File_open(MPI_COMM_WORLD, fin, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_in);
    MPI_File_open(MPI_COMM_WORLD, fout, MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh_out);
    for(it=1;it<=nt;it++){
        view_offset=(it-1)*n1*sizeof(float);
        view_offset2=(nt-it)*n1*sizeof(float);
        MPI_File_read_at_all(fh_in, view_offset, &data[1], n1, MPI_FLOAT, &status);
        MPI_File_write_at_all(fh_out,view_offset2,&data[1],n1, MPI_FLOAT, &status);
    }
    MPI_File_close(&fh_in);
    MPI_File_close(&fh_out);

    MPI_Finalize();
    return 0;
    
}