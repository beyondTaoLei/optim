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
@brief Merge the model or snapshot because of model decomposition strategy
@details This program is meant to be run for merging the snapshots 
    for IFOS2D software. which uses the domain decomposition 
    strategy. Each processor saves the snapshots in the local disk.
@file mpi_merge_domain.c
@author Tao Lei (leit@sustech.edu.cn), Wei Zhang
@date 18 Sep 2019
@param[in]  verbose \b 1: print the details; \b 0(\e default): run in quiet pattern;
@param[in]  del     \b 1: delete the files in the subdomain; \b 0(\e default): keep the above files;
@param[in]  NPROCX  the number of processors in X direction;
@param[in]  NPROCY  the number of processors in Y direction;
@param[in]  NX      the global number of grids in X direction;
@param[in]  NY      the global number of grids in Y direction;
@param[in]  fin     the prefix of local file in each processors, the final name is
                    as \e fin.[coordinate_of_processor_in_X].[coordinate_of_processor_in_Y].bin;
                    the list of snapshot files could be like:
                    test.sp.shot1.0.0.bin (Nodel 0, N0)
                    test.sp.shot1.1.0.bin (Nodel 1, N1)
                    test.sp.shot1.0.1.bin (Nodel 2, N2)
                    test.sp.shot1.1.1.bin (Nodel 3, N3)
@param[out] fout the file which contains all wavefields from each processors <em>(optional & default: fin.bin)</em>;
@code 
mpirun -np 8 mpi_merge_domain NPROCX=4 NPROCY=2 NX=200 NY=200 fin=./snap/inc/test.sp.shot1 verbose=1 del=1
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
    
    /* IN and OUT */
    int verbose, del, NPROCX, NPROCY, NX, NY;
    char *fin, *fout;
    
    /* Local variables */
    int NP, MYID, POS[3];
    int nx_loc, ny_loc, size1, num_loc, num_glob, nt, it, view_offset;
 
    int sizes[2], subsizes[2], starts[2];
    char fin_a[256], filename[256], fout_a[256];
    float **mat_loc=NULL;
    FILE *myfp=NULL;
    MPI_File mpi_fh;
    MPI_Datatype new_type;
    MPI_Status status;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &NP);
    MPI_Comm_rank(MPI_COMM_WORLD, &MYID);
    
    /* 1. read formatted input list */
    /* initialize getpar */
    initargs(argc,argv);
    
    if(MYID==0){
        verbose=0; getparint("verbose",&verbose);
        del=0; getparint("del",&del);
        if(!getparint("NPROCX",&NPROCX)) throw_error("NPROCX must be specified!\n");
        if(!getparint("NPROCY",&NPROCY)) throw_error("NPROCY must be specified!\n");
        if(!getparint("NX",&NX)) throw_error("NX must be specified!\n");
        if(!getparint("NY",&NY)) throw_error("NY must be specified!\n");
        if(!getparstring("fin", &fin)) throw_error("fin must be specified!\n");
        if(!getparstring("fout", &fout)) fout=fin;
        strcpy(fin_a, fin);
        sprintf(fout_a,"%s.bin",fout);
        
        sprintf(filename,"%s.0.0.bin",fin);
        size1=get_file_size(filename);
        num_loc= NX*NY/(NPROCX*NPROCY);
        if(size1%(4*num_loc)==0){
            nt=size1/(4*num_loc);
            //printf("There are %d snapshots\n", nt);
        }else{
            throw_error("Make sure the number of snapshot in time direction is correct!\n");
        }
        
    }
    
    /* bcast integers */
    MPI_Bcast(&verbose,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&del,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&NPROCX,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&NPROCY,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&NX,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&NY,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&nt,1,MPI_INT,0,MPI_COMM_WORLD);
    /* bcast strings */
    MPI_Bcast(&fin_a,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(&fout_a,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    
    /*1. preparation work */
    if(NP!=NPROCX*NPROCY)
        throw_error("Make sure the number of processors is correct!\n");
    else if(NX%NPROCX==0&&NY%NPROCY==0){
        POS[1] = MYID%NPROCX;       /*  x coordinate */
	POS[2] = (MYID/NPROCX);     /*  y coordinate */
        nx_loc=NX/NPROCX;
        ny_loc=NY/NPROCY;
        num_loc=nx_loc*ny_loc;
        num_glob=NX*NY;
        
    }else{
        throw_error("Make sure the number of processors in X and Y direction is correct!\n");
    }
    
    if(MYID==1||NP==1){
        //printf("\nMYID=%d, NP=%d\n",MYID, NP);
        if(verbose){
            printf("\n******\n");
            printf("del         \t| %d\n",del);
            printf("NPROCX      \t| %d\n",NPROCX);
            printf("NPROCY      \t| %d\n",NPROCY);
            printf("NX          \t| %d\n",NX);
            printf("NY          \t| %d\n",NY);
            printf("nx_loc      \t| %d\n",nx_loc);
            printf("ny_loc      \t| %d\n",ny_loc);
            printf("nt          \t| %d\n",nt);
            printf("fin_a       \t| %s\n",fin_a);
            printf("fout_a      \t| %s\n",fout_a);
            printf("******\n");
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    /* allocate the memory */
    mat_loc = matrix(1, nx_loc, 1, ny_loc);
    
    /*open the files storing the data in distributed system */
    sprintf(filename,"%s.%i.%i.bin",fin_a,POS[1],POS[2]);
    //sprintf(fout_a,"%s.bin",fin_a);
    
    myfp=fopen(filename,"r");
    
    sizes[0] = NX;
    sizes[1] = NY;
    subsizes[0] = nx_loc;
    subsizes[1] = ny_loc;
    starts[0] = POS[1]*nx_loc;
    starts[1] = POS[2]*ny_loc;
    //printf("myid, filename =%d, %s\n",MYID, filename);
    //printf("myid, fout_a=%d, %s\n",MYID, fout_a);
    //printf("myid:%d %d %d %d %d %d %d\n",MYID, sizes[0], sizes[1],subsizes[0],subsizes[1],starts[0],starts[1]);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_FLOAT, &new_type);
    MPI_Type_commit(&new_type);
    // Writes the array mat_loc using a view
    MPI_File_open(MPI_COMM_WORLD,fout_a,MPI_MODE_WRONLY | MPI_MODE_CREATE,MPI_INFO_NULL, &mpi_fh);
    for(it=1;it<=nt;it++){
        /*read the data */
        fread(&mat_loc[1][1], sizeof(float), num_loc, myfp);
        /* write data into disk */
        view_offset=(it-1)*num_glob*sizeof(float);
        MPI_File_set_view(mpi_fh, view_offset, MPI_FLOAT, new_type, "native", MPI_INFO_NULL);
        MPI_File_write_all(mpi_fh, &mat_loc[1][1], num_loc, MPI_FLOAT, &status);
    }
    
    MPI_File_close(&mpi_fh);
    fclose(myfp);
    
    /*delete the file in the subdomain */
    if(del==1) remove(filename);
    
    /* deallocate the memory */
    free_matrix(mat_loc, 1, nx_loc, 1, ny_loc);
    
    
    
    MPI_Finalize();
    return 0;
    
}